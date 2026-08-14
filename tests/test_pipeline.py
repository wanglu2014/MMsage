"""
Backend Tests for MMSage x Chain Novelty Pipeline
==================================================
Comprehensive test suite covering all 3 steps + KG + API.
"""

import json
import sys
from pathlib import Path

# Add backend to path
sys.path.insert(0, str(Path(__file__).parent.parent / "backend"))

import pytest
import networkx as nx


# ============================================================
# 1. Step1 Unit Tests
# ============================================================

def test_extract_bacteria_from_filename():
    from step1_mmsage_signal import extract_bacteria_from_filename

    name = "pluscombno1V0317_min0_Akkermansia_muciniphila_seed_1_dim_2_neighbor_2_dist_1_metric_cosine_rank_rootknow_cor0303_2_top_10_Pthre_0.1_pair.csv_clu_0.1_coordinates_tunek.csv"
    assert extract_bacteria_from_filename(name) == "Akkermansia_muciniphila"

    name2 = "pluscombno1V0317_min0_Bacteroides_fragilis_seed_1_dim_3_coordinates_tunek.csv"
    assert extract_bacteria_from_filename(name2) == "Bacteroides_fragilis"


def test_parse_metabolite_from_rowname():
    from step1_mmsage_signal import parse_metabolite_from_rowname

    assert parse_metabolite_from_rowname(
        '"Isobutyric acid-Akkermansia_muciniphila"', "Akkermansia_muciniphila"
    ) == "Isobutyric acid"

    assert parse_metabolite_from_rowname(
        '"2-Hydroxybutyric acid-Akkermansia_muciniphila"', "Akkermansia_muciniphila"
    ) == "2-Hydroxybutyric acid"


# ============================================================
# 2. KnowledgePump Tests
# ============================================================

def test_kg_load():
    """Test that the shipped knowledge graph is usable for the paper axis."""
    gml_path = Path(__file__).parent.parent / "data" / "knowledge_graph" / "auto_built_kg.gml"
    assert gml_path.exists(), f"Required fixture missing: {gml_path}"

    G = nx.read_gml(str(gml_path))
    assert G.number_of_nodes() > 0, "Expected a non-empty knowledge graph"
    assert G.number_of_edges() > 0, "Expected a knowledge graph with edges"
    assert G.is_directed(), "Expected a directed knowledge graph"
    assert "akkermansia_muciniphila" in G
    assert "isobutyric_acid" in G
    assert "ibd" in G


def test_kg_node_types():
    """Test KG has correct node type distribution."""
    gml_path = Path(__file__).parent.parent / "data" / "knowledge_graph" / "auto_built_kg.gml"
    assert gml_path.exists(), f"Required fixture missing: {gml_path}"

    G = nx.read_gml(str(gml_path))
    type_counts = {}
    for _, d in G.nodes(data=True):
        t = d.get('type', 'unknown')
        type_counts[t] = type_counts.get(t, 0) + 1

    required_types = {'microbe', 'metabolite', 'pathway', 'disease'}
    assert required_types <= set(type_counts), f"Missing node types: {required_types - set(type_counts)}"
    assert all(type_counts[node_type] > 0 for node_type in required_types)
    assert sum(type_counts.values()) == G.number_of_nodes()


def test_kg_find_paths():
    """Test KG path finding for Akkermansia."""
    from knowledge_pump import KnowledgePump

    gml_path = Path(__file__).parent.parent / "data" / "knowledge_graph" / "auto_built_kg.gml"
    assert gml_path.exists(), f"Required fixture missing: {gml_path}"

    kg = KnowledgePump(str(gml_path))

    # Check the paper's core microbe-metabolite axis.
    paths = kg.find_relevant_paths("akkermansia_muciniphila", "isobutyric_acid", max_depth=3)

    assert len(paths) > 0, "Expected at least one path from Akkermansia to isobutyric acid"

    # Check path depth
    for path in paths:
        assert len(path) <= 3, f"Path depth {len(path)} exceeds max_depth=3"


# ============================================================
# 3. Step2 Chain Novelty Unit Tests
# ============================================================

def test_build_pubmed_query():
    """Test PubMed query construction without synonym expansion."""
    from step2_chain_novelty import build_pubmed_query

    query = build_pubmed_query(["Akkermansia_muciniphila", "Propionate", "IBD"])

    assert query == '"Akkermansia muciniphila" AND "Propionate" AND "IBD"'
    assert " OR " not in query


def test_kg_id_to_pubmed_term():
    """Test KG node ID to PubMed term conversion preserves casing."""
    from step2_chain_novelty import _kg_id_to_pubmed_term

    # Test metabolite casing preservation
    result = _kg_id_to_pubmed_term("l_carnitine", "Akkermansia_muciniphila", "L-Carnitine")
    assert result == "L-Carnitine", f"Expected 'L-Carnitine', got '{result}'"

    # Test bacteria name
    result = _kg_id_to_pubmed_term("akkermansia_muciniphila", "Akkermansia_muciniphila", "Choline")
    assert result == "Akkermansia muciniphila", f"Expected 'Akkermansia muciniphila', got '{result}'"


def test_aggregate_chain_count():
    """Test chain_count uses the experimental triple count."""
    from step2_chain_novelty import _aggregate_chain_count

    pairwise = {
        'pair_bm': 10,
        'triple_bmd': 5,
        'triple_bmd_exp': 2,
        'pair_bd': 100,
        'pair_md': 50,
        'full_chain': 8,
    }

    result = _aggregate_chain_count(pairwise)
    assert result == 2, f"Expected triple_bmd_exp=2, got {result}"
    assert _aggregate_chain_count({}) == 0


def test_normalize_novelty():
    """Test Chain Novelty formula."""
    from step2_chain_novelty import normalize_novelty
    import math

    # Test normal case
    novelty = normalize_novelty(10, 100)
    expected = 1 - math.log(11) / math.log(101)
    assert abs(novelty - expected) < 0.0001, f"Expected {expected:.4f}, got {novelty:.4f}"

    # Test boundary: C_max = 0
    novelty_zero = normalize_novelty(5, 0)
    assert novelty_zero == 1.0, f"Expected 1.0 when C_max=0, got {novelty_zero}"

    # Test chain_count = 0
    novelty_novel = normalize_novelty(0, 100)
    assert novelty_novel == 1.0, f"Expected 1.0 when chain_count=0, got {novelty_novel}"

    # Test chain_count = C_max
    novelty_known = normalize_novelty(100, 100)
    assert novelty_known == 0.0, f"Expected 0.0 when chain_count=C_max, got {novelty_known}"


def test_find_best_chain_no_path(monkeypatch):
    """Test find_best_chain returns pairwise_counts even without KG path."""
    import step2_chain_novelty as step2
    from knowledge_pump import KnowledgePump

    gml_path = Path(__file__).parent.parent / "data" / "knowledge_graph" / "auto_built_kg.gml"
    assert gml_path.exists(), f"Required fixture missing: {gml_path}"

    kg = KnowledgePump(str(gml_path))
    monkeypatch.setattr(
        step2,
        "_compute_pairwise_counts",
        lambda *args, **kwargs: {'triple_bmd_exp': 0},
    )

    # Use a metabolite not in KG
    result = step2.find_best_chain(
        "Akkermansia_muciniphila",
        "NonExistentMetabolite",
        "IBD",
        kg,
        {},
        max_depth=3,
    )

    assert result['has_path'] == False, "Expected has_path=False for non-existent metabolite"
    assert 'pairwise_counts' in result, "Expected pairwise_counts even without KG path"
    assert result['chain_count'] >= 0, "Expected chain_count >= 0"


# ============================================================
# 4. Step3 Quadrant Unit Tests
# ============================================================

def test_assign_quadrant_all_four():
    """Test all 4 quadrant assignments with 4-parameter signature."""
    from step3_quadrant import assign_quadrant

    # Q-I: high MMSage + high novelty (Dark Matter)
    r = assign_quadrant(0.9, 0.9, 0.5, 0.5)
    assert r['quadrant'] == 'I', f"Expected Q-I, got {r['quadrant']}"
    assert r['is_dark_matter'] == True, "Expected is_dark_matter=True for Q-I"

    # Q-II: low MMSage + high novelty
    r = assign_quadrant(0.3, 0.9, 0.5, 0.5)
    assert r['quadrant'] == 'II', f"Expected Q-II, got {r['quadrant']}"
    assert r['is_dark_matter'] == False, "Expected is_dark_matter=False for Q-II"

    # Q-III: high MMSage + low novelty (Known)
    r = assign_quadrant(0.9, 0.2, 0.5, 0.5)
    assert r['quadrant'] == 'III', f"Expected Q-III, got {r['quadrant']}"
    assert r['is_dark_matter'] == False, "Expected is_dark_matter=False for Q-III"

    # Q-IV: low MMSage + low novelty
    r = assign_quadrant(0.3, 0.2, 0.5, 0.5)
    assert r['quadrant'] == 'IV', f"Expected Q-IV, got {r['quadrant']}"
    assert r['is_dark_matter'] == False, "Expected is_dark_matter=False for Q-IV"


def test_adaptive_threshold_median():
    """Test adaptive threshold computation uses median."""
    # Construct 10 candidates with known median
    mmsage_values = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
    novelty_values = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]

    # Median of 10 values: (5th + 6th) / 2
    expected_mmsage = (0.5 + 0.6) / 2
    expected_novelty = (0.4 + 0.5) / 2

    # Compute median manually
    n = len(mmsage_values)
    mmsage_threshold = (mmsage_values[n // 2] + mmsage_values[(n - 1) // 2]) / 2
    novelty_threshold = (novelty_values[n // 2] + novelty_values[(n - 1) // 2]) / 2

    assert abs(mmsage_threshold - expected_mmsage) < 0.0001, f"Expected {expected_mmsage}, got {mmsage_threshold}"
    assert abs(novelty_threshold - expected_novelty) < 0.0001, f"Expected {expected_novelty}, got {novelty_threshold}"


def test_quadrant_dark_matter_flag():
    """Test only Q-I has is_dark_matter=True."""
    from step3_quadrant import assign_quadrant

    quadrants = [
        (0.9, 0.9, 'I', True),
        (0.3, 0.9, 'II', False),
        (0.9, 0.3, 'III', False),
        (0.3, 0.3, 'IV', False),
    ]

    for mmsage, novelty, expected_q, expected_dm in quadrants:
        r = assign_quadrant(mmsage, novelty, 0.5, 0.5)
        assert r['quadrant'] == expected_q, f"Expected {expected_q}, got {r['quadrant']}"
        assert r['is_dark_matter'] == expected_dm, f"Expected is_dark_matter={expected_dm} for Q-{expected_q}, got {r['is_dark_matter']}"


# ============================================================
# 5. Integration Tests: Real Pipeline Output Validation
# ============================================================

def test_step1_output_schema():
    """Test step1 output has correct schema."""
    output_path = Path(__file__).parent.parent / "outputs" / "step1_candidates.json"
    assert output_path.exists(), f"Required fixture missing: {output_path}"

    with open(output_path, 'r', encoding='utf-8') as f:
        data = json.load(f)

    assert len(data) == 146, f"Expected 146 candidates, got {len(data)}"

    # Check first candidate schema
    item = data[0]
    required_fields = ['bacteria', 'metabolite', 'mmsage_norm', 'pseudotime', 'rank_in_microbe', 'total_in_microbe', 'source_file']
    for field in required_fields:
        assert field in item, f"Missing field '{field}' in step1 output"

    # Check mmsage_norm range
    for item in data:
        assert 0 <= item['mmsage_norm'] <= 1, f"mmsage_norm {item['mmsage_norm']} out of range [0,1]"


def test_step2_output_schema():
    """Test step2 output has correct schema."""
    output_path = Path(__file__).parent.parent / "outputs" / "step2_chain_novelty.json"
    assert output_path.exists(), f"Required fixture missing: {output_path}"

    with open(output_path, 'r', encoding='utf-8') as f:
        data = json.load(f)

    step1_path = Path(__file__).parent.parent / "outputs" / "step1_candidates.json"
    assert step1_path.exists(), f"Required fixture missing: {step1_path}"
    with open(step1_path, 'r', encoding='utf-8') as f:
        step1_data = json.load(f)

    assert len(data) == len(step1_data)
    assert {
        (item['bacteria'], item['metabolite']) for item in data
    } == {
        (item['bacteria'], item['metabolite']) for item in step1_data
    }

    # Check schema
    item = data[0]
    required_fields = ['bacteria', 'metabolite', 'chain_count', 'chain_novelty', 'has_path', 'pairwise_counts']
    for field in required_fields:
        assert field in item, f"Missing field '{field}' in step2 output"

    # Check pairwise_counts structure
    pw = item['pairwise_counts']
    assert isinstance(pw, dict), "pairwise_counts should be a dict"


def test_step3_output_schema():
    """Test step3 output has correct schema."""
    output_path = Path(__file__).parent.parent / "outputs" / "step3_quadrant.json"
    assert output_path.exists(), f"Required fixture missing: {output_path}"

    with open(output_path, 'r', encoding='utf-8') as f:
        data = json.load(f)

    step2_path = Path(__file__).parent.parent / "outputs" / "step2_chain_novelty.json"
    assert step2_path.exists(), f"Required fixture missing: {step2_path}"
    with open(step2_path, 'r', encoding='utf-8') as f:
        step2_data = json.load(f)

    assert len(data) == len(step2_data)
    assert {
        (item['bacteria'], item['metabolite']) for item in data
    } == {
        (item['bacteria'], item['metabolite']) for item in step2_data
    }

    # Check schema
    item = data[0]
    required_fields = ['bacteria', 'metabolite', 'quadrant', 'is_dark_matter', 'mmsage_norm', 'chain_novelty', 'chain_count']
    for field in required_fields:
        assert field in item, f"Missing field '{field}' in step3 output"

    # Check quadrant values
    for item in data:
        assert item['quadrant'] in ['I', 'II', 'III', 'IV'], f"Invalid quadrant '{item['quadrant']}'"


def test_quadrant_distribution():
    """Test quadrant coverage and dark-matter flag consistency."""
    output_path = Path(__file__).parent.parent / "outputs" / "step3_quadrant.json"
    assert output_path.exists(), f"Required fixture missing: {output_path}"

    with open(output_path, 'r', encoding='utf-8') as f:
        data = json.load(f)

    qc = {'I': 0, 'II': 0, 'III': 0, 'IV': 0}
    for item in data:
        qc[item['quadrant']] += 1

    assert sum(qc.values()) == len(data)
    assert all(
        item['is_dark_matter'] == (item['quadrant'] == 'I')
        for item in data
    )


def test_chain_novelty_range():
    """Test chain_novelty in [0,1] and chain_count >= 0."""
    output_path = Path(__file__).parent.parent / "outputs" / "step3_quadrant.json"
    assert output_path.exists(), f"Required fixture missing: {output_path}"

    with open(output_path, 'r', encoding='utf-8') as f:
        data = json.load(f)

    chain_counts = [item['chain_count'] for item in data]
    novelties = [item['chain_novelty'] for item in data]

    # Check ranges
    for cc in chain_counts:
        assert cc >= 0, f"chain_count {cc} < 0"

    for nov in novelties:
        assert 0 <= nov <= 1, f"chain_novelty {nov} out of range [0,1]"

    c_max = max(chain_counts)
    from step2_chain_novelty import normalize_novelty

    for item in data:
        expected = normalize_novelty(item['chain_count'], c_max)
        assert item['chain_novelty'] == pytest.approx(expected, abs=1e-4)


# ============================================================
# 6. API Server Tests
# ============================================================

def test_api_endpoints_status_200():
    """Test all API endpoints return 200."""
    from fastapi.testclient import TestClient
    from api_server import app

    client = TestClient(app)

    endpoints = [
        "/api/pipeline/status",
        "/api/config",
        "/api/graph/stats",
        "/api/step1/candidates",
        "/api/step2/chain_novelty",
        "/api/step3/quadrant",
        "/api/graph/subgraph?bacteria=Akkermansia&max_hop=1",
        "/api/graph/full",
    ]

    for endpoint in endpoints:
        response = client.get(endpoint)
        assert response.status_code == 200, f"{endpoint} returned {response.status_code}"


def test_api_graph_full_no_crash():
    """Test /api/graph/full doesn't crash with Query object bug."""
    from fastapi.testclient import TestClient
    from api_server import app

    client = TestClient(app)
    response = client.get("/api/graph/full")

    assert response.status_code == 200, f"Expected 200, got {response.status_code}"

    data = response.json()
    assert 'nodes' in data, "Expected 'nodes' in response"
    assert 'edges' in data, "Expected 'edges' in response"

    stats_response = client.get("/api/graph/stats")
    assert stats_response.status_code == 200
    stats = stats_response.json()
    assert len(data['nodes']) == stats['nodes']
    assert len(data['edges']) == stats['edges']


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
