"""
Step 2: ScoreSP Path Energy Scoring
=====================================
Compute ScoreSP = E_path / R_path for each candidate pair
using knowledge graph path analysis.
"""

import json
import math
from pathlib import Path
from typing import Optional, List, Dict, Any, Tuple

import networkx as nx


CURRENT_YEAR = 2026


def load_knowledge_graph(gml_path: str) -> nx.DiGraph:
    """Load GML knowledge graph."""
    G = nx.read_gml(gml_path)
    print(f"Loaded KG: {G.number_of_nodes()} nodes, {G.number_of_edges()} edges")
    return G


def find_node_id(G: nx.DiGraph, name: str, node_type: Optional[str] = None) -> Optional[str]:
    """
    Find node ID by name matching. Uses two passes:
    1. Exact match (case/space/dash insensitive)
    2. Synonym matching (acid/ate, common name variants)
    Avoids substring matching to prevent false positives.
    """
    def normalize(s: str) -> str:
        return s.lower().strip().replace(' ', '_').replace('-', '_')

    def _acid_ate_variants(n: str) -> list:
        """Generate acid/conjugate-base name variants dynamically."""
        variants = []
        if n.endswith('_acid'):
            stem = n[:-5]  # remove '_acid'
            # propionic_acid -> propionate, butyric_acid -> butyrate, etc.
            if stem.endswith('ic'):
                variants.append(stem[:-2] + 'ate')
            elif stem.endswith('ous'):
                variants.append(stem[:-3] + 'ite')
        elif n.endswith('ate'):
            stem = n[:-3]
            variants.append(stem + 'ic_acid')
        elif n.endswith('ite'):
            stem = n[:-3]
            variants.append(stem + 'ous_acid')
        return variants

    name_norm = normalize(name)

    # Pass 1: exact match on node ID or label
    for node_id, attrs in G.nodes(data=True):
        if node_type and attrs.get('type', '') != node_type:
            continue
        nid_norm = normalize(node_id)
        label_norm = normalize(attrs.get('label', ''))

        if nid_norm == name_norm or (label_norm and label_norm == name_norm):
            return node_id

    # Pass 2: dynamic acid/salt synonym match
    for variant in _acid_ate_variants(name_norm):
        for node_id, attrs in G.nodes(data=True):
            if node_type and attrs.get('type', '') != node_type:
                continue
            nid_norm = normalize(node_id)
            label_norm = normalize(attrs.get('label', ''))
            if nid_norm == variant or (label_norm and label_norm == variant):
                return node_id

    return None


def extract_khop_subgraph(G: nx.DiGraph, seed_nodes: List[str], max_hop: int = 2) -> nx.DiGraph:
    """
    Extract K-hop subgraph from seed nodes.
    """
    visited = set()
    frontier = set(seed_nodes)

    for _ in range(max_hop):
        next_frontier = set()
        for node in frontier:
            if node in visited:
                continue
            visited.add(node)
            if node in G:
                next_frontier.update(G.successors(node))
                next_frontier.update(G.predecessors(node))
        frontier = next_frontier - visited

    visited.update(frontier)
    return G.subgraph(visited).copy()


def find_paths(
    G,
    source: str,
    target: str,
    max_depth: int = 4,
    exclude_types: Optional[List[str]] = None,
) -> List[List[str]]:
    """
    Find simple paths between source and target (undirected search).
    Excludes intermediate nodes of specified types (e.g. 'disease') to avoid
    biologically meaningless paths through hub disease nodes.
    Returns list of unique node-id paths, deduplicated.
    """
    if source not in G or target not in G:
        return []

    if exclude_types is None:
        exclude_types = ['disease']

    # Build filtered subgraph: remove disease-type intermediate nodes
    # but keep source and target even if they are disease nodes
    nodes_to_keep = set()
    for n, d in G.nodes(data=True):
        ntype = d.get('type', '')
        if n == source or n == target:
            nodes_to_keep.add(n)
        elif ntype not in exclude_types:
            nodes_to_keep.add(n)

    G_filtered = G.subgraph(nodes_to_keep)

    # Search on undirected view for reachability
    G_undirected = G_filtered.to_undirected(as_view=True)

    # For MultiGraph, convert to simple Graph to avoid duplicate paths
    if isinstance(G_undirected, nx.MultiGraph):
        G_simple = nx.Graph(G_undirected)
    else:
        G_simple = G_undirected

    try:
        paths = list(nx.all_simple_paths(G_simple, source, target, cutoff=max_depth))
        # Deduplicate (same node sequence)
        seen = set()
        unique_paths = []
        for p in paths:
            key = tuple(p)
            if key not in seen:
                seen.add(key)
                unique_paths.append(p)
        # Sort by length (shortest first)
        unique_paths.sort(key=len)
        return unique_paths[:5]
    except (nx.NetworkXError, nx.NodeNotFound):
        return []
    except (nx.NetworkXError, nx.NodeNotFound):
        return []


def compute_node_weight(G, node_id: str) -> Tuple[float, float, float]:
    """
    Compute node weight w_i = Novelty_i * (1 + Cred_i).

    Returns: (w_i, novelty, cred)
    """
    if node_id not in G:
        return (0.0, 0.0, 0.0)

    # Novelty: 1 if terminal node (out_degree=0), else 0
    # Per spec: only out_degree matters
    out_deg = G.out_degree(node_id)
    novelty = 1.0 if out_deg == 0 else 0.0

    # Credibility: average impact_factor of connected edges
    ifs = []
    # Outgoing edges
    if isinstance(G, (nx.MultiDiGraph, nx.MultiGraph)):
        for _, target, key, edata in G.edges(node_id, data=True, keys=True):
            if_val = edata.get('edge_impact_factor', 0)
            if if_val:
                ifs.append(float(if_val))
        for source, _, key, edata in G.in_edges(node_id, data=True, keys=True):
            if_val = edata.get('edge_impact_factor', 0)
            if if_val:
                ifs.append(float(if_val))
    else:
        for _, _, edata in G.edges(node_id, data=True):
            if_val = edata.get('edge_impact_factor', 0)
            if if_val:
                ifs.append(float(if_val))
        for _, _, edata in G.in_edges(node_id, data=True):
            if_val = edata.get('edge_impact_factor', 0)
            if if_val:
                ifs.append(float(if_val))

    cred = sum(ifs) / len(ifs) if ifs else 1.0

    w_i = novelty * (1 + cred)
    return (w_i, novelty, cred)


def _get_edge_data_list(G, source: str, target: str) -> List[Dict[str, Any]]:
    """Get all edge data dicts between source and target (handles MultiDiGraph)."""
    results = []
    if G.has_edge(source, target):
        if isinstance(G, nx.MultiDiGraph) or isinstance(G, nx.MultiGraph):
            for key in G[source][target]:
                results.append(dict(G[source][target][key]))
        else:
            results.append(dict(G.edges[source, target]))
    if G.has_edge(target, source):
        if isinstance(G, nx.MultiDiGraph) or isinstance(G, nx.MultiGraph):
            for key in G[target][source]:
                results.append(dict(G[target][source][key]))
        else:
            results.append(dict(G.edges[target, source]))
    return results


def compute_edge_resistance(G, source: str, target: str) -> Tuple[float, float, float]:
    """
    Compute edge resistance r_e = 1 / (1 + Support_e).

    Support_e = sum(RefImp) for all references on edge.
    RefImp = IF * citations / (delta_year + 1).

    Returns: (r_e, support_e, ref_imp)
    """
    edge_data_list = _get_edge_data_list(G, source, target)
    if not edge_data_list:
        return (1.0, 0.0, 0.0)

    # Sum RefImp across all parallel edges (multiple references)
    total_ref_imp = 0.0
    for edata in edge_data_list:
        impact_factor = float(edata.get('edge_impact_factor', 0))
        citation_count = float(edata.get('edge_citation_count', 0))
        pub_year = edata.get('publication_year', None)

        if pub_year:
            delta_year = max(CURRENT_YEAR - int(pub_year), 0)
        else:
            delta_year = 5

        if impact_factor and citation_count:
            total_ref_imp += impact_factor * citation_count / (delta_year + 1)

    support_e = total_ref_imp
    r_e = 1.0 / (1.0 + support_e)

    return (r_e, support_e, total_ref_imp)


def compute_scoresp_for_path(
    G: nx.DiGraph,
    path_nodes: List[str],
) -> Dict[str, Any]:
    """
    Compute ScoreSP for a single path.

    Returns dict with score details.
    """
    if len(path_nodes) < 2:
        return {'score_sp': 0.0, 'e_path': 0.0, 'r_path': 0.0, 'synergy': 0.0,
                'nodes_detail': [], 'edges_detail': []}

    # Node weights
    nodes_detail = []
    novel_count = 0
    sum_wi_sq = 0.0
    for nid in path_nodes:
        w_i, novelty, cred = compute_node_weight(G, nid)
        sum_wi_sq += w_i ** 2
        if novelty > 0:
            novel_count += 1
        nodes_detail.append({
            'name': nid,
            'type': G.nodes[nid].get('type', 'unknown') if nid in G else 'unknown',
            'novelty': novelty,
            'cred': round(cred, 2),
            'w_i': round(w_i, 4),
        })

    # Synergy = sqrt(m) where m = number of novel nodes; 0 if no novel nodes
    synergy = math.sqrt(novel_count) if novel_count > 0 else 0.0

    # E_path
    e_path = synergy * sum_wi_sq

    # Edge resistances
    edges_detail = []
    r_path = 0.0
    for i in range(len(path_nodes) - 1):
        src, tgt = path_nodes[i], path_nodes[i + 1]
        r_e, support_e, ref_imp = compute_edge_resistance(G, src, tgt)
        r_path += r_e
        # Collect PMIDs from edge data
        pmids = []
        for edata in _get_edge_data_list(G, src, tgt):
            pmid = edata.get('pmid', '')
            if pmid:
                pmids.append(str(pmid))
        edges_detail.append({
            'source': src,
            'target': tgt,
            'ref_imp': round(ref_imp, 4),
            'support': round(support_e, 4),
            'r_e': round(r_e, 4),
            'pmids': pmids,
        })

    # ScoreSP = E_path / R_path
    # When E_path=0 but path exists (no terminal nodes), use 1/R_path as
    # fallback score so paths with strong evidence still rank above no-path
    if e_path > 0 and r_path > 0:
        score_sp = e_path / r_path
    elif r_path > 0:
        # Path exists but no novel (terminal) nodes: use inverse resistance
        score_sp = 1.0 / r_path
    else:
        score_sp = 0.0

    return {
        'score_sp': round(score_sp, 4),
        'e_path': round(e_path, 4),
        'r_path': round(r_path, 4),
        'synergy': round(synergy, 4),
        'nodes_detail': nodes_detail,
        'edges_detail': edges_detail,
    }


def compute_scoresp_for_candidate(
    G: nx.DiGraph,
    bacteria: str,
    metabolite: str,
    max_hop: int = 2,
    max_path_depth: int = 4,
) -> Dict[str, Any]:
    """
    Compute ScoreSP for a candidate bacteria-metabolite pair.
    """
    # Find node IDs in graph
    bact_id = find_node_id(G, bacteria, 'microbe')
    meta_id = find_node_id(G, metabolite, 'metabolite')

    result = {
        'bacteria_in_kg': bact_id is not None,
        'metabolite_in_kg': meta_id is not None,
        'has_path': False,
        'score_sp': 0.0,
        'e_path': 0.0,
        'r_path': 0.0,
        'synergy': 0.0,
        'path_nodes': [],
        'path_edges': [],
        'path_str': '',
        'all_paths': [],
    }

    if not bact_id or not meta_id:
        return result

    # Extract K-hop subgraph
    subgraph = extract_khop_subgraph(G, [bact_id, meta_id], max_hop)

    # Find paths between bacteria and metabolite only
    paths = find_paths(subgraph, bact_id, meta_id, max_path_depth)

    if not paths:
        return result

    # Compute ScoreSP for each path, keep the best
    best_score = -1
    best_path_result = None
    all_paths_info = []

    for path in paths:
        path_result = compute_scoresp_for_path(G, path)
        all_paths_info.append({
            'path': ' -> '.join(path),
            'score_sp': path_result['score_sp'],
        })
        if path_result['score_sp'] > best_score:
            best_score = path_result['score_sp']
            best_path_result = path_result
            best_path_result['path_str'] = ' -> '.join(path)

    if best_path_result and best_score > 0:
        result.update({
            'has_path': True,
            'score_sp': best_path_result['score_sp'],
            'e_path': best_path_result['e_path'],
            'r_path': best_path_result['r_path'],
            'synergy': best_path_result['synergy'],
            'path_nodes': best_path_result['nodes_detail'],
            'path_edges': best_path_result['edges_detail'],
            'path_str': best_path_result['path_str'],
            'all_paths': all_paths_info,
        })

    return result


def run_step2(
    candidates_path: str,
    gml_path: str,
    output_path: str,
    subgraph_output_path: Optional[str] = None,
    max_hop: int = 2,
) -> list:
    """
    Run step 2: compute ScoreSP for all candidates.

    Args:
        candidates_path: Path to step1 output JSON
        gml_path: Path to knowledge graph GML
        output_path: Path to save step2 output JSON
        subgraph_output_path: Optional path to save merged subgraph GML
        max_hop: K-hop neighborhood size

    Returns:
        List of scored candidate dicts
    """
    print("=" * 60)
    print("  STEP 2: ScoreSP Path Energy Scoring")
    print("=" * 60)

    # Load inputs
    with open(candidates_path, 'r', encoding='utf-8') as f:
        candidates = json.load(f)

    G = load_knowledge_graph(gml_path)

    # Process each candidate
    results = []
    for i, cand in enumerate(candidates):
        scoresp_result = compute_scoresp_for_candidate(
            G, cand['bacteria'], cand['metabolite'], max_hop
        )
        results.append({
            'candidate': cand,
            **scoresp_result,
        })

        if (i + 1) % 50 == 0:
            print(f"  Processed {i + 1}/{len(candidates)} candidates")

    # Summary
    has_path_count = sum(1 for r in results if r['has_path'])
    print(f"\nResults: {has_path_count}/{len(results)} candidates have KG paths")
    print(f"ScoreSP range: {min(r['score_sp'] for r in results):.4f} - "
          f"{max(r['score_sp'] for r in results):.4f}")

    # Save results
    output_file = Path(output_path)
    output_file.parent.mkdir(parents=True, exist_ok=True)
    with open(output_file, 'w', encoding='utf-8') as f:
        json.dump(results, f, indent=2, ensure_ascii=False)

    print(f"Saved {len(results)} scored candidates to {output_path}")

    # Optionally save merged subgraph
    if subgraph_output_path:
        # Collect all nodes mentioned in paths
        all_path_nodes = set()
        for r in results:
            for n in r.get('path_nodes', []):
                all_path_nodes.add(n['name'])

        if all_path_nodes:
            subgraph = G.subgraph(all_path_nodes).copy()
            nx.write_gml(subgraph, subgraph_output_path)
            print(f"Saved subgraph ({subgraph.number_of_nodes()} nodes) to {subgraph_output_path}")

    return results


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='Step 2: ScoreSP Energy Scoring')
    parser.add_argument('--candidates', default='outputs/step1_candidates.json')
    parser.add_argument('--gml', default='data/knowledge_graph/ibd_test_kg.gml')
    parser.add_argument('--output', default='outputs/step2_scoresp.json')
    parser.add_argument('--subgraph-output', default='outputs/step2_subgraph.gml')
    parser.add_argument('--max-hop', type=int, default=2)
    args = parser.parse_args()

    run_step2(args.candidates, args.gml, args.output, args.subgraph_output, args.max_hop)
