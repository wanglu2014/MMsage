import csv
from pathlib import Path
from typing import Dict, List, Optional

try:
    from .sapbert_retriever import retrieve
except ImportError:
    from sapbert_retriever import retrieve

CHEBI_PATH = Path(__file__).parent.parent.parent / "data" / "databases" / "chebi_relationships_standardized.csv"

_chebi_data: Optional[List[Dict]] = None


def _load_chebi() -> List[Dict]:
    """Load and normalize ChEBI CSV records with caching."""
    global _chebi_data
    if _chebi_data is not None:
        return _chebi_data

    _chebi_data = []
    if not CHEBI_PATH.exists():
        return _chebi_data

    with open(CHEBI_PATH, 'r', encoding='utf-8') as f:
        reader = csv.DictReader(f)
        for row in reader:
            _chebi_data.append({
                'source': (row.get('source_name') or '').strip().lower(),
                'target': (row.get('target_name') or '').strip().lower(),
                'edge_type': (row.get('edge_type') or '').strip().lower(),
                'description': row.get('edge_description', ''),
            })
    return _chebi_data


def _fuzzy_match(query: str, target: str) -> bool:
    """Return whether normalized names match by bidirectional substring."""
    q = query.lower().replace('_', ' ').replace('-', ' ').strip()
    t = target.lower().replace('_', ' ').replace('-', ' ').strip()
    if not q or not t:
        return False
    return q in t or t in q


def check_chebi(metabolite: str, disease: str = "") -> Dict:
    """Match by substring first, then use the SapBERT fallback."""
    data = _load_chebi()

    # Stage 1: exact matching
    matches = []
    for row in data:
        if _fuzzy_match(metabolite, row['source']) or _fuzzy_match(metabolite, row['target']):
            matches.append(row)

    match_method = "exact"

    # Stage 2: SapBERT semantic fallback
    if not matches:
        try:
            met_hits = retrieve(metabolite, entity_type="metabolite", top_k=3, threshold=0.7)
            for hit_name, score in met_hits:
                for row in data:
                    if _fuzzy_match(hit_name, row['source']) or _fuzzy_match(hit_name, row['target']):
                        matches.append(row)
                if matches:
                    match_method = f"sapbert({hit_name},{score:.3f})"
                    break
        except Exception:
            pass

    related = list(set(
        m['target'] if _fuzzy_match(metabolite, m['source']) else m['source']
        for m in matches[:20]
    ))

    return {
        "hit": len(matches) > 0,
        "records": len(matches),
        "match_method": match_method,
        "related_metabolites": related[:10],
        "details": [f"{m['source']}--{m['edge_type']}-->{m['target']}" for m in matches[:5]],
    }
