import csv
from pathlib import Path
from typing import Dict, List, Optional

try:
    from .sapbert_retriever import retrieve
except ImportError:
    from sapbert_retriever import retrieve

BACDIVE_PATH = Path(__file__).parent.parent.parent / "data" / "databases" / "bacdive_edges.tsv"

_bacdive_data: Optional[List[Dict]] = None


def _load_bacdive() -> List[Dict]:
    """Load and normalize BacDive TSV records with caching."""
    global _bacdive_data
    if _bacdive_data is not None:
        return _bacdive_data

    _bacdive_data = []
    if not BACDIVE_PATH.exists():
        return _bacdive_data

    with open(BACDIVE_PATH, 'r', encoding='utf-8') as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            _bacdive_data.append({
                'microbe': (row.get('source_name') or '').strip().lower(),
                'microbe_std': (row.get('standard_sourcename') or '').strip().lower(),
                'edge_type': (row.get('edge_type') or '').strip().lower(),
                'target': (row.get('target_name') or '').strip().lower(),
                'target_attr': (row.get('target_attribute') or '').strip().lower(),
            })
    return _bacdive_data


def _fuzzy_match(query: str, target: str) -> bool:
    """Return whether normalized names match by bidirectional substring."""
    q = query.lower().replace('_', ' ').replace('-', ' ').strip()
    t = target.lower().replace('_', ' ').replace('-', ' ').strip()
    if not q or not t:
        return False
    return q in t or t in q


def check_bacdive(bacteria: str, metabolite: str) -> Dict:
    """Match by substring first, then use the SapBERT fallback."""
    data = _load_bacdive()

    # Stage 1: exact matching
    matches = []
    for row in data:
        if _fuzzy_match(bacteria, row['microbe']) or _fuzzy_match(bacteria, row['microbe_std']):
            matches.append(row)

    match_method = "exact"

    # Stage 2: SapBERT semantic fallback
    if not matches:
        try:
            microbe_hits = retrieve(bacteria, entity_type="microbe", top_k=3, threshold=0.7)
            for hit_name, score in microbe_hits:
                for row in data:
                    if _fuzzy_match(hit_name, row['microbe']) or _fuzzy_match(hit_name, row['microbe_std']):
                        matches.append(row)
                if matches:
                    match_method = f"sapbert({hit_name},{score:.3f})"
                    break
        except Exception:
            pass

    environments = list(set(m['target'] for m in matches[:20] if m['target']))

    return {
        "hit": len(matches) > 0,
        "records": len(matches),
        "match_method": match_method,
        "environments": environments[:10],
        "details": [f"{m['microbe']}--{m['edge_type']}-->{m['target']}" for m in matches[:5]],
    }
