import csv
from pathlib import Path
from typing import Dict, List, Optional
from functools import lru_cache

try:
    from .sapbert_retriever import retrieve, retrieve_batch
except ImportError:
    from sapbert_retriever import retrieve, retrieve_batch

DISBIOME_PATH = Path(__file__).parent.parent.parent / "data" / "databases" / "disbiome_edges.tsv"

_disbiome_data: Optional[List[Dict]] = None


def _load_disbiome() -> List[Dict]:
    """Load and normalize Disbiome TSV records with caching."""
    global _disbiome_data
    if _disbiome_data is not None:
        return _disbiome_data

    _disbiome_data = []
    if not DISBIOME_PATH.exists():
        return _disbiome_data

    with open(DISBIOME_PATH, 'r', encoding='utf-8') as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            _disbiome_data.append({
                'microbe': (row.get('source_name') or '').strip().lower(),
                'microbe_std': (row.get('standard_sourcename') or '').strip().lower(),
                'edge_type': (row.get('edge_type') or '').strip().lower(),
                'disease': (row.get('target_name') or '').strip().lower(),
                'disease_std': (row.get('standard_targetname') or '').strip().lower(),
                'pmid': row.get('PMID', ''),
            })
    return _disbiome_data


def _fuzzy_match(query: str, target: str) -> bool:
    """Return whether normalized names match by bidirectional substring."""
    q = query.lower().replace('_', ' ').replace('-', ' ').strip()
    t = target.lower().replace('_', ' ').replace('-', ' ').strip()
    if not q or not t:
        return False
    return q in t or t in q


DISEASE_SYNONYMS = {
    'ibd': ['inflammatory bowel disease', 'crohn', "crohn's disease", 'ulcerative colitis'],
    'uc': ['ulcerative colitis'],
    'cd': ["crohn's disease", 'crohn disease'],
    't2d': ['type 2 diabetes', 'diabetes mellitus type 2'],
    'crc': ['colorectal cancer', 'colon cancer'],
    'nafld': ['non-alcoholic fatty liver', 'nonalcoholic fatty liver'],
    'nash': ['nonalcoholic steatohepatitis', 'non-alcoholic steatohepatitis'],
}


def _disease_match(query: str, target: str) -> bool:
    """Match a disease directly or through the configured synonyms."""
    if _fuzzy_match(query, target):
        return True
    q = query.lower().strip()
    synonyms = DISEASE_SYNONYMS.get(q, [])
    for syn in synonyms:
        if _fuzzy_match(syn, target):
            return True
    return False


def check_disbiome(bacteria: str, metabolite: str) -> Dict:
    """Match by substring first, then use the SapBERT fallback."""
    data = _load_disbiome()

    # Stage 1: substring matching
    matches = []
    for row in data:
        if _fuzzy_match(bacteria, row['microbe']) or _fuzzy_match(bacteria, row['microbe_std']):
            matches.append(row)

    match_method = "exact"

    # Stage 2: SapBERT semantic fallback
    if not matches:
        try:
            sapbert_hits = retrieve(bacteria, entity_type="microbe", top_k=3, threshold=0.7)
            for hit_name, score in sapbert_hits:
                for row in data:
                    if _fuzzy_match(hit_name, row['microbe']) or _fuzzy_match(hit_name, row['microbe_std']):
                        matches.append(row)
                if matches:
                    match_method = f"sapbert({hit_name},{score:.3f})"
                    break
        except Exception:
            pass

    return {
        "hit": len(matches) > 0,
        "records": len(matches),
        "match_method": match_method,
        "details": [f"{m['microbe']}--{m['edge_type']}-->{m['disease']}" for m in matches[:5]],
    }


def check_disbiome_disease(bacteria: str, disease: str) -> Dict:
    """Match directly first, then use the SapBERT fallback."""
    data = _load_disbiome()

    # Stage 1: exact matching
    matches = []
    for row in data:
        microbe_match = _fuzzy_match(bacteria, row['microbe']) or _fuzzy_match(bacteria, row['microbe_std'])
        disease_match = _disease_match(disease, row['disease']) or _disease_match(disease, row['disease_std'])
        if microbe_match and disease_match:
            matches.append(row)

    match_method = "exact"

    # Stage 2: SapBERT semantic fallback
    if not matches:
        try:
            microbe_hits = retrieve(bacteria, entity_type="microbe", top_k=3, threshold=0.7)
            disease_hits = retrieve(disease, entity_type="disease", top_k=3, threshold=0.7)
            # Match again using the semantic neighbor names.
            alt_microbes = [bacteria] + [h[0] for h in microbe_hits]
            alt_diseases = [disease] + [h[0] for h in disease_hits]
            for am in alt_microbes:
                for ad in alt_diseases:
                    for row in data:
                        microbe_match = _fuzzy_match(am, row['microbe']) or _fuzzy_match(am, row['microbe_std'])
                        disease_match = _disease_match(ad, row['disease']) or _disease_match(ad, row['disease_std'])
                        if microbe_match and disease_match:
                            matches.append(row)
                    if matches:
                        match_method = f"sapbert(m={am},d={ad})"
                        break
                if matches:
                    break
        except Exception:
            pass

    edge_types = list(set(m['edge_type'] for m in matches if m['edge_type']))

    return {
        "hit": len(matches) > 0,
        "records": len(matches),
        "match_method": match_method,
        "edge_types": edge_types,
        "details": [f"{m['microbe']}--{m['edge_type']}-->{m['disease']}" for m in matches[:5]],
    }
