import csv
from pathlib import Path
from typing import Dict, List, Optional

try:
    from .sapbert_retriever import retrieve
except ImportError:
    from sapbert_retriever import retrieve

CTD_PATH = Path(__file__).parent.parent.parent / "data" / "databases" / "CTD_chemicals_diseases.tsv"

_ctd_data: Optional[List[Dict]] = None


def _load_ctd() -> List[Dict]:
    """Load and normalize CTD TSV records with caching."""
    global _ctd_data
    if _ctd_data is not None:
        return _ctd_data

    _ctd_data = []
    if not CTD_PATH.exists():
        return _ctd_data

    with open(CTD_PATH, 'r', encoding='utf-8') as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            _ctd_data.append({
                'chemical': (row.get('source_name') or '').strip().lower(),
                'disease': (row.get('target_name') or '').strip().lower(),
                'edge_type': (row.get('edge_type') or '').strip().lower(),
                'inference_score': row.get('edge_description', ''),
                'pmids': row.get('PMID', ''),
            })
    return _ctd_data


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


def check_ctd(metabolite: str, disease: str) -> Dict:
    """Match directly first, then use the SapBERT fallback."""
    data = _load_ctd()

    # Stage 1: exact matching
    matches = []
    for row in data:
        chem_match = _fuzzy_match(metabolite, row['chemical'])
        dis_match = _disease_match(disease, row['disease'])
        if chem_match and dis_match:
            matches.append(row)

    match_method = "exact"

    # Stage 2: SapBERT semantic fallback
    if not matches:
        try:
            met_hits = retrieve(metabolite, entity_type="metabolite", top_k=3, threshold=0.7)
            dis_hits = retrieve(disease, entity_type="disease", top_k=3, threshold=0.7)
            alt_mets = [metabolite] + [h[0] for h in met_hits]
            alt_diss = [disease] + [h[0] for h in dis_hits]
            for am in alt_mets:
                for ad in alt_diss:
                    for row in data:
                        if _fuzzy_match(am, row['chemical']) and _disease_match(ad, row['disease']):
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
        "details": [f"{m['chemical']}--{m['edge_type']}-->{m['disease']}" for m in matches[:5]],
    }
