import csv
from pathlib import Path
from typing import Dict, List, Optional

BINDED_DB_PATH = Path(__file__).parent.parent.parent / "data" / "databases" / "binded_database_anno_named.csv"

_data: Optional[List[Dict]] = None
_sapbert_available: Optional[bool] = None  # None=untested, True/False=cached result


def _load_data():
    """Load all records into memory on first use."""
    global _data
    if _data is not None:
        return
    if not BINDED_DB_PATH.exists():
        print(f"[BindedDB] File not found: {BINDED_DB_PATH}")
        _data = []
        return
    _data = []
    with open(BINDED_DB_PATH, "r", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for row in reader:
            _data.append(row)
    print(f"[BindedDB] Loaded {len(_data)} records")


def _normalize(s: str) -> str:
    """Normalize an entity name for matching."""
    return s.lower().replace("_", " ").replace("-", " ").strip()


def _fuzzy_match(query: str, target: str) -> bool:
    """Return whether either normalized name contains the other."""
    q = _normalize(query)
    t = _normalize(target)
    if not q or not t:
        return False
    return q in t or t in q


def check_binded_microbe_metabolite(microbe: str, metabolite: str) -> Dict:
    """Check whether a microbe-metabolite pair exists in the database."""
    _load_data()

    # Stage 1: substring matching
    matches = []
    for row in _data:
        sp = row.get("species.name", "")
        cn = row.get("Compound_Name", "")
        if not sp or not cn or cn == "NA":
            continue
        if _fuzzy_match(microbe, sp) and _fuzzy_match(metabolite, cn):
            matches.append({
                "species": sp,
                "compound": cn,
                "cid": row.get("cid.ID", ""),
                "source": row.get("sourcename", ""),
                "genus": row.get("genus.name", ""),
            })

    if matches:
        sources = set(m["source"] for m in matches)
        return {
            "hit": True,
            "records": len(matches),
            "sources": list(sources),
            "match_method": "exact",
            "match_details": matches[:10],
        }

    # Stage 2: SapBERT semantic matching
    global _sapbert_available
    if _sapbert_available is False:
        return {"hit": False, "records": 0, "sources": [], "match_method": "none", "match_details": []}
    try:
        from .sapbert_retriever import retrieve
        _sapbert_available = True
        # Retrieve semantic neighbors for both entities.
        microbe_hits = retrieve(microbe, entity_type="microbe", top_k=5, threshold=0.85)
        metabolite_hits = retrieve(metabolite, entity_type="metabolite", top_k=5, threshold=0.85)

        if not microbe_hits or not metabolite_hits:
            return {"hit": False, "records": 0, "sources": [], "match_method": "none", "match_details": []}

        # Match again using the semantic neighbor names.
        microbe_names = [name for name, _ in microbe_hits]
        metabolite_names = [name for name, _ in metabolite_hits]

        semantic_matches = []
        for row in _data:
            sp = row.get("species.name", "")
            cn = row.get("Compound_Name", "")
            if not sp or not cn or cn == "NA":
                continue
            sp_match = any(_fuzzy_match(mn, sp) for mn in microbe_names)
            cn_match = any(_fuzzy_match(mn, cn) for mn in metabolite_names)
            if sp_match and cn_match:
                semantic_matches.append({
                    "species": sp,
                    "compound": cn,
                    "cid": row.get("cid.ID", ""),
                    "source": row.get("sourcename", ""),
                    "genus": row.get("genus.name", ""),
                })

        if semantic_matches:
            sources = set(m["source"] for m in semantic_matches)
            return {
                "hit": True,
                "records": len(semantic_matches),
                "sources": list(sources),
                "match_method": "semantic",
                "match_details": semantic_matches[:10],
            }
    except Exception as e:
        _sapbert_available = False
        print(f"[BindedDB] SapBERT unavailable ({e}), skipping semantic matching for all future calls")

    return {"hit": False, "records": 0, "sources": [], "match_method": "none", "match_details": []}


def check_binded_microbe(microbe: str) -> Dict:
    """Return metabolites associated with a microbe in the database."""
    _load_data()

    compounds = set()
    sources = set()
    count = 0
    for row in _data:
        sp = row.get("species.name", "")
        if not sp:
            continue
        if _fuzzy_match(microbe, sp):
            cn = row.get("Compound_Name", "")
            if cn and cn != "NA":
                compounds.add(cn)
            sources.add(row.get("sourcename", ""))
            count += 1

    return {
        "hit": count > 0,
        "records": count,
        "compounds": sorted(compounds)[:50],
        "sources": list(sources),
    }
