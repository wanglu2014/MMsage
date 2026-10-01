"""
PubMed Query Module
====================
Query PubMed for co-occurrence counts with caching.
Reused from existing step2_chain_novelty.py
"""

import os
import json
import hashlib
import time
from pathlib import Path
from typing import Dict, List, Optional
from functools import lru_cache

try:
    from Bio import Entrez
    HAS_ENTREZ = True
except ImportError:
    HAS_ENTREZ = False

try:
    import sys as _sys
    _sys.path.insert(0, str(Path(__file__).parent.parent))
    from db_checkers import sapbert_retriever
    _HAS_SAPBERT = True
except ImportError:
    _HAS_SAPBERT = False

# Configuration
ENTREZ_EMAIL = os.getenv("ENTREZ_EMAIL", "").strip() or "entrez-not-configured@invalid"
CACHE_DIR = Path(__file__).parent.parent.parent / "cache" / "pubmed"
CACHE_DIR.mkdir(parents=True, exist_ok=True)
_SAPBERT_SYN_CACHE_FILE = CACHE_DIR / "sapbert_agent_synonyms.json"

# SapBERT-driven synonym cache (replaces hardcoded DEFAULT_SYNONYMS)
DEFAULT_SYNONYMS: Dict[str, List[str]] = {}   # backward-compat export
_sapbert_syn_cache: Dict[str, List[str]] = {}
_sapbert_syn_loaded = False


def _load_agent_syn_cache():
    global _sapbert_syn_cache, _sapbert_syn_loaded
    if _sapbert_syn_loaded:
        return
    if _SAPBERT_SYN_CACHE_FILE.exists():
        try:
            with open(_SAPBERT_SYN_CACHE_FILE, "r", encoding="utf-8") as f:
                _sapbert_syn_cache = json.load(f)
        except Exception:
            pass
    _sapbert_syn_loaded = True


def _save_agent_syn_cache():
    try:
        with open(_SAPBERT_SYN_CACHE_FILE, "w", encoding="utf-8") as f:
            json.dump(_sapbert_syn_cache, f, indent=2, ensure_ascii=False)
    except Exception:
        pass


def _get_sapbert_synonyms(entity_name: str, entity_type: str = None) -> List[str]:
    """Get synonyms via SapBERT. Falls back to [entity_name] if unavailable."""
    _load_agent_syn_cache()
    cache_key = f"{entity_name}|{entity_type or 'any'}"
    if cache_key in _sapbert_syn_cache:
        return _sapbert_syn_cache[cache_key]

    syns = [entity_name.replace("_", " ")]
    if _HAS_SAPBERT:
        try:
            hits = sapbert_retriever.retrieve(
                entity_name.replace("_", " "),
                entity_type=entity_type,
                top_k=5, threshold=0.75, exact_first=True,
            )
            for name, _score in hits:
                clean = name.strip()
                if clean and clean.lower() not in {s.lower() for s in syns}:
                    syns.append(clean)
        except Exception:
            pass

    _sapbert_syn_cache[cache_key] = syns
    _save_agent_syn_cache()
    return syns


def _get_cache_key(query: str) -> str:
    """Generate cache key from query."""
    return hashlib.md5(query.encode()).hexdigest()


def _load_cache(key: str) -> Optional[int]:
    """Load count from cache."""
    cache_file = CACHE_DIR / f"{key}.json"
    if cache_file.exists():
        try:
            with open(cache_file, 'r') as f:
                data = json.load(f)
                return data.get('count')
        except:
            pass
    return None


def _save_cache(key: str, count: int, query: str):
    """Save count to cache."""
    cache_file = CACHE_DIR / f"{key}.json"
    try:
        with open(cache_file, 'w') as f:
            json.dump({'count': count, 'query': query, 'timestamp': time.time()}, f)
    except:
        pass


def build_pubmed_query(terms: List[str], synonym_map: Optional[Dict] = None) -> str:
    """
    Build PubMed query with SapBERT-driven synonym expansion.

    Each term is expanded to its synonyms (via SapBERT retrieval or
    optional synonym_map override), OR-joined, then all groups AND-joined.

    Args:
        terms: List of search terms
        synonym_map: Optional synonym overrides (takes priority over SapBERT)

    Returns:
        PubMed query string
    """
    if synonym_map is None:
        synonym_map = {}

    query_parts = []
    for term in terms:
        term_lower = term.lower().replace('_', ' ')

        # Priority 1: explicit synonym_map override
        synonyms = None
        for key, vals in synonym_map.items():
            if key in term_lower or term_lower in key:
                synonyms = vals + [term]
                break
            for val in vals:
                if val.lower() in term_lower:
                    synonyms = vals + [term]
                    break
            if synonyms:
                break

        # Priority 2: SapBERT dynamic synonym expansion
        if not synonyms:
            entity_type = None
            # Generic disease detection: common suffixes/keywords across diseases
            _disease_hints = ("disease", "syndrome", "disorder", "itis",
                              "osis", "emia", "uria", "pathy", "cancer",
                              "tumor", "carcinoma", "lymphoma", "leukemia")
            if any(kw in term_lower for kw in _disease_hints):
                entity_type = "disease"
            synonyms = _get_sapbert_synonyms(term, entity_type=entity_type)

        unique_syns = list(dict.fromkeys(s for s in synonyms if s))[:4]
        if len(unique_syns) > 1:
            syn_query = ' OR '.join(
                f'"{s}"[Title/Abstract]' for s in unique_syns)
            query_parts.append(f"({syn_query})")
        else:
            query_parts.append(f'"{unique_syns[0]}"[Title/Abstract]')

    return ' AND '.join(query_parts)


def query_pubmed_count(query: str, use_cache: bool = True) -> int:
    """
    Query PubMed for article count.

    Args:
        query: PubMed query string
        use_cache: Whether to use caching

    Returns:
        Number of articles matching query
    """
    # Check cache first
    cache_key = _get_cache_key(query)
    if use_cache:
        cached = _load_cache(cache_key)
        if cached is not None:
            return cached

    if not HAS_ENTREZ:
        # Return 0 if Entrez not available
        _save_cache(cache_key, 0, query)
        return 0

    try:
        Entrez.email = ENTREZ_EMAIL
        handle = Entrez.esearch(db="pubmed", term=query, retmax=0)
        record = Entrez.read(handle)
        handle.close()
        count = int(record.get("Count", 0))

        if use_cache:
            _save_cache(cache_key, count, query)

        # Rate limiting
        time.sleep(0.35)
        return count
    except Exception:
        if use_cache:
            _save_cache(cache_key, 0, query)
        return 0


def normalize_novelty(count: int, c_max: int = 500) -> float:
    """
    Normalize count to novelty score using log formula.

    novelty = 1 - log(1 + count) / log(1 + c_max)

    Args:
        count: Co-occurrence count
        c_max: Normalization constant

    Returns:
        Novelty score in [0, 1]
    """
    import math
    if count <= 0:
        return 1.0
    novelty = 1 - math.log(1 + count) / math.log(1 + c_max)
    return max(0.0, min(1.0, novelty))
