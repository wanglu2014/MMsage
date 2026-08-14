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

try:
    from Bio import Entrez
    HAS_ENTREZ = True
except ImportError:
    HAS_ENTREZ = False

# Configuration
ENTREZ_EMAIL = os.getenv("ENTREZ_EMAIL", "").strip() or "entrez-not-configured@invalid"
CACHE_DIR = Path(__file__).parent.parent.parent / "cache" / "pubmed"
CACHE_DIR.mkdir(parents=True, exist_ok=True)


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
    Build an AND-joined PubMed query without expanding input terms.

    Args:
        terms: List of search terms
        synonym_map: Retained for backward compatibility and ignored

    Returns:
        PubMed query string
    """
    normalized_terms = [
        str(term).replace("_", " ").strip()
        for term in terms
        if str(term).strip()
    ]
    return " AND ".join(
        f'"{term}"[Title/Abstract]' for term in normalized_terms
    )


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
