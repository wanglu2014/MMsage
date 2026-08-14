"""Build KG with an exact Akkermansia + IBD query and experiment filters."""
import sys, io, time
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace', line_buffering=True)
sys.stderr = io.TextIOWrapper(sys.stderr.buffer, encoding='utf-8', errors='replace', line_buffering=True)
sys.path.insert(0, '.')

from pathlib import Path
from typing import Set, Dict, List
from build_kg import (
    DeepSeekPool, DEFAULT_KEYS_CSV, CACHE_DIR,
    search_pubmed, fetch_abstracts, expand_entity_names,
    build_kg_from_abstracts, _make_node_id, _print_kg_stats, _save_gml,
)
import networkx as nx

QUERY = (
    '"Akkermansia muciniphila"'
    ' AND "IBD"'
    ' AND ("in vivo" OR "in vitro" OR "cell line" OR "clinical trial"'
    ' OR "experiment" OR "mouse" OR "mice" OR "rat" OR "rats"'
    ' OR "murine" OR "animal model" OR "patient" OR "patients"'
    ' OR "cohort" OR "randomized" OR "biopsy" OR "fecal"'
    ' OR "stool" OR "culture" OR "fermentation")'
)

PROJECT_ROOT = Path(__file__).resolve().parent.parent
OUTPUT = str(PROJECT_ROOT / "data" / "knowledge_graph" / "auto_built_kg.gml")

print("=" * 65)
print("  BUILD KG: Exact Akkermansia + IBD query with experiment filters")
print("=" * 65)
print(f"  Query: {QUERY[:120]}...")

pool = DeepSeekPool(Path(str(DEFAULT_KEYS_CSV)))

# Phase 1: aliases used only for post-retrieval entity recognition.
print("\nPhase 1: Preparing entity aliases for abstract recognition...")
bacteria_set: Set[str] = {"Akkermansia_muciniphila"}
metabolite_set: Set[str] = set()
disease = "IBD"
bacteria_names, metabolite_names, disease_terms = expand_entity_names(
    bacteria_set, metabolite_set, disease)

# Phase 2: single precise query
print(f"\nPhase 2: PubMed search...")
pmids = search_pubmed(QUERY, max_results=500)
print(f"  Found {len(pmids)} PMIDs")

# Phase 3: fetch abstracts
print(f"\nPhase 3: Fetching abstracts...")
articles = fetch_abstracts(pmids)
print(f"  Got {len(articles)} articles with abstracts")

# Phase 4: LLM extraction (cached results reused)
print(f"\nPhase 4: Extracting relations via DeepSeek...")
cached = len([p for p in pmids if (CACHE_DIR / "extractions" / f"{p}.json").exists()])
print(f"  {cached} already cached, {len(articles) - cached} new to process")

G = build_kg_from_abstracts(
    articles, bacteria_names, metabolite_names,
    disease_terms, disease, pool)

disease_id = _make_node_id(disease)
if disease_id not in G:
    G.add_node(disease_id, label=disease, type="disease", node_type="disease")

_print_kg_stats(G)
_save_gml(G, OUTPUT)
