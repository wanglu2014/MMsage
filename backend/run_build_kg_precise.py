"""
Build a KG with a precise PubMed query and OpenAlex enrichment.

Query experimental Akkermansia and IBD literature, extract triples with
DeepSeek while reusing PMID caches, enrich edge metrics with OpenAlex, and
save the graph to data/knowledge_graph/auto_built_kg.gml.

Run from the backend directory with: python run_build_kg_precise.py
"""
import os
import sys
import io
import time
import re
import json
import urllib.request
import threading
from pathlib import Path
from typing import Set, Dict, List
from datetime import datetime

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace', line_buffering=True)
sys.stderr = io.TextIOWrapper(sys.stderr.buffer, encoding='utf-8', errors='replace', line_buffering=True)
sys.path.insert(0, '.')

import networkx as nx
from build_kg import (
    DeepSeekPool, DEFAULT_KEYS_CSV, CACHE_DIR,
    search_pubmed, fetch_abstracts, expand_entity_names,
    build_kg_from_abstracts, _make_node_id, _print_kg_stats,
)

# ================================================================
# Config
# ================================================================

QUERY = (
    '"Akkermansia muciniphila"'
    ' AND "IBD"'
    ' AND ("in vivo" OR "in vitro" OR "cell line" OR "clinical trial"'
    ' OR "experiment" OR "mouse" OR "mice" OR "rat" OR "rats"'
    ' OR "murine" OR "animal model" OR "patient" OR "patients"'
    ' OR "cohort" OR "randomized" OR "biopsy" OR "fecal"'
    ' OR "stool" OR "culture" OR "fermentation")'
)

_PROJECT_ROOT = Path(__file__).resolve().parent.parent
OUTPUT_GML = str(_PROJECT_ROOT / "data" / "knowledge_graph" / "auto_built_kg.gml")
_OPENALEX_UA_MAIL = os.getenv("ENTREZ_EMAIL", "").strip() or "anonymous"
_OPENALEX_UA = f"MMSage-KG/1.0 (mailto:{_OPENALEX_UA_MAIL})"


# ================================================================
# Step 1-2: PubMed search + DeepSeek LLM extraction
# ================================================================

def build_kg_step():
    print("=" * 65)
    print("  STEP 1-2: PubMed + DeepSeek LLM extraction")
    print("=" * 65)
    print(f"  Query: {QUERY[:100]}...")
    print(f"  Output: {OUTPUT_GML}")

    pool = DeepSeekPool(Path(str(DEFAULT_KEYS_CSV)))

    print("\nPhase 1: Preparing entity aliases for abstract recognition...")
    bacteria_set: Set[str] = {"Akkermansia_muciniphila"}
    metabolite_set: Set[str] = set()
    disease = "IBD"
    bacteria_names, metabolite_names, disease_terms = expand_entity_names(
        bacteria_set, metabolite_set, disease)

    print(f"\nPhase 2: PubMed search...")
    pmids = search_pubmed(QUERY, max_results=500)
    print(f"  Found {len(pmids)} PMIDs")

    print(f"\nPhase 3: Fetching abstracts...")
    articles = fetch_abstracts(pmids)
    print(f"  Got {len(articles)} articles with abstracts")

    print(f"\nPhase 4: DeepSeek LLM extraction...")
    cached = len([p for p in pmids if (CACHE_DIR / "extractions" / f"{p}.json").exists()])
    print(f"  {cached} cached, {len(articles) - cached} new")

    G = build_kg_from_abstracts(
        articles, bacteria_names, metabolite_names,
        disease_terms, disease, pool)

    disease_id = _make_node_id(disease)
    if disease_id not in G:
        G.add_node(disease_id, label=disease, type="disease", node_type="disease")

    _print_kg_stats(G)

    out = Path(OUTPUT_GML)
    out.parent.mkdir(parents=True, exist_ok=True)
    for nid in G.nodes():
        if "label" not in G.nodes[nid]:
            G.nodes[nid]["label"] = str(nid)
    nx.write_gml(G, str(out))
    print(f"\nSaved raw GML to {OUTPUT_GML}")
    return G


# ================================================================
# Step 3: OpenAlex enrichment
# (adapted from rep1221_clin/add_metrics_to_gml.py,
#  using edge_ prefixes and adding publication_year)
# ================================================================

class OpenAlexEnricher:
    def __init__(self, delay=0.4):
        self.delay = delay
        self.cache: Dict[str, dict] = {}
        self.failed: List[str] = []
        self._lock = threading.Lock()
        self._last_req = 0.0

    def _wait(self):
        with self._lock:
            elapsed = time.time() - self._last_req
            if elapsed < self.delay:
                time.sleep(self.delay - elapsed)
            self._last_req = time.time()

    def fetch(self, pmid: str) -> dict:
        clean = re.search(r'\d+', str(pmid))
        if not clean:
            return {}
        pmid = clean.group()
        if pmid in self.cache:
            return self.cache[pmid]

        self._wait()
        try:
            opener = urllib.request.build_opener(urllib.request.ProxyHandler({}))
            url = f"https://api.openalex.org/works/pmid:{pmid}"
            req = urllib.request.Request(url, headers={
                'User-Agent': _OPENALEX_UA})
            with opener.open(req, timeout=30) as resp:
                data = json.loads(resp.read().decode('utf-8'))

            result = {
                'edge_citation_count': data.get('cited_by_count', 0),
                'publication_year': data.get('publication_year', 0),
                'edge_impact_factor': 0,
            }

            loc = data.get('primary_location', {}) or {}
            src = loc.get('source', {}) or {}
            if src.get('type') != 'repository' and src.get('id'):
                src_id = src['id'].split('/')[-1]
                self._wait()
                url2 = f"https://api.openalex.org/sources/{src_id}"
                req2 = urllib.request.Request(url2, headers={
                    'User-Agent': _OPENALEX_UA})
                try:
                    with opener.open(req2, timeout=30) as resp2:
                        sdata = json.loads(resp2.read().decode('utf-8'))
                        result['edge_impact_factor'] = (
                            sdata.get('summary_stats', {}).get('2yr_mean_citedness', 0) or 0)
                except Exception:
                    pass

            self.cache[pmid] = result
            return result
        except Exception as e:
            self.failed.append(pmid)
            return {}

    def enrich(self, gml_path: str):
        print("\n" + "=" * 65)
        print("  STEP 3: OpenAlex enrichment")
        print("=" * 65)
        print(f"  Start: {datetime.now().strftime('%H:%M:%S')}")

        G = nx.read_gml(gml_path)
        print(f"  Loaded: {G.number_of_nodes()} nodes, {G.number_of_edges()} edges")

        unique_pmids: Set[str] = set()
        for _, _, d in G.edges(data=True):
            raw = d.get('pmid', '')
            for p in str(raw).split(','):
                p = p.strip()
                if p and re.search(r'\d+', p):
                    unique_pmids.add(re.search(r'\d+', p).group())

        print(f"  Unique PMIDs: {len(unique_pmids)}")
        est_min = len(unique_pmids) * self.delay * 2 / 60
        print(f"  Estimated time: {est_min:.1f} min")

        t0 = time.time()
        for i, pmid in enumerate(unique_pmids, 1):
            self.fetch(pmid)
            if i % 50 == 0:
                elapsed = time.time() - t0
                speed = i / elapsed if elapsed > 0 else 0
                remain = (len(unique_pmids) - i) / speed if speed > 0 else 0
                print(f"  {i}/{len(unique_pmids)} ({i*100//len(unique_pmids)}%) "
                      f"- {speed:.1f}/s - {remain/60:.1f}min left")

        print(f"\n  Fetched: {len(self.cache)} OK, {len(self.failed)} failed "
              f"({time.time()-t0:.0f}s)")

        enriched = 0
        for u, v, d in G.edges(data=True):
            raw = d.get('pmid', '')
            first_pmid = str(raw).split(',')[0].strip()
            m = re.search(r'\d+', first_pmid)
            if m and m.group() in self.cache:
                metrics = self.cache[m.group()]
                d['edge_impact_factor'] = metrics.get('edge_impact_factor', 0)
                d['edge_citation_count'] = metrics.get('edge_citation_count', 0)
                d['publication_year'] = metrics.get('publication_year', 0)
                enriched += 1

        print(f"  Enriched {enriched}/{G.number_of_edges()} edges")

        nx.write_gml(G, gml_path)
        print(f"  Saved: {gml_path}")
        print(f"  Done: {datetime.now().strftime('%H:%M:%S')}")


# ================================================================
# Main
# ================================================================

if __name__ == "__main__":
    G = build_kg_step()
    enricher = OpenAlexEnricher(delay=0.4)
    enricher.enrich(OUTPUT_GML)
    print("\n" + "=" * 65)
    print("  ALL DONE")
    print("=" * 65)
