"""
Build Knowledge Graph from PubMed literature (v3)
==================================================
DeepSeek API (multi-key pool from CSV, openai SDK) + SapBERT NER.
Prompt/parse aligned with rep1212 ===BEGIN NODES/EDGES=== CSV format.

Two entry points:
  run_build_kg_direct(species, disease)  -- standalone, no candidates needed
  run_build_kg(candidates_path, ...)     -- full pipeline integration

Pipeline:
1. Build PubMed queries from species+disease (or candidates)
2. Fetch abstracts with experimental evidence filter
3. SapBERT-expanded entity recognition
4. DeepSeek relation extraction (===BEGIN NODES/EDGES=== CSV format)
5. Build NetworkX graph and save as GML
"""

import os
import json
import sys
import io
import csv
import time
import re
import random
import threading
from pathlib import Path
from typing import List, Dict, Any, Set, Tuple, Optional

import networkx as nx
from openai import OpenAI

if sys.stdout.encoding and sys.stdout.encoding.lower() != 'utf-8':
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')

try:
    from Bio import Entrez, Medline
    HAS_ENTREZ = True
except ImportError:
    HAS_ENTREZ = False

try:
    from db_checkers import sapbert_retriever
    HAS_SAPBERT = True
except ImportError:
    HAS_SAPBERT = False


# ---- Config ----

ENTREZ_EMAIL = os.getenv("ENTREZ_EMAIL", "").strip() or "entrez-not-configured@invalid"
CACHE_DIR = Path(__file__).parent.parent / "cache" / "kg_build"
SYNONYM_CACHE_FILE = CACHE_DIR / "sapbert_synonym_cache.json"
DEFAULT_DISEASE = "IBD"
DEFAULT_KEYS_CSV = Path(__file__).resolve().parent.parent / "key" / "keys_0605.csv"

BATCH_SIZE = 50
MAX_RETRY = 5
RETRY_DELAY = 2
QUERY_DELAY_S = 0.35

DEEPSEEK_BASE_URL = "https://api.deepseek.com"
DEEPSEEK_MODEL = "deepseek-chat"

EXPERIMENTAL_FILTER = (
    ' AND ("in vivo" OR "in vitro"'
    ' OR "cell line" OR "clinical trial" OR "experiment")'
)

VALID_RELATIONS = frozenset({
    "produces", "metabolizes", "inhibits", "promotes", "protects",
    "associated_with", "degrades", "transports", "modulates",
})


# ================================================================
# DeepSeek API Client Pool (adapted from rep1212 APIClientPool)
# ================================================================

def _load_keys_from_csv(keys_csv: Path) -> List[str]:
    """Load DeepSeek API keys from CSV file (rep1212 format)."""
    keys = []
    with open(keys_csv, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            key = line.split(",")[0].strip().strip('"').strip("'")
            if key.startswith("sk-"):
                keys.append(key)
    return keys


class DeepSeekPool:
    """Simplified multi-key client pool (from rep1212 APIClientPool)."""

    def __init__(self, keys_csv: Path = DEFAULT_KEYS_CSV, max_keys: int = 10):
        all_keys = _load_keys_from_csv(keys_csv)
        if not all_keys:
            raise RuntimeError(f"No API keys found in {keys_csv}")
        use_keys = all_keys[:max_keys]
        print(f"[DeepSeek] Loaded {len(use_keys)} keys from {keys_csv}")

        self._clients = []
        for i, key in enumerate(use_keys):
            client = OpenAI(
                api_key=key,
                base_url=DEEPSEEK_BASE_URL,
                timeout=120,
                max_retries=1,
            )
            self._clients.append({
                "client": client,
                "key_tail": key[-8:],
                "errors": 0,
                "last_used": 0.0,
            })
        self._lock = threading.Lock()

    def call(self, system_prompt: str, user_prompt: str,
             max_retries: int = 3) -> Optional[str]:
        """Call DeepSeek API with key rotation and retry."""
        for attempt in range(max_retries):
            info = self._pick_client()
            try:
                resp = info["client"].chat.completions.create(
                    model=DEEPSEEK_MODEL,
                    messages=[
                        {"role": "system", "content": system_prompt},
                        {"role": "user", "content": user_prompt},
                    ],
                    temperature=0.1,
                    max_tokens=2000,
                    timeout=90,
                )
                text = resp.choices[0].message.content or ""
                # strip <think>...</think> tags (deepseek-reasoner)
                text = re.sub(r"<think>.*?</think>", "", text, flags=re.DOTALL)
                with self._lock:
                    info["errors"] = max(0, info["errors"] - 1)
                return text.strip()
            except Exception as e:
                with self._lock:
                    info["errors"] += 1
                delay = RETRY_DELAY * (2 ** attempt) * random.uniform(0.8, 1.2)
                print(f"  [WARN] DeepSeek call failed (key ...{info['key_tail']}, "
                      f"attempt {attempt+1}): {e}")
                if attempt < max_retries - 1:
                    time.sleep(delay)
        return None

    def _pick_client(self) -> Dict:
        with self._lock:
            available = [c for c in self._clients if c["errors"] < 5]
            if not available:
                for c in self._clients:
                    c["errors"] = 0
                available = self._clients
            pick = min(available, key=lambda c: (c["errors"], c["last_used"]))
            pick["last_used"] = time.time()
            return pick


_pool: Optional[DeepSeekPool] = None


def _get_pool(keys_csv: Path = DEFAULT_KEYS_CSV) -> DeepSeekPool:
    global _pool
    if _pool is None:
        _pool = DeepSeekPool(keys_csv)
    return _pool


# ================================================================
# SapBERT aliases for post-retrieval entity recognition (cached)
# ================================================================

_synonym_cache: Dict[str, List[str]] = {}
_synonym_cache_loaded = False


def _load_synonym_cache():
    global _synonym_cache, _synonym_cache_loaded
    if _synonym_cache_loaded:
        return
    CACHE_DIR.mkdir(parents=True, exist_ok=True)
    if SYNONYM_CACHE_FILE.exists():
        try:
            with open(SYNONYM_CACHE_FILE, "r", encoding="utf-8") as f:
                _synonym_cache = json.load(f)
        except Exception:
            _synonym_cache = {}
    _synonym_cache_loaded = True


def _save_synonym_cache():
    CACHE_DIR.mkdir(parents=True, exist_ok=True)
    with open(SYNONYM_CACHE_FILE, "w", encoding="utf-8") as f:
        json.dump(_synonym_cache, f, indent=2, ensure_ascii=False)


def get_sapbert_synonyms(
    entity_name: str, entity_type: str = None,
    top_k: int = 5, threshold: float = 0.75,
) -> List[str]:
    """Get aliases for post-retrieval entity recognition. Cached."""
    _load_synonym_cache()
    cache_key = f"{entity_name}|{entity_type or 'any'}"
    if cache_key in _synonym_cache:
        return _synonym_cache[cache_key]

    synonyms = [entity_name.replace("_", " ")]
    if HAS_SAPBERT:
        try:
            hits = sapbert_retriever.retrieve(
                entity_name.replace("_", " "),
                entity_type=entity_type, top_k=top_k,
                threshold=threshold, exact_first=True,
            )
            for name, _ in hits:
                clean = name.strip()
                if clean and clean.lower() not in {s.lower() for s in synonyms}:
                    synonyms.append(clean)
        except Exception as e:
            print(f"  [WARN] SapBERT lookup failed for '{entity_name}': {e}")

    _synonym_cache[cache_key] = synonyms
    return synonyms


def expand_entity_names(
    bacteria_set: Set[str], metabolite_set: Set[str], disease: str,
) -> Tuple[Dict[str, str], Dict[str, str], List[str]]:
    """Build aliases for abstract recognition, not PubMed query expansion."""
    bacteria_names: Dict[str, str] = {}
    for bact in bacteria_set:
        for s in get_sapbert_synonyms(bact, entity_type="microbe"):
            bacteria_names[s.lower()] = bact
        parts = bact.replace("_", " ").split()
        if parts:
            bacteria_names[parts[0].lower()] = bact

    metabolite_names: Dict[str, str] = {}
    for met in metabolite_set:
        for s in get_sapbert_synonyms(met, entity_type="metabolite"):
            metabolite_names[s.lower()] = met

    disease_terms = get_sapbert_synonyms(disease, entity_type="disease")
    _DISEASE_FALLBACK_SYNONYMS = {
        "IBD": ["IBD", "inflammatory bowel disease", "Crohn's disease",
                "ulcerative colitis", "colitis"],
    }
    for fallback in _DISEASE_FALLBACK_SYNONYMS.get(disease.upper(), []):
        if fallback.lower() not in {t.lower() for t in disease_terms}:
            disease_terms.append(fallback)

    print(f"  Entity expansion: {len(bacteria_set)} bacteria -> "
          f"{len(bacteria_names)} name variants")
    print(f"  Entity expansion: {len(metabolite_set)} metabolites -> "
          f"{len(metabolite_names)} name variants")
    print(f"  Disease synonyms: {disease_terms}")
    _save_synonym_cache()
    return bacteria_names, metabolite_names, disease_terms


# ================================================================
# PubMed Fetching (robust, from rep1212)
# ================================================================

def search_pubmed(query: str, max_results: int = 40) -> List[str]:
    if not HAS_ENTREZ:
        return []
    Entrez.email = ENTREZ_EMAIL
    record = None
    for attempt in range(MAX_RETRY):
        try:
            handle = Entrez.esearch(
                db="pubmed", term=query, retmax=0,
                usehistory="y", sort="relevance")
            record = Entrez.read(handle)
            handle.close()
            break
        except Exception as e:
            wait = RETRY_DELAY * (attempt + 1)
            print(f"  [WARN] PubMed search attempt {attempt+1} failed: {e}")
            time.sleep(wait)
            if attempt == MAX_RETRY - 1:
                return []

    total = int(record.get("Count", 0)) if record else 0
    if total == 0:
        return []
    actual_max = min(total, max_results)
    pmids: List[str] = []
    for start in range(0, actual_max, BATCH_SIZE):
        batch_size = min(BATCH_SIZE, actual_max - start)
        for attempt in range(MAX_RETRY):
            try:
                time.sleep(QUERY_DELAY_S)
                handle = Entrez.esearch(
                    db="pubmed", term=query,
                    retstart=start, retmax=batch_size, sort="relevance")
                result = Entrez.read(handle)
                handle.close()
                pmids.extend(result.get("IdList", []))
                break
            except Exception:
                if attempt < MAX_RETRY - 1:
                    time.sleep(RETRY_DELAY * (attempt + 1))
    return pmids[:actual_max]


def fetch_abstracts(pmids: List[str]) -> List[Dict[str, str]]:
    if not pmids or not HAS_ENTREZ:
        return []
    Entrez.email = ENTREZ_EMAIL
    articles: List[Dict[str, str]] = []
    for i in range(0, len(pmids), BATCH_SIZE):
        batch = pmids[i:i + BATCH_SIZE]
        for attempt in range(MAX_RETRY):
            try:
                time.sleep(QUERY_DELAY_S)
                handle = Entrez.efetch(
                    db="pubmed", id=",".join(batch),
                    rettype="medline", retmode="text")
                for rec in Medline.parse(handle):
                    ab = rec.get("AB", "")
                    if ab:
                        articles.append({
                            "pmid": rec.get("PMID", ""),
                            "title": rec.get("TI", ""),
                            "abstract": ab,
                            "journal": rec.get("JT", ""),
                            "year": rec.get("DP", "")[:4] if rec.get("DP") else "",
                        })
                handle.close()
                break
            except Exception as e:
                if attempt < MAX_RETRY - 1:
                    time.sleep(RETRY_DELAY * (attempt + 1))
    return articles


# ================================================================
# Entity Recognition
# ================================================================

def _normalize(text: str) -> str:
    return text.lower().strip()

def _term_in_text(term: str, text_lower: str) -> bool:
    t = _normalize(term)
    if len(t) <= 3:
        return bool(re.search(r'\b' + re.escape(t) + r'\b', text_lower))
    return t in text_lower

def _make_node_id(name: str) -> str:
    return re.sub(r'[^a-z0-9]+', '_', name.lower()).strip('_')

def recognize_entities(
    text: str,
    bacteria_names: Dict[str, str],
    metabolite_names: Dict[str, str],
    disease_terms: List[str],
) -> Dict[str, Any]:
    text_lower = text.lower()
    found_bacteria: Set[str] = set()
    found_metabolites: Set[str] = set()
    found_disease = False

    for norm_name, orig_name in bacteria_names.items():
        search_terms = [norm_name]
        parts = orig_name.replace("_", " ").split()
        if parts:
            search_terms.append(parts[0].lower())
        for term in search_terms:
            if _term_in_text(term, text_lower):
                found_bacteria.add(orig_name)
                break

    for norm_name, orig_name in metabolite_names.items():
        if _term_in_text(norm_name, text_lower):
            found_metabolites.add(orig_name)

    for dt in disease_terms:
        if _term_in_text(dt, text_lower):
            found_disease = True
            break

    return {"bacteria": found_bacteria, "metabolites": found_metabolites,
            "disease": found_disease}


# ================================================================
# LLM Relation Extraction (rep1212 ===BEGIN NODES/EDGES=== format)
# ================================================================

_RE_SYSTEM_PROMPT = "You are a biomedical knowledge graph extraction expert."

_RE_USER_TEMPLATE = """Extract biological entities and their relationships from the following PubMed abstract.

===BEGIN CONTEXT===
Title: {title}
PMID: {pmid}
{abstract}
===END CONTEXT===

Entities already identified in this abstract:
{entity_list}

Output strictly in this format:

===BEGIN NODES===
id,type,description
[node_id],[node_type],[description without commas]
===END NODES===

===BEGIN EDGES===
source,target,relation,description
[source_id],[target_id],[relation_type],[evidence sentence without commas]
===END EDGES===

Requirements:
1. Node type must be one of: microbe / metabolite / pathway / disease
2. Node IDs must use lowercase_underscores (e.g. akkermansia_muciniphila)
3. Relation type must be one of: produces / metabolizes / inhibits / promotes / protects / associated_with / degrades / transports / modulates
4. Only extract relations explicitly supported by the abstract text
5. Descriptions must NOT contain commas (use semicolons instead)
6. If a metabolic pathway is mentioned in context include it as a pathway node
7. If no relations found output empty sections
"""


def _parse_csv_sections(response: str) -> Optional[Dict[str, Any]]:
    """Parse rep1212-style ===BEGIN NODES/EDGES=== CSV response."""
    if not response:
        return None

    # strip markdown code fences and <think> tags
    text = re.sub(r"<think>.*?</think>", "", response, flags=re.DOTALL)
    text = text.replace("```text", "```")
    if "```" in text:
        blocks = re.findall(r"```(?:\w*)\n(.*?)```", text, re.DOTALL)
        for b in blocks:
            if "===BEGIN NODES===" in b:
                text = b
                break

    nodes_match = re.search(
        r"===BEGIN NODES===\s*\n(.*?)===END NODES===", text, re.DOTALL)
    edges_match = re.search(
        r"===BEGIN EDGES===\s*\n(.*?)===END EDGES===", text, re.DOTALL)

    if not nodes_match and not edges_match:
        return None

    nodes: List[Dict[str, str]] = []
    edges: List[Dict[str, str]] = []

    if nodes_match:
        for line in nodes_match.group(1).strip().split("\n"):
            line = line.strip()
            if not line or line.startswith("id,") or line.startswith("["):
                continue
            parts = [p.strip() for p in line.split(",", 2)]
            if len(parts) >= 2:
                nid = _make_node_id(parts[0])
                ntype = parts[1].strip().lower()
                desc = parts[2] if len(parts) > 2 else ""
                label = parts[0].replace("_", " ")
                if ntype not in ("microbe", "metabolite", "pathway", "disease"):
                    ntype = "metabolite"
                nodes.append({"id": nid, "type": ntype,
                              "label": label, "description": desc})

    if edges_match:
        for line in edges_match.group(1).strip().split("\n"):
            line = line.strip()
            if not line or line.startswith("source,") or line.startswith("["):
                continue
            parts = [p.strip() for p in line.split(",", 3)]
            if len(parts) >= 3:
                src = _make_node_id(parts[0])
                tgt = _make_node_id(parts[1])
                rel = parts[2].strip().lower()
                desc = parts[3] if len(parts) > 3 else ""
                if rel not in VALID_RELATIONS:
                    rel = "associated_with"
                edges.append({"source": src, "target": tgt,
                              "relation": rel, "description": desc})

    if not nodes and not edges:
        return None
    return {"nodes": nodes, "edges": edges}


def extract_relations_llm(
    abstract: str, title: str, pmid: str,
    found_entities: Dict[str, Any], disease: str,
    pool: DeepSeekPool = None,
) -> Optional[Dict[str, Any]]:
    """Extract relations using DeepSeek API (rep1212 CSV format)."""
    entity_lines = []
    for bact in found_entities["bacteria"]:
        entity_lines.append(f"- {bact.replace('_', ' ')} (microbe)")
    for met in found_entities["metabolites"]:
        entity_lines.append(f"- {met} (metabolite)")
    if found_entities["disease"]:
        entity_lines.append(f"- {disease} (disease)")

    if len(entity_lines) < 1:
        return None

    user_prompt = _RE_USER_TEMPLATE.format(
        title=title, pmid=pmid, abstract=abstract,
        entity_list="\n".join(entity_lines))

    if pool is None:
        pool = _get_pool()

    response = pool.call(_RE_SYSTEM_PROMPT, user_prompt)
    return _parse_csv_sections(response)


# ================================================================
# KG Building
# ================================================================

def _add_to_graph(G: nx.DiGraph, result: Dict, article: Dict):
    pmid = article.get("pmid", "")
    for node in result.get("nodes", []):
        nid = node["id"]
        if nid not in G:
            G.add_node(nid, label=node.get("label", nid),
                       type=node.get("type", "unknown"),
                       node_type=node.get("type", "unknown"),
                       description=node.get("description", ""))
        elif not G.nodes[nid].get("description") and node.get("description"):
            G.nodes[nid]["description"] = node["description"]
    for edge in result.get("edges", []):
        src, tgt = edge["source"], edge["target"]
        if src not in G:
            G.add_node(src, label=src, type="unknown", node_type="unknown")
        if tgt not in G:
            G.add_node(tgt, label=tgt, type="unknown", node_type="unknown")
        if G.has_edge(src, tgt):
            old = G[src][tgt].get("pmid", "")
            if pmid and pmid not in old:
                G[src][tgt]["pmid"] = f"{old},{pmid}" if old else pmid
        else:
            G.add_edge(src, tgt,
                       relation=edge.get("relation", "associated_with"),
                       description=edge.get("description", ""),
                       pmid=pmid,
                       journal=article.get("journal", ""),
                       year=article.get("year", ""))


def build_kg_from_abstracts(
    articles: List[Dict[str, Any]],
    bacteria_names: Dict[str, str],
    metabolite_names: Dict[str, str],
    disease_terms: List[str],
    disease: str,
    pool: DeepSeekPool = None,
) -> nx.DiGraph:
    G = nx.DiGraph()
    llm_ok = llm_fail = skip_no_ent = 0

    for i, article in enumerate(articles):
        pmid = article["pmid"]
        cache_file = CACHE_DIR / "extractions" / f"{pmid}.json"
        if cache_file.exists():
            with open(cache_file, "r", encoding="utf-8") as f:
                result = json.load(f)
            if result.get("nodes") or result.get("edges"):
                _add_to_graph(G, result, article)
                llm_ok += 1
                continue

        text = f"{article.get('title', '')} {article['abstract']}"
        found = recognize_entities(
            text, bacteria_names, metabolite_names, disease_terms)
        if not found["bacteria"] and not found["metabolites"]:
            skip_no_ent += 1
            continue

        result = extract_relations_llm(
            article["abstract"], article.get("title", ""),
            pmid, found, disease, pool)

        if result is None:
            llm_fail += 1
            continue

        llm_ok += 1
        cache_file.parent.mkdir(parents=True, exist_ok=True)
        with open(cache_file, "w", encoding="utf-8") as f:
            json.dump(result, f, indent=2, ensure_ascii=False)
        _add_to_graph(G, result, article)

        if (i + 1) % 10 == 0:
            print(f"  {i+1}/{len(articles)} abstracts | "
                  f"graph: {G.number_of_nodes()} nodes, "
                  f"{G.number_of_edges()} edges | "
                  f"LLM OK:{llm_ok} fail:{llm_fail} skip:{skip_no_ent}")

    print(f"  Summary: LLM OK={llm_ok}, fail={llm_fail}, "
          f"skip(no entities)={skip_no_ent}")
    if llm_fail > 0 and llm_ok == 0:
        print("  [ERROR] All LLM calls failed. Check API keys.")
    return G


def _print_kg_stats(G: nx.DiGraph):
    print(f"\nFinal KG: {G.number_of_nodes()} nodes, {G.number_of_edges()} edges")
    type_counts: Dict[str, int] = {}
    for _, d in G.nodes(data=True):
        t = d.get("type", "unknown")
        type_counts[t] = type_counts.get(t, 0) + 1
    for t, c in sorted(type_counts.items(), key=lambda x: -x[1]):
        print(f"  {t}: {c}")
    rel_counts: Dict[str, int] = {}
    for _, _, d in G.edges(data=True):
        r = d.get("relation", "unknown")
        rel_counts[r] = rel_counts.get(r, 0) + 1
    if rel_counts:
        print("\nEdge relations:")
        for r, c in sorted(rel_counts.items(), key=lambda x: -x[1]):
            print(f"  {r}: {c}")


def _save_gml(G: nx.DiGraph, output_gml_path: str):
    out = Path(output_gml_path)
    out.parent.mkdir(parents=True, exist_ok=True)
    for nid in G.nodes():
        if "label" not in G.nodes[nid]:
            G.nodes[nid]["label"] = str(nid)
    nx.write_gml(G, str(out))
    print(f"\nSaved KG to {output_gml_path}")


# ================================================================
# Direct Entry: species + disease (no candidates needed)
# ================================================================

def run_build_kg_direct(
    species: str,
    disease: str,
    output_gml_path: str,
    keys_csv: str = str(DEFAULT_KEYS_CSV),
    max_articles: int = 40,
    max_queries: int = 10,
) -> nx.DiGraph:
    """
    Build KG directly from a species name + disease.
    No candidates JSON needed.

    Example:
        run_build_kg_direct("Akkermansia muciniphila", "IBD", "output.gml")
    """
    if not HAS_ENTREZ:
        raise RuntimeError("Biopython required: pip install biopython")

    print("=" * 60)
    print("  BUILD KG: species + disease -> PubMed -> DeepSeek -> GML")
    print("=" * 60)
    print(f"  Species: {species}")
    print(f"  Disease: {disease}")

    pool = DeepSeekPool(Path(keys_csv))

    # Phase 1: prepare aliases used only after articles have been fetched.
    print("\nPhase 1: Preparing entity aliases for abstract recognition...")
    bacteria_set = {species.replace(" ", "_")}
    metabolite_set: Set[str] = set()
    bacteria_names, metabolite_names, disease_terms = expand_entity_names(
        bacteria_set, metabolite_set, disease)

    # Phase 2: build queries directly from species + disease
    print("\nPhase 2: Generating PubMed queries...")

    def _build_term(name: str, etype: str = None) -> str:
        del etype  # Retained in the local signature for call-site compatibility.
        normalized = name.replace("_", " ").strip()
        return f'"{normalized}"'

    sp_term = _build_term(species, "microbe")
    dis_term = _build_term(disease, "disease")
    exp = EXPERIMENTAL_FILTER

    queries = [
        f'{sp_term} AND {dis_term} AND (metabolite OR metabolome OR microbiome){exp}',
        f'{sp_term} AND {dis_term} AND (gut OR intestinal){exp}',
        f'{sp_term} AND (butyrate OR propionate OR acetate OR "short-chain fatty acid" OR SCFA){exp}',
        f'{sp_term} AND (tryptophan OR indole OR bile acid){exp}',
        f'{sp_term} AND {dis_term}{exp}',
        f'{sp_term} AND (immune OR inflammation OR barrier OR mucin){exp}',
    ]
    queries = queries[:max_queries]
    print(f"  Generated {len(queries)} queries")

    # Phase 3: fetch
    print(f"\nPhase 3: Fetching PubMed abstracts (max {max_articles}/query)...")
    all_articles: List[Dict[str, str]] = []
    seen_pmids: Set[str] = set()
    for i, q in enumerate(queries):
        pmids = search_pubmed(q, max_results=max_articles)
        new = [p for p in pmids if p not in seen_pmids]
        if new:
            arts = fetch_abstracts(new)
            all_articles.extend(arts)
            seen_pmids.update(new)
        print(f"  Query {i+1}/{len(queries)}: +{len(new)} PMIDs, "
              f"total {len(all_articles)} articles")

    print(f"  Total: {len(all_articles)} unique articles")

    if not all_articles:
        print("  [ERROR] No articles fetched. Check PubMed connectivity.")
        G = nx.DiGraph()
        _save_gml(G, output_gml_path)
        return G

    # For NER: any metabolite mentioned in abstracts will be detected
    # via SapBERT index (if available). We start with empty metabolite_names
    # and let the LLM discover them.

    # Phase 4: extract relations
    print(f"\nPhase 4: Extracting relations via DeepSeek...")
    G = build_kg_from_abstracts(
        all_articles, bacteria_names, metabolite_names,
        disease_terms, disease, pool)

    # ---------------------------------------------------------
    # Merge disease aliases into the canonical disease node in direct mode.
    # ---------------------------------------------------------
    disease_id = _make_node_id(disease)
    if disease_id not in G:
        G.add_node(disease_id, label=disease, type="disease", node_type="disease")

    synonym_ids = set([_make_node_id(dt) for dt in disease_terms])
    
    rogues_to_merge = [node for node in list(G.nodes()) if node in synonym_ids and node != disease_id]

    for rogue in rogues_to_merge:
        for src, tgt, data in list(G.out_edges(rogue, data=True)):
            if tgt != disease_id: 
                G.add_edge(disease_id, tgt, **data)
        
        for src, tgt, data in list(G.in_edges(rogue, data=True)):
            if src != disease_id:
                G.add_edge(src, disease_id, **data)
                
        G.remove_node(rogue)

    G.nodes[disease_id]['label'] = disease
    G.nodes[disease_id]['type'] = 'disease'
    G.nodes[disease_id]['node_type'] = 'disease'
    # ---------------------------------------------------------

    _print_kg_stats(G)
    _save_gml(G, output_gml_path)
    return G


# ================================================================
# Pipeline Entry: candidates JSON (existing interface)
# ================================================================

def generate_search_queries(
    candidates: List[Dict], disease: str,
    max_queries: int = 100, experimental_only: bool = True,
) -> List[Dict[str, str]]:
    queries: List[Dict[str, str]] = []
    seen: Set[str] = set()
    bacteria_set = set(c["bacteria"] for c in candidates)
    metabolite_set = set(c["metabolite"] for c in candidates)

    def _bt(name, etype=None):
        del etype  # Retained in the local signature for call-site compatibility.
        normalized = name.replace("_", " ").strip()
        return f'"{normalized}"'

    dt = _bt(disease, "disease")
    ef = EXPERIMENTAL_FILTER if experimental_only else ""

    for bact in bacteria_set:
        k = f"{bact}|dis"
        if k not in seen:
            seen.add(k)
            queries.append({"query": f'{_bt(bact,"microbe")} AND {dt} AND '
                            f'(metabolite OR metabolome OR microbiome){ef}',
                            "type": "bacteria_disease", "bacteria": bact})

    for c in candidates[:60]:
        k = f"{c['bacteria']}|{c['metabolite']}"
        if k not in seen:
            seen.add(k)
            queries.append({"query": f'{_bt(c["bacteria"],"microbe")} AND '
                            f'{_bt(c["metabolite"],"metabolite")}{ef}',
                            "type": "pair", "bacteria": c["bacteria"],
                            "metabolite": c["metabolite"]})

    seen_m: Set[str] = set()
    for c in candidates[:40]:
        m = c["metabolite"]
        if m not in seen_m:
            seen_m.add(m)
            k = f"{m}|{disease}"
            if k not in seen:
                seen.add(k)
                queries.append({"query": f'{_bt(m,"metabolite")} AND {dt} AND '
                                f'(gut OR microbiome OR microbiota){ef}',
                                "type": "metabolite_disease", "metabolite": m})

    for bact in bacteria_set:
        k = f"{bact}|pw"
        if k not in seen:
            seen.add(k)
            queries.append({"query": f'{_bt(bact,"microbe")} AND (metabolite OR '
                            f'"short-chain fatty acid" OR SCFA){ef}',
                            "type": "bacteria_pathway", "bacteria": bact})

    return queries[:max_queries]


def run_build_kg(
    candidates_path: str, output_gml_path: str,
    disease: str = DEFAULT_DISEASE,
    max_articles_per_query: int = 40,
    max_queries: int = 100,
    experimental_only: bool = True,
    keys_csv: str = str(DEFAULT_KEYS_CSV),
) -> nx.DiGraph:
    """Build KG from step1 candidates (full pipeline integration)."""
    if not HAS_ENTREZ:
        raise RuntimeError("Biopython required: pip install biopython")

    print("=" * 60)
    print("  BUILD KNOWLEDGE GRAPH FROM LITERATURE (v3)")
    print("  DeepSeek API + SapBERT NER (rep1212 aligned)")
    print("=" * 60)

    pool = DeepSeekPool(Path(keys_csv))

    with open(candidates_path, "r", encoding="utf-8") as f:
        candidates = json.load(f)
    print(f"Loaded {len(candidates)} candidates")

    bacteria_set = set(c["bacteria"] for c in candidates)
    metabolite_set = set(c["metabolite"] for c in candidates)

    print("\nPhase 1: Preparing entity aliases for abstract recognition...")
    bacteria_names, metabolite_names, disease_terms = expand_entity_names(
        bacteria_set, metabolite_set, disease)

    print("\nPhase 2: PubMed queries...")
    queries = generate_search_queries(
        candidates, disease, max_queries, experimental_only)
    print(f"  Generated {len(queries)} queries")

    print(f"\nPhase 2b: Fetching (max {max_articles_per_query}/query)...")
    all_articles: List[Dict[str, str]] = []
    seen_pmids: Set[str] = set()
    for i, q in enumerate(queries):
        pmids = search_pubmed(q["query"], max_results=max_articles_per_query)
        new = [p for p in pmids if p not in seen_pmids]
        if new:
            all_articles.extend(fetch_abstracts(new))
            seen_pmids.update(new)
        if (i + 1) % 10 == 0:
            print(f"  {i+1}/{len(queries)} queries, "
                  f"{len(all_articles)} articles")

    print(f"  Total: {len(all_articles)} articles")

    print("\nPhase 3: DeepSeek relation extraction...")
    G = build_kg_from_abstracts(
        all_articles, bacteria_names, metabolite_names,
        disease_terms, disease, pool)

    # ---------------------------------------------------------
    # Merge disease aliases into the canonical disease node.
    # ---------------------------------------------------------
    did = _make_node_id(disease)
    if did not in G:
        G.add_node(did, label=disease, type="disease", node_type="disease")

    # Normalize SapBERT aliases to graph node IDs.
    synonym_ids = set([_make_node_id(dt) for dt in disease_terms])
    
    # Identify alias nodes already present in the graph.
    rogues_to_merge = [node for node in list(G.nodes()) if node in synonym_ids and node != did]

    # Transfer alias edges to the canonical disease node.
    for rogue in rogues_to_merge:
        # Transfer outgoing edges.
        for src, tgt, data in list(G.out_edges(rogue, data=True)):
            if tgt != did: # Avoid self-loops.
                G.add_edge(did, tgt, **data)
        
        # Transfer incoming edges.
        for src, tgt, data in list(G.in_edges(rogue, data=True)):
            if src != did:
                G.add_edge(src, did, **data)
                
        # Remove the merged alias node.
        G.remove_node(rogue)

    # Normalize the canonical disease node attributes.
    G.nodes[did]['label'] = disease
    G.nodes[did]['type'] = 'disease'
    G.nodes[did]['node_type'] = 'disease'
    # ---------------------------------------------------------

    _print_kg_stats(G)
    _save_gml(G, output_gml_path)
    return G


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(
        description="Build KG (v3: DeepSeek + SapBERT, rep1212 aligned)")
    sub = parser.add_subparsers(dest="mode")

    # Mode 1: direct (species + disease)
    p1 = sub.add_parser("direct", help="Species + disease -> GML")
    p1.add_argument("--species", required=True)
    p1.add_argument("--disease", default="IBD")
    p1.add_argument("--output", default=str(
        Path(__file__).parent.parent / "data" / "knowledge_graph" / "direct_kg.gml"))
    p1.add_argument("--keys-csv", default=str(DEFAULT_KEYS_CSV))
    p1.add_argument("--max-articles", type=int, default=40)
    p1.add_argument("--max-queries", type=int, default=10)

    # Mode 2: candidates (pipeline integration)
    p2 = sub.add_parser("candidates", help="Candidates JSON -> GML")
    p2.add_argument("--candidates", default="outputs/step1_candidates.json")
    p2.add_argument("--output", default=str(
        Path(__file__).parent.parent / "data" / "knowledge_graph" / "auto_built_kg.gml"))
    p2.add_argument("--disease", default="IBD")
    p2.add_argument("--keys-csv", default=str(DEFAULT_KEYS_CSV))
    p2.add_argument("--max-articles", type=int, default=40)
    p2.add_argument("--max-queries", type=int, default=100)
    p2.add_argument("--no-experimental-filter", action="store_true")

    args = parser.parse_args()

    if args.mode == "direct":
        run_build_kg_direct(
            args.species, args.disease, args.output,
            keys_csv=args.keys_csv,
            max_articles=args.max_articles,
            max_queries=args.max_queries)
    elif args.mode == "candidates":
        run_build_kg(
            args.candidates, args.output,
            disease=args.disease, keys_csv=args.keys_csv,
            max_articles_per_query=args.max_articles,
            max_queries=args.max_queries,
            experimental_only=not args.no_experimental_filter)
    else:
        parser.print_help()
