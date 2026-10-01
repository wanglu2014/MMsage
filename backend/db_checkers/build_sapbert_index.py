#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
SapBERT Embedding Index Builder
================================
从 5 个数据源（Disbiome, CTD, ChEBI, BacDive, binded_database）中提取所有
唯一实体名，用 SapBERT 编码后存为 .npy + mapping JSON，供 RAG 检索使用。

用法:
    C:/Users/WLPC/.conda/envs/mindmap_env/python.exe build_sapbert_index.py

输出:
    sapbert_index/
        entity_embeddings.npy     (N x 768 float32)
        entity_mapping.json       (index -> {name, source, type})
"""

import csv
import json
import sys
import numpy as np
from pathlib import Path
from typing import Dict, Set, Tuple

# ── 数据源路径 ──────────────────────────────────────────────────
_DATA_DIR = Path(__file__).parent.parent.parent / "data" / "databases"

SOURCES = {
    "disbiome": _DATA_DIR / "disbiome_edges.tsv",
    "ctd": _DATA_DIR / "CTD_chemicals_diseases.tsv",
    "chebi": _DATA_DIR / "chebi_relationships_standardized.csv",
    "bacdive": _DATA_DIR / "bacdive_edges.tsv",
    "binded": _DATA_DIR / "binded_database_anno_named.csv",
}

OUTPUT_DIR = Path(__file__).parent / "sapbert_index"


def extract_entities() -> Dict[str, Set[str]]:
    """从各数据源提取唯一实体名，按类型分组。"""
    microbes = set()
    metabolites = set()
    diseases = set()

    # ── Disbiome: 微生物 + 疾病 ──
    p = SOURCES["disbiome"]
    if p.exists():
        with open(p, "r", encoding="utf-8") as f:
            for row in csv.DictReader(f, delimiter="\t"):
                for col in ("source_name", "standard_sourcename"):
                    v = (row.get(col) or "").strip()
                    if v:
                        microbes.add(v)
                for col in ("target_name", "standard_targetname"):
                    v = (row.get(col) or "").strip()
                    if v:
                        diseases.add(v)
        print(f"[Disbiome] microbes={len(microbes)}, diseases={len(diseases)}")

    # ── CTD: 化学物质 + 疾病 ──
    p = SOURCES["ctd"]
    if p.exists():
        with open(p, "r", encoding="utf-8") as f:
            for row in csv.DictReader(f, delimiter="\t"):
                v = (row.get("source_name") or "").strip()
                if v:
                    metabolites.add(v)
                v = (row.get("target_name") or "").strip()
                if v:
                    diseases.add(v)
        print(f"[CTD] metabolites={len(metabolites)}, diseases={len(diseases)}")

    # ── ChEBI: 代谢物 ──
    p = SOURCES["chebi"]
    if p.exists():
        with open(p, "r", encoding="utf-8") as f:
            for row in csv.DictReader(f):
                for col in ("source_name", "target_name"):
                    v = (row.get(col) or "").strip()
                    if v:
                        metabolites.add(v)
        print(f"[ChEBI] metabolites={len(metabolites)}")

    # ── BacDive: 微生物 ──
    p = SOURCES["bacdive"]
    if p.exists():
        with open(p, "r", encoding="utf-8") as f:
            for row in csv.DictReader(f, delimiter="\t"):
                for col in ("source_name", "standard_sourcename"):
                    v = (row.get(col) or "").strip()
                    if v:
                        microbes.add(v)
        print(f"[BacDive] microbes={len(microbes)}")

    # ── binded_database: 微生物 + 代谢物 ──
    p = SOURCES["binded"]
    if p.exists():
        with open(p, "r", encoding="utf-8") as f:
            for row in csv.DictReader(f):
                sp = (row.get("species.name") or "").strip()
                if sp:
                    microbes.add(sp)
                cn = (row.get("Compound_Name") or "").strip()
                if cn and cn != "NA":
                    metabolites.add(cn)
        print(f"[binded] microbes={len(microbes)}, metabolites={len(metabolites)}")

    print(f"\n=== Total unique entities ===")
    print(f"  microbes:    {len(microbes)}")
    print(f"  metabolites: {len(metabolites)}")
    print(f"  diseases:    {len(diseases)}")
    print(f"  TOTAL:       {len(microbes) + len(metabolites) + len(diseases)}")

    return {"microbe": microbes, "metabolite": metabolites, "disease": diseases}


def build_index(entities: Dict[str, Set[str]], batch_size: int = 256):
    """用 SapBERT 编码所有实体，存为 .npy + mapping JSON。"""
    from sentence_transformers import SentenceTransformer

    print("\nLoading SapBERT model...")
    model = SentenceTransformer("cambridgeltl/SapBERT-from-PubMedBERT-fulltext")
    print(f"Model loaded. Embedding dim: {model.get_sentence_embedding_dimension()}")

    # 构建 (name, type) 列表
    all_entities = []
    for etype, names in entities.items():
        for name in sorted(names):
            all_entities.append((name, etype))

    total = len(all_entities)
    print(f"\nEncoding {total} entities in batches of {batch_size}...")

    names = [e[0] for e in all_entities]
    embeddings = model.encode(names, batch_size=batch_size, show_progress_bar=True,
                              normalize_embeddings=True)  # L2 归一化，cosine sim = dot product

    # 保存
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    np.save(str(OUTPUT_DIR / "entity_embeddings.npy"), embeddings.astype(np.float32))

    mapping = {}
    for i, (name, etype) in enumerate(all_entities):
        mapping[str(i)] = {"name": name, "type": etype}

    with open(OUTPUT_DIR / "entity_mapping.json", "w", encoding="utf-8") as f:
        json.dump(mapping, f, ensure_ascii=False, indent=1)

    print(f"\nSaved to {OUTPUT_DIR}/")
    print(f"  entity_embeddings.npy: {embeddings.shape}")
    print(f"  entity_mapping.json: {len(mapping)} entries")


if __name__ == "__main__":
    entities = extract_entities()
    build_index(entities)
