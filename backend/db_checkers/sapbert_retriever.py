# ═══════════════════════════════════════════════════════════════════
# 【模块：SapBERT 语义检索器】
# ═══════════════════════════════════════════════════════════════════
#
# 【承上启下】被各 db_checker 调用，替代原有的模糊子串匹配。
#             使用预计算的 SapBERT embedding 索引，通过 cosine similarity
#             实现跨数据库的实体名语义匹配。
#
# 【核心思想】生物医学实体名变体极多（butyrate / butyric acid / n-butyrate），
#             SapBERT 在 UMLS 上预训练，能将同义实体映射到相近的向量空间位置。
#             预计算所有数据库实体的 embedding，查询时只需编码查询词 + dot product。
#
# 【参考】 AgentConc/rep1212/sapbert_matcher.py 的两阶段匹配策略
#
# 【输入】 query (实体名), entity_type (microbe/metabolite/disease), top_k
# 【输出】 List[(name, similarity_score)]
#
# ───────────────────────────────────────────────────────────────────

import json
import numpy as np
from pathlib import Path
from typing import List, Tuple, Optional, Dict

# ═══════════════════════════════════════════════════════════════════
# 【模块 1：索引加载与缓存】
# ═══════════════════════════════════════════════════════════════════

INDEX_DIR = Path(__file__).parent / "sapbert_index"

_embeddings: Optional[np.ndarray] = None  # ← (N, 768) float32, L2 归一化
_mapping: Optional[Dict] = None           # ← {index: {name, type}}
_type_indices: Optional[Dict[str, np.ndarray]] = None  # ← {type: array of indices}


def _load_index():
    """加载预计算的 SapBERT embedding 索引（带缓存）。"""
    global _embeddings, _mapping, _type_indices
    if _embeddings is not None:
        return

    emb_path = INDEX_DIR / "entity_embeddings.npy"
    map_path = INDEX_DIR / "entity_mapping.json"

    if not emb_path.exists() or not map_path.exists():
        print(f"[SapBERT] Index not found at {INDEX_DIR}. Run build_sapbert_index.py first.")
        _embeddings = np.array([])
        _mapping = {}
        _type_indices = {}
        return

    _embeddings = np.load(str(emb_path))  # ← (N, 768)
    with open(map_path, "r", encoding="utf-8") as f:
        _mapping = json.load(f)

    # 按实体类型建立索引，加速按类型过滤的检索
    _type_indices = {}
    for idx_str, info in _mapping.items():
        etype = info.get("type", "unknown")
        if etype not in _type_indices:
            _type_indices[etype] = []
        _type_indices[etype].append(int(idx_str))

    for etype in _type_indices:
        _type_indices[etype] = np.array(_type_indices[etype], dtype=np.int64)

    print(f"[SapBERT] Loaded index: {_embeddings.shape[0]} entities, "
          f"types: { {k: len(v) for k, v in _type_indices.items()} }")


# ═══════════════════════════════════════════════════════════════════
# 【模块 2：SapBERT 模型加载（惰性）】
# ═══════════════════════════════════════════════════════════════════

_model = None


def _get_model():
    """惰性加载 SapBERT 模型（首次查询时加载，后续复用）。"""
    global _model
    if _model is None:
        from sentence_transformers import SentenceTransformer
        _model = SentenceTransformer("cambridgeltl/SapBERT-from-PubMedBERT-fulltext")
        print("[SapBERT] Model loaded.")
    return _model


def encode_query(query: str) -> np.ndarray:
    """以查询字符串为输入，返回 L2 归一化的 embedding 向量 (768,)。"""
    model = _get_model()
    emb = model.encode([query], normalize_embeddings=True)
    return emb[0]


# ═══════════════════════════════════════════════════════════════════
# 【模块 3：两阶段检索】
# ═══════════════════════════════════════════════════════════════════
#
# 【承上启下】参考 AgentConc/rep1212/step2_node_matcher.py 的两阶段策略：
#   Stage 1: 精确字符串匹配（快速，零误差）
#   Stage 2: SapBERT 语义匹配（处理同义词/变体）
#
# ───────────────────────────────────────────────────────────────────

def _exact_match(query: str, entity_type: str = None) -> List[str]:
    """Stage 1: 精确子串匹配（归一化后）。"""
    _load_index()
    if not _mapping:
        return []

    q = query.lower().replace("_", " ").replace("-", " ").strip()
    results = []

    indices = range(len(_mapping)) if entity_type is None else _type_indices.get(entity_type, [])
    for idx in indices:
        info = _mapping[str(int(idx))]
        name = info["name"].lower().replace("_", " ").replace("-", " ").strip()
        if q == name or q in name or name in q:
            results.append(info["name"])

    return results


def retrieve(
    query: str,
    entity_type: str = None,
    top_k: int = 10,
    threshold: float = 0.7,
    exact_first: bool = True,
) -> List[Tuple[str, float]]:
    """
    两阶段检索：先精确匹配，再 SapBERT 语义匹配。

    【输入】
      query: 实体名（如 "butyrate", "Akkermansia muciniphila"）
      entity_type: 过滤类型（"microbe" / "metabolite" / "disease" / None=全部）
      top_k: 返回前 K 个结果
      threshold: 最低 cosine similarity 阈值
      exact_first: 是否优先返回精确匹配结果

    【输出】 List[(entity_name, similarity_score)]，按相似度降序
    """
    _load_index()
    if _embeddings is None or len(_embeddings) == 0:
        return []

    results = []

    # Stage 1: 精确匹配
    if exact_first:
        exact_hits = _exact_match(query, entity_type)
        for name in exact_hits[:top_k]:
            results.append((name, 1.0))
        if len(results) >= top_k:
            return results[:top_k]

    # Stage 2: SapBERT 语义匹配
    query_emb = encode_query(query)  # ← (768,)

    # 按类型过滤
    if entity_type and entity_type in _type_indices:
        indices = _type_indices[entity_type]
        sub_emb = _embeddings[indices]  # ← (M, 768)
    else:
        indices = np.arange(len(_embeddings))
        sub_emb = _embeddings

    # cosine similarity = dot product（因为已 L2 归一化）
    scores = sub_emb @ query_emb  # ← (M,)

    # Top-K
    top_indices = np.argsort(scores)[::-1][:top_k * 2]  # 多取一些，去重后截断

    exact_names = {r[0] for r in results}
    for idx in top_indices:
        score = float(scores[idx])
        if score < threshold:
            break
        real_idx = int(indices[idx])
        name = _mapping[str(real_idx)]["name"]
        if name not in exact_names:
            results.append((name, score))
            exact_names.add(name)
        if len(results) >= top_k:
            break

    return results[:top_k]


# ═══════════════════════════════════════════════════════════════════
# 【模块 4：批量检索（供 db_checker 使用）】
# ═══════════════════════════════════════════════════════════════════

def retrieve_batch(
    queries: List[str],
    entity_type: str = None,
    top_k: int = 5,
    threshold: float = 0.7,
) -> Dict[str, List[Tuple[str, float]]]:
    """批量检索，对多个查询词一次性编码后并行计算相似度。"""
    _load_index()
    if _embeddings is None or len(_embeddings) == 0:
        return {q: [] for q in queries}

    model = _get_model()
    query_embs = model.encode(queries, normalize_embeddings=True)  # ← (Q, 768)

    if entity_type and entity_type in _type_indices:
        indices = _type_indices[entity_type]
        sub_emb = _embeddings[indices]
    else:
        indices = np.arange(len(_embeddings))
        sub_emb = _embeddings

    # 批量 dot product
    all_scores = query_embs @ sub_emb.T  # ← (Q, M)

    results = {}
    for qi, query in enumerate(queries):
        scores = all_scores[qi]
        top_idx = np.argsort(scores)[::-1][:top_k]
        hits = []
        for idx in top_idx:
            score = float(scores[idx])
            if score < threshold:
                break
            real_idx = int(indices[idx])
            hits.append((_mapping[str(real_idx)]["name"], score))
        results[query] = hits

    return results
