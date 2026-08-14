import json
import numpy as np
from pathlib import Path
from typing import List, Tuple, Optional, Dict

INDEX_DIR = Path(__file__).parent / "sapbert_index"

_embeddings: Optional[np.ndarray] = None  # Shape: (N, 768), float32, L2-normalized
_mapping: Optional[Dict] = None  # Maps each index to its name and type.
_type_indices: Optional[Dict[str, np.ndarray]] = None  # Maps types to index arrays.


def _load_index():
    """Load and cache the precomputed SapBERT embedding index."""
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

    _embeddings = np.load(str(emb_path))  # Shape: (N, 768)
    with open(map_path, "r", encoding="utf-8") as f:
        _mapping = json.load(f)

    # Build index arrays by entity type.
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


_model = None


def _get_model():
    """Load the SapBERT model lazily and reuse it."""
    global _model
    if _model is None:
        from sentence_transformers import SentenceTransformer
        _model = SentenceTransformer("cambridgeltl/SapBERT-from-PubMedBERT-fulltext")
        print("[SapBERT] Model loaded.")
    return _model


def encode_query(query: str) -> np.ndarray:
    """Return an L2-normalized embedding for a query."""
    model = _get_model()
    emb = model.encode([query], normalize_embeddings=True)
    return emb[0]


def _exact_match(query: str, entity_type: str = None) -> List[str]:
    """Return exact normalized name matches."""
    _load_index()
    if not _mapping:
        return []

    q = query.lower().replace("_", " ").replace("-", " ").strip()
    results = []

    indices = range(len(_mapping)) if entity_type is None else _type_indices.get(entity_type, [])
    for idx in indices:
        info = _mapping[str(int(idx))]
        name = info["name"].lower().replace("_", " ").replace("-", " ").strip()
        if q == name:
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
    Retrieve exact normalized matches first, then SapBERT semantic matches.

    Results are returned as ``(entity_name, similarity_score)`` pairs in
    descending similarity order.
    """
    _load_index()
    if _embeddings is None or len(_embeddings) == 0:
        return []

    results = []

    # Stage 1: exact matching
    if exact_first:
        exact_hits = _exact_match(query, entity_type)
        if exact_hits:
            return [(name, 1.0) for name in exact_hits[:top_k]]

    # Stage 2: SapBERT semantic matching
    query_emb = encode_query(query)  # Shape: (768,)

    # Filter by entity type.
    if entity_type and entity_type in _type_indices:
        indices = _type_indices[entity_type]
        sub_emb = _embeddings[indices]  # Shape: (M, 768)
    else:
        indices = np.arange(len(_embeddings))
        sub_emb = _embeddings

    # Cosine similarity equals the dot product after L2 normalization.
    scores = sub_emb @ query_emb  # Shape: (M,)

    # Top-K
    top_indices = np.argsort(scores)[::-1][:top_k * 2]  # Overfetch before deduplication.

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


def retrieve_batch(
    queries: List[str],
    entity_type: str = None,
    top_k: int = 5,
    threshold: float = 0.7,
) -> Dict[str, List[Tuple[str, float]]]:
    """Retrieve matches for a batch of query terms."""
    _load_index()
    if _embeddings is None or len(_embeddings) == 0:
        return {q: [] for q in queries}

    model = _get_model()
    query_embs = model.encode(queries, normalize_embeddings=True)  # Shape: (Q, 768)

    if entity_type and entity_type in _type_indices:
        indices = _type_indices[entity_type]
        sub_emb = _embeddings[indices]
    else:
        indices = np.arange(len(_embeddings))
        sub_emb = _embeddings

    # Compute all pairwise dot products in one batch.
    all_scores = query_embs @ sub_emb.T  # Shape: (Q, M)

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
