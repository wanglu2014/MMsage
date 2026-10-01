# ═══════════════════════════════════════════════════════════════════
# 【模块：ChEBI 代谢物关系查询器】
# ═══════════════════════════════════════════════════════════════════
#
# 【承上启下】被 MetaboliteDiseaseAgent 调用。ChEBI（Chemical Entities of
#             Biological Interest）是代谢物化学关系的权威数据库，记录了
#             代谢物之间的 is_a（分类）、enantiomer_of（对映体）、
#             has_functional_parent（功能母体）等化学关系。
#
# 【核心思想】ChEBI 不直接记录代谢物-疾病关系，但如果一个代谢物在 ChEBI
#             中有丰富的化学关系网络，说明它是一个被充分特征化的化合物。
#             这是间接证据——被充分特征化的代谢物更可能有可靠的生物学功能注释。
#
# 【数据源】 H:\...\3.chebi\chebi_relationships_standardized.csv (~180K 条记录)
# 【列结构】 source_name(代谢物A) | edge_type(contain/association) | target_name(代谢物B)
#
# 【输入】 metabolite (代谢物名), disease (未使用，ChEBI 无疾病信息)
# 【输出】 {"hit": bool, "records": int, "related_metabolites": list, "details": list}
#
# ───────────────────────────────────────────────────────────────────

import csv
from pathlib import Path
from typing import Dict, List, Optional

try:
    from .sapbert_retriever import retrieve
except ImportError:
    from sapbert_retriever import retrieve

# ═══════════════════════════════════════════════════════════════════
# 【模块 1：数据加载与缓存】
# ═══════════════════════════════════════════════════════════════════

CHEBI_PATH = Path(__file__).parent.parent.parent / "data" / "databases" / "chebi_relationships_standardized.csv"

_chebi_data: Optional[List[Dict]] = None  # ← 模块级缓存


def _load_chebi() -> List[Dict]:
    """以 CHEBI_PATH 为输入，解析 CSV 文件，生成标准化的记录列表（带缓存）。"""
    global _chebi_data
    if _chebi_data is not None:
        return _chebi_data

    _chebi_data = []
    if not CHEBI_PATH.exists():
        return _chebi_data

    with open(CHEBI_PATH, 'r', encoding='utf-8') as f:
        reader = csv.DictReader(f)  # ← CSV 格式（非 TSV），自动识别列名
        for row in reader:
            _chebi_data.append({
                'source': (row.get('source_name') or '').strip().lower(),     # ← 源代谢物名
                'target': (row.get('target_name') or '').strip().lower(),     # ← 目标代谢物名
                'edge_type': (row.get('edge_type') or '').strip().lower(),    # ← 关系类型：contain/association
                'description': row.get('edge_description', ''),               # ← 关系描述（含原始 ChEBI 关系名）
            })
    return _chebi_data


# ═══════════════════════════════════════════════════════════════════
# 【模块 2：名称匹配工具函数】
# ═══════════════════════════════════════════════════════════════════

def _fuzzy_match(query: str, target: str) -> bool:
    """以 query 和 target 为输入，归一化（含连字符处理）后执行双向子串匹配。"""
    q = query.lower().replace('_', ' ').replace('-', ' ').strip()
    t = target.lower().replace('_', ' ').replace('-', ' ').strip()
    if not q or not t:
        return False
    return q in t or t in q


# ═══════════════════════════════════════════════════════════════════
# 【模块 3：ChEBI 代谢物关系查询】
# ═══════════════════════════════════════════════════════════════════
#
# 【承上启下】在 ChEBI 的代谢物-代谢物关系网络中查找目标代谢物。
#             匹配 source 或 target 列，提取关联的其他代谢物列表。
#             例如 butyrate 在 ChEBI 中有 348 条关系记录，关联到
#             3-hydroxybutyrate、poly(3-hydroxybutyrate) 等化合物。
#
# 【输入】 metabolite (代谢物名), disease (未使用)
# 【输出】 {"hit": bool, "records": int, "related_metabolites": list, "details": list}
#
# ───────────────────────────────────────────────────────────────────

def check_chebi(metabolite: str, disease: str = "") -> Dict:
    """两阶段匹配：先精确子串匹配，未命中则用 SapBERT 语义匹配。"""
    data = _load_chebi()

    # Stage 1: 精确匹配
    matches = []
    for row in data:
        if _fuzzy_match(metabolite, row['source']) or _fuzzy_match(metabolite, row['target']):
            matches.append(row)

    match_method = "exact"

    # Stage 2: SapBERT 语义匹配（仅在 Stage 1 未命中时触发）
    if not matches:
        try:
            met_hits = retrieve(metabolite, entity_type="metabolite", top_k=3, threshold=0.7)
            for hit_name, score in met_hits:
                for row in data:
                    if _fuzzy_match(hit_name, row['source']) or _fuzzy_match(hit_name, row['target']):
                        matches.append(row)
                if matches:
                    match_method = f"sapbert({hit_name},{score:.3f})"
                    break
        except Exception:
            pass

    related = list(set(
        m['target'] if _fuzzy_match(metabolite, m['source']) else m['source']
        for m in matches[:20]
    ))

    return {
        "hit": len(matches) > 0,
        "records": len(matches),
        "match_method": match_method,
        "related_metabolites": related[:10],
        "details": [f"{m['source']}--{m['edge_type']}-->{m['target']}" for m in matches[:5]],
    }
