# ═══════════════════════════════════════════════════════════════════
# 【模块：BacDive 微生物特征化查询器】
# ═══════════════════════════════════════════════════════════════════
#
# 【承上启下】被 MicrobeMetaboliteAgent 调用。BacDive（Bacterial Diversity
#             Metadatabase）是微生物表型特征化的权威数据库，记录了微生物的
#             生长环境（inhabit）、表达酶（express）、培养基（medium）等信息。
#
# 【核心思想】BacDive 不直接记录微生物-代谢物关系，但如果一个微生物在
#             BacDive 中有丰富的表型记录，说明它是一个被充分研究的菌种。
#             这是间接证据——被充分特征化的微生物更可能有可靠的代谢功能注释。
#             例如 Akkermansia muciniphila 在 BacDive 中有 44 条记录，
#             包括其生长环境（#gastrointestinal tract）和表达的酶
#             （beta-galactosidase）。
#
# 【数据源】 H:\...\microbe\bacdive_edges.tsv
# 【列结构】 source_name(微生物) | edge_type(inhabit/express) | target_name(环境/酶/培养基)
#
# 【输入】 bacteria (微生物名), metabolite (未使用，保持接口一致)
# 【输出】 {"hit": bool, "records": int, "environments": list, "details": list}
#
# ───────────────────────────────────────────────────────────────────

import csv  # ← 导入 CSV 解析器，用于读取 TSV 文件
from pathlib import Path
from typing import Dict, List, Optional

try:
    from .sapbert_retriever import retrieve
except ImportError:
    from sapbert_retriever import retrieve

# ═══════════════════════════════════════════════════════════════════
# 【模块 1：数据加载与缓存】
# ═══════════════════════════════════════════════════════════════════
#
# 【核心思想】首次调用时从磁盘加载全部 BacDive 记录到内存，后续调用
#             直接返回缓存。所有字段统一转为小写，便于后续模糊匹配。
#
# 【输入】 BACDIVE_PATH (TSV 文件路径)
# 【输出】 _bacdive_data (List[Dict], 全部记录的内存缓存)
#
# ───────────────────────────────────────────────────────────────────

BACDIVE_PATH = Path(__file__).parent.parent.parent / "data" / "databases" / "bacdive_edges.tsv"

_bacdive_data: Optional[List[Dict]] = None  # ← 模块级缓存，避免重复读取磁盘


def _load_bacdive() -> List[Dict]:
    """以 BACDIVE_PATH 为输入，解析 TSV 文件，生成标准化的记录列表（带缓存）。"""
    global _bacdive_data
    if _bacdive_data is not None:
        return _bacdive_data  # ← 缓存命中，直接返回

    _bacdive_data = []
    if not BACDIVE_PATH.exists():
        return _bacdive_data  # ← 文件不存在时返回空列表，不报错（优雅降级）

    with open(BACDIVE_PATH, 'r', encoding='utf-8') as f:
        reader = csv.DictReader(f, delimiter='\t')  # ← 以 TSV 格式解析
        for row in reader:
            _bacdive_data.append({
                'microbe': (row.get('source_name') or '').strip().lower(),       # ← 微生物原始名（小写化）
                'microbe_std': (row.get('standard_sourcename') or '').strip().lower(),  # ← 微生物标准名
                'edge_type': (row.get('edge_type') or '').strip().lower(),       # ← 关系类型：inhabit/express
                'target': (row.get('target_name') or '').strip().lower(),        # ← 目标实体（环境/酶/培养基）
                'target_attr': (row.get('target_attribute') or '').strip().lower(),  # ← 目标属性
            })
    return _bacdive_data


# ═══════════════════════════════════════════════════════════════════
# 【模块 2：名称匹配工具函数】
# ═══════════════════════════════════════════════════════════════════
#
# 【核心思想】与其他 checker 共享相同的模糊匹配策略：
#             归一化（下划线→空格）后做双向子串包含匹配。
#
# ───────────────────────────────────────────────────────────────────

def _fuzzy_match(query: str, target: str) -> bool:
    """以 query 和 target 为输入，归一化后执行双向子串包含匹配，返回是否匹配。"""
    q = query.lower().replace('_', ' ').replace('-', ' ').strip()
    t = target.lower().replace('_', ' ').replace('-', ' ').strip()
    if not q or not t:
        return False
    return q in t or t in q  # ← 双向匹配


# ═══════════════════════════════════════════════════════════════════
# 【模块 3：BacDive 微生物特征化查询】
# ═══════════════════════════════════════════════════════════════════
#
# 【承上启下】在 BacDive 中查找目标微生物的所有表型记录。
#             返回匹配记录数和关联的环境/酶/培养基列表。
#             例如查询 "Akkermansia" 会返回 44 条记录，包括：
#             - inhabit → #gastrointestinal tract, #feces (stool)
#             - express → beta-galactosidase, n-acetyl-beta-glucosaminidase
#             - inhabit → medium 187, columbia blood medium 等培养基
#
# 【输入】 bacteria (微生物名), metabolite (未使用，保持接口一致)
# 【输出】 {"hit": bool, "records": int, "environments": list, "details": list}
#
# ───────────────────────────────────────────────────────────────────

def check_bacdive(bacteria: str, metabolite: str) -> Dict:
    """两阶段匹配：先精确子串匹配，未命中则用 SapBERT 语义匹配。"""
    data = _load_bacdive()

    # Stage 1: 精确匹配
    matches = []
    for row in data:
        if _fuzzy_match(bacteria, row['microbe']) or _fuzzy_match(bacteria, row['microbe_std']):
            matches.append(row)

    match_method = "exact"

    # Stage 2: SapBERT 语义匹配（仅在 Stage 1 未命中时触发）
    if not matches:
        try:
            microbe_hits = retrieve(bacteria, entity_type="microbe", top_k=3, threshold=0.7)
            for hit_name, score in microbe_hits:
                for row in data:
                    if _fuzzy_match(hit_name, row['microbe']) or _fuzzy_match(hit_name, row['microbe_std']):
                        matches.append(row)
                if matches:
                    match_method = f"sapbert({hit_name},{score:.3f})"
                    break
        except Exception:
            pass

    environments = list(set(m['target'] for m in matches[:20] if m['target']))

    return {
        "hit": len(matches) > 0,
        "records": len(matches),
        "match_method": match_method,
        "environments": environments[:10],
        "details": [f"{m['microbe']}--{m['edge_type']}-->{m['target']}" for m in matches[:5]],
    }
