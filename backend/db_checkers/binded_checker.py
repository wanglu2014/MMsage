# ═══════════════════════════════════════════════════════════════════
# 【模块：Binded Database 微生物-代谢物关联查询器】
# ═══════════════════════════════════════════════════════════════════
#
# 【承上启下】作为第一跳Agent（MicrobeMetaboliteAgent）的核心数据源，
#             本模块查询 binded_database_anno_named.csv（41.8万条记录），
#             该文件整合了 AGORA、Mambo、Pamet、WOM 四个微生物代谢模型数据库，
#             记录了微生物能产生/代谢哪些化合物（通过 PubChem CID 关联）。
#
# 【核心思想】两阶段匹配：
#   Stage 1: 精确子串匹配（快速，确定性）
#   Stage 2: SapBERT 语义匹配（处理名称变体，如 butyrate ↔ butyric acid）
#
# 【输入】 microbe (微生物名), metabolite (代谢物名, 可选)
# 【输出】 dict: {hit, records, sources, compounds/microbes, match_details}
#
# ───────────────────────────────────────────────────────────────────

import csv
from pathlib import Path
from typing import Dict, List, Optional

# ═══════════════════════════════════════════════════════════════════
# 【子模块 1：数据加载与缓存】
# ═══════════════════════════════════════════════════════════════════

BINDED_DB_PATH = Path(__file__).parent.parent.parent / "data" / "databases" / "binded_database_anno_named.csv"

_data: Optional[List[Dict]] = None
_sapbert_available: Optional[bool] = None  # None=untested, True/False=cached result


def _load_data():
    """首次调用时加载全量数据到内存。"""
    global _data
    if _data is not None:
        return
    if not BINDED_DB_PATH.exists():
        print(f"[BindedDB] File not found: {BINDED_DB_PATH}")
        _data = []
        return
    _data = []
    with open(BINDED_DB_PATH, "r", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for row in reader:
            _data.append(row)
    print(f"[BindedDB] Loaded {len(_data)} records")


# ═══════════════════════════════════════════════════════════════════
# 【子模块 2：模糊匹配工具】
# ═══════════════════════════════════════════════════════════════════

def _normalize(s: str) -> str:
    """归一化实体名：小写 + 下划线/连字符转空格 + 去首尾空格。"""
    return s.lower().replace("_", " ").replace("-", " ").strip()


def _fuzzy_match(query: str, target: str) -> bool:
    """双向子串包含匹配。"""
    q = _normalize(query)
    t = _normalize(target)
    if not q or not t:
        return False
    return q in t or t in q


# ═══════════════════════════════════════════════════════════════════
# 【子模块 3：查询接口】
# ═══════════════════════════════════════════════════════════════════

def check_binded_microbe_metabolite(microbe: str, metabolite: str) -> Dict:
    """
    查询微生物-代谢物是否在 binded_database 中有记录。

    【两阶段匹配】
    Stage 1: 精确子串匹配（species.name + Compound_Name）
    Stage 2: 如果 Stage 1 无结果，尝试 SapBERT 语义匹配

    【输入】 microbe (微生物名), metabolite (代谢物名)
    【输出】 dict: {hit, records, sources, match_type, match_details}
    """
    _load_data()

    # ── Stage 1: 精确子串匹配 ──
    matches = []
    for row in _data:
        sp = row.get("species.name", "")
        cn = row.get("Compound_Name", "")
        if not sp or not cn or cn == "NA":
            continue
        if _fuzzy_match(microbe, sp) and _fuzzy_match(metabolite, cn):
            matches.append({
                "species": sp,
                "compound": cn,
                "cid": row.get("cid.ID", ""),
                "source": row.get("sourcename", ""),
                "genus": row.get("genus.name", ""),
            })

    if matches:
        sources = set(m["source"] for m in matches)
        return {
            "hit": True,
            "records": len(matches),
            "sources": list(sources),
            "match_method": "exact",
            "match_details": matches[:10],
        }

    # ── Stage 2: SapBERT 语义匹配 ──
    global _sapbert_available
    if _sapbert_available is False:
        return {"hit": False, "records": 0, "sources": [], "match_method": "none", "match_details": []}
    try:
        from .sapbert_retriever import retrieve
        _sapbert_available = True
        # 分别检索微生物和代谢物的语义近邻
        microbe_hits = retrieve(microbe, entity_type="microbe", top_k=5, threshold=0.85)
        metabolite_hits = retrieve(metabolite, entity_type="metabolite", top_k=5, threshold=0.85)

        if not microbe_hits or not metabolite_hits:
            return {"hit": False, "records": 0, "sources": [], "match_method": "none", "match_details": []}

        # 用语义近邻的名称重新做精确匹配
        microbe_names = [name for name, _ in microbe_hits]
        metabolite_names = [name for name, _ in metabolite_hits]

        semantic_matches = []
        for row in _data:
            sp = row.get("species.name", "")
            cn = row.get("Compound_Name", "")
            if not sp or not cn or cn == "NA":
                continue
            sp_match = any(_fuzzy_match(mn, sp) for mn in microbe_names)
            cn_match = any(_fuzzy_match(mn, cn) for mn in metabolite_names)
            if sp_match and cn_match:
                semantic_matches.append({
                    "species": sp,
                    "compound": cn,
                    "cid": row.get("cid.ID", ""),
                    "source": row.get("sourcename", ""),
                    "genus": row.get("genus.name", ""),
                })

        if semantic_matches:
            sources = set(m["source"] for m in semantic_matches)
            return {
                "hit": True,
                "records": len(semantic_matches),
                "sources": list(sources),
                "match_method": "semantic",
                "match_details": semantic_matches[:10],
            }
    except Exception as e:
        _sapbert_available = False
        print(f"[BindedDB] SapBERT unavailable ({e}), skipping semantic matching for all future calls")

    return {"hit": False, "records": 0, "sources": [], "match_method": "none", "match_details": []}


def check_binded_microbe(microbe: str) -> Dict:
    """
    查询微生物在 binded_database 中关联了哪些代谢物。

    【输入】 microbe (微生物名)
    【输出】 dict: {hit, records, compounds, sources}
    """
    _load_data()

    compounds = set()
    sources = set()
    count = 0
    for row in _data:
        sp = row.get("species.name", "")
        if not sp:
            continue
        if _fuzzy_match(microbe, sp):
            cn = row.get("Compound_Name", "")
            if cn and cn != "NA":
                compounds.add(cn)
            sources.add(row.get("sourcename", ""))
            count += 1

    return {
        "hit": count > 0,
        "records": count,
        "compounds": sorted(compounds)[:50],  # 最多返回 50 个
        "sources": list(sources),
    }
