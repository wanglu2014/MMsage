# ═══════════════════════════════════════════════════════════════════
# 【模块：CTD 化学物质-疾病关联查询器】
# ═══════════════════════════════════════════════════════════════════
#
# 【承上启下】被 MetaboliteDiseaseAgent 调用。CTD（Comparative Toxicogenomics
#             Database，比较毒理基因组学数据库）是化学物质-疾病关联的权威数据源，
#             包含 treat（治疗）、marker（标志物）等关系类型，以及推断分数和
#             文献 PMID。命中 CTD 是代谢物-疾病关系的强证据。
#
# 【核心思想】通过化学物质名模糊匹配 + 疾病同义词扩展，在 CTD TSV 文件中
#             查找代谢物与疾病的已知关联。
#
# 【数据源】 H:\...\27.ctdchemdis\CTD_chemicals_diseases.tsv (~750K 条记录)
# 【列结构】 source_name(化学物质) | target_name(疾病) | edge_type | PMID
#
# 【输入】 metabolite (代谢物/化学物质名), disease (疾病名/缩写)
# 【输出】 {"hit": bool, "records": int, "edge_types": list, "details": list}
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
#
# 【核心思想】首次调用时从磁盘加载全部 CTD 记录到内存，后续调用直接返回缓存。
#
# ───────────────────────────────────────────────────────────────────

CTD_PATH = Path(__file__).parent.parent.parent / "data" / "databases" / "CTD_chemicals_diseases.tsv"

_ctd_data: Optional[List[Dict]] = None  # ← 模块级缓存


def _load_ctd() -> List[Dict]:
    """以 CTD_PATH 为输入，解析 TSV 文件，生成标准化的记录列表（带缓存）。"""
    global _ctd_data
    if _ctd_data is not None:
        return _ctd_data  # ← 缓存命中

    _ctd_data = []
    if not CTD_PATH.exists():
        return _ctd_data  # ← 文件不存在时优雅降级

    with open(CTD_PATH, 'r', encoding='utf-8') as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            _ctd_data.append({
                'chemical': (row.get('source_name') or '').strip().lower(),  # ← 化学物质名（小写化）
                'disease': (row.get('target_name') or '').strip().lower(),   # ← 疾病名（小写化）
                'edge_type': (row.get('edge_type') or '').strip().lower(),   # ← 关系类型：treat/marker 等
                'inference_score': row.get('edge_description', ''),          # ← 推断分数（CTD 特有）
                'pmids': row.get('PMID', ''),  # ← 支持文献的 PMID 列表
            })
    return _ctd_data


# ═══════════════════════════════════════════════════════════════════
# 【模块 2：名称匹配工具函数】
# ═══════════════════════════════════════════════════════════════════
#
# 【核心思想】与 disbiome_checker 共享相同的匹配策略：
#             模糊子串匹配 + 疾病同义词扩展。
#
# ───────────────────────────────────────────────────────────────────

def _fuzzy_match(query: str, target: str) -> bool:
    """以 query 和 target 为输入，归一化后执行双向子串包含匹配。"""
    q = query.lower().replace('_', ' ').replace('-', ' ').strip()
    t = target.lower().replace('_', ' ').replace('-', ' ').strip()
    if not q or not t:
        return False
    return q in t or t in q


DISEASE_SYNONYMS = {
    'ibd': ['inflammatory bowel disease', 'crohn', "crohn's disease", 'ulcerative colitis'],
    'uc': ['ulcerative colitis'],
    'cd': ["crohn's disease", 'crohn disease'],
    't2d': ['type 2 diabetes', 'diabetes mellitus type 2'],
    'crc': ['colorectal cancer', 'colon cancer'],
    'nafld': ['non-alcoholic fatty liver', 'nonalcoholic fatty liver'],
    'nash': ['nonalcoholic steatohepatitis', 'non-alcoholic steatohepatitis'],
}


def _disease_match(query: str, target: str) -> bool:
    """以疾病查询词为输入，先直接匹配，再同义词扩展匹配。"""
    if _fuzzy_match(query, target):
        return True
    q = query.lower().strip()
    synonyms = DISEASE_SYNONYMS.get(q, [])
    for syn in synonyms:
        if _fuzzy_match(syn, target):
            return True
    return False


# ═══════════════════════════════════════════════════════════════════
# 【模块 3：CTD 化学物质-疾病关联查询】
# ═══════════════════════════════════════════════════════════════════
#
# 【承上启下】被 MetaboliteDiseaseAgent 调用。同时匹配化学物质名和疾病名，
#             返回所有匹配的 CTD 记录。CTD 的 edge_type 字段（如 treat、
#             marker）提供了关系的具体类型，是高质量的结构化证据。
#
# 【输入】 metabolite (代谢物名), disease (疾病名/缩写)
# 【输出】 {"hit": bool, "records": int, "edge_types": list, "details": list}
#
# ───────────────────────────────────────────────────────────────────

def check_ctd(metabolite: str, disease: str) -> Dict:
    """两阶段匹配：先精确匹配，未命中则用 SapBERT 语义匹配。"""
    data = _load_ctd()

    # Stage 1: 精确匹配
    matches = []
    for row in data:
        chem_match = _fuzzy_match(metabolite, row['chemical'])
        dis_match = _disease_match(disease, row['disease'])
        if chem_match and dis_match:
            matches.append(row)

    match_method = "exact"

    # Stage 2: SapBERT 语义匹配（仅在 Stage 1 未命中时触发）
    if not matches:
        try:
            met_hits = retrieve(metabolite, entity_type="metabolite", top_k=3, threshold=0.7)
            dis_hits = retrieve(disease, entity_type="disease", top_k=3, threshold=0.7)
            alt_mets = [metabolite] + [h[0] for h in met_hits]
            alt_diss = [disease] + [h[0] for h in dis_hits]
            for am in alt_mets:
                for ad in alt_diss:
                    for row in data:
                        if _fuzzy_match(am, row['chemical']) and _disease_match(ad, row['disease']):
                            matches.append(row)
                    if matches:
                        match_method = f"sapbert(m={am},d={ad})"
                        break
                if matches:
                    break
        except Exception:
            pass

    edge_types = list(set(m['edge_type'] for m in matches if m['edge_type']))

    return {
        "hit": len(matches) > 0,
        "records": len(matches),
        "match_method": match_method,
        "edge_types": edge_types,
        "details": [f"{m['chemical']}--{m['edge_type']}-->{m['disease']}" for m in matches[:5]],
    }
