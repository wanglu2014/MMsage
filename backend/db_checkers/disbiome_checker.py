# ═══════════════════════════════════════════════════════════════════
# 【模块：Disbiome 数据库查询器】
# ═══════════════════════════════════════════════════════════════════
#
# 【承上启下】作为 db_checkers 的核心模块之一，本查询器负责从 Disbiome
#             数据库（微生物-疾病关联的权威数据源，10,648 条记录）中检索
#             微生物的已知关联。它被两个Agent调用：
#             - MicrobeMetaboliteAgent 调用 check_disbiome() 检查微生物是否有记录（间接证据）
#             - MicrobeDiseaseAgent 调用 check_disbiome_disease() 检查微生物-疾病直接关联
#
# 【核心思想】通过模糊名称匹配 + 疾病同义词扩展，在 TSV 文件中查找匹配记录。
#             疾病同义词扩展是关键设计——Disbiome 中存储的是具体疾病名
#             （如 "ulcerative colitis"），而用户输入的往往是缩写（如 "IBD"）。
#
# 【数据源】 H:\...\microbe\disbiome_edges.tsv
# 【列结构】 source_name(微生物) | edge_type(increase/decrease) | target_name(疾病)
#
# ───────────────────────────────────────────────────────────────────

import csv  # ← 导入 CSV 解析器，用于读取 TSV 文件
from pathlib import Path
from typing import Dict, List, Optional
from functools import lru_cache

try:
    from .sapbert_retriever import retrieve, retrieve_batch  # ← SapBERT 语义检索（包内导入）
except ImportError:
    from sapbert_retriever import retrieve, retrieve_batch  # ← 直接运行时的回退导入

# ═══════════════════════════════════════════════════════════════════
# 【模块 1：数据加载与缓存】
# ═══════════════════════════════════════════════════════════════════
#
# 【核心思想】首次调用时从磁盘加载全部记录到内存，后续调用直接返回缓存。
#             所有字段统一转为小写，便于后续模糊匹配。
#
# 【输入】 DISBIOME_PATH (TSV 文件路径)
# 【输出】 _disbiome_data (List[Dict], 全部记录的内存缓存)
#
# ───────────────────────────────────────────────────────────────────

DISBIOME_PATH = Path(__file__).parent.parent.parent / "data" / "databases" / "disbiome_edges.tsv"

_disbiome_data: Optional[List[Dict]] = None  # ← 模块级缓存，避免重复读取磁盘


def _load_disbiome() -> List[Dict]:
    """以 DISBIOME_PATH 为输入，解析 TSV 文件，生成标准化的记录列表（带缓存）。"""
    global _disbiome_data
    if _disbiome_data is not None:
        return _disbiome_data  # ← 缓存命中，直接返回

    _disbiome_data = []
    if not DISBIOME_PATH.exists():
        return _disbiome_data  # ← 文件不存在时返回空列表，不报错（优雅降级）

    with open(DISBIOME_PATH, 'r', encoding='utf-8') as f:
        reader = csv.DictReader(f, delimiter='\t')  # ← 以 TSV 格式解析，自动识别列名
        for row in reader:
            _disbiome_data.append({
                'microbe': (row.get('source_name') or '').strip().lower(),       # ← 微生物原始名（小写化）
                'microbe_std': (row.get('standard_sourcename') or '').strip().lower(),  # ← 微生物标准名
                'edge_type': (row.get('edge_type') or '').strip().lower(),       # ← 关系类型：increase/decrease
                'disease': (row.get('target_name') or '').strip().lower(),       # ← 疾病原始名
                'disease_std': (row.get('standard_targetname') or '').strip().lower(),  # ← 疾病标准名
                'pmid': row.get('PMID', ''),  # ← 文献 PMID（可能为空）
            })
    return _disbiome_data


# ═══════════════════════════════════════════════════════════════════
# 【模块 2：名称匹配工具函数】
# ═══════════════════════════════════════════════════════════════════
#
# 【核心思想】生物医学实体名称格式多样（下划线/空格/大小写混用），
#             需要归一化后做子串包含匹配。疾病名还需要同义词扩展。
#
# ───────────────────────────────────────────────────────────────────

def _fuzzy_match(query: str, target: str) -> bool:
    """以 query 和 target 为输入，归一化后执行双向子串包含匹配，返回是否匹配。"""
    q = query.lower().replace('_', ' ').replace('-', ' ').strip()
    t = target.lower().replace('_', ' ').replace('-', ' ').strip()
    if not q or not t:
        return False
    return q in t or t in q  # ← 双向匹配：query 包含在 target 中，或 target 包含在 query 中


# 疾病同义词扩展表：将常用缩写映射到数据库中实际使用的疾病全称
# 这是 Disbiome 查询的关键——数据库中存的是 "ulcerative colitis"，用户输入的是 "IBD"
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
    """以疾病查询词为输入，先尝试直接匹配，再尝试同义词扩展匹配，返回是否命中。"""
    if _fuzzy_match(query, target):
        return True  # ← 直接匹配成功
    q = query.lower().strip()
    synonyms = DISEASE_SYNONYMS.get(q, [])  # ← 查找同义词列表
    for syn in synonyms:
        if _fuzzy_match(syn, target):
            return True  # ← 同义词匹配成功（如 "IBD" → "crohn's disease" 匹配到数据库记录）
    return False


# ═══════════════════════════════════════════════════════════════════
# 【模块 3：微生物特征化查询（check_disbiome）】
# ═══════════════════════════════════════════════════════════════════
#
# 【承上启下】被 MicrobeMetaboliteAgent 调用。Disbiome 本身是微生物-疾病
#             数据库，不直接记录微生物-代谢物关系。但如果一个微生物在
#             Disbiome 中有记录，说明它是一个被研究过的微生物——这本身
#             就是间接证据，增加该微生物相关发现的可信度。
#
# 【输入】 bacteria (微生物名), metabolite (未使用，保持接口一致)
# 【输出】 {"hit": bool, "records": int, "details": list}
#
# ───────────────────────────────────────────────────────────────────

def check_disbiome(bacteria: str, metabolite: str) -> Dict:
    """
    两阶段匹配：先精确子串匹配，未命中则用 SapBERT 语义匹配找相近微生物名再查。
    """
    data = _load_disbiome()

    # Stage 1: 精确子串匹配
    matches = []
    for row in data:
        if _fuzzy_match(bacteria, row['microbe']) or _fuzzy_match(bacteria, row['microbe_std']):
            matches.append(row)

    match_method = "exact"

    # Stage 2: SapBERT 语义匹配（仅在 Stage 1 未命中时触发）
    if not matches:
        try:
            sapbert_hits = retrieve(bacteria, entity_type="microbe", top_k=3, threshold=0.7)
            for hit_name, score in sapbert_hits:
                for row in data:
                    if _fuzzy_match(hit_name, row['microbe']) or _fuzzy_match(hit_name, row['microbe_std']):
                        matches.append(row)
                if matches:
                    match_method = f"sapbert({hit_name},{score:.3f})"
                    break
        except Exception:
            pass

    return {
        "hit": len(matches) > 0,
        "records": len(matches),
        "match_method": match_method,
        "details": [f"{m['microbe']}--{m['edge_type']}-->{m['disease']}" for m in matches[:5]],
    }


# ═══════════════════════════════════════════════════════════════════
# 【模块 4：微生物-疾病直接关联查询（check_disbiome_disease）】
# ═══════════════════════════════════════════════════════════════════
#
# 【承上启下】被 MicrobeDiseaseAgent 调用。这是 Disbiome 的核心用法——
#             同时匹配微生物名和疾病名，查找直接的微生物-疾病关联。
#             使用疾病同义词扩展确保 "IBD" 能匹配到 "Crohn's disease"
#             和 "ulcerative colitis" 等具体疾病名。
#
# 【输入】 bacteria (微生物名), disease (疾病名/缩写)
# 【输出】 {"hit": bool, "records": int, "edge_types": list, "details": list}
#
# ───────────────────────────────────────────────────────────────────

def check_disbiome_disease(bacteria: str, disease: str) -> Dict:
    """两阶段匹配：先精确匹配，未命中则用 SapBERT 语义匹配。"""
    data = _load_disbiome()

    # Stage 1: 精确匹配
    matches = []
    for row in data:
        microbe_match = _fuzzy_match(bacteria, row['microbe']) or _fuzzy_match(bacteria, row['microbe_std'])
        disease_match = _disease_match(disease, row['disease']) or _disease_match(disease, row['disease_std'])
        if microbe_match and disease_match:
            matches.append(row)

    match_method = "exact"

    # Stage 2: SapBERT 语义匹配（仅在 Stage 1 未命中时触发）
    if not matches:
        try:
            microbe_hits = retrieve(bacteria, entity_type="microbe", top_k=3, threshold=0.7)
            disease_hits = retrieve(disease, entity_type="disease", top_k=3, threshold=0.7)
            # 用语义匹配到的名称重新查
            alt_microbes = [bacteria] + [h[0] for h in microbe_hits]
            alt_diseases = [disease] + [h[0] for h in disease_hits]
            for am in alt_microbes:
                for ad in alt_diseases:
                    for row in data:
                        microbe_match = _fuzzy_match(am, row['microbe']) or _fuzzy_match(am, row['microbe_std'])
                        disease_match = _disease_match(ad, row['disease']) or _disease_match(ad, row['disease_std'])
                        if microbe_match and disease_match:
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
        "details": [f"{m['microbe']}--{m['edge_type']}-->{m['disease']}" for m in matches[:5]],
    }
