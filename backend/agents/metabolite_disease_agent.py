# ═══════════════════════════════════════════════════════════════════
# 【模块：代谢物-疾病 Agent (Metabolite-Disease Agent)】
# ═══════════════════════════════════════════════════════════════════
#
# 【承上启下】作为三跳推理链的第二跳，本Agent负责收集"代谢物↔疾病"
#             这一环节的所有证据。它查询 PubMed 获取文献共现数，查询
#             CTD（比较毒理基因组学数据库）确认化学物质-疾病的已知关联，
#             查询 ChEBI 确认代谢物是否有已知的化学关系网络。
#             返回的 HopEvidence 将传递给 EvidenceAggregator 进行聚合。
#
# 【核心思想】从文献和化学/毒理学数据库两个维度量化代谢物与疾病之间的
#             已知关联程度。CTD 提供直接的化学物质-疾病证据，ChEBI 提供
#             代谢物的化学特征化程度（间接证据）。
#
# 【输入】 metabolite (代谢物名), disease (疾病名)
# 【输出】 HopEvidence (含 pubmed_count, db_hits, sources)
#
# ───────────────────────────────────────────────────────────────────

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent))  # ← 将 backend/ 加入搜索路径

from .hop_evidence import HopEvidence  # ← 导入统一的证据数据结构
from .pubmed_query import build_pubmed_query, query_pubmed_count, DEFAULT_SYNONYMS  # ← 导入 PubMed 查询工具
from db_checkers.ctd_checker import check_ctd  # ← 导入 CTD 化学物质-疾病关联查询器
from db_checkers.chebi_checker import check_chebi  # ← 导入 ChEBI 代谢物关系查询器


class MetaboliteDiseaseAgent:
    """代谢物-疾病关系证据收集Agent。负责链条第二跳的证据查询。"""

    def __init__(self, synonym_map=None):
        self.synonym_map = synonym_map or DEFAULT_SYNONYMS  # ← 初始化同义词扩展词典

    def run(self, metabolite: str, disease: str) -> HopEvidence:
        """
        执行代谢物-疾病证据查询。

        【输入】 metabolite (代谢物名), disease (疾病名)
        【输出】 HopEvidence (第二跳的完整证据包)
        """

        # ═══════════════════════════════════════════════════════════════
        # 【子模块 1：PubMed 文献共现查询】
        # ═══════════════════════════════════════════════════════════════
        #
        # 【核心思想】量化文献中对该代谢物-疾病关系的研究程度。
        #             例如 butyrate+IBD 有 612 篇文献，说明该关系已被充分研究。
        #
        query = build_pubmed_query([metabolite, disease], self.synonym_map)  # ← 以代谢物名和疾病名为输入，构建带同义词扩展的 PubMed 查询
        pubmed_count = query_pubmed_count(query)  # ← 以查询字符串为输入，调用 Entrez API，获取匹配文章数

        # ═══════════════════════════════════════════════════════════════
        # 【子模块 2：CTD 化学物质-疾病关联查询】
        # ═══════════════════════════════════════════════════════════════
        #
        # 【核心思想】CTD 是权威的化学物质-疾病关联数据库，包含 treat（治疗）、
        #             marker（标志物）等关系类型。命中 CTD 是强证据。
        #             内置疾病同义词扩展（IBD → Crohn's/UC/colitis）。
        #
        ctd_result = check_ctd(metabolite, disease)  # ← 以代谢物名和疾病名为输入，在 CTD_chemicals_diseases.tsv 中模糊匹配

        # ═══════════════════════════════════════════════════════════════
        # 【子模块 3：ChEBI 代谢物特征化查询】
        # ═══════════════════════════════════════════════════════════════
        #
        # 【核心思想】ChEBI 记录代谢物之间的化学关系（is_a、enantiomer_of 等）。
        #             如果代谢物在 ChEBI 中有丰富的关系网络，说明它是一个
        #             被充分特征化的化合物（间接证据，增加可信度）。
        #
        chebi_result = check_chebi(metabolite, disease)  # ← 以代谢物名为输入，在 chebi_relationships_standardized.csv 中模糊匹配

        # ═══════════════════════════════════════════════════════════════
        # 【子模块 4：证据来源汇总】
        # ═══════════════════════════════════════════════════════════════
        sources = []
        if pubmed_count > 0:
            sources.append(f"PubMed: {pubmed_count} articles for [{metabolite}] AND [{disease}]")
        if ctd_result['hit']:
            sources.append(f"CTD: {ctd_result['records']} records ({', '.join(ctd_result['edge_types'][:3])})")
        if chebi_result['hit']:
            sources.append(f"ChEBI: {chebi_result['records']} relationships for {metabolite}")

        return HopEvidence(  # ← 封装为标准 HopEvidence 结构返回
            hop_type="metabolite_disease",
            pubmed_count=pubmed_count,
            db_hits={
                "ctd": ctd_result['hit'],
                "chebi": chebi_result['hit'],
            },
            db_details={
                "ctd": ctd_result,
                "chebi": chebi_result,
            },
            sources=sources,
            query_used=query,
        )
