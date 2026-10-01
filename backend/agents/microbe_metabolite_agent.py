# ═══════════════════════════════════════════════════════════════════
# 【模块：微生物-代谢物 Agent (Microbe-Metabolite Agent)】
# ═══════════════════════════════════════════════════════════════════
#
# 【承上启下】作为三跳推理链的第一跳，本Agent负责收集"微生物↔代谢物"
#             这一环节的所有证据。它查询 PubMed 获取文献共现数，查询
#             Disbiome 确认微生物是否有已知的疾病关联记录（间接表征），
#             查询 BacDive 确认微生物是否有培养/代谢特征数据。
#             返回的 HopEvidence 将传递给 EvidenceAggregator 进行聚合。
#
# 【核心思想】从文献和数据库两个维度量化微生物与代谢物之间的已知关联程度。
#
# 【输入】 bacteria (微生物名), metabolite (代谢物名)
# 【输出】 HopEvidence (含 pubmed_count, db_hits, sources)
#
# ───────────────────────────────────────────────────────────────────

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent))  # ← 将 backend/ 加入搜索路径，确保 db_checkers 可导入

from .hop_evidence import HopEvidence  # ← 导入统一的证据数据结构
from .pubmed_query import build_pubmed_query, query_pubmed_count, DEFAULT_SYNONYMS  # ← 导入 PubMed 查询工具和同义词表
from db_checkers.disbiome_checker import check_disbiome  # ← 导入 Disbiome 微生物特征查询器
from db_checkers.bacdive_checker import check_bacdive  # ← 导入 BacDive 微生物环境/代谢查询器
from db_checkers.binded_checker import check_binded_microbe_metabolite  # ← 导入 Binded Database 微生物-代谢物直接关联查询器（AGORA/Mambo/Pamet/WOM）


class MicrobeMetaboliteAgent:
    """微生物-代谢物关系证据收集Agent。负责链条第一跳的证据查询。"""

    def __init__(self, synonym_map=None):
        self.synonym_map = synonym_map or DEFAULT_SYNONYMS  # ← 以用户自定义或默认同义词表为输入，初始化 PubMed 查询扩展词典

    def run(self, bacteria: str, metabolite: str) -> HopEvidence:
        """
        执行微生物-代谢物证据查询。

        【输入】 bacteria (微生物名), metabolite (代谢物名)
        【输出】 HopEvidence (第一跳的完整证据包)
        """

        # ═══════════════════════════════════════════════════════════════
        # 【子模块 1：PubMed 文献共现查询】
        # ═══════════════════════════════════════════════════════════════
        #
        # 【核心思想】通过 PubMed 两两共现计数，量化文献中对该微生物-代谢物
        #             关系的研究程度。共现数越高，说明该关系越被充分研究。
        #
        query = build_pubmed_query([bacteria, metabolite], self.synonym_map)  # ← 以微生物名和代谢物名为输入，应用同义词扩展，构建 PubMed AND 查询字符串
        pubmed_count = query_pubmed_count(query)  # ← 以查询字符串为输入，调用 Entrez API（带缓存），获取匹配文章数

        # ═══════════════════════════════════════════════════════════════
        # 【子模块 2：Disbiome 数据库查询】
        # ═══════════════════════════════════════════════════════════════
        #
        # 【核心思想】Disbiome 记录了微生物-疾病的直接关联。如果该微生物在
        #             Disbiome 中有记录，说明它是一个被研究过的微生物（间接证据）。
        #
        disbiome_result = check_disbiome(bacteria, metabolite)  # ← 以微生物名为输入，在 disbiome_edges.tsv 中模糊匹配，返回命中情况和记录数

        # ═══════════════════════════════════════════════════════════════
        # 【子模块 3：BacDive 数据库查询】
        # ═══════════════════════════════════════════════════════════════
        #
        # 【核心思想】BacDive 记录了微生物的培养条件和环境特征。如果该微生物
        #             在 BacDive 中有记录，说明它有实验室培养数据（间接证据）。
        #
        bacdive_result = check_bacdive(bacteria, metabolite)  # ← 以微生物名为输入，在 bacdive_edges.tsv 中模糊匹配，返回命中情况和环境列表

        # ═══════════════════════════════════════════════════════════════
        # 【子模块 3b：Binded Database 查询（AGORA/Mambo/Pamet/WOM）】
        # ═══════════════════════════════════════════════════════════════
        #
        # 【核心思想】binded_database 整合了 4 个微生物代谢模型数据库（41.8万条），
        #             记录了微生物能产生/代谢哪些化合物。这是第一跳最直接的证据——
        #             如果微生物-代谢物对在代谢模型中有记录，说明存在已知的代谢关系。
        #             使用 SapBERT 语义匹配处理名称变体（butyrate ↔ butyric acid）。
        #
        binded_result = check_binded_microbe_metabolite(bacteria, metabolite)  # ← 两阶段匹配：精确子串 + SapBERT 语义

        # ═══════════════════════════════════════════════════════════════
        # 【子模块 4：证据来源汇总】
        # ═══════════════════════════════════════════════════════════════
        #
        # 【核心思想】将各数据源的查询结果整理为人类可读的描述列表，
        #             供前端展示和报告生成使用。
        #
        sources = []
        if pubmed_count > 0:
            sources.append(f"PubMed: {pubmed_count} articles for [{bacteria}] AND [{metabolite}]")  # ← 记录 PubMed 命中描述
        if disbiome_result['hit']:
            sources.append(f"Disbiome: {disbiome_result['records']} records for {bacteria}")  # ← 记录 Disbiome 命中描述
        if bacdive_result['hit']:
            sources.append(f"BacDive: {bacdive_result['records']} records for {bacteria}")  # ← 记录 BacDive 命中描述
        if binded_result['hit']:
            sources.append(f"BindedDB: {binded_result['records']} records ({','.join(binded_result.get('sources',[])[:3])})")  # ← 记录 Binded Database 命中描述

        return HopEvidence(  # ← 以上述所有查询结果为输入，封装为标准 HopEvidence 结构，返回给 MasterAgent
            hop_type="microbe_metabolite",
            pubmed_count=pubmed_count,
            db_hits={
                "disbiome": disbiome_result['hit'],
                "bacdive": bacdive_result['hit'],
                "binded": binded_result['hit'],
            },
            db_details={
                "disbiome": disbiome_result,
                "bacdive": bacdive_result,
                "binded": binded_result,
            },
            sources=sources,
            query_used=query,
        )
