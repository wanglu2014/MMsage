# ═══════════════════════════════════════════════════════════════════
# 【模块：微生物-疾病 Agent (Microbe-Disease Agent)】
# ═══════════════════════════════════════════════════════════════════
#
# 【承上启下】作为三跳推理链的辅助验证跳，本Agent负责收集"微生物↔疾病"
#             的直接关联证据。它不参与 chain_novelty 的主计算（主计算只用
#             第一跳和第二跳的 min 值），但提供重要的交叉验证信息：
#             如果微生物-疾病共现很高但 chain_novelty 也很高，说明中间的
#             代谢物环节是新发现——这正是我们要找的"暗物质"。
#
# 【核心思想】通过微生物-疾病的直接证据，为链条的新颖性判断提供参照基线。
#
# 【输入】 bacteria (微生物名), disease (疾病名)
# 【输出】 HopEvidence (含 pubmed_count, db_hits, sources)
#
# ───────────────────────────────────────────────────────────────────

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent))  # ← 将 backend/ 加入搜索路径

from .hop_evidence import HopEvidence  # ← 导入统一的证据数据结构
from .pubmed_query import build_pubmed_query, query_pubmed_count, DEFAULT_SYNONYMS  # ← 导入 PubMed 查询工具
from db_checkers.disbiome_checker import check_disbiome_disease  # ← 导入 Disbiome 微生物-疾病直接关联查询器


class MicrobeDiseaseAgent:
    """微生物-疾病关系证据收集Agent。负责辅助验证跳的证据查询。"""

    def __init__(self, synonym_map=None):
        self.synonym_map = synonym_map or DEFAULT_SYNONYMS  # ← 初始化同义词扩展词典

    def run(self, bacteria: str, disease: str) -> HopEvidence:
        """
        执行微生物-疾病证据查询。

        【输入】 bacteria (微生物名), disease (疾病名)
        【输出】 HopEvidence (辅助验证跳的完整证据包)
        """

        # ═══════════════════════════════════════════════════════════════
        # 【子模块 1：PubMed 文献共现查询】
        # ═══════════════════════════════════════════════════════════════
        #
        # 【核心思想】量化文献中对该微生物-疾病直接关系的研究程度。
        #             例如 Akkermansia+IBD 有 158 篇文献，说明该微生物
        #             与该疾病的关系已被广泛研究。
        #
        query = build_pubmed_query([bacteria, disease], self.synonym_map)  # ← 以微生物名和疾病名为输入，构建带同义词扩展的 PubMed 查询
        pubmed_count = query_pubmed_count(query)  # ← 以查询字符串为输入，调用 Entrez API，获取匹配文章数

        # ═══════════════════════════════════════════════════════════════
        # 【子模块 2：Disbiome 微生物-疾病直接关联查询】
        # ═══════════════════════════════════════════════════════════════
        #
        # 【核心思想】Disbiome 是专门记录微生物-疾病关联的数据库，包含
        #             increase（在疾病中增多）和 decrease（在疾病中减少）
        #             两种关系类型。内置疾病同义词扩展（IBD → Crohn's/UC）。
        #
        disbiome_result = check_disbiome_disease(bacteria, disease)  # ← 以微生物名和疾病名为输入，在 disbiome_edges.tsv 中双重模糊匹配

        # ═══════════════════════════════════════════════════════════════
        # 【子模块 3：证据来源汇总】
        # ═══════════════════════════════════════════════════════════════
        sources = []
        if pubmed_count > 0:
            sources.append(f"PubMed: {pubmed_count} articles for [{bacteria}] AND [{disease}]")
        if disbiome_result['hit']:
            sources.append(f"Disbiome: {disbiome_result['records']} records ({', '.join(disbiome_result['edge_types'][:3])})")

        return HopEvidence(  # ← 封装为标准 HopEvidence 结构返回
            hop_type="microbe_disease",
            pubmed_count=pubmed_count,
            db_hits={
                "disbiome": disbiome_result['hit'],
            },
            db_details={
                "disbiome": disbiome_result,
            },
            sources=sources,
            query_used=query,
        )
