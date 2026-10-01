# ═══════════════════════════════════════════════════════════════════
# 【模块：证据聚合器 (Evidence Aggregator)】
# ═══════════════════════════════════════════════════════════════════
#
# 【承上启下】继三个子Agent分别返回各自的 HopEvidence 后，本模块负责
#             将三跳证据聚合为一个综合的 chain_novelty 分数。这是整个
#             多Agent系统的核心计算模块，其输出将直接用于四象限分类。
#
# 【核心思想】一条推理链的新颖度取决于最薄弱的一环（min策略），
#             数据库命中提供额外的可信度加分（db_bonus）。
#
# 【核心公式】
#   chain_count = min(count_mm, count_md)     ← 取两跳中的最小值
#   db_bonus = 10 × (命中数据库数), 上限50   ← 数据库加分（5个数据库）
#   total_count = chain_count + db_bonus
#   chain_novelty = 1 - log(1 + total_count) / log(1 + C_max)
#
# 【输入】 hop_results: List[HopEvidence] (三个子Agent的返回值)
# 【输出】 dict (含 chain_novelty, bottleneck, hop_counts, recommendation 等)
#
# ───────────────────────────────────────────────────────────────────

import math  # ← 导入数学库，用于对数计算
from typing import List, Dict
from .hop_evidence import HopEvidence  # ← 导入统一的证据数据结构


class EvidenceAggregator:
    """证据聚合器：将三跳证据合并为综合 chain_novelty 分数。"""

    def __init__(self, c_max: int = 500):
        self.c_max = c_max  # ← 归一化上限常数，控制 novelty 公式的压缩程度。500 表示共现数达到 500 时 novelty 趋近于 0

    def aggregate(self, hop_results: List[HopEvidence]) -> Dict:
        """
        聚合三跳证据，计算综合 chain_novelty。

        【输入】 hop_results (三个 HopEvidence 实例)
        【输出】 dict (完整的聚合结果，含 chain_novelty, bottleneck 等)
        """

        # ═══════════════════════════════════════════════════════════════
        # 【子模块 1：提取各跳 PubMed 共现数】
        # ═══════════════════════════════════════════════════════════════
        #
        # 【核心思想】从三个 HopEvidence 中提取 pubmed_count，按跳类型索引。
        #
        counts = {h.hop_type: h.pubmed_count for h in hop_results}  # ← 以 hop_results 为输入，构建 {跳类型: 共现数} 映射字典
        mm_count = counts.get("microbe_metabolite", 0)  # ← 提取微生物-代谢物共现数，缺失时默认为 0
        md_count = counts.get("metabolite_disease", 0)   # ← 提取代谢物-疾病共现数
        bd_count = counts.get("microbe_disease", 0)       # ← 提取微生物-疾病共现数（辅助，不参与主计算）

        # ═══════════════════════════════════════════════════════════════
        # 【子模块 2：计算 chain_count（最薄弱环节策略）】
        # ═══════════════════════════════════════════════════════════════
        #
        # 【核心思想】取两跳中的最小值。如果微生物→代谢物有 500 篇文献，
        #             但代谢物→疾病只有 2 篇，整条链就是新颖的（瓶颈在第二跳）。
        #
        chain_count = min(mm_count, md_count)  # ← 以两跳共现数为输入，取最小值，得到链条的原始共现强度

        # ═══════════════════════════════════════════════════════════════
        # 【子模块 3：计算数据库命中加分】
        # ═══════════════════════════════════════════════════════════════
        #
        # 【核心思想】每命中一个数据库加 10 分（固定加分，非加权），上限 50。
        #             数据库命中表示该关系有结构化数据支持，增加可信度。
        #             v7: 新增 binded_database（AGORA/Mambo/Pamet/WOM），共 5 个数据库。
        #
        db_bonus = 0
        for h in hop_results:  # ← 遍历三个 HopEvidence
            for db_name, hit in h.db_hits.items():  # ← 遍历每个Agent查询的数据库
                if hit:
                    db_bonus += 10  # ← 每命中一个数据库 +10
        db_bonus = min(db_bonus, 50)  # ← 5 个数据库上限 50

        # ═══════════════════════════════════════════════════════════════
        # 【子模块 4：计算最终 chain_novelty】
        # ═══════════════════════════════════════════════════════════════
        #
        # 【核心思想】使用对数归一化公式，将共现数压缩到 [0, 1] 区间。
        #             公式与原始 step2 的 Chain Novelty 完全一致，保持兼容性。
        #
        total_count = chain_count + db_bonus  # ← 以原始共现数和数据库加分为输入，求和得到最终共现强度
        if self.c_max <= 0:
            chain_novelty = 1.0  # ← 退化情况：C_max=0 时所有候选对 novelty 均为 1.0
        else:
            chain_novelty = 1.0 - math.log(1 + total_count) / math.log(1 + self.c_max)  # ← 以 total_count 和 C_max 为输入，应用对数归一化公式
            chain_novelty = max(0.0, min(1.0, chain_novelty))  # ← 裁剪到 [0, 1] 区间，防止数值溢出

        # ═══════════════════════════════════════════════════════════════
        # 【子模块 5：识别瓶颈环节】
        # ═══════════════════════════════════════════════════════════════
        #
        # 【核心思想】共现数更低的那一跳就是瓶颈——实验验证应优先针对该环节。
        #
        if mm_count <= md_count:
            bottleneck = "microbe_metabolite"  # ← 微生物-代谢物环节证据更少，是瓶颈
        else:
            bottleneck = "metabolite_disease"  # ← 代谢物-疾病环节证据更少，是瓶颈

        # ═══════════════════════════════════════════════════════════════
        # 【子模块 6：汇总所有证据来源】
        # ═══════════════════════════════════════════════════════════════
        all_sources = []
        all_db_details = {}
        for h in hop_results:  # ← 遍历三个 HopEvidence，合并 sources 和 db_details
            all_sources.extend(h.sources)
            all_db_details.update(h.db_details)

        return {  # ← 以上述所有计算结果为输入，组装为完整的聚合结果字典
            "chain_count": total_count,
            "chain_count_raw": chain_count,
            "chain_novelty": round(chain_novelty, 4),
            "bottleneck": bottleneck,
            "hop_counts": {
                "microbe_metabolite": mm_count,
                "metabolite_disease": md_count,
                "microbe_disease": bd_count,
            },
            "db_bonus": db_bonus,
            "db_hits": {h.hop_type: h.db_hits for h in hop_results},
            "sources": all_sources,
            "db_details": all_db_details,
            "recommendation": self._recommend(chain_novelty, bottleneck),
        }

    def _recommend(self, novelty: float, bottleneck: str) -> str:
        """
        根据 novelty 分数和瓶颈位置生成研究建议。

        【输入】 novelty (0~1), bottleneck (瓶颈跳类型)
        【输出】 str (人类可读的研究建议)
        """
        bn_label = {
            "microbe_metabolite": "microbe-metabolite",
            "metabolite_disease": "metabolite-disease",
        }.get(bottleneck, bottleneck)  # ← 将内部标识转换为可读标签

        if novelty > 0.7:
            return f"High novelty. Bottleneck: {bn_label} link. Worth investigating."
        elif novelty > 0.3:
            return f"Moderate novelty. {bn_label} link has limited evidence."
        else:
            return "Low novelty. This chain is well-studied."
