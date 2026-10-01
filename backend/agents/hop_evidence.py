# ═══════════════════════════════════════════════════════════════════
# 【模块：HopEvidence 数据结构定义】
# ═══════════════════════════════════════════════════════════════════
#
# 【承上启下】作为多Agent协作系统的基础数据契约，本模块定义了所有子Agent
#             的统一返回格式 HopEvidence，确保三个子Agent（微生物-代谢物、
#             代谢物-疾病、微生物-疾病）的输出可以被 EvidenceAggregator
#             无差别地解析和聚合。
#
# 【核心思想】通过 dataclass 定义标准化的证据容器，实现Agent间的松耦合通信。
#
# 【输入】 各子Agent的查询结果（PubMed计数、数据库命中等）
# 【输出】 HopEvidence 实例，可序列化为 dict 供 JSON 输出
#
# ───────────────────────────────────────────────────────────────────

from dataclasses import dataclass, field, asdict  # ← 导入 dataclass 装饰器和序列化工具，获得声明式数据类定义能力
from typing import Dict, List  # ← 导入类型注解，提升代码可读性和IDE支持


@dataclass
class HopEvidence:
    """每个子Agent的标准返回格式——一跳（hop）的证据包。"""

    hop_type: str           # ← 跳类型标识："microbe_metabolite" / "metabolite_disease" / "microbe_disease"
    pubmed_count: int       # ← PubMed 两两共现文章数，是 novelty 计算的核心输入

    db_hits: Dict[str, bool] = field(default_factory=dict)    # ← 各数据库命中情况，如 {"disbiome": True, "ctd": False}
    db_details: Dict[str, dict] = field(default_factory=dict)  # ← 各数据库的完整查询结果，供前端展示详情
    sources: List[str] = field(default_factory=list)           # ← 人类可读的证据来源描述列表，用于报告生成
    query_used: str = ""    # ← 实际使用的 PubMed 查询字符串，便于结果复现和调试

    def to_dict(self) -> dict:
        """以 self 为输入，调用 asdict 递归转换，生成纯 dict 结构供 JSON 序列化。"""
        return asdict(self)
