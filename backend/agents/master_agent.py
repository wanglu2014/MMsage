# ═══════════════════════════════════════════════════════════════════
# 【模块：Master Coordinator Agent（主协调器）】
# ═══════════════════════════════════════════════════════════════════
#
# 【承上启下】作为多Agent协作系统的入口和调度中心，本模块接收一个候选三元组
#             (bacteria, metabolite, disease)，并行调度三个子Agent分别查询
#             各自负责的跳（hop），收集返回的 HopEvidence 后交给
#             EvidenceAggregator 计算综合 chain_novelty。
#             它是 run_pipeline.py 中 step2b 调用的核心对象。
#
# 【核心思想】通过并行调度实现高效的多源证据收集，通过统一聚合实现
#             可解释的新颖性评分。
#
# 【输入】 bacteria (微生物名), metabolite (代谢物名), disease (疾病名)
# 【输出】 dict (聚合结果，含 chain_novelty, hop_counts, bottleneck,
#               hop_evidence, sources, recommendation 等)
#
# ───────────────────────────────────────────────────────────────────

import sys
from pathlib import Path
from typing import Dict, List, Optional
from concurrent.futures import ThreadPoolExecutor, as_completed  # ← 导入线程池，用于三个Agent的并行执行

sys.path.insert(0, str(Path(__file__).parent.parent))  # ← 将 backend/ 加入搜索路径

from .hop_evidence import HopEvidence  # ← 导入统一的证据数据结构
from .microbe_metabolite_agent import MicrobeMetaboliteAgent  # ← 导入第一跳Agent
from .metabolite_disease_agent import MetaboliteDiseaseAgent  # ← 导入第二跳Agent
from .microbe_disease_agent import MicrobeDiseaseAgent  # ← 导入辅助验证Agent
from .evidence_aggregator import EvidenceAggregator  # ← 导入证据聚合器


class MasterAgent:
    """
    主协调器：调度三个子Agent并聚合结果。

    【使用方式】
        master = MasterAgent(c_max=500)
        result = master.run("Akkermansia_muciniphila", "butyrate", "IBD")
    """

    def __init__(self, synonym_map=None, c_max: int = 500, parallel: bool = True):
        self.mm_agent = MicrobeMetaboliteAgent(synonym_map)  # ← 初始化第一跳Agent（微生物↔代谢物）
        self.md_agent = MetaboliteDiseaseAgent(synonym_map)   # ← 初始化第二跳Agent（代谢物↔疾病）
        self.bd_agent = MicrobeDiseaseAgent(synonym_map)       # ← 初始化辅助Agent（微生物↔疾病）
        self.aggregator = EvidenceAggregator(c_max=c_max)      # ← 初始化证据聚合器，设定归一化上限
        self.parallel = parallel  # ← 是否启用并行执行（默认 True，使用线程池）

    def run(self, bacteria: str, metabolite: str, disease: str) -> Dict:
        """
        对单个候选三元组执行完整的多Agent证据查询和聚合。

        【输入】 bacteria, metabolite, disease (候选三元组)
        【输出】 dict (完整聚合结果 + 原始 hop_evidence)
        """

        # ═══════════════════════════════════════════════════════════════
        # 【子模块 1：调度三个子Agent】
        # ═══════════════════════════════════════════════════════════════
        #
        # 【核心思想】根据 parallel 标志选择并行或串行执行。并行模式下
        #             三个Agent同时查询 PubMed 和数据库，显著减少总耗时。
        #
        if self.parallel:
            hop_results = self._run_parallel(bacteria, metabolite, disease)  # ← 并行调度三个Agent
        else:
            hop_results = self._run_sequential(bacteria, metabolite, disease)  # ← 串行调度三个Agent

        # ═══════════════════════════════════════════════════════════════
        # 【子模块 2：聚合证据并附加元数据】
        # ═══════════════════════════════════════════════════════════════
        #
        # 【核心思想】调用 EvidenceAggregator 计算 chain_novelty，
        #             然后附加候选三元组信息和原始 hop_evidence 供前端展示。
        #
        result = self.aggregator.aggregate(hop_results)  # ← 以三个 HopEvidence 为输入，执行 min+bonus 聚合，生成 chain_novelty 等
        result['bacteria'] = bacteria      # ← 附加微生物名到结果
        result['metabolite'] = metabolite  # ← 附加代谢物名到结果
        result['disease'] = disease        # ← 附加疾病名到结果
        result['hop_evidence'] = [h.to_dict() for h in hop_results]  # ← 将三个 HopEvidence 序列化为 dict 列表，保留完整的原始证据
        return result

    def _run_parallel(self, bacteria: str, metabolite: str, disease: str) -> List[HopEvidence]:
        """
        并行执行三个子Agent。

        【核心思想】使用 ThreadPoolExecutor(max_workers=3) 同时提交三个任务，
                    通过 as_completed 收集结果。任何Agent失败时返回空证据，
                    不影响其他Agent的结果。
        """
        results = []
        with ThreadPoolExecutor(max_workers=3) as executor:  # ← 创建3线程的线程池
            futures = {
                executor.submit(self.mm_agent.run, bacteria, metabolite): "mm",  # ← 提交第一跳任务
                executor.submit(self.md_agent.run, metabolite, disease): "md",    # ← 提交第二跳任务
                executor.submit(self.bd_agent.run, bacteria, disease): "bd",      # ← 提交辅助跳任务
            }
            for future in as_completed(futures):  # ← 按完成顺序收集结果
                try:
                    results.append(future.result())  # ← 获取Agent返回的 HopEvidence
                except Exception:
                    hop_type = {  # ← Agent执行失败时，构造空的 HopEvidence 作为降级处理
                        "mm": "microbe_metabolite",
                        "md": "metabolite_disease",
                        "bd": "microbe_disease",
                    }[futures[future]]
                    results.append(HopEvidence(hop_type=hop_type, pubmed_count=0))
        return results

    def _run_sequential(self, bacteria: str, metabolite: str, disease: str) -> List[HopEvidence]:
        """串行执行三个子Agent（用于调试或无并行需求的场景）。"""
        return [
            self.mm_agent.run(bacteria, metabolite),  # ← 串行执行第一跳
            self.md_agent.run(metabolite, disease),    # ← 串行执行第二跳
            self.bd_agent.run(bacteria, disease),      # ← 串行执行辅助跳
        ]

    def run_batch(self, candidates: List[Dict], disease: str = "IBD") -> List[Dict]:
        """
        批量执行：对一组候选对逐一运行多Agent查询。

        【输入】 candidates (含 bacteria/metabolite 键的字典列表), disease (疾病上下文)
        【输出】 List[Dict] (每个候选对的聚合结果列表)
        """
        results = []
        for i, cand in enumerate(candidates):
            bacteria = cand.get('bacteria', '')
            metabolite = cand.get('metabolite', '')
            if not bacteria or not metabolite:
                continue  # ← 跳过缺失关键字段的候选对
            result = self.run(bacteria, metabolite, disease)  # ← 对单个候选对执行完整查询
            result['candidate'] = cand  # ← 附加原始候选对信息
            results.append(result)
        return results
