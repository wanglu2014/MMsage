import sys
from pathlib import Path
from typing import Dict, List, Optional
from concurrent.futures import ThreadPoolExecutor, as_completed

sys.path.insert(0, str(Path(__file__).parent.parent))

from .hop_evidence import HopEvidence
from .microbe_metabolite_agent import MicrobeMetaboliteAgent
from .metabolite_disease_agent import MetaboliteDiseaseAgent
from .microbe_disease_agent import MicrobeDiseaseAgent
from .evidence_aggregator import EvidenceAggregator


class MasterAgent:
    """
    Coordinate the three evidence agents and aggregate their results.

    Example:
        master = MasterAgent(c_max=500)
        result = master.run("Akkermansia_muciniphila", "butyrate", "IBD")
    """

    def __init__(self, synonym_map=None, c_max: int = 500, parallel: bool = True):
        self.mm_agent = MicrobeMetaboliteAgent(synonym_map)
        self.md_agent = MetaboliteDiseaseAgent(synonym_map)
        self.bd_agent = MicrobeDiseaseAgent(synonym_map)
        self.aggregator = EvidenceAggregator(c_max=c_max)
        self.parallel = parallel

    def run(self, bacteria: str, metabolite: str, disease: str) -> Dict:
        """Query and aggregate evidence for one candidate triplet."""

        if self.parallel:
            hop_results = self._run_parallel(bacteria, metabolite, disease)
        else:
            hop_results = self._run_sequential(bacteria, metabolite, disease)

        result = self.aggregator.aggregate(hop_results)
        result['bacteria'] = bacteria
        result['metabolite'] = metabolite
        result['disease'] = disease
        result['hop_evidence'] = [h.to_dict() for h in hop_results]
        return result

    def _run_parallel(self, bacteria: str, metabolite: str, disease: str) -> List[HopEvidence]:
        """Run the three evidence agents concurrently."""
        results = []
        with ThreadPoolExecutor(max_workers=3) as executor:
            futures = {
                executor.submit(self.mm_agent.run, bacteria, metabolite): "mm",
                executor.submit(self.md_agent.run, metabolite, disease): "md",
                executor.submit(self.bd_agent.run, bacteria, disease): "bd",
            }
            for future in as_completed(futures):
                try:
                    results.append(future.result())
                except Exception:
                    hop_type = {
                        "mm": "microbe_metabolite",
                        "md": "metabolite_disease",
                        "bd": "microbe_disease",
                    }[futures[future]]
                    results.append(HopEvidence(hop_type=hop_type, pubmed_count=0))
        return results

    def _run_sequential(self, bacteria: str, metabolite: str, disease: str) -> List[HopEvidence]:
        """Run the three evidence agents sequentially."""
        return [
            self.mm_agent.run(bacteria, metabolite),
            self.md_agent.run(metabolite, disease),
            self.bd_agent.run(bacteria, disease),
        ]

    def run_batch(self, candidates: List[Dict], disease: str = "IBD") -> List[Dict]:
        """Run evidence queries for a list of candidate pairs."""
        results = []
        for i, cand in enumerate(candidates):
            bacteria = cand.get('bacteria', '')
            metabolite = cand.get('metabolite', '')
            if not bacteria or not metabolite:
                continue
            result = self.run(bacteria, metabolite, disease)
            result['candidate'] = cand
            results.append(result)
        return results
