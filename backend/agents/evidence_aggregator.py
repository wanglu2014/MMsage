import math
from typing import List, Dict
from .hop_evidence import HopEvidence


class EvidenceAggregator:
    """Combine hop evidence into a chain-novelty score."""

    def __init__(self, c_max: int = 500):
        self.c_max = c_max

    def aggregate(self, hop_results: List[HopEvidence]) -> Dict:
        """Aggregate three hop results and calculate chain novelty."""

        counts = {h.hop_type: h.pubmed_count for h in hop_results}
        mm_count = counts.get("microbe_metabolite", 0)
        md_count = counts.get("metabolite_disease", 0)
        bd_count = counts.get("microbe_disease", 0)

        chain_count = min(mm_count, md_count)

        db_bonus = 0
        for h in hop_results:
            for db_name, hit in h.db_hits.items():
                if hit:
                    db_bonus += 10
        db_bonus = min(db_bonus, 50)

        total_count = chain_count + db_bonus
        if self.c_max <= 0:
            chain_novelty = 1.0
        else:
            chain_novelty = 1.0 - math.log(1 + total_count) / math.log(1 + self.c_max)
            chain_novelty = max(0.0, min(1.0, chain_novelty))

        if mm_count <= md_count:
            bottleneck = "microbe_metabolite"
        else:
            bottleneck = "metabolite_disease"

        all_sources = []
        all_db_details = {}
        for h in hop_results:
            all_sources.extend(h.sources)
            all_db_details.update(h.db_details)

        return {
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
        """Generate a recommendation from novelty and bottleneck values."""
        bn_label = {
            "microbe_metabolite": "microbe-metabolite",
            "metabolite_disease": "metabolite-disease",
        }.get(bottleneck, bottleneck)

        if novelty > 0.7:
            return f"High novelty. Bottleneck: {bn_label} link. Worth investigating."
        elif novelty > 0.3:
            return f"Moderate novelty. {bn_label} link has limited evidence."
        else:
            return "Low novelty. This chain is well-studied."
