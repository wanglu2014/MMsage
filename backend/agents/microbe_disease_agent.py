import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent))

from .hop_evidence import HopEvidence
from .pubmed_query import build_pubmed_query, query_pubmed_count
from db_checkers.disbiome_checker import check_disbiome_disease


class MicrobeDiseaseAgent:
    """Collect evidence for the microbe-disease validation hop."""

    def __init__(self, synonym_map=None):
        self.synonym_map = synonym_map  # Deprecated compatibility input; query terms are not expanded.

    def run(self, bacteria: str, disease: str) -> HopEvidence:
        """Query evidence for a microbe-disease pair."""

        query = build_pubmed_query([bacteria, disease], self.synonym_map)
        pubmed_count = query_pubmed_count(query)

        disbiome_result = check_disbiome_disease(bacteria, disease)

        sources = []
        if pubmed_count > 0:
            sources.append(f"PubMed: {pubmed_count} articles for [{bacteria}] AND [{disease}]")
        if disbiome_result['hit']:
            sources.append(f"Disbiome: {disbiome_result['records']} records ({', '.join(disbiome_result['edge_types'][:3])})")

        return HopEvidence(
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
