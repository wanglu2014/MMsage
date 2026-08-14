import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent))

from .hop_evidence import HopEvidence
from .pubmed_query import build_pubmed_query, query_pubmed_count
from db_checkers.disbiome_checker import check_disbiome
from db_checkers.bacdive_checker import check_bacdive
from db_checkers.binded_checker import check_binded_microbe_metabolite


class MicrobeMetaboliteAgent:
    """Collect evidence for the microbe-metabolite hop."""

    def __init__(self, synonym_map=None):
        self.synonym_map = synonym_map  # Deprecated compatibility input; query terms are not expanded.

    def run(self, bacteria: str, metabolite: str) -> HopEvidence:
        """Query evidence for a microbe-metabolite pair."""

        query = build_pubmed_query([bacteria, metabolite], self.synonym_map)
        pubmed_count = query_pubmed_count(query)

        disbiome_result = check_disbiome(bacteria, metabolite)

        bacdive_result = check_bacdive(bacteria, metabolite)

        binded_result = check_binded_microbe_metabolite(bacteria, metabolite)

        sources = []
        if pubmed_count > 0:
            sources.append(f"PubMed: {pubmed_count} articles for [{bacteria}] AND [{metabolite}]")
        if disbiome_result['hit']:
            sources.append(f"Disbiome: {disbiome_result['records']} records for {bacteria}")
        if bacdive_result['hit']:
            sources.append(f"BacDive: {bacdive_result['records']} records for {bacteria}")
        if binded_result['hit']:
            sources.append(f"BindedDB: {binded_result['records']} records ({','.join(binded_result.get('sources',[])[:3])})")

        return HopEvidence(
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
