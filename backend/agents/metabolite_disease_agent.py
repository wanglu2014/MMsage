import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent))

from .hop_evidence import HopEvidence
from .pubmed_query import build_pubmed_query, query_pubmed_count
from db_checkers.ctd_checker import check_ctd
from db_checkers.chebi_checker import check_chebi


class MetaboliteDiseaseAgent:
    """Collect evidence for the metabolite-disease hop."""

    def __init__(self, synonym_map=None):
        self.synonym_map = synonym_map  # Deprecated compatibility input; query terms are not expanded.

    def run(self, metabolite: str, disease: str) -> HopEvidence:
        """Query evidence for a metabolite-disease pair."""

        query = build_pubmed_query([metabolite, disease], self.synonym_map)
        pubmed_count = query_pubmed_count(query)

        ctd_result = check_ctd(metabolite, disease)

        chebi_result = check_chebi(metabolite, disease)

        sources = []
        if pubmed_count > 0:
            sources.append(f"PubMed: {pubmed_count} articles for [{metabolite}] AND [{disease}]")
        if ctd_result['hit']:
            sources.append(f"CTD: {ctd_result['records']} records ({', '.join(ctd_result['edge_types'][:3])})")
        if chebi_result['hit']:
            sources.append(f"ChEBI: {chebi_result['records']} relationships for {metabolite}")

        return HopEvidence(
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
