from dataclasses import dataclass, field, asdict
from typing import Dict, List


@dataclass
class HopEvidence:
    """Standard evidence payload returned by each hop agent."""

    hop_type: str
    pubmed_count: int

    db_hits: Dict[str, bool] = field(default_factory=dict)
    db_details: Dict[str, dict] = field(default_factory=dict)
    sources: List[str] = field(default_factory=list)
    query_used: str = ""

    def to_dict(self) -> dict:
        """Convert the evidence payload to a JSON-serializable dictionary."""
        return asdict(self)
