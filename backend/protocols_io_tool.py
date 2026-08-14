"""
Read-only protocols.io retrieval for Step 3 experimental design.
"""

from __future__ import annotations

import json
import os
import re
from http.client import RemoteDisconnected
from typing import Any, Dict, List
from urllib.error import HTTPError, URLError
from urllib.parse import urlencode
from urllib.request import Request, urlopen


API_BASE = "https://www.protocols.io/api"
REQUEST_TIMEOUT = 20
MAX_SEARCH_RESULTS = 5
MAX_PROTOCOL_DETAILS = 3

VITRO_TERMS = (
    "in vitro",
    "cell",
    "culture",
    "organoid",
    "epithelial",
    "macrophage",
    "fermentation",
    "metabolite",
)
VIVO_TERMS = (
    "in vivo",
    "mouse",
    "mice",
    "murine",
    "animal",
    "colitis",
    "dss",
    "gavage",
)


def _clean_text(value: Any, limit: int = 1200) -> str:
    text = re.sub(r"<[^>]+>", " ", str(value or ""))
    text = re.sub(r"\s+", " ", text).strip()
    return text[:limit]


def _dedupe(items: List[str]) -> List[str]:
    seen = set()
    result = []
    for item in items:
        value = _clean_text(item)
        key = value.lower()
        if value and key not in seen:
            seen.add(key)
            result.append(value)
    return result


def _request_json(path: str, token: str, params: Dict[str, Any] | None = None) -> dict:
    url = f"{API_BASE}{path}"
    if params:
        url = f"{url}?{urlencode(params)}"
    request = Request(
        url,
        headers={
            "Authorization": f"Bearer {token}",
            "Accept": "application/json",
            "User-Agent": "MMSage-Step3/1.0",
        },
    )
    with urlopen(request, timeout=REQUEST_TIMEOUT) as response:
        return json.loads(response.read().decode("utf-8", errors="replace"))


def _flatten_component_text(value: Any, output: List[str]) -> None:
    if isinstance(value, dict):
        for key, nested in value.items():
            if key.lower() in {
                "title",
                "name",
                "description",
                "text",
                "content",
                "value",
                "step",
                "amount",
                "duration",
                "temperature",
            }:
                if isinstance(nested, (str, int, float)):
                    output.append(str(nested))
            _flatten_component_text(nested, output)
    elif isinstance(value, list):
        for nested in value:
            _flatten_component_text(nested, output)


def _step_summaries(payload: dict) -> List[str]:
    summaries = []
    for index, step in enumerate(payload.get("steps") or [], start=1):
        parts: List[str] = []
        _flatten_component_text(step, parts)
        text = _clean_text(" | ".join(_dedupe(parts)), 900)
        if text:
            summaries.append(f"Step {index}: {text}")
    return summaries[:20]


def _material_summaries(payload: dict) -> List[str]:
    materials = []
    for item in payload.get("materials") or []:
        if not isinstance(item, dict):
            continue
        parts = [
            item.get("name"),
            f"vendor: {(item.get('vendor') or {}).get('name')}" if isinstance(item.get("vendor"), dict) else "",
            f"SKU: {item.get('sku')}" if item.get("sku") else "",
            item.get("url"),
        ]
        text = _clean_text(" | ".join(str(part) for part in parts if part), 500)
        if text:
            materials.append(text)
    return _dedupe(materials)[:20]


def _protocol_url(item: dict) -> str:
    uri = str(item.get("uri") or "").strip()
    return f"https://www.protocols.io/view/{uri}" if uri else ""


def _citation(item: dict) -> str:
    doi = str(item.get("doi") or "").strip()
    if doi:
        doi = doi.replace("https://", "").replace("http://", "")
        return f"protocols.io DOI: {doi}"
    uri = str(item.get("uri") or "").strip()
    return f"protocols.io: {uri}" if uri else f"protocols.io ID: {item.get('id')}"


def _classify_protocol(text: str) -> List[str]:
    lowered = text.lower()
    types = []
    if any(term in lowered for term in VITRO_TERMS):
        types.append("in_vitro")
    if any(term in lowered for term in VIVO_TERMS):
        types.append("in_vivo")
    return types or ["unspecified"]


def _score_protocol(item: dict, query_terms: List[str]) -> int:
    haystack = " ".join(
        [
            str(item.get("title") or ""),
            str(item.get("description") or ""),
            str(item.get("before_start") or ""),
        ]
    ).lower()
    score = sum(2 for term in query_terms if term and term.lower() in haystack)
    if item.get("doi"):
        score += 1
    if item.get("peer_reviewed"):
        score += 2
    return score


def _is_relevant_protocol(item: dict, required_terms: List[str]) -> bool:
    haystack = " ".join(
        [
            str(item.get("title") or ""),
            str(item.get("description") or ""),
            str(item.get("before_start") or ""),
        ]
    ).lower()
    return any(term.lower() in haystack for term in required_terms if len(term.strip()) >= 3)


def _search_public_protocols(query: str, token: str) -> List[dict]:
    payload = _request_json(
        "/v3/protocols",
        token,
        {
            "filter": "public",
            "key": query,
            "order_field": "activity",
            "order_dir": "desc",
            "page_size": MAX_SEARCH_RESULTS,
            "page_id": 1,
        },
    )
    return [item for item in payload.get("items") or [] if isinstance(item, dict)]


def _fetch_protocol(item: dict, token: str) -> dict:
    identifier = item.get("id") or item.get("uri")
    detail_payload = _request_json(
        f"/v4/protocols/{identifier}",
        token,
        {"last_version": 1, "content_format": "markdown"},
    )
    detail = detail_payload.get("payload") or detail_payload.get("protocol") or detail_payload
    if not isinstance(detail, dict):
        detail = {}

    steps = _request_json(
        f"/v4/protocols/{identifier}/steps",
        token,
        {"last_version": 1, "content_format": "markdown"},
    )
    materials = _request_json(f"/v3/protocols/{identifier}/materials", token)
    merged = {**item, **detail}
    step_summaries = _step_summaries(steps)
    material_summaries = _material_summaries(materials)
    searchable = " ".join(
        [
            str(merged.get("title") or ""),
            str(merged.get("description") or ""),
            str(merged.get("before_start") or ""),
            " ".join(step_summaries),
            " ".join(material_summaries),
        ]
    )
    return {
        "id": merged.get("id") or identifier,
        "title": _clean_text(merged.get("title"), 300),
        "description": _clean_text(merged.get("description"), 700),
        "before_start": _clean_text(merged.get("before_start"), 700),
        "citation": _citation(merged),
        "doi": _clean_text(merged.get("doi"), 300),
        "uri": _clean_text(merged.get("uri"), 300),
        "url": _protocol_url(merged),
        "peer_reviewed": bool(merged.get("peer_reviewed")),
        "authors": _dedupe(
            [
                author.get("name")
                for author in merged.get("authors") or []
                if isinstance(author, dict) and author.get("name")
            ]
        ),
        "steps": step_summaries,
        "materials": material_summaries,
        "model_types": _classify_protocol(searchable),
    }


def collect_protocol_evidence(
    bacteria: str,
    metabolite: str,
    disease: str,
    research_question: str = "",
) -> dict:
    token = os.getenv("PROTOCOLS_IO_ACCESS_TOKEN", "").strip()
    result = {
        "enabled": bool(token),
        "status": "disabled",
        "message": "PROTOCOLS_IO_ACCESS_TOKEN is not configured.",
        "queries": [],
        "in_vitro": [],
        "in_vivo": [],
        "all": [],
    }
    if not token:
        return result

    bacteria_parts = bacteria.replace("_", " ").split()
    bacteria_short = " ".join(bacteria_parts[:2]) if bacteria_parts else bacteria
    disease_lower = disease.lower()
    queries = _dedupe(
        [
            bacteria_short,
            metabolite,
            f"{bacteria} culture {metabolite}",
            f"{metabolite} cell culture assay",
            f"{disease} mouse model",
            "DSS" if "ibd" in disease_lower else "",
            "DSS colitis" if "ibd" in disease_lower else f"{disease} animal model",
            "mouse colitis" if "ibd" in disease_lower else "",
            research_question,
        ]
    )[:8]
    result["queries"] = queries
    candidates: Dict[str, dict] = {}
    errors = []
    try:
        for query in queries:
            try:
                query_results = _search_public_protocols(query, token)
            except (HTTPError, URLError, TimeoutError, RemoteDisconnected, OSError, json.JSONDecodeError) as exc:
                errors.append(f"{query}: {exc}")
                continue
            for item in query_results:
                key = str(item.get("id") or item.get("uri") or "")
                if key:
                    candidates[key] = item

        query_terms = _dedupe(
            [
                bacteria,
                bacteria_parts[0] if bacteria_parts else "",
                metabolite,
                disease,
                "DSS",
                "colitis",
            ]
            + research_question.split()
        )[:20]
        candidates = {
            key: item
            for key, item in candidates.items()
            if _is_relevant_protocol(item, query_terms)
        }
        ranked = sorted(
            candidates.values(),
            key=lambda item: _score_protocol(item, query_terms),
            reverse=True,
        )
        for item in ranked[:MAX_PROTOCOL_DETAILS]:
            try:
                protocol = _fetch_protocol(item, token)
            except (HTTPError, URLError, TimeoutError, RemoteDisconnected, OSError, json.JSONDecodeError) as exc:
                errors.append(f"protocol {item.get('id') or item.get('uri')}: {exc}")
                continue
            result["all"].append(protocol)
            for model_type in protocol["model_types"]:
                if model_type in ("in_vitro", "in_vivo"):
                    result[model_type].append(protocol)

        result["status"] = "completed"
        result["message"] = (
            f"Retrieved {len(result['all'])} protocols.io operational references"
            f" with {len(errors)} recoverable request errors."
        )
    except (HTTPError, URLError, TimeoutError, RemoteDisconnected, OSError, json.JSONDecodeError, ValueError) as exc:
        result["status"] = "error"
        result["message"] = f"protocols.io retrieval failed: {exc}"
    return result
