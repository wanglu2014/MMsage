"""
Step 3 protocol refinement layer.

This module keeps storage file-based and does not modify database/catalog logic.
"""

from __future__ import annotations

import json
import re
import time
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, List, Optional
from uuid import uuid4

try:
    import requests
except ImportError:  # pragma: no cover - runtime guard
    requests = None

from agents.pubmed_query import build_pubmed_query
from llm_reasoning import OPENAI_API_BASE, OPENAI_API_KEY, OPENAI_MODEL

BACKEND_DIR = Path(__file__).parent
PROJECT_DIR = BACKEND_DIR.parent
OUTPUTS_DIR = PROJECT_DIR / "outputs"
DATA_DIR = PROJECT_DIR / "data"
SESSIONS_DIR = DATA_DIR / "protocol_refiner_sessions"

MODULE_ORDER = [
    "Hypothesis",
    "Experimental Design",
    "Mechanism Validation",
    "Controls & Readouts",
    "Feasibility & Timeline",
]

MODULE_SPECS = {
    "Hypothesis": {
        "module_id": "hypothesis",
        "goal": "Define the core bacteria-metabolite-disease hypothesis and the main causal chain to test.",
    },
    "Experimental Design": {
        "module_id": "experimental_design",
        "goal": "Define the overall validation strategy across in vitro, in vivo, and clinical directions as applicable.",
    },
    "Mechanism Validation": {
        "module_id": "mechanism_validation",
        "goal": "Define the key pathway, receptor, enzyme, or downstream mechanism experiments.",
    },
    "Controls & Readouts": {
        "module_id": "controls_readouts",
        "goal": "Define controls, endpoints, readouts, and interpretation rules.",
    },
    "Feasibility & Timeline": {
        "module_id": "feasibility_timeline",
        "goal": "Define feasibility, stage-gated execution, resource needs, and timeline.",
    },
}

REVIEW_AGENTS = {
    "design": "Evaluate experimental logic, controls, readouts, and causal validity. Do not discuss budget or literature volume unless they directly affect design validity.",
    "evidence": "Evaluate whether the module's claims are supported by the current candidate evidence and literature. Do not focus on budget or operations.",
    "feasibility": "Evaluate time, budget, operational complexity, sample access, and alignment with user requirements. Do not optimize for ideal scientific completeness at any cost.",
}


def _now_iso() -> str:
    return datetime.now(timezone.utc).isoformat()


def _slug(text: str, max_len: int = 30) -> str:
    clean = re.sub(r"[^\w]+", "_", (text or "session").strip()).strip("_").lower()
    return clean[:max_len] or "session"


def _save_json(path: Path, data: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8") as f:
        json.dump(data, f, indent=2, ensure_ascii=False)


def _load_json(path: Path) -> Any:
    if not path.exists():
        return None
    with open(path, "r", encoding="utf-8") as f:
        return json.load(f)


def _session_dir(session_id: str) -> Path:
    return SESSIONS_DIR / session_id


def _read_output_json(filename: str) -> list:
    path = OUTPUTS_DIR / filename
    if not path.exists():
        return []
    with open(path, "r", encoding="utf-8") as f:
        return json.load(f)


def _load_pair_from_outputs(bacteria: str, metabolite: str) -> Optional[dict]:
    quadrants = _read_output_json("step3_quadrant.json")
    target = None
    for row in quadrants:
        if row.get("bacteria") == bacteria and row.get("metabolite") == metabolite:
            target = dict(row)
            break
    if not target:
        return None

    step2 = _read_output_json("step2_chain_novelty.json")
    for row in step2:
        if row.get("bacteria") == bacteria and row.get("metabolite") == metabolite:
            for key in (
                "pairwise_counts",
                "edge_cooccurrences",
                "chain_path_str",
                "chain_path",
                "bottleneck_edge",
                "agent_evidence",
                "agent_recommendation",
                "chain_query",
                "query_terms",
                "has_path",
                "chain_count",
                "chain_novelty",
                "all_paths_info",
            ):
                if row.get(key) is not None:
                    target[key] = row.get(key)
            break

    if not target.get("agent_evidence"):
        step2b = _read_output_json("step2b_agent_evidence.json")
        for row in step2b:
            if row.get("bacteria") == bacteria and row.get("metabolite") == metabolite:
                target["agent_evidence"] = row
                break
    return target


def load_candidate_record(run_id: Optional[str], bacteria: str, metabolite: str) -> Optional[dict]:
    if run_id:
        try:
            from export_catalog import get_record

            rec = get_record(run_id, bacteria, metabolite)
            if rec:
                return rec
        except Exception:
            pass
    return _load_pair_from_outputs(bacteria, metabolite)


def _candidate_metrics(rec: Optional[dict], disease: str) -> Dict[str, Any]:
    rec = rec or {}
    ae = rec.get("agent_evidence") or {}
    step3 = rec.get("step3_detail") or {}
    return {
        "bacteria": rec.get("bacteria", ""),
        "metabolite": rec.get("metabolite", ""),
        "disease": rec.get("disease") or disease,
        "mmsage_norm": rec.get("mmsage_norm", 0),
        "pair_bm_exp": rec.get("pair_bm_exp", 0),
        "pair_md_exp": rec.get("pair_md_exp", 0),
        "chain_count": rec.get("chain_count", 0),
        "chain_novelty": ae.get("chain_novelty", rec.get("chain_novelty", 1.0)),
        "quadrant": rec.get("quadrant", ""),
        "quadrant_label": rec.get("quadrant_label", ""),
        "chain_path_str": rec.get("chain_path_str", ""),
        "bottleneck_edge": rec.get("bottleneck_edge"),
        "pairwise_counts": rec.get("pairwise_counts") or {},
        "agent_recommendation": ae.get("recommendation") or rec.get("agent_recommendation"),
        "hop_evidence": step3.get("hop_evidence") or [],
    }


def _format_candidate_context(metrics: Dict[str, Any], mechanism_summary: str = "") -> str:
    lines = [
        f"Candidate: {metrics.get('bacteria', '').replace('_', ' ')} x {metrics.get('metabolite', '')}",
        f"Disease: {metrics.get('disease', '')}",
        f"MMSage signal: {metrics.get('mmsage_norm', 0):.3f}",
        f"Microbe-Metabolite papers (bm_exp): {metrics.get('pair_bm_exp', 0)}",
        f"Metabolite-Disease papers (md_exp): {metrics.get('pair_md_exp', 0)}",
        f"Chain count: {metrics.get('chain_count', 0)}",
        f"Chain novelty: {metrics.get('chain_novelty', 1.0):.3f}",
        f"Quadrant: {metrics.get('quadrant_label') or metrics.get('quadrant') or 'Unknown'}",
    ]
    if metrics.get("chain_path_str"):
        lines.append(f"Chain path: {metrics['chain_path_str']}")
    if metrics.get("bottleneck_edge"):
        edge = metrics["bottleneck_edge"]
        lines.append(
            f"Bottleneck edge: {edge.get('source', '?')} -> {edge.get('target', '?')} ({edge.get('cooccurrence', 0)})"
        )
    if metrics.get("agent_recommendation"):
        lines.append(f"Agent recommendation: {metrics['agent_recommendation']}")
    if mechanism_summary:
        lines.append(f"Additional user mechanism summary: {mechanism_summary}")
    return "\n".join(lines)


def _fetch_supporting_literature(bacteria: str, metabolite: str, disease: str) -> List[dict]:
    try:
        from build_kg import fetch_abstracts, search_pubmed

        bm_query = build_pubmed_query([bacteria, metabolite])
        md_query = build_pubmed_query([metabolite, disease])
        bm_articles = fetch_abstracts(search_pubmed(bm_query, max_results=3))
        md_articles = fetch_abstracts(search_pubmed(md_query, max_results=3))
        articles = []
        for article in bm_articles:
            a = dict(article)
            a["topic"] = "microbe_metabolite"
            articles.append(a)
        for article in md_articles:
            a = dict(article)
            a["topic"] = "metabolite_disease"
            articles.append(a)
        return articles
    except Exception:
        return []


def _format_literature_context(articles: List[dict]) -> str:
    if not articles:
        return "No additional literature abstracts were retrieved."
    chunks = []
    for article in articles[:6]:
        chunks.append(
            f"[{article.get('topic', 'literature')}] PMID {article.get('pmid', '')}: "
            f"{article.get('title', '')} ({article.get('journal', '')}, {article.get('year', '')})\n"
            f"{(article.get('abstract', '') or '')[:700]}"
        )
    return "\n\n".join(chunks)


def _prepare_protocol_support(
    bacteria: str,
    metabolite: str,
    disease: str = "IBD",
    run_id: Optional[str] = None,
    mechanism_summary: str = "",
) -> dict:
    record = load_candidate_record(run_id, bacteria, metabolite)
    metrics = _candidate_metrics(record, disease)
    literature = _fetch_supporting_literature(bacteria, metabolite, disease)
    candidate_context = _format_candidate_context(metrics, mechanism_summary)
    literature_context = _format_literature_context(literature)
    return {
        "candidate_metrics": metrics,
        "literature": literature,
        "candidate_context": candidate_context,
        "literature_context": literature_context,
    }


def _chat_json(
    system_prompt: str,
    user_prompt: str,
    max_tokens: int = 2000,
    temperature: float = 0.3,
    request_timeout: float = 240,
) -> dict:
    if not OPENAI_API_KEY or requests is None:
        raise RuntimeError("LLM is not configured for protocol refinement.")

    json_mode_instruction = (
        "\n\nOutput exactly one valid json object with no Markdown wrapper or surrounding text."
    )
    token_budgets = [max_tokens, min(max(max_tokens * 2, max_tokens + 2000), 8000)]
    last_error: Optional[Exception] = None
    for attempt, token_budget in enumerate(token_budgets, start=1):
        retry_instruction = ""
        if attempt > 1:
            retry_instruction = (
                "\n\nIMPORTANT RETRY: The previous response was truncated or invalid JSON. "
                "Return one complete, compact JSON object only. Keep every string concise, "
                "avoid repeated evidence, and close all arrays and objects."
            )
        response = None
        for network_attempt in range(1, 3):
            try:
                response = requests.post(
                    f"{OPENAI_API_BASE}/chat/completions",
                    headers={
                        "Content-Type": "application/json",
                        "Authorization": f"Bearer {OPENAI_API_KEY}",
                    },
                    json={
                        "model": OPENAI_MODEL,
                        "messages": [
                            {"role": "system", "content": system_prompt + retry_instruction},
                            {"role": "user", "content": user_prompt + json_mode_instruction},
                        ],
                        "temperature": temperature if attempt == 1 else min(temperature, 0.1),
                        "max_tokens": token_budget,
                        "response_format": {"type": "json_object"},
                    },
                    timeout=request_timeout,
                )
                break
            except requests.RequestException as exc:
                last_error = exc
                if network_attempt == 2:
                    raise RuntimeError(
                        f"LLM network request failed after {network_attempt} attempts: {exc}"
                    ) from exc
                time.sleep(2)
        if response is None:
            raise RuntimeError(f"LLM network request failed: {last_error}")
        if response.status_code != 200:
            raise RuntimeError(f"LLM API error: {response.status_code} - {response.text[:1000]}")

        payload = response.json()
        choices = payload.get("choices") or []
        if not choices:
            raise RuntimeError("LLM API returned no choices.")
        choice = choices[0]
        content = str((choice.get("message") or {}).get("content") or "").strip()
        finish_reason = str(choice.get("finish_reason") or "")
        try:
            parsed = json.loads(content)
            if not isinstance(parsed, dict):
                raise ValueError("LLM JSON response is not an object.")
            return parsed
        except (json.JSONDecodeError, ValueError) as exc:
            last_error = exc
            if attempt == len(token_budgets):
                raise RuntimeError(
                    f"LLM returned invalid JSON after {attempt} attempts "
                    f"(finish_reason={finish_reason or 'unknown'}, chars={len(content)}): {exc}"
                ) from exc

    raise RuntimeError(f"LLM JSON generation failed: {last_error}")


def generate_protocol_draft(
    bacteria: str,
    metabolite: str,
    disease: str = "IBD",
    run_id: Optional[str] = None,
    mechanism_summary: str = "",
) -> dict:
    support = _prepare_protocol_support(
        bacteria=bacteria,
        metabolite=metabolite,
        disease=disease,
        run_id=run_id,
        mechanism_summary=mechanism_summary,
    )
    metrics = support["candidate_metrics"]
    literature = support["literature"]
    candidate_context = support["candidate_context"]
    literature_context = support["literature_context"]

    system_prompt = """You are an expert microbiome-metabolomics experimental strategist.

Generate an initial validation protocol draft for a bacteria-metabolite-disease candidate.
The draft should be concrete, scientifically cautious, and actionable.

Return JSON with:
- title
- executive_summary
- global_hypothesis
- protocol_full_text
- major_risks (array of strings)

Requirements:
- Cover hypothesis, experimental design, mechanism validation, controls/readouts, and feasibility/timeline.
- Use evidence strength to calibrate claims.
- Cite PMIDs inline when literature supports a claim.
- Evaluate each PMID independently. Direct microbial production requires a same-paper controlled culture experiment with direct metabolite measurement; never assemble direct evidence from separate indirect papers.
- Preserve competing direct-production, indirect-production, and parallel-effect hypotheses until branch-specific experiments distinguish them.
- Order the study as culture/cell gates, conditional animal experiments, and human validation only after an interpretable animal result.
- Do not output markdown fences."""

    user_prompt = f"""Candidate context:
{candidate_context}

Literature context:
{literature_context}

Generate the initial protocol draft."""

    result = _chat_json(system_prompt, user_prompt, max_tokens=2200, temperature=0.5)
    protocol = result.get("protocol_full_text", "").strip()
    if not protocol:
        raise RuntimeError("LLM returned an empty protocol draft.")
    return {
        "protocol": protocol,
        "generated_by": "ai",
        "model": OPENAI_MODEL,
        "title": result.get("title", f"{bacteria} x {metabolite} protocol"),
        "executive_summary": result.get("executive_summary", ""),
        "global_hypothesis": result.get("global_hypothesis", ""),
        "major_risks": result.get("major_risks", []),
        "candidate_metrics": metrics,
        "literature": articles_to_export(literature),
        "candidate_context": candidate_context,
        "literature_context": literature_context,
    }


def build_protocol_payload_from_text(
    protocol_text: str,
    bacteria: str,
    metabolite: str,
    disease: str = "IBD",
    run_id: Optional[str] = None,
    mechanism_summary: str = "",
) -> dict:
    support = _prepare_protocol_support(
        bacteria=bacteria,
        metabolite=metabolite,
        disease=disease,
        run_id=run_id,
        mechanism_summary=mechanism_summary,
    )
    clean_protocol = (protocol_text or "").strip()
    if not clean_protocol:
        raise RuntimeError("Protocol text is empty.")
    return {
        "protocol": clean_protocol,
        "generated_by": "dashboard_existing_protocol",
        "model": OPENAI_MODEL,
        "title": f"{bacteria} x {metabolite} protocol",
        "executive_summary": "",
        "global_hypothesis": "",
        "major_risks": [],
        "candidate_metrics": support["candidate_metrics"],
        "literature": articles_to_export(support["literature"]),
        "candidate_context": support["candidate_context"],
        "literature_context": support["literature_context"],
    }


def articles_to_export(articles: List[dict]) -> List[dict]:
    exported = []
    for article in articles:
        exported.append(
            {
                "topic": article.get("topic"),
                "pmid": article.get("pmid"),
                "title": article.get("title"),
                "journal": article.get("journal"),
                "year": article.get("year"),
                "abstract": article.get("abstract"),
            }
        )
    return exported


def _fallback_decompose(protocol_text: str) -> List[dict]:
    return [
        {
            "module_id": MODULE_SPECS[name]["module_id"],
            "module_name": name,
            "goal": MODULE_SPECS[name]["goal"],
            "summary": f"Fallback module for {name}.",
            "content": protocol_text,
            "enabled": True,
            "priority": "high" if idx < 3 else "medium",
            "user_notes": "",
        }
        for idx, name in enumerate(MODULE_ORDER)
    ]


def _normalize_modules(modules: List[dict], protocol_text: str) -> List[dict]:
    by_name = {m.get("module_name"): m for m in modules if m.get("module_name") in MODULE_SPECS}
    normalized = []
    for idx, name in enumerate(MODULE_ORDER):
        src = by_name.get(name, {})
        normalized.append(
            {
                "module_id": src.get("module_id") or MODULE_SPECS[name]["module_id"],
                "module_name": name,
                "goal": src.get("goal") or MODULE_SPECS[name]["goal"],
                "summary": src.get("summary") or f"{name} section derived from the protocol draft.",
                "content": src.get("content") or protocol_text,
                "enabled": src.get("enabled", True),
                "priority": src.get("priority") or ("high" if idx < 3 else "medium"),
                "user_notes": src.get("user_notes", ""),
            }
        )
    return normalized


def decompose_protocol(protocol_text: str, candidate_context: str) -> List[dict]:
    system_prompt = """You decompose a validation protocol into exactly five fixed modules.

Allowed module names only:
1. Hypothesis
2. Experimental Design
3. Mechanism Validation
4. Controls & Readouts
5. Feasibility & Timeline

Return JSON with a top-level key "modules".
Each module must include:
- module_id
- module_name
- goal
- summary
- content
- enabled
- priority
- user_notes

Do not invent extra module types."""

    user_prompt = f"""Candidate context:
{candidate_context}

Protocol draft:
{protocol_text}

Decompose the protocol into the five fixed modules."""
    try:
        result = _chat_json(system_prompt, user_prompt, max_tokens=2200, temperature=0.2)
        return _normalize_modules(result.get("modules") or [], protocol_text)
    except Exception:
        return _fallback_decompose(protocol_text)


def _build_review_prompt(
    agent_type: str,
    module: dict,
    candidate_context: str,
    literature_context: str,
    extra_requirements: str,
) -> tuple[str, str]:
    system_prompt = f"""You are the {agent_type.title()} Agent in a Step 3 protocol refinement workflow.

Role:
{REVIEW_AGENTS[agent_type]}

Return JSON with exactly these fields:
- agent_type
- module_id
- overall_judgment
- strengths
- major_issues
- minor_issues
- blocking_issues
- recommended_revisions
- optional_improvements
- confidence
- confidence_reason

Allowed overall_judgment values:
- strong
- acceptable
- weak
- blocked

All issue and revision fields must be arrays of short, specific strings."""

    user_prompt = f"""Candidate context:
{candidate_context}

Literature context:
{literature_context}

User extra requirements:
{extra_requirements or 'None'}

Module under review:
Name: {module.get('module_name')}
Goal: {module.get('goal')}
Summary: {module.get('summary')}
Priority: {module.get('priority', 'medium')}
User notes: {module.get('user_notes') or 'None'}
Content:
{module.get('content')}

Review this module from your assigned perspective only."""
    return system_prompt, user_prompt


def review_module(agent_type: str, module: dict, candidate_context: str, literature_context: str, extra_requirements: str) -> dict:
    system_prompt, user_prompt = _build_review_prompt(
        agent_type, module, candidate_context, literature_context, extra_requirements
    )
    result = _chat_json(system_prompt, user_prompt, max_tokens=1600, temperature=0.2)
    result["agent_type"] = agent_type
    result["module_id"] = module["module_id"]
    return result


def synthesize_module(
    module: dict,
    reviews: List[dict],
    candidate_context: str,
    extra_requirements: str,
) -> dict:
    system_prompt = """You are the Judge/Synthesizer Agent in a Step 3 protocol refinement workflow.

Your job is to merge review feedback, resolve conflicts conservatively, and produce a revised module.

Rules:
- Prioritize blocking issues.
- Respect explicit user requirements.
- Keep the same module scope.
- If evidence is weak, downgrade claims rather than deleting the whole module when possible.

Return JSON with these fields:
- module_id
- overall_decision
- agreed_issues
- conflicting_issues
- blocking_issues
- applied_revisions
- deferred_revisions
- rejected_revisions
- revised_content
- change_log
- remaining_risks

Allowed overall_decision values:
- accept_with_minor_edits
- revise
- blocked"""

    user_prompt = f"""Candidate context:
{candidate_context}

User extra requirements:
{extra_requirements or 'None'}

Original module:
Name: {module.get('module_name')}
Goal: {module.get('goal')}
Priority: {module.get('priority', 'medium')}
User notes: {module.get('user_notes') or 'None'}
Content:
{module.get('content')}

Reviewer outputs:
{json.dumps(reviews, ensure_ascii=False, indent=2)}

Produce the revised module."""

    result = _chat_json(system_prompt, user_prompt, max_tokens=2000, temperature=0.2)
    result["module_id"] = module["module_id"]
    return result


def _needs_second_round(reviews: List[dict], synthesis: dict) -> bool:
    if synthesis.get("blocking_issues"):
        return True
    for review in reviews:
        if review.get("overall_judgment") in {"weak", "blocked"}:
            return True
    return False


def _build_final_protocol(modules: List[dict], candidate_context: str, extra_requirements: str, debate_points: List[str]) -> dict:
    joined_modules = "\n\n".join(
        [f"## {m['module_name']}\n{m['content']}" for m in modules if m.get("enabled", True)]
    )
    try:
        system_prompt = """You are assembling a final enhanced experimental protocol after multi-agent refinement.

Return JSON with:
- protocol_full_text
- change_summary
- remaining_risks

The final text should be polished but still concrete and experimentally oriented."""
        user_prompt = f"""Candidate context:
{candidate_context}

User extra requirements:
{extra_requirements or 'None'}

Refined modules:
{joined_modules}

Key debate points:
{json.dumps(debate_points, ensure_ascii=False)}

Assemble the final enhanced protocol."""
        result = _chat_json(system_prompt, user_prompt, max_tokens=2200, temperature=0.3)
        if result.get("protocol_full_text"):
            return result
    except Exception:
        pass

    return {
        "protocol_full_text": joined_modules,
        "change_summary": debate_points[:10],
        "remaining_risks": [],
    }


def start_refinement_session(
    bacteria: str,
    metabolite: str,
    disease: str = "IBD",
    run_id: Optional[str] = None,
    mechanism_summary: str = "",
    protocol_text: str = "",
) -> dict:
    if (protocol_text or "").strip():
        draft = build_protocol_payload_from_text(
            protocol_text=protocol_text,
            bacteria=bacteria,
            metabolite=metabolite,
            disease=disease,
            run_id=run_id,
            mechanism_summary=mechanism_summary,
        )
    else:
        draft = generate_protocol_draft(
            bacteria=bacteria,
            metabolite=metabolite,
            disease=disease,
            run_id=run_id,
            mechanism_summary=mechanism_summary,
        )
    session_id = f"step3ref_{datetime.now().strftime('%Y%m%d_%H%M%S')}_{_slug(bacteria, 12)}_{uuid4().hex[:6]}"
    session = {
        "session_id": session_id,
        "status": "draft_ready",
        "created_at": _now_iso(),
        "updated_at": _now_iso(),
        "candidate": {
            "run_id": run_id,
            "bacteria": bacteria,
            "metabolite": metabolite,
            "disease": disease,
        },
        "title": draft.get("title"),
        "executive_summary": draft.get("executive_summary"),
        "global_hypothesis": draft.get("global_hypothesis"),
        "generated_by": draft.get("generated_by"),
        "model": draft.get("model"),
        "candidate_metrics": draft.get("candidate_metrics") or {},
        "extra_requirements": "",
    }

    modules = decompose_protocol(draft["protocol"], draft["candidate_context"])
    session["status"] = "modules_ready"

    session_dir = _session_dir(session_id)
    session_dir.mkdir(parents=True, exist_ok=True)
    _save_json(session_dir / "session.json", session)
    _save_json(
        session_dir / "protocol_v0.json",
        {
            "protocol": draft["protocol"],
            "candidate_context": draft["candidate_context"],
            "literature_context": draft["literature_context"],
            "literature": draft["literature"],
            "major_risks": draft.get("major_risks", []),
        },
    )
    _save_json(session_dir / "modules_v0.json", {"modules": modules})
    return {
        **session,
        "protocol_v0": draft["protocol"],
        "major_risks": draft.get("major_risks", []),
        "modules": modules,
    }


def confirm_refinement_session(session_id: str, modules: List[dict], extra_requirements: str = "") -> dict:
    session_dir = _session_dir(session_id)
    session = _load_json(session_dir / "session.json")
    if not session:
        raise FileNotFoundError("Refinement session not found.")

    existing = _load_json(session_dir / "modules_v0.json") or {"modules": []}
    by_id = {m.get("module_id"): dict(m) for m in existing.get("modules") or []}
    for incoming in modules:
        mid = incoming.get("module_id")
        if mid in by_id:
            by_id[mid]["enabled"] = incoming.get("enabled", by_id[mid].get("enabled", True))
            by_id[mid]["priority"] = incoming.get("priority", by_id[mid].get("priority", "medium"))
            by_id[mid]["user_notes"] = incoming.get("user_notes", by_id[mid].get("user_notes", ""))
    merged = [by_id[mid] for mid in [MODULE_SPECS[name]["module_id"] for name in MODULE_ORDER] if mid in by_id]

    session["status"] = "confirmed"
    session["updated_at"] = _now_iso()
    session["extra_requirements"] = extra_requirements or ""
    _save_json(session_dir / "session.json", session)
    _save_json(session_dir / "modules_confirmed.json", {"modules": merged})
    return {"session_id": session_id, "status": session["status"], "modules": merged, "extra_requirements": session["extra_requirements"]}


def run_refinement_session(session_id: str) -> dict:
    session_dir = _session_dir(session_id)
    session = _load_json(session_dir / "session.json")
    if not session:
        raise FileNotFoundError("Refinement session not found.")

    protocol_v0 = _load_json(session_dir / "protocol_v0.json") or {}
    modules_blob = _load_json(session_dir / "modules_confirmed.json") or _load_json(session_dir / "modules_v0.json") or {"modules": []}
    modules = modules_blob.get("modules") or []
    candidate_context = protocol_v0.get("candidate_context", "")
    literature_context = protocol_v0.get("literature_context", "")
    extra_requirements = session.get("extra_requirements", "")

    session["status"] = "review_running"
    session["updated_at"] = _now_iso()
    _save_json(session_dir / "session.json", session)

    refined_modules = []
    module_results = []
    debate_points: List[str] = []

    for module in modules:
        if not module.get("enabled", True):
            refined_modules.append(module)
            module_results.append(
                {
                    "module_id": module["module_id"],
                    "module_name": module["module_name"],
                    "status": "disabled",
                    "rounds": [],
                    "final_module": module,
                }
            )
            continue

        current_module = dict(module)
        rounds = []
        for round_no in (1, 2):
            reviews = [
                review_module(agent_type, current_module, candidate_context, literature_context, extra_requirements)
                for agent_type in REVIEW_AGENTS
            ]
            synthesis = synthesize_module(current_module, reviews, candidate_context, extra_requirements)
            rounds.append(
                {
                    "round": round_no,
                    "reviews": reviews,
                    "synthesis": synthesis,
                }
            )
            debate_points.extend((synthesis.get("agreed_issues") or [])[:3])
            debate_points.extend((synthesis.get("conflicting_issues") or [])[:2])
            revised_content = synthesis.get("revised_content") or current_module.get("content")
            current_module = {
                **current_module,
                "content": revised_content,
                "summary": synthesis.get("review_summary") or current_module.get("summary", ""),
                "remaining_risks": synthesis.get("remaining_risks") or [],
                "change_log": synthesis.get("change_log") or [],
            }
            if round_no == 1 and not _needs_second_round(reviews, synthesis):
                break

        refined_modules.append(current_module)
        module_results.append(
            {
                "module_id": module["module_id"],
                "module_name": module["module_name"],
                "status": "reviewed",
                "rounds": rounds,
                "final_module": current_module,
            }
        )

    final_protocol = _build_final_protocol(refined_modules, candidate_context, extra_requirements, debate_points)
    result = {
        "session_id": session_id,
        "status": "completed",
        "candidate": session.get("candidate") or {},
        "protocol_v0": protocol_v0.get("protocol", ""),
        "extra_requirements": extra_requirements,
        "modules": refined_modules,
        "module_results": module_results,
        "final_protocol": final_protocol.get("protocol_full_text", ""),
        "change_summary": final_protocol.get("change_summary", []),
        "remaining_risks": final_protocol.get("remaining_risks", []),
        "updated_at": _now_iso(),
    }

    session["status"] = "completed"
    session["updated_at"] = result["updated_at"]
    _save_json(session_dir / "session.json", session)
    _save_json(session_dir / "result.json", result)
    return result


def get_refinement_session_result(session_id: str) -> Optional[dict]:
    session_dir = _session_dir(session_id)
    result = _load_json(session_dir / "result.json")
    if result:
        return result
    session = _load_json(session_dir / "session.json")
    if not session:
        return None
    protocol_v0 = _load_json(session_dir / "protocol_v0.json") or {}
    modules = (_load_json(session_dir / "modules_confirmed.json") or _load_json(session_dir / "modules_v0.json") or {"modules": []}).get("modules") or []
    return {
        "session_id": session_id,
        "status": session.get("status"),
        "candidate": session.get("candidate") or {},
        "protocol_v0": protocol_v0.get("protocol", ""),
        "modules": modules,
        "extra_requirements": session.get("extra_requirements", ""),
        "title": session.get("title"),
        "executive_summary": session.get("executive_summary"),
        "global_hypothesis": session.get("global_hypothesis"),
        "candidate_metrics": session.get("candidate_metrics") or {},
    }
