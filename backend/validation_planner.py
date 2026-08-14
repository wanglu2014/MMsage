"""
Evidence-driven Step 3 validation planning.
"""

from __future__ import annotations

import json
import re
from typing import Any, Callable, Dict, List, Optional
from urllib.error import HTTPError, URLError
from urllib.parse import quote, urlencode
from urllib.request import Request, urlopen

from agents.pubmed_query import build_pubmed_query
from build_kg import fetch_abstracts, search_pubmed
from protocol_refiner import _chat_json, _prepare_protocol_support, articles_to_export
from protocols_io_tool import collect_protocol_evidence

MAX_QUERY_RESULTS = 6
QUESTION_MAX_QUERY_RESULTS = 12
MAX_QUERY_ROUNDS = 3
FULLTEXT_MAX_ARTICLES = 6
FULLTEXT_REQUEST_TIMEOUT = 20
PMC_BIOC_URL = "https://www.ncbi.nlm.nih.gov/research/bionlp/RESTful/pmcoa.cgi/BioC_json"
UNIVERSAL_EXPERIMENT_DESIGN_PROMPT_RULES = """Universal experimental-design reasoning contract:
- Start with the core scientific question and the most consequential unresolved uncertainty. Build the smallest non-duplicative set of experiments needed to distinguish the competing explanations; never add modules to meet a numerical quota.
- Keep the response at research-strategy level. State what is studied, why it is needed, the study object, explicit groups, collected samples and timing, primary and secondary indicators, key controls, and the positive/negative decision branches. Do not write an SOP, equipment list, numbered procedure, exhaustive dose schedule, or full statistical analysis plan.
- Return structured group objects with group_name, exposure_or_condition, and control_purpose. Do not return group_count; it is derived from the groups array.
- Use one scientific question and one primary indicator per module. Additional indicators must be clearly secondary.
- Distinguish route priority from evidence state. route_status describes where a module belongs in the route; result_status records not_run, positive, negative, or inconclusive; execution_status records whether its activation gate is currently open.
- A downstream experiment is conditional until its structured activation_gate is satisfied by actual upstream result_status values. A module being ready to run is not evidence that it has produced a positive result.
- Use only these canonical module IDs and roles in a microbe-metabolite-disease workflow: V1=direct_production_gate, V2=indirect_ecology_gate, V3=host_response_gate, A1=effect_interaction_gate, A2=direct_causal_rescue, A3=indirect_ecology_causal, and H1=human_translation. Never invent M1/M2/M3 aliases.
- For a microbe-metabolite-disease workflow, every live competing hypothesis must be linked to at least one named experimental module. Use hypothesis_ids to make that link explicit; a conditional module still needs a complete decision-level experiment card and is not a placeholder.
- unlock_rule means the upstream condition that permits this module to start; it must never depend on the module's own result. Put consequences of this module's result only in branch_if_positive and branch_if_negative. The application will canonicalize activation gates and unlock rules.
- State an explicit claim_boundary for every module. Association, host activity, efficacy, interaction, microbial production, and causal mediation are different claims and must not be substituted for one another.
- Exact strains, media, doses, routes, concentrations, and time values are optional in this high-level design. Include them only when strongly candidate-relevant evidence supports them; otherwise state that protocol optimization follows after the decision gate is selected.
- Cite only evidence that directly supports the stated edge of the evidence chain. For each cited claim state what the paper supports and what it does not support. Do not combine separate papers into an untested complete causal chain."""
MICROBE_METABOLITE_IN_VITRO_PROMPT_RULES = """Universal microbe-metabolite in vitro and cell rules:
- Separate three explanations when they are relevant: direct microbial production, indirect ecological regulation, and independent or interactive host effects. Never assume the metabolite is the microbe's only or principal mediator.
- V1 is the direct-production source gate. At minimum distinguish live candidate culture, uninoculated medium, and a biologically justified substrate or precursor control; include baseline and later culture samples; measure net formation of the target metabolite together with microbial growth or abundance. A time-dependent concentration change alone is insufficient if carryover, medium background, substrate conversion, or normalization to microbial biomass remains unresolved.
- V2 is the indirect-ecology source gate. If direct production is negative or inconclusive, first select one ecological system rather than listing conditioned medium, co-culture, and a defined community as interchangeable groups. Confirm any proposed producer independently, then compare matched candidate-present and candidate-absent conditions within that same system while measuring producer abundance or activity and target-metabolite output. Correlation alone is not production evidence.
- A positive direct-production result does not exclude ecological amplification. Keep indirect ecology deferred unless direct output is insufficient or ecological amplification is itself a study objective.
- V3 may run in parallel with V1. Choose one disease-relevant host model, one disease-relevant challenge, and a challenge-condition microbe-by-metabolite 2 x 2 core with an unchallenged baseline reference. Prespecify one directional primary host indicator and require the effect to occur without unacceptable toxicity; do not return a menu of cell models or a primary-indicator menu.
- Map H1 only to V1 direct-source testing and, if later warranted, A2 direct causal follow-up; map H2 only to V2 ecological-source testing and, if later warranted, A3 ecological causal follow-up; map H3 to V3 host activity and A1 effect/interaction testing. Do not substitute module IDs for hypothesis IDs.
- A cell response demonstrates host activity only. It does not demonstrate microbial production, disease efficacy, or metabolite mediation.
- Candidate-specific substrates, media, tracers, cell lines, and analytical techniques must come from the current candidate and evidence. Never transfer a candidate-specific detail into an unrelated plan."""
MICROBE_METABOLITE_IN_VIVO_PROMPT_RULES = """Universal microbe-metabolite animal rules:
- The first animal module tests disease effects and the microbe-by-metabolite interaction without claiming mediation. Use a healthy reference plus the disease-condition 2 x 2 core: disease control, disease plus microbe, disease plus metabolite, and disease plus both, unless a different grouping is scientifically required and explained.
- In A1, state one animal object; the healthy reference plus all four disease-condition factorial groups; baseline and endpoint samples; direct confirmation of microbial exposure or colonization and target-metabolite exposure; one primary disease indicator; and a small set of secondary disease or mechanism indicators. Define separate decision contrasts for the microbial main effect, metabolite main effect, and their interaction rather than calling improvement in any arm an interaction.
- The effect-and-interaction module may proceed after positive host activity or strong prior host-effect evidence; it does not require direct or indirect microbial production to be positive.
- Direct mediation follow-up requires positive direct-production and animal-effect results, then a function-loss or otherwise specific necessity test plus metabolite rescue when feasible.
- Ecological mediation follow-up requires positive indirect-ecology and animal-effect results, then a controlled producer/community necessity test plus rescue when feasible.
- If the animal effect is positive but the causal follow-up fails, interpret the microbe and metabolite as parallel or partly overlapping effects, not as demonstrated mediation.
- A1 is the H3 effect-and-interaction experiment. A2 is the conditional H1 causal follow-up and A3 is the conditional H2 causal follow-up; retain those hypothesis_ids even when their execution status is conditional_future.
- Return brief conditional blueprints for A2 and A3 whenever their corresponding V1 or V2 source branch and A1 are present. Their result_status remains not_run and their code-owned gates keep them conditional; including a blueprint does not authorize execution. Do not invent mutants, producers, doses, disease models, or mechanisms."""
MICROBE_METABOLITE_HUMAN_PROMPT_RULES = """Universal conditional human-study rules:
- Preserve one concise observational human blueprint for a microbe-metabolite-disease workflow unless the user explicitly excludes human research.
- The blueprint must list the study population, groups or disease strata, biospecimens and timing, primary and secondary indicators, major confounder domains, and its unlock rule. Do not provide intervention doses or recruitment procedures.
- Keep the module conditional until a relevant animal branch is positive and independently replicated and ethics plus the analysis plan are ready.
- Initial human work tests association, direction, and temporal compatibility only. It cannot prove microbial production, causal mediation, or disease causation.
- Disease subtype and strata must follow the user's question or strong disease-specific evidence. Do not silently replace a broad disease with a preferred subtype."""
METHOD_SECTION_TERMS = ("method", "methods", "materials", "experimental", "materials|methods")
CELL_TERMS = (
    "in vitro",
    "cell",
    "cell line",
    "epithelial",
    "organoid",
    "macrophage",
    "monocyte",
    "co-culture",
)
ANIMAL_TERMS = (
    "mouse",
    "mice",
    "murine",
    "animal",
    "rat",
    "colitis",
    "gnotobiotic",
    "germ-free",
)
CONTRADICT_TERMS = (
    "no significant",
    "not significant",
    "contradict",
    "negative result",
    "null finding",
    "not associated",
    "inconsistent",
    "failed to",
)
DIRECT_CULTURE_TERMS = (
    "monoculture",
    "mono-culture",
    "pure culture",
    "axenic culture",
    "defined medium",
    "culture medium",
    "culture supernatant",
    "cell-free supernatant",
    "bacterial culture",
    "fermentation culture",
)
STRICT_DIRECT_CULTURE_TERMS = (
    "monoculture",
    "mono-culture",
    "pure culture",
    "axenic culture",
    "axenic",
    "isolated strain",
    "single strain",
    "was cultured",
    "were cultured",
    "was grown",
    "were grown",
    "cultured anaerobically",
    "grown anaerobically",
)
METABOLITE_PRODUCTION_TERMS = (
    "produce",
    "produces",
    "produced",
    "production",
    "biosynthesis",
    "biosynthetic",
    "synthesize",
    "synthesizes",
    "secreted",
    "secretion",
    "accumulated",
    "accumulation",
    "formed",
    "formation",
    "fermentation product",
)
METABOLITE_MEASUREMENT_TERMS = (
    "gc-ms",
    "gc ms",
    "lc-ms",
    "lc ms",
    "mass spectrometry",
    "targeted metabolomics",
    "quantified",
    "quantification",
    "measured",
    "concentration",
    "time course",
    "isotope tracing",
    "stable isotope",
    "gas chromatography",
    "liquid chromatography",
    "chromatography",
    "hplc",
    "nmr",
    "detected",
    "main metabolites",
    "metabolites were",
    "mg/l",
    "mmol/l",
    "μmol/l",
    "umol/l",
)
INDIRECT_PRODUCTION_TERMS = (
    "co-culture",
    "coculture",
    "cross-feeding",
    "cross feeding",
    "conditioned medium",
    "defined community",
    "microbial consortium",
    "community metabolism",
    "substrate supply",
    "mucin degradation",
    "metabolic interaction",
)
REVIEW_LIKE_TERMS = (
    "systematic review",
    "scoping review",
    "narrative review",
    "review article",
    "a review",
    "meta-analysis",
    "this review",
    "we review",
    "review of",
    "overview of",
    "an overview",
    "perspective article",
    "commentary article",
)
STRONG_CITATION_EVIDENCE_TYPES = {
    "direct_monoculture_production",
    "direct_monoculture_nonproduction",
    "candidate_microbe_disease_intervention",
    "candidate_metabolite_disease_intervention",
    "indirect_ecological_evidence",
}
DISEASE_SYNONYMS = {
    "ibd": [
        "IBD",
        "inflammatory bowel disease",
        "Crohn disease",
        "Crohn's disease",
        "ulcerative colitis",
        "colitis",
    ],
    "type 2 diabetes": ["type 2 diabetes", "type II diabetes", "T2D", "T2DM"],
    "chronic kidney disease": ["chronic kidney disease", "CKD"],
    "nonalcoholic fatty liver disease": ["nonalcoholic fatty liver disease", "NAFLD"],
    "metabolic dysfunction-associated steatotic liver disease": [
        "metabolic dysfunction-associated steatotic liver disease",
        "MASLD",
    ],
    "colorectal cancer": ["colorectal cancer", "CRC"],
}
CELL_MODEL_PATTERN = re.compile(
    r"\b(?:Caco-2|HT-29|RAW ?264\.7|THP-1|HCT116|IEC-6|MODE-K|LoVo|SW480|MC38|HIEC|CCD[- ]841|T84|LS174T)\b",
    re.IGNORECASE,
)
ANIMAL_MODEL_PATTERN = re.compile(
    r"\b(?:C57BL/6(?:J)?|BALB/c|Sprague[- ]Dawley|Wistar|germ-free|specific pathogen-free|SPF|DSS(?:-induced)? colitis|TNBS(?:-induced)? colitis|AOM/DSS|Il10-/-|IL-10 knockout|NOD/SCID)\b",
    re.IGNORECASE,
)
CONDITION_PATTERN = re.compile(
    r"\b(?:\d+(?:\.\d+)?(?:\s?[x×]\s?10\^?\d+)?)\s?(?:uM|μM|mM|nM|pM|mg/kg(?:/day)?|g/kg|mg/mL|ug/mL|μg/mL|ng/mL|g/L|CFU(?:/mL)?|CFU|MOI|% DSS)\b",
    re.IGNORECASE,
)
TIME_PATTERN = re.compile(
    r"\b\d+(?:\.\d+)?\s?(?:h|hr|hrs|hour|hours|day|days|week|weeks)\b",
    re.IGNORECASE,
)


def _report_progress(
    progress_callback: Optional[Callable[[str, str, int, Optional[Dict[str, Any]]], None]],
    stage: str,
    message: str,
    percent: int,
    extra: Optional[Dict[str, Any]] = None,
) -> None:
    if progress_callback is None:
        return
    progress_callback(stage, message, percent, extra or {})


def _as_str_list(value: Any) -> List[str]:
    if not isinstance(value, list):
        return []
    return [str(item).strip() for item in value if str(item).strip()]


def _dedupe_preserve(items: List[str]) -> List[str]:
    seen = set()
    ordered = []
    for item in items:
        clean = str(item).strip()
        key = clean.lower()
        if clean and key not in seen:
            seen.add(key)
            ordered.append(clean)
    return ordered


def _truncate(text: str, limit: int = 520) -> str:
    clean = re.sub(r"\s+", " ", str(text or "")).strip()
    if len(clean) <= limit:
        return clean
    return clean[: limit - 3].rstrip() + "..."


def _contains_any(text: str, terms: List[str] | tuple[str, ...]) -> bool:
    lowered = str(text or "").lower()
    return any(term.lower() in lowered for term in terms if term)


def _contains_any_complete_term(text: str, terms: List[str] | tuple[str, ...]) -> bool:
    lowered = str(text or "").lower()
    for term in terms:
        clean_term = str(term or "").lower().strip()
        if not clean_term:
            continue
        if re.search(r"(?<![a-z0-9])" + re.escape(clean_term) + r"(?![a-z0-9])", lowered):
            return True
    return False


def _entity_aliases(entity: str) -> List[str]:
    value = str(entity or "").strip()
    if not value:
        return []
    aliases = [value, value.replace("_", " ")]
    parts = value.replace("_", " ").split()
    if parts:
        aliases.append(parts[0])
    lowered = value.replace("_", " ").lower()
    if "akkermansia muciniphila" in lowered:
        aliases.extend(["A. muciniphila", "Akkermansia"])
    if lowered in {
        "isobutyric acid",
        "isobutyrate",
        "iso-butyrate",
        "2-methylpropanoic acid",
        "sodium isobutyrate",
    }:
        aliases.extend(
            [
                "isobutyric acid",
                "isobutyrate",
                "iso-butyrate",
                "2-methylpropanoic acid",
                "sodium isobutyrate",
            ]
        )
    return _dedupe_preserve(aliases)


def _citation_entity_aliases(entity: str, is_microbe: bool = False) -> List[str]:
    """Return conservative aliases for deciding whether a paper may be cited.

    Literature retrieval may use broad genus or first-token variants, but those
    variants are unsafe for citation adjudication because they can match a
    different species or a chemically distinct multi-word metabolite.
    """
    value = str(entity or "").strip()
    if not value:
        return []
    normalized = re.sub(r"\s+", " ", value.replace("_", " ")).strip()
    aliases = [value, normalized]
    parts = normalized.split()
    if is_microbe and len(parts) >= 2 and parts[0] and parts[1]:
        aliases.extend(
            [
                f"{parts[0][0]}. {parts[1]}",
                f"{parts[0][0]}.{parts[1]}",
            ]
        )
    lowered = normalized.lower()
    if lowered in {
        "isobutyric acid",
        "isobutyrate",
        "iso-butyrate",
        "2-methylpropanoic acid",
        "sodium isobutyrate",
    }:
        aliases.extend(
            [
                "isobutyric acid",
                "isobutyrate",
                "iso-butyrate",
                "2-methylpropanoic acid",
                "sodium isobutyrate",
            ]
        )
    return _dedupe_preserve(aliases)


def _disease_aliases(disease: str) -> List[str]:
    aliases = _entity_aliases(disease)
    lowered = str(disease or "").strip().lower()
    for key, group in DISEASE_SYNONYMS.items():
        if lowered == key or any(lowered == str(item).strip().lower() for item in group):
            aliases.extend(group)
    return _dedupe_preserve(aliases)


def _citation_disease_aliases(disease: str) -> List[str]:
    value = str(disease or "").strip()
    if not value:
        return []
    aliases = [value, value.replace("_", " ")]
    lowered = value.replace("_", " ").lower()
    for key, group in DISEASE_SYNONYMS.items():
        if lowered == key or any(lowered == str(item).strip().lower() for item in group):
            aliases.extend(group)
    tokens = re.findall(r"[A-Za-z]+|\d+", value.replace("_", " "))
    acronym = "".join(token if token.isdigit() else token[0] for token in tokens if token.lower() not in {"of", "and", "the"})
    if 2 <= len(acronym) <= 8 and sum(char.isalpha() for char in acronym) >= 2:
        aliases.append(acronym.upper())
    return _dedupe_preserve(aliases)


QUESTION_QUERY_STOPWORDS = {
    "about", "after", "against", "among", "and", "animal", "animals", "basic",
    "because", "between", "brief", "cell", "cells", "clinical", "cohort",
    "common", "complete", "constraint", "constraints", "control", "controls",
    "data", "decision", "design", "disease", "dose", "doses", "effect",
    "evidence", "experiment", "experimental", "experiments", "focus", "for",
    "from", "group", "groups", "host", "human", "humans", "in", "into",
    "keep", "less", "materials", "mechanism", "methods", "model", "models",
    "mouse", "mice", "need", "not", "operational", "or", "patient", "patients",
    "plan", "priority", "procedure", "protocol", "protocols", "question",
    "readout", "readouts", "relationship", "response", "rule", "settings",
    "simple", "step", "study", "summary", "than", "that", "the", "their",
    "them", "then", "these", "this", "time", "timing", "use", "using",
    "validation", "what", "when", "where", "which", "with", "without",
}


def _extract_question_keywords(text: str, limit: int = 8) -> List[str]:
    tokens = re.findall(r"[A-Za-z][A-Za-z0-9_\-/+]{2,}", str(text or ""))
    picked: List[str] = []
    for token in tokens:
        clean = token.strip(" _-/+").lower()
        if not clean or clean in QUESTION_QUERY_STOPWORDS:
            continue
        if re.fullmatch(r"\d+", clean):
            continue
        picked.append(token)
    return _dedupe_preserve(picked)[:limit]


def _append_query_spec(
    specs: List[dict],
    seen_queries: set[str],
    query: str,
    label: str,
    rationale: str,
    round_id: int = 1,
) -> None:
    clean_query = re.sub(r"\s+", " ", str(query or "")).strip()
    if not clean_query:
        return
    key = clean_query.lower()
    if key in seen_queries:
        return
    seen_queries.add(key)
    specs.append(
        {
            "round": round_id,
            "label": label,
            "rationale": rationale,
            "query": clean_query[:600],
        }
    )


def _build_question_keyword_query_specs(
    research_question: str,
    prompt_constraints: str,
    disease: str = "",
    max_specs: int = 8,
) -> List[dict]:
    source_text = " ".join(
        part for part in [research_question, prompt_constraints] if str(part or "").strip()
    )
    keywords = _extract_question_keywords(source_text, limit=24)
    disease_block = build_pubmed_query([disease]) if str(disease or "").strip() else ""
    specs: List[dict] = []
    seen_queries: set[str] = set()

    if not keywords:
        _append_query_spec(
            specs,
            seen_queries,
            _fallback_question_query(research_question, disease),
            "Question keyword fallback",
            "Fallback keyword search when no stable English keywords can be extracted from the question.",
        )
        return specs

    windows: List[List[str]] = []
    windows.append(keywords[:6])
    windows.append(keywords[:4])
    if len(keywords) > 4:
        windows.append(keywords[2:6])
    if len(keywords) > 6:
        windows.append(keywords[4:8])
    if len(keywords) > 8:
        windows.append(keywords[6:10])
    if len(keywords) > 10:
        windows.append(keywords[8:12])

    chunk_size = 3 if len(keywords) >= 9 else 2
    for idx in range(0, len(keywords), chunk_size):
        windows.append(keywords[idx : idx + chunk_size])

    for idx, terms in enumerate(windows, start=1):
        cleaned_terms = _dedupe_preserve([str(term).strip() for term in terms if str(term).strip()])
        if not cleaned_terms:
            continue
        keyword_block = build_pubmed_query(cleaned_terms[:6])
        query = f"{keyword_block} AND {disease_block}" if disease_block else keyword_block
        _append_query_spec(
            specs,
            seen_queries,
            query,
            f"Question keyword set {idx}",
            "Keyword-derived search from the standalone question to broaden literature retrieval.",
        )
        if len(specs) >= max_specs:
            break

    if not specs:
        _append_query_spec(
            specs,
            seen_queries,
            _fallback_question_query(research_question, disease),
            "Question keyword fallback",
            "Fallback keyword search when keyword chunking did not yield a valid query.",
        )
    return specs[:max_specs]


def _entity_pattern(alias: str) -> re.Pattern[str]:
    term = re.sub(r"\s+", " ", str(alias or "").strip().lower())
    escaped = re.escape(term).replace(r"\ ", r"\s+")
    return re.compile(r"(?<![a-z0-9])" + escaped + r"(?![a-z0-9])", re.IGNORECASE)


def _contains_entity(text: str, aliases: List[str]) -> bool:
    clean = re.sub(r"\s+", " ", str(text or "")).strip()
    return any(_entity_pattern(alias).search(clean) for alias in aliases if str(alias or "").strip())


def _contains_nonnegated_entity(text: str, aliases: List[str]) -> bool:
    lowered = re.sub(r"\s+", " ", str(text or "")).strip().lower()
    negation_pattern = re.compile(r"\b(?:no|not|without|neither|nor|did not|was not|were not|is not|are not|failed to|unable to)\b")
    for alias in aliases:
        if not str(alias or "").strip():
            continue
        for match in _entity_pattern(alias).finditer(lowered):
            prefix = lowered[max(0, match.start() - 55) : match.start()]
            if not negation_pattern.search(prefix):
                return True
    return False


def _iter_local_windows(text: str, max_sentences: int = 2) -> List[str]:
    clean = re.sub(r"\s+", " ", str(text or "")).strip()
    sentences = [part.strip() for part in re.split(r"(?<=[.!?])\s+", clean) if part.strip()]
    if not sentences:
        return []
    windows: List[str] = []
    for index in range(len(sentences)):
        for width in range(1, max_sentences + 1):
            window = " ".join(sentences[index : index + width]).strip()
            if window:
                windows.append(window)
    return _dedupe_preserve(windows)


def _match_is_negated(text: str, start: int, end: int) -> bool:
    lowered = str(text or "").lower()
    context = lowered[max(0, start - 65) : min(len(lowered), end + 25)]
    return bool(
        re.search(
            r"\b(?:no|not|without|neither|nor|did not|does not|was not|were not|failed to|unable to|lack(?:ed|s|ing)?)\b",
            context,
        )
    )


def _candidate_linked_culture_signal(window: str, bacteria_aliases: List[str]) -> bool:
    lowered = str(window or "").lower()
    strict_single = _contains_any_complete_term(window, STRICT_DIRECT_CULTURE_TERMS)
    ecology_terms = ("co-culture", "coculture", "mixed culture", "community", "consortium", "together with")
    if _contains_any(lowered, ecology_terms) and not _contains_any_complete_term(
        window, ("monoculture", "mono-culture", "pure culture", "axenic culture", "single strain")
    ):
        return False
    for alias in bacteria_aliases:
        if not str(alias or "").strip():
            continue
        entity = _entity_pattern(alias).pattern
        patterns = (
            rf"{entity}.{{0,100}}\b(?:was|were|is|are)?\s*(?:culture|cultures|cultured|grown|incubated|fermented)\b",
            rf"\b(?:culture|cultures|cultured|grown|incubated|fermented)\b.{{0,100}}{entity}",
        )
        if any(re.search(pattern, lowered, flags=re.IGNORECASE) for pattern in patterns):
            return strict_single or not _contains_any(lowered, ecology_terms)
    return False


def _positive_metabolite_production_signal(window: str, metabolite_aliases: List[str]) -> bool:
    lowered = str(window or "").lower()
    production = r"(?:produc(?:e|es|ed|tion)|synthesi[sz](?:e|es|ed|is)|secret(?:e|es|ed|ion)|form(?:ed|ation)|accumulat(?:e|es|ed|ion)|generated?)"
    for alias in metabolite_aliases:
        if not str(alias or "").strip():
            continue
        entity = _entity_pattern(alias).pattern
        patterns = (
            rf"\b{production}\b.{{0,140}}{entity}",
            rf"{entity}.{{0,100}}\b(?:was|were|is|are|showed|shows)?\s*{production}\b",
            rf"\b(?:production|formation|accumulation|secretion)\s+of\s+{entity}",
        )
        for pattern in patterns:
            for match in re.finditer(pattern, lowered, flags=re.IGNORECASE):
                if not _match_is_negated(lowered, match.start(), match.end()):
                    return True
    return False


def _negative_metabolite_production_signal(window: str, metabolite_aliases: List[str]) -> bool:
    lowered = str(window or "").lower()
    for alias in metabolite_aliases:
        if not str(alias or "").strip():
            continue
        entity = _entity_pattern(alias).pattern
        patterns = (
            rf"\b(?:no|without)\b.{{0,60}}\b(?:production|formation|accumulation|secretion)\b.{{0,80}}{entity}",
            rf"{entity}.{{0,45}}\b(?:was|were|is|are)?\s*(?:not|never)\s+(?:produced|formed|detected|accumulated|increased|secreted)\b",
            rf"\b(?:failed|unable)\s+to\s+(?:produce|form|generate|secrete)\b.{{0,80}}{entity}",
            rf"\b(?:no|not)\s+(?:detectable|significant)?\s*{entity}.{{0,50}}\b(?:production|formation|accumulation|increase)\b",
        )
        if any(re.search(pattern, lowered, flags=re.IGNORECASE) for pattern in patterns):
            return True
    return False


def _production_attributed_to_other_microbe(
    window: str,
    bacteria_aliases: List[str],
    metabolite_aliases: List[str],
) -> bool:
    candidate_text = " ".join(str(alias or "").lower() for alias in bacteria_aliases)
    relation = r"(?:produced|produces|formed|forms|generated|generates|secreted|secretes)"
    binomial = r"([A-Z][a-z]{2,}\s+[a-z][a-z0-9-]{2,})"
    for metabolite_alias in metabolite_aliases:
        if not str(metabolite_alias or "").strip():
            continue
        metabolite = rf"(?i:{_entity_pattern(metabolite_alias).pattern})"
        patterns = (
            rf"{binomial}\s+{relation}.{{0,100}}{metabolite}",
            rf"{metabolite}.{{0,100}}{relation}\s+by\s+{binomial}",
        )
        for pattern in patterns:
            for match in re.finditer(pattern, str(window or "")):
                named_microbe = str(match.group(1) or "").lower().strip()
                if named_microbe and named_microbe not in candidate_text:
                    return True
    return False


def _candidate_direct_culture_signal(
    text: str,
    bacteria_aliases: List[str],
    metabolite_aliases: List[str],
) -> bool:
    """Conservative screen for same-paper, candidate-specific direct culture evidence."""
    clean = re.sub(r"\s+", " ", str(text or "")).strip()
    if not clean or not _contains_entity(clean, bacteria_aliases) or not _contains_entity(clean, metabolite_aliases):
        return False
    if _contains_any(clean, REVIEW_LIKE_TERMS):
        return False
    for window in _iter_local_windows(clean, max_sentences=3):
        if not (_contains_entity(window, bacteria_aliases) and _contains_entity(window, metabolite_aliases)):
            continue
        if (
            _candidate_linked_culture_signal(window, bacteria_aliases)
            and _positive_metabolite_production_signal(window, metabolite_aliases)
            and not _negative_metabolite_production_signal(window, metabolite_aliases)
            and not _production_attributed_to_other_microbe(window, bacteria_aliases, metabolite_aliases)
            and _contains_any_complete_term(window, METABOLITE_MEASUREMENT_TERMS)
        ):
            return True
    return False


def _candidate_direct_culture_nonproduction_signal(
    text: str,
    bacteria_aliases: List[str],
    metabolite_aliases: List[str],
) -> bool:
    clean = re.sub(r"\s+", " ", str(text or "")).strip()
    if not clean or _contains_any(clean, REVIEW_LIKE_TERMS):
        return False
    for window in _iter_local_windows(clean, max_sentences=3):
        if not (_contains_entity(window, bacteria_aliases) and _contains_entity(window, metabolite_aliases)):
            continue
        if (
            _candidate_linked_culture_signal(window, bacteria_aliases)
            and _negative_metabolite_production_signal(window, metabolite_aliases)
            and _contains_any_complete_term(window, METABOLITE_MEASUREMENT_TERMS)
        ):
            return True
    return False


def _candidate_intervention_signal(
    text: str,
    entity_aliases: List[str],
    disease_aliases: List[str],
) -> bool:
    """Require wording that the named entity itself, not a cocktail or diet, was administered."""
    clean = re.sub(r"\s+", " ", str(text or "")).strip()
    if not clean or not _contains_nonnegated_entity(clean, entity_aliases) or not _contains_nonnegated_entity(clean, disease_aliases):
        return False
    if _contains_any(clean, REVIEW_LIKE_TERMS):
        return False
    action = r"(?:administered|treated|supplemented|gavaged|dosed|exposed|incubated|challenged)"
    noun = r"(?:administration|supplementation|gavage|treatment|dosing|exposure|addition)"
    cocktail_terms = (
        "cocktail", "consortium", "multi-strain", "multistrain", "compound probiotic",
        "mixture", "blend", "synbiotic", "formula", "as part of", "containing strains",
    )
    for window in _iter_local_windows(clean, max_sentences=2):
        if not (
            _contains_nonnegated_entity(window, entity_aliases)
            and _contains_nonnegated_entity(window, disease_aliases)
        ):
            continue
        lowered = window.lower()
        if _contains_any(lowered, cocktail_terms) and not _contains_any(
            lowered, ("alone", "only group", "single-strain", "monotherapy", "individual strain")
        ):
            continue
        for alias in entity_aliases:
            if not str(alias or "").strip():
                continue
            entity = _entity_pattern(alias).pattern
            patterns = (
                rf"\b{action}\b\s+(?:(?:cells|animals|mice|rats|organoids)\s+)?(?:(?:with|to)\s+)?(?:live\s+|pasteurized\s+|heat-killed\s+)?{entity}",
                rf"\b{noun}\b\s+(?:of|with)\s+(?:live\s+|pasteurized\s+|heat-killed\s+)?{entity}",
                rf"{entity}\s+(?:itself\s+)?(?:was|were)\s+{action}\b",
                rf"{entity}\s+{noun}\b",
            )
            for pattern in patterns:
                for match in re.finditer(pattern, lowered, flags=re.IGNORECASE):
                    if not _match_is_negated(lowered, match.start(), match.end()):
                        return True
    return False


def _candidate_effect_signal(
    text: str,
    entity_aliases: List[str],
    disease_aliases: List[str],
) -> bool:
    """Accept an explicit candidate effect while excluding observational level changes."""
    clean = re.sub(r"\s+", " ", str(text or "")).strip()
    if not clean or _contains_any(clean, REVIEW_LIKE_TERMS):
        return False
    observational_terms = (
        "associated with", "correlated with", "correlation", "predictive of", "higher fecal",
        "lower fecal", "fecal level", "fecal concentration", "plasma level", "serum level",
        "circulating level", "abundance", "observed level", "measured level",
    )
    effect = r"(?:ameliorat(?:e|es|ed)|alleviat(?:e|es|ed)|attenuat(?:e|es|ed)|suppress(?:es|ed)?|inhibit(?:s|ed)?|improv(?:e|es|ed)|protect(?:s|ed)?|restor(?:e|es|ed)|reduc(?:e|es|ed))"
    model_terms = (
        "in vivo", "in vitro", "mouse", "mice", "murine", "animal", "cell", "cells",
        "organoid", "epithelial", "macrophage", "disease model",
    )
    for window in _iter_local_windows(clean, max_sentences=2):
        if not (
            _contains_nonnegated_entity(window, entity_aliases)
            and _contains_nonnegated_entity(window, disease_aliases)
            and _contains_any(window, model_terms)
        ):
            continue
        if _contains_any(window, observational_terms):
            continue
        lowered = window.lower()
        for alias in entity_aliases:
            if not str(alias or "").strip():
                continue
            entity = _entity_pattern(alias).pattern
            patterns = (
                rf"{entity}\s+(?:(?:itself|treatment|exposure|supplementation)\s+)?(?:was|were|is|are)?\s*(?:found\s+to\s+)?{effect}\b",
                rf"\b{effect}\b.{{0,45}}(?:by|with|after)\s+{entity}",
            )
            for pattern in patterns:
                for match in re.finditer(pattern, lowered, flags=re.IGNORECASE):
                    if not _match_is_negated(lowered, match.start(), match.end()):
                        return True
    return False


def _candidate_entity_manipulation_signal(text: str, entity_aliases: List[str]) -> bool:
    """Require the named entity itself to be added, removed, inoculated, cultured, or administered."""
    clean = re.sub(r"\s+", " ", str(text or "")).strip()
    if not clean or not _contains_entity(clean, entity_aliases) or _contains_any(clean, REVIEW_LIKE_TERMS):
        return False
    action = r"(?:inoculated|cultured|grown|incubated|added|removed|depleted|administered|supplemented|gavaged|dosed)"
    noun = r"(?:inoculation|culture|addition|removal|depletion|administration|supplementation|gavage|dosing)"
    for alias in entity_aliases:
        if not str(alias or "").strip():
            continue
        entity = _entity_pattern(alias).pattern
        patterns = (
            rf"\b{action}\b\s+(?:with\s+|of\s+)?{entity}",
            rf"\b{noun}\b\s+(?:of|with)\s+{entity}",
            rf"{entity}.{{0,60}}\b{action}\b",
            rf"{entity}.{{0,60}}\b{noun}\b",
        )
        for pattern in patterns:
            for match in re.finditer(pattern, clean.lower(), flags=re.IGNORECASE):
                if not _match_is_negated(clean, match.start(), match.end()):
                    return True
    return False


def _candidate_indirect_ecology_signal(
    text: str,
    bacteria_aliases: List[str],
    metabolite_aliases: List[str],
) -> bool:
    clean = re.sub(r"\s+", " ", str(text or "")).strip()
    if not clean or _contains_any(clean, REVIEW_LIKE_TERMS):
        return False
    for window in _iter_local_windows(clean, max_sentences=2):
        if not (
            _contains_nonnegated_entity(window, bacteria_aliases)
            and _contains_entity(window, metabolite_aliases)
            and _contains_any(window, INDIRECT_PRODUCTION_TERMS)
            and _candidate_entity_manipulation_signal(window, bacteria_aliases)
            and _contains_any_complete_term(window, METABOLITE_MEASUREMENT_TERMS)
        ):
            continue
        return True
    return False


def _classify_candidate_article_evidence(
    text: str,
    bacteria_aliases: List[str],
    metabolite_aliases: List[str],
    disease_aliases: List[str],
) -> tuple[str, str, bool]:
    """Return evidence type, candidate relevance, and whether the PMID may be cited."""
    has_b = _contains_entity(text, bacteria_aliases)
    has_m = _contains_entity(text, metabolite_aliases)
    has_d = _contains_entity(text, disease_aliases)
    direct_signal = _candidate_direct_culture_signal(text, bacteria_aliases, metabolite_aliases)
    direct_negative_signal = _candidate_direct_culture_nonproduction_signal(
        text, bacteria_aliases, metabolite_aliases
    )
    microbe_intervention = has_b and has_d and _candidate_intervention_signal(
        text, bacteria_aliases, disease_aliases
    )
    metabolite_intervention = has_m and has_d and _candidate_intervention_signal(
        text, metabolite_aliases, disease_aliases
    )
    indirect_signal = has_b and has_m and _candidate_indirect_ecology_signal(
        text, bacteria_aliases, metabolite_aliases
    )
    if direct_signal:
        return "direct_monoculture_production", "direct", True
    if direct_negative_signal:
        return "direct_monoculture_nonproduction", "direct", True
    if microbe_intervention:
        return "candidate_microbe_disease_intervention", "direct", True
    if metabolite_intervention:
        return "candidate_metabolite_disease_intervention", "direct", True
    if indirect_signal:
        return "indirect_ecological_evidence", "direct", True
    if has_b and has_m:
        return "candidate_pair_association", "partial", False
    if (has_b or has_m) and _contains_any(text, DIRECT_CULTURE_TERMS):
        return "analogous_method_only", "analogous", False
    if has_b or has_m or has_d:
        return "background_only", "partial", False
    return "unrelated", "unrelated", False


def _http_get_json(url: str, params: Optional[Dict[str, Any]] = None) -> Any:
    full_url = url
    if params:
        full_url = f"{url}?{urlencode(params)}"
    request = Request(
        full_url,
        headers={
            "User-Agent": "MMSage-ValidationPlan/1.0",
            "Accept": "application/json",
        },
    )
    with urlopen(request, timeout=FULLTEXT_REQUEST_TIMEOUT) as response:
        return json.loads(response.read().decode("utf-8", errors="replace"))


def _fetch_bioc_document(identifier: str) -> Optional[dict]:
    if not identifier:
        return None
    encoded = quote(identifier, safe="")
    try:
        payload = _http_get_json(f"{PMC_BIOC_URL}/{encoded}/unicode")
    except (HTTPError, URLError, TimeoutError, json.JSONDecodeError):
        return None
    if isinstance(payload, dict):
        documents = payload.get("documents")
        if isinstance(documents, list) and documents:
            return documents[0]
    return None


def _is_methods_passage(passage: dict) -> bool:
    infons = passage.get("infons") or {}
    surface = " ".join(str(value) for value in infons.values() if value).lower()
    return any(term in surface for term in METHOD_SECTION_TERMS)


def _normalize_phrase(value: str) -> str:
    return re.sub(r"\s+", " ", str(value or "")).strip(" .;,:")


def _extract_condition_markers(text: str) -> List[str]:
    markers = []
    for pattern in (CONDITION_PATTERN, TIME_PATTERN, CELL_MODEL_PATTERN, ANIMAL_MODEL_PATTERN):
        markers.extend(_normalize_phrase(match.group(0)) for match in pattern.finditer(text or ""))
    return _dedupe_preserve(markers)


def _classify_method_snippet(text: str) -> str:
    lowered = str(text or "").lower()
    has_cell = bool(CELL_MODEL_PATTERN.search(lowered)) or _contains_any(lowered, CELL_TERMS)
    has_animal = bool(ANIMAL_MODEL_PATTERN.search(lowered)) or _contains_any(lowered, ANIMAL_TERMS)
    if has_cell and not has_animal:
        return "in_vitro"
    if has_animal:
        return "in_vivo"
    return "unspecified"


def _summarize_method_snippet(text: str) -> str:
    clean = _normalize_phrase(text)
    return _truncate(clean, 280)


def _extract_method_passages_from_bioc(document: dict) -> List[dict]:
    passages = []
    for passage in document.get("passages", []) if isinstance(document, dict) else []:
        if not isinstance(passage, dict) or not _is_methods_passage(passage):
            continue
        text = str(passage.get("text") or "").strip()
        if not text:
            continue
        markers = _extract_condition_markers(text)
        if not markers and len(text) < 80:
            continue
        passages.append(
            {
                "text": text,
                "section": str((passage.get("infons") or {}).get("section_type") or "").strip(),
                "markers": markers[:8],
                "model_type": _classify_method_snippet(text),
            }
        )
    return passages


def _build_method_evidence_item(article: dict, document: dict, passage: dict) -> dict:
    infons = document.get("infons") or {}
    pmid = str(infons.get("article-id_pmid") or article.get("pmid") or "").strip()
    pmcid = str(infons.get("article-id_pmc") or "").strip()
    title = str(article.get("title") or "").strip()
    citation = f"PMID {pmid}"
    if pmcid:
        citation = f"{citation}; PMCID {pmcid}"
    return {
        "pmid": pmid,
        "pmcid": pmcid,
        "title": title,
        "citation": citation,
        "model_type": passage.get("model_type") or "unspecified",
        "reported_conditions": _as_str_list(passage.get("markers"))[:8],
        "support_summary": f"{citation}: {_summarize_method_snippet(passage.get('text') or '')}",
    }


def collect_fulltext_method_evidence(
    articles: List[dict],
    bacteria: str,
    metabolite: str,
    disease: str,
) -> dict:
    del bacteria, metabolite, disease
    evidence = {"in_vitro": [], "in_vivo": [], "all": []}
    for article in articles[:FULLTEXT_MAX_ARTICLES]:
        pmid = str(article.get("pmid") or "").strip()
        if not pmid:
            continue
        document = _fetch_bioc_document(pmid)
        if not document:
            continue
        for passage in _extract_method_passages_from_bioc(document):
            item = _build_method_evidence_item(article, document, passage)
            evidence["all"].append(item)
            if item["model_type"] == "in_vitro":
                evidence["in_vitro"].append(item)
            elif item["model_type"] == "in_vivo":
                evidence["in_vivo"].append(item)

    for key in ("all", "in_vitro", "in_vivo"):
        deduped = []
        seen = set()
        for item in evidence[key]:
            marker_key = (item.get("citation"), item.get("support_summary"))
            if marker_key in seen:
                continue
            seen.add(marker_key)
            deduped.append(item)
        evidence[key] = deduped[:8]
    return evidence


def _format_method_evidence_context(items: List[dict], limit: int = 6) -> str:
    if not items:
        return "No open-access full-text methods evidence was retrieved."
    lines = []
    for item in items[:limit]:
        lines.append(
            f"- {item.get('citation')}: {item.get('title')}\n"
            f"  Conditions: {', '.join(item.get('reported_conditions') or ['not explicitly parsed'])}\n"
            f"  Methods summary: {_truncate(item.get('support_summary') or '', 180)}"
        )
    return "\n".join(lines)


def _method_evidence_strings(items: List[dict], field: str, limit: int = 6) -> List[str]:
    values = []
    for item in items[:limit]:
        value = item.get(field)
        if isinstance(value, list):
            values.extend(str(v).strip() for v in value if str(v).strip())
        elif str(value or "").strip():
            values.append(str(value).strip())
    return _dedupe_preserve(values)


def _backfill_method_fields(items: List[dict], method_evidence: List[dict]) -> List[dict]:
    support = _method_evidence_strings(method_evidence, "support_summary", limit=4)
    conditions = _method_evidence_strings(method_evidence, "reported_conditions", limit=8)
    citations = _method_evidence_strings(method_evidence, "citation", limit=6)
    for item in items:
        if not isinstance(item, dict):
            continue
        if not _as_str_list(item.get("fulltext_method_support")) and support:
            item["fulltext_method_support"] = support[:3]
        if not _as_str_list(item.get("reported_conditions")) and conditions:
            item["reported_conditions"] = conditions[:6]
        if not _as_str_list(item.get("source_citations")) and citations:
            item["source_citations"] = citations[:4]
    return items


def _format_protocol_evidence_context(items: List[dict], limit: int = 4) -> str:
    if not items:
        return "No protocols.io operational reference was retrieved."
    chunks = []
    for item in items[:limit]:
        chunks.append(
            f"- {item.get('citation')}: {item.get('title')}\n"
            f"  URL: {item.get('url') or 'not available'}\n"
            f"  Materials: {'; '.join((_dedupe_preserve(item.get('materials') or []))[:5]) or 'not listed'}\n"
            f"  Procedure: {'; '.join((_dedupe_preserve(item.get('steps') or []))[:6]) or 'not listed'}"
        )
    return "\n".join(chunks)


def _backfill_protocol_fields(items: List[dict], protocol_evidence: List[dict]) -> List[dict]:
    citations = _dedupe_preserve(
        [str(item.get("citation") or "").strip() for item in protocol_evidence if item.get("citation")]
    )
    urls = _dedupe_preserve(
        [str(item.get("url") or "").strip() for item in protocol_evidence if item.get("url")]
    )
    materials = _dedupe_preserve(
        [value for item in protocol_evidence for value in (item.get("materials") or [])]
    )
    steps = _dedupe_preserve(
        [value for item in protocol_evidence for value in (item.get("steps") or [])]
    )
    for item in items:
        if not isinstance(item, dict):
            continue
        if not _as_str_list(item.get("protocols_io_support")) and steps:
            item["protocols_io_support"] = steps[:8]
        if not _as_str_list(item.get("protocols_io_materials")) and materials:
            item["protocols_io_materials"] = materials[:8]
        if not _as_str_list(item.get("protocols_io_citations")) and citations:
            item["protocols_io_citations"] = citations[:4]
        if not _as_str_list(item.get("protocols_io_urls")) and urls:
            item["protocols_io_urls"] = urls[:4]
    return items


def _preserve_experiment_source_fields(revised: Any, original: List[dict]) -> List[dict]:
    if not isinstance(revised, list):
        return original

    revised_items = [dict(item) for item in revised if isinstance(item, dict)]
    used_revised: set[int] = set()
    derived_fields = {
        "activation_gate",
        "activation_gate_reasons",
        "activation_gate_status",
        "completion_issues",
        "design_issues",
        "design_status",
        "execution_status",
        "generation_status",
        "group_definitions",
        "sampling",
        "endpoints",
        "analysis_plan",
        "decision",
        "parameters",
        "operationalization_issues",
        "prerequisite_module_ids",
        "unlock_rule",
    }

    def _match_revised(source: dict, source_index: int) -> Optional[int]:
        module_id = str(source.get("module_id") or "").strip()
        if module_id:
            for revised_index, candidate in enumerate(revised_items):
                if revised_index in used_revised:
                    continue
                if str(candidate.get("module_id") or "").strip() == module_id:
                    return revised_index
        role = str(source.get("experiment_role") or "").strip()
        if role:
            for revised_index, candidate in enumerate(revised_items):
                if revised_index in used_revised:
                    continue
                candidate_id = str(candidate.get("module_id") or "").strip()
                if module_id and candidate_id:
                    continue
                if str(candidate.get("experiment_role") or "").strip() == role:
                    return revised_index
        if source_index < len(revised_items) and source_index not in used_revised:
            candidate_id = str(revised_items[source_index].get("module_id") or "").strip()
            if not module_id or not candidate_id:
                return source_index
        return None

    def _merge_nonempty(source: dict, candidate: dict) -> dict:
        merged = {key: value for key, value in source.items() if key not in derived_fields}
        for key, value in candidate.items():
            if key in derived_fields:
                continue
            if isinstance(value, str):
                if value.strip():
                    merged[key] = value
            elif isinstance(value, (list, dict)):
                if value:
                    merged[key] = value
            elif value is not None:
                merged[key] = value
        return merged

    merged_items: List[dict] = []
    for source_index, source in enumerate(original):
        if not isinstance(source, dict):
            continue
        revised_index = _match_revised(source, source_index)
        if revised_index is None:
            merged_items.append(_merge_nonempty(source, {}))
            continue
        used_revised.add(revised_index)
        merged_items.append(_merge_nonempty(source, revised_items[revised_index]))

    for revised_index, candidate in enumerate(revised_items):
        if revised_index not in used_revised:
            merged_items.append(candidate)
    return merged_items


def _normalize_audit_claims(items: Any, include_reason: bool = False) -> List[dict]:
    normalized = []
    if not isinstance(items, list):
        return normalized
    for item in items[:12]:
        if not isinstance(item, dict):
            continue
        normalized_item = {
            "claim": str(item.get("claim") or "").strip(),
            "support_level": str(item.get("support_level") or "speculative").strip().lower(),
            "pmids": _as_str_list(item.get("pmids")),
        }
        if include_reason:
            normalized_item["reason"] = str(item.get("reason") or "").strip()
        normalized.append(normalized_item)
    return normalized


def _normalize_evidence_items(items: Any) -> List[dict]:
    normalized = []
    if not isinstance(items, list):
        return normalized
    for item in items[:14]:
        if not isinstance(item, dict):
            continue
        normalized.append(
            {
                "claim_type": str(item.get("claim_type") or "unspecified").strip(),
                "claim": str(item.get("claim") or "").strip(),
                "support_level": str(item.get("support_level") or "speculative").strip().lower(),
                "pmids": _as_str_list(item.get("pmids")),
                "model_type": str(item.get("model_type") or "unspecified").strip(),
                "evidence_summary": str(item.get("evidence_summary") or "").strip(),
            }
        )
    return normalized


def _normalize_production_evidence_item(item: Any, direct: bool = False) -> dict:
    if not isinstance(item, dict):
        return {}
    pmid = str(item.get("pmid") or "").strip()
    evidence_type = str(item.get("evidence_type") or ("direct_monoculture_production" if direct else "uncertain")).strip().lower()
    allowed_evidence_types = {
        "direct_monoculture_production",
        "direct_monoculture_nonproduction",
        "candidate_microbe_disease_intervention",
        "candidate_metabolite_disease_intervention",
        "candidate_pair_association",
        "indirect_ecological_evidence",
        "analogous_method_only",
        "background_only",
        "unrelated",
        "uncertain",
    }
    if evidence_type not in allowed_evidence_types:
        evidence_type = "uncertain"
    candidate_relevance = str(item.get("candidate_relevance") or "uncertain").strip().lower()
    if candidate_relevance not in {"direct", "partial", "analogous", "unrelated", "uncertain"}:
        candidate_relevance = "uncertain"
    normalized = {
        "pmid": pmid,
        "title": str(item.get("title") or "").strip(),
        "claim": str(item.get("claim") or item.get("finding") or "").strip(),
        "evidence_type": evidence_type,
        "candidate_relevance": candidate_relevance,
        "citation_eligible": bool(item.get("citation_eligible", False)),
        "claim_scope": str(item.get("claim_scope") or "").strip(),
        "model_system": str(item.get("model_system") or item.get("culture_system") or "").strip(),
        "measured_output": str(item.get("measured_output") or "").strip(),
        "why_this_support_level": str(
            item.get("why_this_support_level") or item.get("why_direct") or item.get("reason") or ""
        ).strip(),
    }
    return normalized if pmid or normalized["claim"] else {}


def _normalize_production_evidence_assessment(value: Any) -> dict:
    if not isinstance(value, dict):
        value = {}
    direct_evidence = [
        normalized
        for normalized in (
            _normalize_production_evidence_item(item, direct=True)
            for item in (value.get("direct_evidence") or [])[:8]
        )
        if normalized
        and normalized.get("pmid")
        and normalized.get("evidence_type") == "direct_monoculture_production"
        and normalized.get("candidate_relevance") == "direct"
    ]
    paper_findings = [
        normalized
        for normalized in (
            _normalize_production_evidence_item(item)
            for item in (value.get("paper_findings") or [])[:16]
        )
        if normalized and normalized.get("pmid")
    ]
    status = str(value.get("status") or "not_assessed").strip().lower()
    allowed_statuses = {"direct_supported", "indirect_only", "not_found", "conflicting", "not_assessed"}
    if status not in allowed_statuses:
        status = "not_assessed"
    if status == "direct_supported" and not direct_evidence:
        status = "indirect_only" if paper_findings else "not_found"
    return {
        "status": status,
        "conclusion": str(value.get("conclusion") or "").strip(),
        "direct_evidence": direct_evidence,
        "paper_findings": paper_findings,
        "evidence_limitations": _as_str_list(value.get("evidence_limitations")),
    }


def _paper_reported_model_scope(text: str) -> str:
    """Describe only the model explicitly visible in the retrieved paper text."""
    lowered = str(text or "").lower()
    dss = bool(re.search(r"\b(?:dss|dextran sodium sulfate)\b", lowered))
    mouse = bool(re.search(r"\b(?:mouse|mice|murine)\b", lowered))
    if dss and mouse:
        return "a DSS-induced mouse colitis model"
    if dss:
        return "a DSS-induced colitis model"
    if mouse and any(term in lowered for term in ("colitis", "inflammatory bowel", "ibd")):
        return "a paper-reported mouse colitis model"
    if any(term in lowered for term in ("organoid", "epithelial cell", "macrophage", "cell line", "in vitro")):
        return "the paper-reported in vitro host model"
    if any(term in lowered for term in ("patients", "participants", "cohort", "clinical")):
        return "the paper-reported human population"
    if mouse:
        return "the paper-reported mouse model"
    if any(term in lowered for term in ("rat", "animal", "in vivo")):
        return "the paper-reported preclinical model"
    return "the paper-reported experimental model"


def _ground_production_assessment_to_articles(
    assessment: dict,
    articles: List[dict],
    bacteria: str,
    metabolite: str,
    disease: str,
    fulltext_method_evidence: Optional[Dict[str, List[dict]]] = None,
) -> dict:
    grounded = _normalize_production_evidence_assessment(assessment)
    article_map = {
        str(article.get("pmid") or "").strip(): article
        for article in articles
        if str(article.get("pmid") or "").strip()
    }
    bacteria_aliases = _citation_entity_aliases(bacteria, is_microbe=True)
    metabolite_aliases = _citation_entity_aliases(metabolite)
    disease_aliases = _citation_disease_aliases(disease)
    fulltext_by_pmid: Dict[str, List[str]] = {}
    for method_item in ((fulltext_method_evidence or {}).get("all") or []):
        if not isinstance(method_item, dict):
            continue
        pmid = str(method_item.get("pmid") or "").strip()
        if not pmid:
            continue
        method_text = " ".join(
            part
            for part in [
                str(method_item.get("support_summary") or "").strip(),
                " ".join(_as_str_list(method_item.get("reported_conditions"))),
            ]
            if part
        )
        if method_text:
            fulltext_by_pmid.setdefault(pmid, []).append(method_text)

    raw_items_by_pmid: Dict[str, dict] = {}
    for source_item in (grounded.get("direct_evidence") or []) + (grounded.get("paper_findings") or []):
        if not isinstance(source_item, dict):
            continue
        pmid = str(source_item.get("pmid") or "").strip()
        if pmid and pmid in article_map and pmid not in raw_items_by_pmid:
            raw_items_by_pmid[pmid] = dict(source_item)

    def grounded_claim(evidence_type: str, paper_scope: str) -> tuple[str, str, str]:
        if evidence_type == "direct_monoculture_production":
            return (
                f"This paper reports candidate-specific controlled-culture evidence consistent with direct {metabolite} production by {bacteria}.",
                "Limited to the exact strain, medium, culture system, analytical assay, and conditions reported in this paper; disease mediation is not established.",
                "The retrieved title/abstract and any matched full-text method excerpt contain the named candidate pair, strict single-strain culture language, production language, and direct metabolite measurement language.",
            )
        if evidence_type == "direct_monoculture_nonproduction":
            return (
                f"This paper directly tested controlled-culture {metabolite} production by {bacteria} but reported that production was not detected or not supported under the tested conditions.",
                "This is strong candidate-relevant negative evidence limited to the reported strain, medium, analytical sensitivity, time points, and culture conditions; it does not prove production is impossible under every condition.",
                "One retrieved title/abstract or one self-contained full-text passage contains the named candidate pair, controlled single-strain culture, direct measurement, and explicit nonproduction language.",
            )
        if evidence_type == "candidate_microbe_disease_intervention":
            return (
                f"This paper directly evaluates {bacteria} as an intervention in {paper_scope}.",
                f"Limited to {paper_scope}; it does not establish that {bacteria} produces {metabolite} or that {metabolite} mediates the microbial effect.",
                "The named microbe, disease context, and intervention language occur in the retrieved paper.",
            )
        if evidence_type == "candidate_metabolite_disease_intervention":
            return (
                f"This paper directly evaluates {metabolite} as an intervention in {paper_scope}.",
                f"Limited to {paper_scope}; it does not establish {bacteria} as the direct or indirect source of {metabolite}, or {metabolite} as a mediator of {bacteria}'s effect.",
                "The named metabolite, disease context, and intervention language occur in the retrieved paper.",
            )
        if evidence_type == "candidate_pair_association":
            return (
                f"This paper reports {bacteria} and {metabolite} in the same study, supporting a candidate-pair association only.",
                "Does not establish direct production, temporal order, ecological mediation, or a disease-causal pathway.",
                "Both named candidate entities occur in the retrieved paper, but the strict direct-culture screen is not met.",
            )
        if evidence_type == "indirect_ecological_evidence":
            return (
                f"This paper provides candidate-linked ecological evidence involving {bacteria} and {metabolite}.",
                "Supports an indirect ecological or cross-feeding hypothesis only; it does not establish direct production by the named microbe.",
                "Both candidate entities and ecological interaction language occur in the retrieved paper without strict direct-culture evidence.",
            )
        if evidence_type == "analogous_method_only":
            return (
                "This paper provides an analogous culture or metabolite-measurement method only.",
                "Operational analogy only; it is not candidate-specific biological evidence.",
                "A relevant method context is present, but the named candidate pair is not studied together.",
            )
        if evidence_type == "background_only":
            return (
                "This paper provides background evidence for only part of the candidate chain.",
                "Background only; it cannot be used as direct evidence for the complete candidate relationship.",
                "Only a subset of the named candidate entities or contexts occurs in the retrieved paper.",
            )
        return (
            "This paper does not provide candidate-specific evidence for the proposed relationship.",
            "Unrelated to the candidate chain for evidence-adjudication purposes.",
            "The retrieved title and abstract do not contain the candidate evidence pattern required for a stronger category.",
        )

    def ground_item(item: dict) -> Optional[dict]:
        pmid = str(item.get("pmid") or "").strip()
        article = article_map.get(pmid)
        if not article:
            return None
        abstract_text = " ".join(
            part
            for part in [
                str(article.get("title") or "").strip(),
                str(article.get("abstract") or "").strip(),
            ]
            if part
        )
        candidates = [abstract_text] + fulltext_by_pmid.get(pmid, [])
        classified = [
            _classify_candidate_article_evidence(
                passage,
                bacteria_aliases,
                metabolite_aliases,
                disease_aliases,
            )
            for passage in candidates
            if str(passage or "").strip()
        ]
        classification_priority = {
            "direct_monoculture_production": 0,
            "direct_monoculture_nonproduction": 1,
            "candidate_microbe_disease_intervention": 2,
            "candidate_metabolite_disease_intervention": 3,
            "indirect_ecological_evidence": 4,
            "candidate_pair_association": 5,
            "analogous_method_only": 6,
            "background_only": 7,
            "unrelated": 8,
        }
        evidence_type, relevance, citation_eligible = min(
            classified or [("unrelated", "unrelated", False)],
            key=lambda value: classification_priority.get(value[0], 99),
        )

        paper_scope = _paper_reported_model_scope(abstract_text)
        claim, claim_scope, reason = grounded_claim(evidence_type, paper_scope)
        item["evidence_type"] = evidence_type
        item["candidate_relevance"] = relevance
        item["citation_eligible"] = citation_eligible
        item["title"] = str(article.get("title") or "").strip()
        item["claim"] = claim
        item["claim_scope"] = claim_scope
        item["why_this_support_level"] = reason
        if evidence_type in {"direct_monoculture_production", "direct_monoculture_nonproduction"}:
            item["model_system"] = "Candidate-specific controlled single-strain culture reported in the retrieved paper"
            item["measured_output"] = (
                f"Direct measurement of {metabolite} reported in the retrieved paper"
                if evidence_type == "direct_monoculture_production"
                else f"Direct measurement reported no supported {metabolite} production under the tested conditions"
            )
        elif evidence_type in {
            "candidate_microbe_disease_intervention",
            "candidate_metabolite_disease_intervention",
        }:
            item["model_system"] = paper_scope
            item["measured_output"] = "Paper-reported disease or host-response outcome"
        else:
            item["model_system"] = ""
            item["measured_output"] = ""
        return item

    paper_items: List[dict] = []
    for pmid, article in article_map.items():
        seed = dict(raw_items_by_pmid.get(pmid) or {"pmid": pmid})
        grounded_item = ground_item(seed)
        if grounded_item is not None:
            paper_items.append(grounded_item)
    evidence_priority = {
        "direct_monoculture_production": 0,
        "direct_monoculture_nonproduction": 1,
        "candidate_microbe_disease_intervention": 2,
        "candidate_metabolite_disease_intervention": 3,
        "indirect_ecological_evidence": 4,
        "candidate_pair_association": 5,
        "analogous_method_only": 6,
        "background_only": 7,
        "unrelated": 8,
    }
    paper_items.sort(key=lambda item: evidence_priority.get(str(item.get("evidence_type")), 99))
    strong_paper_items = [item for item in paper_items if bool(item.get("citation_eligible"))]
    excluded_count = len(paper_items) - len(strong_paper_items)
    direct_items = [
        dict(item)
        for item in strong_paper_items
        if item.get("evidence_type") == "direct_monoculture_production"
    ]
    direct_negative_items = [
        dict(item)
        for item in strong_paper_items
        if item.get("evidence_type") == "direct_monoculture_nonproduction"
    ]
    candidate_non_direct = {
        "candidate_microbe_disease_intervention",
        "candidate_metabolite_disease_intervention",
        "indirect_ecological_evidence",
    }
    if direct_items and direct_negative_items:
        status = "conflicting"
    elif direct_items:
        status = "conflicting" if grounded.get("status") == "conflicting" else "direct_supported"
    elif any(item.get("evidence_type") in candidate_non_direct for item in strong_paper_items):
        status = "indirect_only"
    else:
        status = "not_found"

    if status == "direct_supported":
        conclusion = (
            f"At least one retrieved PMID met the conservative same-paper screen for candidate-specific controlled-culture measurement of {metabolite} by {bacteria}; interpretation remains limited to the reported strain and conditions."
        )
    elif status == "conflicting":
        conclusion = (
            f"The retrieved set contains candidate-specific direct-culture evidence but also conflicting interpretation; verify the full methods and results paper by paper before treating direct production as established."
        )
    elif status == "indirect_only":
        conclusion = (
            f"No retrieved PMID met the strict same-paper direct-production threshold for {bacteria} and {metabolite}; the retrieved candidate-specific evidence is non-direct only."
        )
    elif direct_negative_items:
        conclusion = (
            f"At least one retrieved PMID directly tested controlled-culture {metabolite} production by {bacteria} but did not support production under the reported conditions; direct production remains unconfirmed rather than impossible."
        )
    else:
        conclusion = (
            f"No retrieved PMID met the strict same-paper direct-production threshold for {bacteria} and {metabolite} in the current search."
        )
    limitations = list(grounded.get("evidence_limitations") or [])
    if excluded_count:
        limitations.append(
            f"{excluded_count} weak co-occurrence, analogous-method, background, or unrelated search hits were excluded from citations because they did not directly test a candidate experimental link."
        )
    return {
        "status": status,
        "conclusion": conclusion,
        "direct_evidence": direct_items,
        "paper_findings": strong_paper_items,
        "evidence_limitations": _dedupe_preserve(limitations),
    }


def _normalize_hypothesis_branches(items: Any) -> List[dict]:
    normalized: List[dict] = []
    if not isinstance(items, list):
        return normalized
    for index, item in enumerate(items[:6], start=1):
        if not isinstance(item, dict):
            continue
        normalized.append(
            {
                "hypothesis_id": str(item.get("hypothesis_id") or f"H{index}").strip(),
                "statement": str(item.get("statement") or "").strip(),
                "current_evidence_status": str(item.get("current_evidence_status") or "unresolved").strip().lower(),
                "evidence_basis": _as_str_list(item.get("evidence_basis")),
                "discriminating_prediction": str(item.get("discriminating_prediction") or "").strip(),
                "in_vitro_gate": str(item.get("in_vitro_gate") or "").strip(),
                "animal_gate": str(item.get("animal_gate") or "").strip(),
                "human_gate": str(item.get("human_gate") or "").strip(),
                "falsification_or_redirection": str(item.get("falsification_or_redirection") or "").strip(),
            }
        )
    return normalized


def _default_candidate_hypothesis_branches(
    bacteria: str,
    metabolite: str,
    disease: str,
    production_status: str = "not_assessed",
) -> List[dict]:
    direct_status = "supported_by_current_culture_evidence" if production_status == "direct_supported" else "unresolved"
    return _normalize_hypothesis_branches(
        [
            {
                "hypothesis_id": "H1",
                "statement": f"{bacteria} directly produces {metabolite} under controlled culture conditions.",
                "current_evidence_status": direct_status,
                "evidence_basis": [],
                "discriminating_prediction": f"Live {bacteria} will generate a time-dependent increase in newly formed {metabolite} relative to uninoculated and inactivated controls.",
                "in_vitro_gate": "Controlled culture with direct metabolite quantification and, when appropriate, stable-isotope tracing.",
                "animal_gate": "Only after the metabolic function is established, compare function-intact and function-deficient microbial interventions with metabolite rescue at comparable exposure or colonization.",
                "human_gate": "Only after branch-specific animal causality is reproducible, test whether microbial and metabolite changes are directionally coupled to disease activity in humans.",
                "falsification_or_redirection": "If newly formed metabolite is not detected, stop claiming direct production and move to the indirect ecological branch.",
            },
            {
                "hypothesis_id": "H2",
                "statement": f"{bacteria} indirectly increases {metabolite} by changing substrates, conditioned media, cross-feeding partners, or community metabolism.",
                "current_evidence_status": "plausible_but_unresolved",
                "evidence_basis": [],
                "discriminating_prediction": f"{metabolite} will increase only in conditioned-medium, co-culture, or defined-community settings and will track a partner organism or substrate-transfer process.",
                "in_vitro_gate": "Conditioned-medium, co-culture, or defined-community experiment with substrate tracing and partner controls.",
                "animal_gate": "Use a defined community or co-colonization design that measures substrate transfer, candidate producer activity, metabolite output, and disease phenotype.",
                "human_gate": "If the ecological branch passes in animals, include community composition or pathway covariates rather than modeling the microbe as the direct producer.",
                "falsification_or_redirection": "If neither controlled culture nor ecological systems increase the metabolite, treat the original pair as association-only and test independent disease effects.",
            },
            {
                "hypothesis_id": "H3",
                "statement": f"{bacteria} and {metabolite} independently or interactively influence {disease}, but {metabolite} is not the sole mediator of the microbial effect.",
                "current_evidence_status": "plausible_alternative",
                "evidence_basis": [],
                "discriminating_prediction": "Microbial protection will persist partly or fully when the candidate metabolite pathway is absent, or inactivated microbial material will retain an effect inconsistent with a purely metabolic mechanism.",
                "in_vitro_gate": "Compare live microbe, inactivated microbe or structural material when appropriate, metabolite alone, and combined exposure in a host-cell model.",
                "animal_gate": "Separate an effect-and-interaction study from any later mediation study and include controls for non-metabolic microbial effects.",
                "human_gate": "Model the microbe and metabolite as separate exposures plus interaction unless preclinical mediation criteria were met.",
                "falsification_or_redirection": "If microbial effects disappear only when metabolite production is lost and return with metabolite rescue, downgrade this alternative in favor of mediation.",
            },
        ]
    )




def _proposed_experiment_source_fields(note: str) -> dict:
    return {
        "evidence_basis": [note],
        "query_round_support": [],
        "fulltext_method_support": [],
        "reported_conditions": [note],
        "source_citations": [],
        "protocols_io_support": [],
        "protocols_io_materials": [],
        "protocols_io_citations": [],
        "protocols_io_urls": [],
    }


def _generic_candidate_in_vitro_plan(bacteria: str, metabolite: str, disease: str) -> List[dict]:
    """Compatibility hook; fixed candidate module templates are intentionally disabled."""
    return []

def _generic_candidate_in_vivo_plan(bacteria: str, metabolite: str, disease: str) -> List[dict]:
    """Compatibility hook; fixed candidate module templates are intentionally disabled."""
    return []


def _candidate_role_semantic_text(item: dict) -> str:
    """Collect decision-level content used to verify a model-supplied role."""
    values: List[str] = []
    for field in (
        "scientific_question",
        "why_needed",
        "study_object",
        "primary_indicator",
        "positive_gate",
        "aim",
        "biological_question",
        "model_system",
        "model",
        "experimental_material",
        "design",
        "group_logic",
        "primary_endpoint",
    ):
        value = item.get(field)
        if value:
            values.append(str(value))
    for field in (
        "hypothesis_ids",
        "hypothesis_tested",
        "secondary_indicators",
        "key_controls",
        "controls",
        "readouts",
        "primary_endpoints",
        "secondary_endpoints",
        "mechanistic_endpoints",
    ):
        values.extend(_as_str_list(item.get(field)))
    for group in item.get("groups") or item.get("group_definitions") or []:
        if isinstance(group, dict):
            values.extend(
                str(group.get(field) or "")
                for field in ("group_name", "exposure_or_condition", "control_purpose")
            )
        elif group:
            values.append(str(group))
    return " ".join(values).lower()


def _infer_candidate_experiment_role(item: dict, is_in_vivo: bool = False) -> str:
    declared_role = _canonicalize_experiment_role(item.get("experiment_role"))
    valid = {
        "effect_interaction_gate", "direct_causal_rescue", "indirect_ecology_causal"
    } if is_in_vivo else {
        "direct_production_gate", "indirect_ecology_gate", "host_response_gate"
    }
    text = _candidate_role_semantic_text(item)
    if is_in_vivo:
        if any(
            term in text
            for term in (
                "defined community",
                "gnotobiotic",
                "validated producer",
                "producer-depleted",
                "producer depleted",
                "community necessity",
                "cross-feeding",
                "cross feeding",
            )
        ):
            return "indirect_ecology_causal"
        if any(
            term in text
            for term in (
                "function-deficient",
                "function deficient",
                "function-loss",
                "function loss",
                "production-deficient",
                "production deficient",
                "knockout",
                "genetic complementation",
                "source-function necessity",
                "source function necessity",
            )
        ):
            return "direct_causal_rescue"
        if declared_role in valid:
            return declared_role
        return "effect_interaction_gate"
    if any(
        term in text
        for term in (
            "host cell",
            "host-cell",
            "organoid",
            "epithelial",
            "macrophage",
            "barrier integrity",
            "host response",
            "host-response",
            "cell viability",
            "cytotoxicity",
        )
    ):
        return "host_response_gate"
    if any(
        term in text
        for term in (
            "indirectly increase",
            "indirect ecological",
            "conditioned medium",
            "cross-feeding",
            "cross feeding",
            "defined community",
            "co-culture",
            "coculture",
            "candidate-present",
            "candidate-absent",
            "candidate present",
            "candidate absent",
        )
    ):
        return "indirect_ecology_gate"
    if any(
        term in text
        for term in (
            "direct production",
            "directly produce",
            "monoculture",
            "mono-culture",
            "pure culture",
            "uninoculated medium",
            "culture supernatant",
        )
    ):
        return "direct_production_gate"
    if declared_role in valid:
        return declared_role
    return "direct_production_gate"


def _merge_complete_experiment_item(default: dict, generated: Optional[dict], is_in_vivo: bool = False) -> dict:
    merged = dict(default)
    if not isinstance(generated, dict):
        return merged
    list_fields = {
        "prerequisite_module_ids", "key_materials_equipment", "procedure_steps", "controls", "groups", "readouts",
        "primary_endpoints", "secondary_endpoints", "mechanistic_endpoints", "key_confounders",
        "query_round_support", "fulltext_method_support", "reported_conditions",
        "protocols_io_support", "protocols_io_materials", "protocols_io_citations", "protocols_io_urls",
    }
    enrichment_string_fields = {
        "module_id", "execution_status", "experimental_unit", "replication_and_sampling", "sample_size_basis",
        "randomization_and_blinding", "safety_and_stopping_rules",
        "model_system", "model", "why_this_model", "experimental_material", "design", "group_logic",
        "intervention", "intervention_route", "dose_timing_logic", "timeline", "primary_endpoint",
        "data_analysis",
    }
    for field in list_fields:
        merged[field] = _dedupe_preserve(_as_str_list(default.get(field)) + _as_str_list(generated.get(field)))
    for field in enrichment_string_fields:
        candidate_value = str(generated.get(field) or "").strip()
        if candidate_value:
            merged[field] = candidate_value
    if generated.get("parameter_provenance"):
        merged["parameter_provenance"] = _normalize_parameter_provenance(generated.get("parameter_provenance"))
    merged["experiment_role"] = str(default.get("experiment_role") or "")
    return merged


def _complete_candidate_experiment_plan(
    generated_items: List[dict],
    default_items: List[dict],
    is_in_vivo: bool = False,
) -> List[dict]:
    normalized_generated = _normalize_experiment_items(generated_items, is_in_vivo=is_in_vivo)
    by_role: Dict[str, dict] = {}
    for item in normalized_generated:
        role = _infer_candidate_experiment_role(item, is_in_vivo=is_in_vivo)
        by_role.setdefault(role, item)
    completed = [
        _merge_complete_experiment_item(default, by_role.get(str(default.get("experiment_role"))), is_in_vivo=is_in_vivo)
        for default in default_items
    ]
    return _normalize_experiment_items(completed, is_in_vivo=is_in_vivo)


CANDIDATE_ROLE_HYPOTHESIS_IDS = {
    "direct_production_gate": ["H1"],
    "indirect_ecology_gate": ["H2"],
    "host_response_gate": ["H3"],
    "effect_interaction_gate": ["H3"],
    "direct_causal_rescue": ["H1"],
    "indirect_ecology_causal": ["H2"],
}


def _candidate_module_completeness_score(item: dict) -> int:
    """Prefer the most decision-complete draft when a model repeats one role."""
    scalar_fields = (
        "scientific_question",
        "why_needed",
        "study_object",
        "samples_and_timing",
        "primary_indicator",
        "positive_gate",
        "branch_if_positive",
        "branch_if_negative",
        "claim_boundary",
    )
    list_fields = ("secondary_indicators", "key_controls")
    score = sum(3 for field in scalar_fields if str(item.get(field) or "").strip())
    score += sum(2 for field in list_fields if _as_str_list(item.get(field)))
    groups = item.get("groups") if isinstance(item.get("groups"), list) else []
    score += min(len(groups), 6)
    score += sum(
        1
        for group in groups
        if isinstance(group, dict)
        and str(group.get("group_name") or "").strip()
        and str(group.get("exposure_or_condition") or "").strip()
        and str(group.get("control_purpose") or "").strip()
    )
    return score


def _dedupe_candidate_modules_by_role(items: List[dict]) -> List[dict]:
    """Expose one canonical module per candidate role, never V1.2-style aliases."""
    role_order: List[str] = []
    selected: Dict[str, tuple[int, dict]] = {}
    for item in items:
        role = _canonicalize_experiment_role(item.get("experiment_role"))
        score = _candidate_module_completeness_score(item)
        if role not in selected:
            role_order.append(role)
            selected[role] = (score, item)
            continue
        if score > selected[role][0]:
            selected[role] = (score, item)
    return [selected[role][1] for role in role_order]


def _assign_candidate_experiment_roles(
    generated_items: List[dict],
    is_in_vivo: bool = False,
) -> List[dict]:
    normalized = _normalize_experiment_items(generated_items, is_in_vivo=is_in_vivo)
    role_order = (
        ("effect_interaction_gate", "direct_causal_rescue", "indirect_ecology_causal")
        if is_in_vivo
        else ("direct_production_gate", "indirect_ecology_gate", "host_response_gate")
    )
    valid_roles = set(role_order)
    for item in normalized:
        inferred_role = _infer_candidate_experiment_role(item, is_in_vivo=is_in_vivo)
        if inferred_role in valid_roles:
            item["experiment_role"] = inferred_role
            item["hypothesis_ids"] = list(CANDIDATE_ROLE_HYPOTHESIS_IDS[inferred_role])
    deduped = _dedupe_candidate_modules_by_role(normalized)
    order_index = {role: index for index, role in enumerate(role_order)}
    return sorted(
        deduped,
        key=lambda item: order_index.get(
            _canonicalize_experiment_role(item.get("experiment_role")), len(role_order)
        ),
    )


def _strong_candidate_evidence_items(assessment: Any, evidence_types: Optional[set[str]] = None) -> List[dict]:
    normalized = _normalize_production_evidence_assessment(assessment)
    items = (normalized.get("direct_evidence") or []) + (normalized.get("paper_findings") or [])
    selected: List[dict] = []
    seen_pmids = set()
    for item in items:
        if not isinstance(item, dict) or not bool(item.get("citation_eligible")):
            continue
        evidence_type = str(item.get("evidence_type") or "")
        if evidence_type not in STRONG_CITATION_EVIDENCE_TYPES:
            continue
        if evidence_types is not None and evidence_type not in evidence_types:
            continue
        pmid = str(item.get("pmid") or "").strip()
        if not pmid or pmid in seen_pmids:
            continue
        seen_pmids.add(pmid)
        selected.append(item)
    return selected


def _attach_strong_evidence_to_plan(plan: dict) -> dict:
    assessment = plan.get("direct_production_evidence_assessment") or {}
    role_types = {
        "direct_production_gate": {"direct_monoculture_production", "direct_monoculture_nonproduction"},
        "indirect_ecology_gate": {"indirect_ecological_evidence"},
        "host_response_gate": {"candidate_metabolite_disease_intervention"},
        "effect_interaction_gate": {"candidate_microbe_disease_intervention", "candidate_metabolite_disease_intervention"},
        "direct_causal_rescue": {"direct_monoculture_production", "candidate_microbe_disease_intervention", "candidate_metabolite_disease_intervention"},
        "indirect_ecology_causal": {"indirect_ecological_evidence"},
    }
    for key in ("in_vitro_plan", "in_vivo_plan"):
        for item in plan.get(key) or []:
            if not isinstance(item, dict):
                continue
            evidence = _strong_candidate_evidence_items(
                assessment,
                role_types.get(str(item.get("experiment_role") or ""), set()),
            )
            item["source_citations"] = []
            if not evidence:
                item["evidence_basis"] = [
                    "No strongly candidate-relevant paper was attached to this experimental role; all operational values remain proposed and require candidate-specific optimization."
                ]
                item["fulltext_method_support"] = []
                continue
            citations = [f"PMID {entry.get('pmid')}" for entry in evidence]
            basis = [
                f"PMID {entry.get('pmid')} ({entry.get('title') or 'title unavailable'}): {entry.get('claim')}"
                for entry in evidence
            ]
            item["source_citations"] = _dedupe_preserve(citations)
            item["evidence_basis"] = _dedupe_preserve(basis)
            evidence_pmids = {str(entry.get("pmid") or "").strip() for entry in evidence}
            method_bucket = (
                (plan.get("fulltext_method_evidence") or {}).get("in_vivo", [])
                if key == "in_vivo_plan"
                else (plan.get("fulltext_method_evidence") or {}).get("in_vitro", [])
            )
            matched_methods = [
                method
                for method in method_bucket
                if isinstance(method, dict) and str(method.get("pmid") or "").strip() in evidence_pmids
            ]
            item["fulltext_method_support"] = _dedupe_preserve(
                [str(method.get("support_summary") or "").strip() for method in matched_methods]
            )
            item["reported_conditions"] = _dedupe_preserve(
                [
                    condition
                    for method in matched_methods
                    for condition in _as_str_list(method.get("reported_conditions"))
                ]
            ) or _as_str_list(item.get("reported_conditions"))

    branch_types = {
        "H1": {"direct_monoculture_production", "direct_monoculture_nonproduction"},
        "H2": {"indirect_ecological_evidence"},
        "H3": {"candidate_microbe_disease_intervention", "candidate_metabolite_disease_intervention"},
    }
    for branch in plan.get("hypothesis_branches") or []:
        if not isinstance(branch, dict):
            continue
        evidence = _strong_candidate_evidence_items(
            assessment,
            branch_types.get(str(branch.get("hypothesis_id") or ""), set()),
        )
        branch["evidence_basis"] = _dedupe_preserve(
            _as_str_list(branch.get("evidence_basis")) + [f"PMID {entry.get('pmid')}" for entry in evidence]
        )
    return plan


def _sanitize_nested_pmid_references(value: Any, allowed_pmids: set[str]) -> Any:
    if isinstance(value, dict):
        cleaned: Dict[str, Any] = {}
        for key, item in value.items():
            if key in {"pmids", "representative_pmids"} and isinstance(item, list):
                cleaned[key] = [pmid for pmid in _as_str_list(item) if pmid in allowed_pmids]
            else:
                cleaned[key] = _sanitize_nested_pmid_references(item, allowed_pmids)
        return cleaned
    if isinstance(value, list):
        return [_sanitize_nested_pmid_references(item, allowed_pmids) for item in value]
    if isinstance(value, str):
        def replace(match: re.Match[str]) -> str:
            return match.group(0) if match.group(1) in allowed_pmids else ""

        cleaned = re.sub(r"\bPMID\s*:?\s*(\d+)\b", replace, value, flags=re.IGNORECASE)
        cleaned = re.sub(r"\[\s*\]", "", cleaned)
        return re.sub(r"\s+", " ", cleaned).strip(" ;,|")
    return value


def _sanitize_candidate_plan_references(plan: dict) -> dict:
    """Keep candidate-plan PubMed references on one grounded strong-evidence whitelist."""
    assessment = _normalize_production_evidence_assessment(
        plan.get("direct_production_evidence_assessment")
    )
    strong_items = _strong_candidate_evidence_items(assessment)
    allowed_pmids = {
        str(item.get("pmid") or "").strip()
        for item in strong_items
        if str(item.get("pmid") or "").strip()
    }

    plan["evidence_basis"] = [
        {
            "claim_type": str(item.get("evidence_type") or "candidate_specific_evidence"),
            "claim": str(item.get("claim") or "").strip(),
            "support_level": "direct",
            "pmids": [str(item.get("pmid"))],
            "model_type": str(item.get("model_system") or "candidate-specific experimental evidence"),
            "evidence_summary": str(item.get("claim_scope") or "").strip(),
        }
        for item in strong_items
    ]

    audit = plan.get("protocol_audit") if isinstance(plan.get("protocol_audit"), dict) else {}
    audit["supported_claims"] = [
        {
            "claim": str(item.get("claim") or "").strip(),
            "support_level": "direct",
            "pmids": [str(item.get("pmid"))],
        }
        for item in strong_items
    ]
    for bucket in ("inference_only_claims", "unsupported_or_overstated_claims"):
        for item in audit.get(bucket) or []:
            if isinstance(item, dict):
                item["pmids"] = []
    plan["protocol_audit"] = audit

    candidate = plan.get("candidate") or {}
    plan["working_hypothesis"] = {
        "statement": (
            f"Distinguish direct production, indirect ecological production, and independent or interactive "
            f"effects of {candidate.get('bacteria') or 'the candidate microbe'} and "
            f"{candidate.get('metabolite') or 'the candidate metabolite'} on "
            f"{candidate.get('disease') or 'the disease phenotype'}."
        ),
        "direct_support": [
            f"PMID {item.get('pmid')}: {item.get('claim')}" for item in strong_items
        ],
        "inference_only": [
            "No collection of separate papers is treated as proof of the complete microbe-metabolite-disease chain."
        ],
    }

    plan["retrieved_literature"] = [
        item
        for item in (plan.get("retrieved_literature") or [])
        if str(item.get("pmid") or "").strip() in allowed_pmids
    ]
    method_bundle = plan.get("fulltext_method_evidence")
    if isinstance(method_bundle, dict):
        for bucket in ("in_vitro", "in_vivo", "all"):
            if isinstance(method_bundle.get(bucket), list):
                method_bundle[bucket] = [
                    item
                    for item in method_bundle.get(bucket) or []
                    if str(item.get("pmid") or "").strip() in allowed_pmids
                ]

    for item in plan.get("evidence_strength_map") or []:
        if not isinstance(item, dict):
            continue
        item["representative_pmids"] = [
            pmid for pmid in _as_str_list(item.get("representative_pmids")) if pmid in allowed_pmids
        ]
    for item in plan.get("self_reflection") or []:
        if not isinstance(item, dict):
            continue
        checked = []
        for value in _as_str_list(item.get("evidence_checked")):
            mentioned = re.findall(r"\bPMID\s*:?\s*(\d+)\b", value, flags=re.IGNORECASE)
            if not mentioned or all(pmid in allowed_pmids for pmid in mentioned):
                checked.append(value)
        item["evidence_checked"] = checked
    for field in ("protocol_audit", "self_reflection", "overall_risk_flags"):
        plan[field] = _sanitize_nested_pmid_references(plan.get(field), allowed_pmids)
    return plan


def _question_is_microbe_metabolite_study(research_question: str, prompt_constraints: str = "") -> bool:
    raw_text = f"{research_question} {prompt_constraints}"
    text = raw_text.lower()
    microbe_signal = bool(
        re.search(r"\b(?:microbe|microbial|microbiome|bacter(?:ia|ium|ial)?|probiotic|strain|colonization)\b", text)
        or re.search(r"(?<![a-z0-9])akk(?:ermansia(?: muciniphila)?)?(?=$|[^a-z0-9])", text)
        or re.search(r"\b[A-Z][a-z]{2,}[_\s][a-z][a-z0-9-]{2,}\b", raw_text)
        or any(
            term in raw_text
            for term in (
                "\u5fae\u751f\u7269",
                "\u7ec6\u83cc",
                "\u83cc\u682a",
                "\u76ca\u751f\u83cc",
                "\u5b9a\u690d",
                "\u80a0\u9053\u83cc\u7fa4",
            )
        )
    )
    metabolite_signal = bool(
        re.search(r"\b(?:metabolite|metabolic product|short-chain fatty acid|scfa|bcfa|butyrate|propionate|acetate|succinate|isobutyrate|indole|bile acid)\b", text)
        or " acid" in text
        or any(
            term in raw_text
            for term in (
                "\u4ee3\u8c22\u7269",
                "\u4ee3\u8c22\u4ea7\u7269",
                "\u77ed\u94fe\u8102\u80aa\u9178",
                "\u4e01\u9178",
                "\u5f02\u4e01\u9178",
                "\u4e19\u9178",
                "\u4e59\u9178",
                "\u7425\u73c0\u9178",
                "\u5432\u54da",
                "\u80c6\u6c41\u9178",
            )
        )
    )
    return microbe_signal and metabolite_signal


def _filter_question_relevant_articles(
    articles: List[dict],
    research_question: str,
    disease: str = "",
) -> List[dict]:
    keywords = [keyword.lower() for keyword in _extract_question_keywords(research_question, limit=18)]
    disease_aliases = _citation_disease_aliases(disease) if disease else []
    experimental_terms = (
        "culture", "cultured", "intervention", "administered", "treated", "supplemented",
        "exposed", "cell", "organoid", "mouse", "mice", "animal", "randomized", "cohort",
        "patients", "measured", "quantified", "assay", "metabolomics", "mass spectrometry",
    )
    selected: List[dict] = []
    for article in articles:
        text = f"{article.get('title', '')} {article.get('abstract', '')}"
        if _contains_any(text, REVIEW_LIKE_TERMS) or not _contains_any(text, experimental_terms):
            continue
        matched = sum(1 for keyword in keywords if _contains_any_complete_term(text, (keyword,)))
        threshold = 2 if len(keywords) >= 3 else 1
        disease_match = _contains_entity(text, disease_aliases) if disease_aliases else False
        if matched >= threshold or (disease_match and matched >= 1):
            selected.append(article)
    return selected


def _sanitize_question_plan_references(plan: dict) -> dict:
    allowed_pmids = {
        str(item.get("pmid") or "").strip()
        for item in (plan.get("retrieved_literature") or [])
        if str(item.get("pmid") or "").strip()
    }
    plan = _sanitize_nested_pmid_references(plan, allowed_pmids)
    audit = plan.get("protocol_audit") if isinstance(plan.get("protocol_audit"), dict) else {}
    audit["supported_claims"] = [
        item
        for item in (audit.get("supported_claims") or [])
        if isinstance(item, dict) and _as_str_list(item.get("pmids"))
    ]
    plan["protocol_audit"] = audit
    plan["evidence_basis"] = [
        item
        for item in (plan.get("evidence_basis") or [])
        if isinstance(item, dict) and _as_str_list(item.get("pmids"))
    ]
    for key in ("in_vitro_plan", "in_vivo_plan", "human_plan"):
        for item in plan.get(key) or []:
            if not isinstance(item, dict):
                continue
            item["source_citations"] = [
                citation
                for citation in _as_str_list(item.get("source_citations"))
                if any(
                    pmid in allowed_pmids
                    for pmid in re.findall(r"\bPMID\s*:?\s*(\d+)\b", citation, flags=re.IGNORECASE)
                )
            ]
    return plan


def _normalize_parameter_provenance(items: Any) -> List[dict]:
    normalized: List[dict] = []
    if not isinstance(items, list):
        return normalized
    allowed_statuses = {"reported", "adapted", "proposed_pilot", "unresolved"}
    for item in items[:24]:
        if isinstance(item, str):
            text = item.strip()
            if text:
                normalized.append(
                    {
                        "parameter": "unspecified parameter",
                        "value": text,
                        "status": "unresolved",
                        "source": "NONE",
                        "source_context": "",
                        "transfer_rationale": "",
                        "pilot_check": "",
                    }
                )
            continue
        if not isinstance(item, dict):
            continue
        status = str(item.get("status") or "unresolved").strip().lower()
        if status not in allowed_statuses:
            status = "unresolved"
        normalized.append(
            {
                "parameter": str(item.get("parameter") or "unspecified parameter").strip(),
                "value": str(item.get("value") or "").strip(),
                "status": status,
                "source": str(item.get("source") or "NONE").strip(),
                "source_context": str(item.get("source_context") or "").strip(),
                "transfer_rationale": str(item.get("transfer_rationale") or "").strip(),
                "pilot_check": str(item.get("pilot_check") or "").strip(),
            }
        )
    return normalized


CANONICAL_EXPERIMENT_ROLES = {
    "direct_production_gate",
    "indirect_ecology_gate",
    "host_response_gate",
    "effect_interaction_gate",
    "direct_causal_rescue",
    "indirect_ecology_causal",
    "human_translation",
}
EXPERIMENT_ROLE_ALIASES = {
    "direct_production": "direct_production_gate",
    "direct_production_assay": "direct_production_gate",
    "indirect_ecology": "indirect_ecology_gate",
    "indirect_production": "indirect_ecology_gate",
    "host_response": "host_response_gate",
    "host_activity": "host_response_gate",
    "host_effect": "host_response_gate",
    "effect_and_interaction": "effect_interaction_gate",
    "effect_interaction": "effect_interaction_gate",
    "animal_study_gate": "effect_interaction_gate",
    "animal_effect": "effect_interaction_gate",
    "direct_causal": "direct_causal_rescue",
    "indirect_causal": "indirect_ecology_causal",
    "human": "human_translation",
    "human_observation": "human_translation",
    "human_observational": "human_translation",
    "human_observational_validation": "human_translation",
    "human_validation": "human_translation",
    "human_translation_gate": "human_translation",
}
VALID_EXECUTION_STATUSES = {"ready_now", "conditional_future", "excluded"}
VALID_RESULT_STATUSES = {"not_run", "positive", "negative", "inconclusive"}
VALID_DESIGN_STATUSES = {"execution_complete", "decision_complete", "needs_resolution"}
VALID_ACTIVATION_GATE_STATUSES = {"met", "unmet", "unknown", "not_evaluated"}
PUBLIC_EXPERIMENT_MODULE_FIELDS = (
    "module_id",
    "experiment_role",
    "hypothesis_ids",
    "route_status",
    "execution_status",
    "result_status",
    "design_status",
    "activation_gate",
    "scientific_question",
    "why_needed",
    "study_object",
    "groups",
    "group_count",
    "samples_and_timing",
    "primary_indicator",
    "secondary_indicators",
    "key_controls",
    "positive_gate",
    "branch_if_positive",
    "branch_if_negative",
    "unlock_rule",
    "claim_boundary",
    "completion_issues",
)


def _canonicalize_experiment_role(value: Any) -> str:
    role = re.sub(r"[^a-z0-9]+", "_", str(value or "").strip().lower()).strip("_")
    if not role:
        return "unspecified"
    return EXPERIMENT_ROLE_ALIASES.get(role, role)


def _normalize_result_status(value: Any) -> str:
    status = str(value or "not_run").strip().lower()
    aliases = {
        "pending": "not_run",
        "planned": "not_run",
        "not_started": "not_run",
        "running": "not_run",
        "in_progress": "not_run",
        "complete_positive": "positive",
        "passed": "positive",
        "pass": "positive",
        "complete_negative": "negative",
        "failed": "negative",
        "fail": "negative",
        "equivocal": "inconclusive",
        "ambiguous": "inconclusive",
        "na": "not_run",
        "n_a": "not_run",
        "not_applicable": "not_run",
    }
    status = aliases.get(status, status)
    return status if status in VALID_RESULT_STATUSES else "not_run"


def _normalize_design_status(value: Any) -> str:
    status = str(value or "").strip().lower()
    aliases = {
        "complete": "execution_complete",
        "ready": "execution_complete",
        "protocol_complete": "execution_complete",
        "outline_only": "decision_complete",
        "blueprint": "decision_complete",
        "degraded": "needs_resolution",
        "incomplete": "needs_resolution",
    }
    status = aliases.get(status, status)
    return status if status in VALID_DESIGN_STATUSES else ""


def _first_nonempty_text(*values: Any) -> str:
    for value in values:
        if isinstance(value, str) and value.strip():
            return value.strip()
    return ""


def _normalize_group_definitions(value: Any) -> List[dict]:
    raw_items: List[Any] = []
    if isinstance(value, list):
        raw_items = value
    elif isinstance(value, dict):
        for group_type, entries in value.items():
            if isinstance(entries, list):
                for entry in entries:
                    if isinstance(entry, dict):
                        raw_items.append({**entry, "group_type": entry.get("group_type") or group_type})
                    else:
                        raw_items.append({"label": entry, "group_type": group_type})
            elif entries:
                raw_items.append({"label": entries, "group_type": group_type})

    normalized: List[dict] = []
    for index, entry in enumerate(raw_items[:24], start=1):
        if isinstance(entry, str):
            label = entry.strip()
            if label:
                normalized.append(
                    {
                        "group_name": label,
                        "exposure_or_condition": label,
                        "control_purpose": "Prespecified comparison for this module",
                    }
                )
            continue
        if not isinstance(entry, dict):
            continue
        label = _first_nonempty_text(
            entry.get("group_name"), entry.get("label"), entry.get("name"), entry.get("description")
        )
        factor_levels = entry.get("factor_levels") if isinstance(entry.get("factor_levels"), dict) else entry.get("factors")
        exposure = _first_nonempty_text(
            entry.get("exposure_or_condition"),
            entry.get("condition"),
            entry.get("intervention"),
        )
        if not exposure and isinstance(factor_levels, dict):
            exposure = ", ".join(f"{key}={value}" for key, value in factor_levels.items())
        exposure = exposure or label
        normalized.append(
            {
                "group_name": label or f"Group {index}",
                "exposure_or_condition": exposure or "Prespecified condition",
                "control_purpose": _first_nonempty_text(
                    entry.get("control_purpose"), entry.get("purpose")
                ) or "Prespecified comparison for this module",
            }
        )
    return normalized


def _group_definition_labels(items: List[dict]) -> List[str]:
    labels: List[str] = []
    for item in items:
        label = str(item.get("group_name") or item.get("exposure_or_condition") or "").strip()
        if label:
            labels.append(label)
    return _dedupe_preserve(labels)


def _normalize_primary_endpoint(value: Any, fallback_threshold: str = "") -> dict:
    if isinstance(value, dict):
        return {
            "measure": _first_nonempty_text(value.get("measure"), value.get("name"), value.get("endpoint")),
            "timepoint_or_summary": _first_nonempty_text(
                value.get("timepoint_or_summary"), value.get("timepoint"), value.get("summary_measure")
            ),
            "estimand": str(value.get("estimand") or "").strip(),
            "success_criterion": _first_nonempty_text(
                value.get("success_criterion"), value.get("threshold"), fallback_threshold
            ),
        }
    return {
        "measure": str(value or "").strip(),
        "timepoint_or_summary": "",
        "estimand": "",
        "success_criterion": str(fallback_threshold or "").strip(),
    }


def _normalize_experiment_items(
    items: Any,
    is_in_vivo: bool = False,
    is_human: bool = False,
) -> List[dict]:
    normalized = []
    if not isinstance(items, list):
        return normalized
    for item in items[:8]:
        if not isinstance(item, dict):
            continue
        role = _canonicalize_experiment_role(item.get("experiment_role"))
        if is_human:
            role = "human_translation"

        objective = _first_nonempty_text(
            item.get("objective"),
            item.get("scientific_question"),
            item.get("aim"),
            item.get("biological_question"),
        )
        model_spec_input = item.get("model_spec")
        if isinstance(model_spec_input, dict):
            model_name = _first_nonempty_text(
                model_spec_input.get("name"),
                model_spec_input.get("system"),
                model_spec_input.get("description"),
            )
            model_rationale = _first_nonempty_text(
                model_spec_input.get("rationale"), model_spec_input.get("why_this_model")
            )
        else:
            model_name = _first_nonempty_text(model_spec_input, item.get("study_object"))
            model_rationale = ""
        if is_in_vivo:
            model_name = _first_nonempty_text(model_name, item.get("model"), item.get("model_system"))
            model_rationale = _first_nonempty_text(model_rationale, item.get("why_this_model"))
        else:
            model_name = _first_nonempty_text(model_name, item.get("model_system"), item.get("model"))

        group_definitions = _normalize_group_definitions(
            item.get("group_definitions") if item.get("group_definitions") is not None else item.get("groups")
        )
        if not group_definitions:
            group_definitions = _normalize_group_definitions(
                item.get("key_controls") if item.get("key_controls") is not None else item.get("controls")
            )

        sampling_input = item.get("sampling") if isinstance(item.get("sampling"), dict) else {}
        sampling = {
            "experimental_unit": _first_nonempty_text(
                sampling_input.get("experimental_unit"), item.get("experimental_unit")
            ),
            "replication_and_sampling": _first_nonempty_text(
                sampling_input.get("replication_and_sampling"),
                sampling_input.get("replication"),
                item.get("replication_and_sampling"),
            ),
            "timepoints": _as_str_list(sampling_input.get("timepoints")),
            "sampling_mode": str(sampling_input.get("sampling_mode") or "").strip(),
            "sample_size_basis": _first_nonempty_text(
                sampling_input.get("sample_size_basis"), item.get("sample_size_basis")
            ),
            "allocation_and_blinding": _first_nonempty_text(
                sampling_input.get("allocation_and_blinding"),
                sampling_input.get("randomization_and_blinding"),
                item.get("randomization_and_blinding"),
            ),
        }

        endpoints_input = item.get("endpoints") if isinstance(item.get("endpoints"), dict) else {}
        legacy_primary_endpoints = _as_str_list(item.get("primary_endpoints"))
        primary_endpoint_input = endpoints_input.get("primary")
        if primary_endpoint_input is None:
            primary_endpoint_input = item.get("primary_indicator") or item.get("primary_endpoint") or (
                legacy_primary_endpoints[0] if legacy_primary_endpoints else ""
            )
        primary_endpoint = _normalize_primary_endpoint(
            primary_endpoint_input,
            _first_nonempty_text(
                item.get("positive_gate"), item.get("success_threshold"), item.get("go_no_go_threshold")
            ),
        )
        secondary_endpoints = (
            _as_str_list(item.get("secondary_indicators"))
            or _as_str_list(endpoints_input.get("secondary"))
            or _as_str_list(item.get("secondary_endpoints"))
        )
        mechanistic_endpoints = _as_str_list(
            endpoints_input.get("mechanistic_or_exposure")
        ) or _as_str_list(item.get("mechanistic_endpoints"))
        safety_endpoints = _as_str_list(endpoints_input.get("safety"))
        endpoints = {
            "primary": primary_endpoint,
            "secondary": secondary_endpoints,
            "mechanistic_or_exposure": mechanistic_endpoints,
            "safety": safety_endpoints,
        }

        analysis_input = item.get("analysis_plan")
        if analysis_input is None:
            analysis_input = item.get("analysis")
        if isinstance(analysis_input, dict):
            analysis_plan = {
                "model": _first_nonempty_text(analysis_input.get("model"), analysis_input.get("summary")),
                "planned_contrasts": _as_str_list(analysis_input.get("planned_contrasts")),
                "repeated_or_cluster_structure": _first_nonempty_text(
                    analysis_input.get("repeated_or_cluster_structure"),
                    analysis_input.get("repeated_measures"),
                    analysis_input.get("clustering"),
                ),
                "multiplicity": str(analysis_input.get("multiplicity") or "").strip(),
                "missing_data": str(analysis_input.get("missing_data") or "").strip(),
            }
        else:
            analysis_plan = {
                "model": _first_nonempty_text(analysis_input, item.get("data_analysis")),
                "planned_contrasts": [],
                "repeated_or_cluster_structure": "",
                "multiplicity": "",
                "missing_data": "",
            }

        decision_input = item.get("decision") if isinstance(item.get("decision"), dict) else {}
        decision = {
            "success_criterion": _first_nonempty_text(
                decision_input.get("success_criterion"),
                item.get("positive_gate"),
                primary_endpoint.get("success_criterion"),
            ),
            "if_positive": _first_nonempty_text(
                decision_input.get("if_positive"),
                item.get("branch_if_positive"),
                item.get("positive_result_interpretation"),
            ),
            "if_negative": _first_nonempty_text(
                decision_input.get("if_negative"),
                item.get("branch_if_negative"),
                item.get("negative_result_interpretation"),
            ),
            "failure_action": _first_nonempty_text(
                decision_input.get("failure_action"), item.get("failure_action")
            ),
            "go_no_go": _first_nonempty_text(decision_input.get("go_no_go"), item.get("go_no_go")),
        }

        legacy_prerequisites = _as_str_list(item.get("prerequisite_module_ids"))
        activation_gate = _normalize_activation_gate(item.get("activation_gate"))
        if not activation_gate and legacy_prerequisites:
            activation_gate = _activation_gate_from_legacy_prerequisites(legacy_prerequisites, role=role)

        parameter_input = item.get("parameters")
        if parameter_input is None:
            parameter_input = item.get("parameter_provenance")
        parameter_provenance = _normalize_parameter_provenance(parameter_input)
        execution_status = str(item.get("execution_status") or "").strip().lower()
        result_status = _normalize_result_status(item.get("result_status"))
        design_status = _normalize_design_status(item.get("design_status") or item.get("generation_status"))
        activation_gate_status = str(item.get("activation_gate_status") or "not_evaluated").strip().lower()
        if activation_gate_status not in VALID_ACTIVATION_GATE_STATUSES:
            activation_gate_status = "not_evaluated"

        stage = str(item.get("stage") or "").strip().lower()
        if not stage:
            if is_human or role == "human_translation":
                stage = "human"
            elif is_in_vivo:
                stage = "animal"
            elif role == "host_response_gate":
                stage = "host_cell"
            else:
                stage = "source_validation"

        hypothesis_ids = _as_str_list(item.get("hypothesis_ids"))
        if not hypothesis_ids:
            for hypothesis in _as_str_list(item.get("hypothesis_tested")):
                matches = re.findall(r"\bH\d+\b", hypothesis, flags=re.IGNORECASE)
                hypothesis_ids.extend(match.upper() for match in matches)
        hypothesis_ids = _dedupe_preserve(hypothesis_ids)
        route_status = _normalize_route_status(item.get("route_status"), role)
        scientific_question = _first_nonempty_text(
            item.get("scientific_question"), item.get("biological_question"), objective
        )
        why_needed = _first_nonempty_text(
            item.get("why_needed"),
            item.get("priority_rationale"),
            item.get("gap_addressed"),
            item.get("design"),
        )
        study_object = _first_nonempty_text(item.get("study_object"), model_name)
        samples_and_timing = _first_nonempty_text(
            item.get("samples_and_timing"),
            item.get("timeline"),
            item.get("dose_timing_logic"),
            sampling.get("replication_and_sampling"),
        )
        primary_indicator = _first_nonempty_text(
            item.get("primary_indicator"), primary_endpoint.get("measure")
        )
        secondary_indicators = _as_str_list(item.get("secondary_indicators")) or secondary_endpoints
        if not secondary_indicators:
            secondary_indicators = [
                value
                for value in _as_str_list(item.get("readouts"))
                if value.strip().lower() != primary_indicator.strip().lower()
            ]
        key_controls = _as_str_list(item.get("key_controls")) or _as_str_list(item.get("controls"))
        positive_gate = _first_nonempty_text(
            item.get("positive_gate"), decision.get("success_criterion")
        )
        unlock_rule = _first_nonempty_text(
            item.get("unlock_rule"), item.get("stage_gate"), item.get("prerequisite_result")
        )
        claim_boundary = _first_nonempty_text(
            item.get("claim_boundary"), item.get("claim_scope"), DEFAULT_CLAIM_SCOPES.get(role)
        )

        normalized_item = {
            "module_id": str(item.get("module_id") or "").strip(),
            "stage": stage,
            "experiment_role": role,
            "hypothesis_ids": hypothesis_ids,
            "route_status": route_status,
            "scientific_question": scientific_question,
            "why_needed": why_needed,
            "study_object": study_object,
            "groups": group_definitions,
            "group_count": len(group_definitions),
            "samples_and_timing": samples_and_timing,
            "primary_indicator": primary_indicator,
            "secondary_indicators": secondary_indicators,
            "key_controls": key_controls,
            "positive_gate": positive_gate,
            "unlock_rule": unlock_rule,
            "claim_boundary": claim_boundary,
            "objective": objective,
            "claim_scope": claim_boundary,
            "execution_status": execution_status,
            "result_status": result_status,
            "design_status": design_status,
            "activation_gate": activation_gate,
            "activation_gate_status": activation_gate_status,
            "activation_gate_reasons": _as_str_list(item.get("activation_gate_reasons")),
            "prerequisite_module_ids": legacy_prerequisites,
            "model_spec": {"name": model_name, "rationale": model_rationale},
            "group_definitions": group_definitions,
            "sampling": sampling,
            "endpoints": endpoints,
            "analysis_plan": analysis_plan,
            "decision": decision,
            "parameters": parameter_provenance,
            # Legacy projections remain until every caller consumes the compact schema.
            "aim": _first_nonempty_text(item.get("aim"), objective),
            "hypothesis_tested": _as_str_list(item.get("hypothesis_tested")) or hypothesis_ids,
            "prerequisite_result": str(item.get("prerequisite_result") or "").strip(),
            "stage_gate": _first_nonempty_text(item.get("stage_gate"), unlock_rule),
            "branch_if_positive": _first_nonempty_text(item.get("branch_if_positive"), decision.get("if_positive")),
            "branch_if_negative": _first_nonempty_text(item.get("branch_if_negative"), decision.get("if_negative")),
            "biological_question": _first_nonempty_text(item.get("biological_question"), objective),
            "priority_rationale": _first_nonempty_text(item.get("priority_rationale"), why_needed),
            "gap_addressed": str(item.get("gap_addressed") or "").strip(),
            "model_system": _first_nonempty_text(item.get("model_system"), model_name),
            "experimental_material": str(item.get("experimental_material") or "").strip(),
            "key_materials_equipment": _as_str_list(item.get("key_materials_equipment")),
            "design": str(item.get("design") or "").strip(),
            "group_logic": str(item.get("group_logic") or "").strip(),
            "procedure_steps": _as_str_list(item.get("procedure_steps")),
            "intervention": str(item.get("intervention") or "").strip(),
            "dose_timing_logic": str(item.get("dose_timing_logic") or "").strip(),
            "controls": key_controls,
            "readouts": _as_str_list(item.get("readouts")) or _dedupe_preserve(
                ([primary_indicator] if primary_indicator else []) + secondary_indicators
            ),
            "readout_rationale": str(item.get("readout_rationale") or "").strip(),
            "evidence_basis": _as_str_list(item.get("evidence_basis")),
            "query_round_support": _as_str_list(item.get("query_round_support")),
            "fulltext_method_support": _as_str_list(item.get("fulltext_method_support")),
            "reported_conditions": _as_str_list(item.get("reported_conditions")),
            "source_citations": _as_str_list(item.get("source_citations")),
            "protocols_io_support": _as_str_list(item.get("protocols_io_support")),
            "protocols_io_materials": _as_str_list(item.get("protocols_io_materials")),
            "protocols_io_citations": _as_str_list(item.get("protocols_io_citations")),
            "protocols_io_urls": _as_str_list(item.get("protocols_io_urls")),
            "key_confounders": _as_str_list(item.get("key_confounders")),
            "primary_endpoint": _first_nonempty_text(item.get("primary_endpoint"), primary_endpoint.get("measure")),
            "success_threshold": _first_nonempty_text(item.get("success_threshold"), decision.get("success_criterion")),
            "failure_action": _first_nonempty_text(item.get("failure_action"), decision.get("failure_action")),
            "support_result": str(item.get("support_result") or "").strip(),
            "positive_result_interpretation": str(item.get("positive_result_interpretation") or "").strip(),
            "negative_result_interpretation": str(item.get("negative_result_interpretation") or "").strip(),
            "go_no_go": str(item.get("go_no_go") or "").strip(),
            "mechanism_support_criterion": str(item.get("mechanism_support_criterion") or "").strip(),
            "phenomenology_support_criterion": str(item.get("phenomenology_support_criterion") or "").strip(),
            "decision_impact": str(item.get("decision_impact") or "").strip(),
            "data_analysis": _first_nonempty_text(item.get("data_analysis"), analysis_plan.get("model")),
            "experimental_unit": sampling.get("experimental_unit") or "",
            "replication_and_sampling": sampling.get("replication_and_sampling") or "",
            "sample_size_basis": sampling.get("sample_size_basis") or "",
            "randomization_and_blinding": sampling.get("allocation_and_blinding") or "",
            "parameter_provenance": parameter_provenance,
            "safety_and_stopping_rules": str(item.get("safety_and_stopping_rules") or "").strip(),
            "completion_issues": _as_str_list(item.get("completion_issues")),
            "generation_status": str(item.get("generation_status") or "").strip().lower(),
        }
        if is_in_vivo:
            normalized_item["model"] = _first_nonempty_text(item.get("model"), model_name)
            normalized_item["why_this_model"] = _first_nonempty_text(item.get("why_this_model"), model_rationale)
            normalized_item["intervention_route"] = str(item.get("intervention_route") or "").strip()
            normalized_item["timeline"] = str(item.get("timeline") or "").strip()
            normalized_item["primary_endpoints"] = legacy_primary_endpoints or (
                [primary_endpoint.get("measure")] if primary_endpoint.get("measure") else []
            )
            normalized_item["secondary_endpoints"] = secondary_endpoints
            normalized_item["mechanistic_endpoints"] = mechanistic_endpoints
            normalized_item["go_no_go_threshold"] = _first_nonempty_text(
                item.get("go_no_go_threshold"), decision.get("success_criterion")
            )
        normalized.append(normalized_item)
    return normalized


EXPERIMENT_ROLE_MODULE_BASE = {
    "direct_production_gate": "V1",
    "indirect_ecology_gate": "V2",
    "host_response_gate": "V3",
    "effect_interaction_gate": "A1",
    "direct_causal_rescue": "A2",
    "indirect_ecology_causal": "A3",
    "human_translation": "H1",
}
EXPERIMENT_ROLE_PREREQUISITES = {
    "direct_production_gate": [],
    "indirect_ecology_gate": ["V1"],
    "host_response_gate": [],
    "effect_interaction_gate": ["V3", "STRONG_HOST_EFFECT_EVIDENCE"],
    "direct_causal_rescue": ["V1", "A1"],
    "indirect_ecology_causal": ["V2", "A1"],
    "human_translation": ["A1_OR_A2_OR_A3", "INDEPENDENT_REPLICATION", "ETHICS_AND_ANALYSIS_PLAN"],
}

EXPERIMENT_ROLE_ACTIVATION_GATES = {
    "direct_production_gate": {"operator": "all_of", "clauses": []},
    "indirect_ecology_gate": {
        "operator": "any_of",
        "clauses": [
            {"module_id": "V1", "accepted_results": ["negative", "inconclusive"]},
            {
                "operator": "all_of",
                "clauses": [
                    {"module_id": "V1", "accepted_results": ["positive"]},
                    {"external_gate": "direct_output_insufficient", "accepted_results": ["positive"]},
                ],
            },
            {
                "external_gate": "ecological_amplification_objective",
                "accepted_results": ["positive"],
            },
        ],
    },
    "host_response_gate": {"operator": "all_of", "clauses": []},
    "effect_interaction_gate": {
        "operator": "any_of",
        "clauses": [
            {"module_id": "V3", "accepted_results": ["positive"]},
            {"external_gate": "strong_host_effect_evidence", "accepted_results": ["positive"]},
        ],
    },
    "direct_causal_rescue": {
        "operator": "all_of",
        "clauses": [
            {"module_id": "V1", "accepted_results": ["positive"]},
            {"module_id": "A1", "accepted_results": ["positive"]},
        ],
    },
    "indirect_ecology_causal": {
        "operator": "all_of",
        "clauses": [
            {"module_id": "V2", "accepted_results": ["positive"]},
            {"module_id": "A1", "accepted_results": ["positive"]},
        ],
    },
    "human_translation": {
        "operator": "all_of",
        "clauses": [
            {
                "operator": "any_of",
                "clauses": [
                    {"module_id": "A1", "accepted_results": ["positive"]},
                    {"module_id": "A2", "accepted_results": ["positive"]},
                    {"module_id": "A3", "accepted_results": ["positive"]},
                ],
            },
            {"external_gate": "independent_replication", "accepted_results": ["positive"]},
            {"external_gate": "ethics_and_analysis_plan", "accepted_results": ["positive"]},
        ],
    },
}
DEFAULT_CLAIM_SCOPES = {
    "direct_production_gate": "candidate-specific source attribution under the tested culture conditions",
    "indirect_ecology_gate": "candidate-dependent ecological modulation under the tested community conditions",
    "host_response_gate": "host activity without microbial-source or mediation inference",
    "effect_interaction_gate": "animal main effects and interaction without mediation inference",
    "direct_causal_rescue": "direct-branch causal contribution when loss, complementation, and rescue align",
    "indirect_ecology_causal": "defined-community ecological causal contribution",
    "human_translation": "association and temporal compatibility only; no production, mediation, or disease-causation claim",
}
DEFAULT_ROUTE_STATUSES = {
    "direct_production_gate": "first_priority",
    "indirect_ecology_gate": "next_if_negative",
    "host_response_gate": "parallel_support",
    "effect_interaction_gate": "next_if_positive",
    "direct_causal_rescue": "next_if_positive",
    "indirect_ecology_causal": "next_if_positive",
    "human_translation": "deferred",
}
VALID_ROUTE_STATUSES = {
    "first_priority",
    "parallel_support",
    "next_if_positive",
    "next_if_negative",
    "alternative",
    "deferred",
}


def _normalize_route_status(value: Any, role: str) -> str:
    status = str(value or "").strip().lower()
    aliases = {
        "first_gate": "first_priority",
        "first": "first_priority",
        "parallel": "parallel_support",
        "conditional_future": DEFAULT_ROUTE_STATUSES.get(role, "deferred"),
        "conditional": DEFAULT_ROUTE_STATUSES.get(role, "deferred"),
        "future": "deferred",
    }
    status = aliases.get(status, status)
    if status not in VALID_ROUTE_STATUSES:
        status = DEFAULT_ROUTE_STATUSES.get(role, "deferred")
    return status


def _normalize_gate_result_values(value: Any) -> List[str]:
    raw_values = value if isinstance(value, list) else [value]
    normalized: List[str] = []
    for raw in raw_values:
        if raw is None or str(raw).strip() == "":
            continue
        status = _normalize_result_status(raw)
        # Activation gates must be driven by observed outcomes. Planning or
        # execution labels normalize to not_run and must never unlock a branch.
        if status == "not_run":
            continue
        if status not in normalized:
            normalized.append(status)
    return normalized or ["positive"]


def _normalize_activation_gate(value: Any) -> dict:
    if value is None:
        return {}
    if isinstance(value, list):
        return {
            "operator": "all_of",
            "clauses": [gate for gate in (_normalize_activation_gate(item) for item in value) if gate],
        }
    if isinstance(value, str):
        return _activation_gate_from_legacy_prerequisites([value])
    if not isinstance(value, dict):
        return {}

    if value.get("module_id"):
        return {
            "module_id": str(value.get("module_id") or "").strip(),
            "accepted_results": _normalize_gate_result_values(
                value.get("accepted_results")
                if value.get("accepted_results") is not None
                else value.get("required_results") or value.get("required_result")
            ),
        }
    if value.get("external_gate"):
        return {
            "external_gate": str(value.get("external_gate") or "").strip().lower(),
            "accepted_results": _normalize_gate_result_values(
                value.get("accepted_results")
                if value.get("accepted_results") is not None
                else value.get("required_results") or value.get("required_result")
            ),
        }

    operator = str(value.get("operator") or "").strip().lower()
    clauses = value.get("clauses")
    if operator not in {"all_of", "any_of"}:
        if isinstance(value.get("all_of"), list):
            operator = "all_of"
            clauses = value.get("all_of")
        elif isinstance(value.get("any_of"), list):
            operator = "any_of"
            clauses = value.get("any_of")
    if operator not in {"all_of", "any_of"}:
        return {}
    if not isinstance(clauses, list):
        clauses = []
    return {
        "operator": operator,
        "clauses": [gate for gate in (_normalize_activation_gate(item) for item in clauses) if gate],
    }


def _legacy_prerequisite_leaf(value: str, role: str = "") -> dict:
    clean = str(value or "").strip()
    upper = clean.upper()
    if upper in {
        "INDEPENDENT_REPLICATION",
        "ETHICS_AND_ANALYSIS_PLAN",
        "STRONG_HOST_EFFECT_EVIDENCE",
    }:
        return {
            "external_gate": upper.lower(),
            "accepted_results": ["positive"],
        }
    accepted_results = ["positive"]
    if role == "indirect_ecology_gate" and upper == "V1":
        accepted_results = ["negative", "inconclusive"]
    return {"module_id": clean, "accepted_results": accepted_results}


def _activation_gate_from_legacy_prerequisites(values: Any, role: str = "") -> dict:
    clauses: List[dict] = []
    for value in _as_str_list(values):
        alternatives = [part.strip() for part in value.split("_OR_") if part.strip()]
        if len(alternatives) > 1:
            clauses.append(
                {
                    "operator": "any_of",
                    "clauses": [_legacy_prerequisite_leaf(part, role=role) for part in alternatives],
                }
            )
        elif alternatives:
            clauses.append(_legacy_prerequisite_leaf(alternatives[0], role=role))
    return {"operator": "all_of", "clauses": clauses}


def _default_activation_gate(role: str) -> dict:
    return _normalize_activation_gate(EXPERIMENT_ROLE_ACTIVATION_GATES.get(role)) or {
        "operator": "all_of",
        "clauses": [],
    }


def _activation_gate_has_constraints(gate: Any) -> bool:
    normalized = _normalize_activation_gate(gate)
    if not normalized:
        return False
    if normalized.get("module_id") or normalized.get("external_gate"):
        return True
    return any(_activation_gate_has_constraints(clause) for clause in normalized.get("clauses") or [])


def _enforce_role_activation_gate(role: str, provided_gate: Any) -> dict:
    """Keep every canonical role gate code-owned, including empty V1/V3 gates."""
    canonical_role = _canonicalize_experiment_role(role)
    required = _default_activation_gate(canonical_role)
    provided = _normalize_activation_gate(provided_gate)
    if canonical_role in CANONICAL_EXPERIMENT_ROLES:
        return required
    return provided or required


def _canonical_gate_contract(role: str, provided_gate: Any = None) -> dict:
    canonical_role = _canonicalize_experiment_role(role)
    gate = _enforce_role_activation_gate(canonical_role, provided_gate)
    summary = _activation_gate_summary(gate)
    if summary:
        unlock_rule = f"Unlock when {summary}."
    elif canonical_role == "direct_production_gate":
        unlock_rule = "No upstream result prerequisite; run as the first source-attribution gate."
    elif canonical_role == "host_response_gate":
        unlock_rule = "No upstream result prerequisite; this host-activity module may run in parallel with V1."
    else:
        unlock_rule = "No upstream result prerequisite; execution still requires a decision-complete design."
    return {
        "activation_gate": gate,
        "prerequisite_module_ids": _activation_gate_module_references(gate),
        "unlock_rule": unlock_rule,
    }


def _canonicalize_module_reference_text(value: Any) -> str:
    """Normalize obsolete route labels without rewriting biological uses such as M1 macrophages."""
    text = str(value or "")
    text = re.sub(
        r"\(\s*(?:module\s+)?M([123])\s*\)",
        lambda match: f"(V{match.group(1)})",
        text,
        flags=re.IGNORECASE,
    )
    text = re.sub(
        r"\bmodule\s+M([123])\b",
        lambda match: f"module V{match.group(1)}",
        text,
        flags=re.IGNORECASE,
    )
    text = re.sub(
        r"\b(from|to|after|before|requires?|requiring)\s+M([123])\b",
        lambda match: f"{match.group(1)} V{match.group(2)}",
        text,
        flags=re.IGNORECASE,
    )
    return re.sub(r"\s+", " ", text).strip()


def _remap_activation_gate_module_ids(gate: Any, module_id_map: Dict[str, str]) -> dict:
    normalized = _normalize_activation_gate(gate)
    if not normalized:
        return {}
    module_id = str(normalized.get("module_id") or "").strip()
    if module_id:
        normalized["module_id"] = module_id_map.get(module_id, module_id)
        return normalized
    clauses = normalized.get("clauses") or []
    normalized["clauses"] = [
        remapped
        for remapped in (
            _remap_activation_gate_module_ids(clause, module_id_map) for clause in clauses
        )
        if remapped
    ]
    return normalized


def _activation_gate_module_references(gate: Any) -> List[str]:
    normalized = _normalize_activation_gate(gate)
    if not normalized:
        return []
    module_id = str(normalized.get("module_id") or "").strip()
    if module_id:
        return [module_id]
    references: List[str] = []
    for clause in normalized.get("clauses") or []:
        references.extend(_activation_gate_module_references(clause))
    return _dedupe_preserve(references)


def _gate_external_result(external_results: Any, gate_name: str) -> str:
    if not isinstance(external_results, dict):
        return "not_run"
    value = external_results.get(gate_name)
    if value is None:
        value = next(
            (
                candidate
                for key, candidate in external_results.items()
                if str(key or "").strip().lower() == gate_name
            ),
            None,
        )
    if isinstance(value, dict):
        value = value.get("result_status") if value.get("result_status") is not None else value.get("status")
    return _normalize_result_status(value)


def _evaluate_activation_gate(
    gate: Any,
    modules_by_id: Dict[str, dict],
    external_results: Optional[Dict[str, Any]] = None,
    _visited_module_ids: Optional[set[str]] = None,
) -> dict:
    normalized = _normalize_activation_gate(gate)
    if not normalized:
        return {"status": "met", "reasons": []}

    visited_module_ids = set(_visited_module_ids or set())
    module_id = str(normalized.get("module_id") or "").strip()
    if module_id:
        if module_id in visited_module_ids:
            return {
                "status": "unknown",
                "reasons": [f"Activation gate contains a dependency cycle at module {module_id}."],
            }
        target = modules_by_id.get(module_id)
        if target is None:
            return {
                "status": "unknown",
                "reasons": [f"Activation gate references missing module {module_id}."],
            }
        actual = _normalize_result_status(target.get("result_status"))
        accepted = _normalize_gate_result_values(normalized.get("accepted_results"))
        if actual not in accepted:
            return {
                "status": "unmet",
                "reasons": [
                    f"Module {module_id} result_status={actual}; accepted_results={','.join(accepted)}."
                ],
            }
        target_design_status = _normalize_design_status(target.get("design_status"))
        if target_design_status == "needs_resolution":
            return {
                "status": "unmet",
                "reasons": [
                    f"Module {module_id} has result_status={actual}, but design_status=needs_resolution makes the result ineligible for downstream activation."
                ],
            }
        target_role = _canonicalize_experiment_role(target.get("experiment_role"))
        target_gate = _enforce_role_activation_gate(target_role, target.get("activation_gate"))
        target_gate_result = _evaluate_activation_gate(
            target_gate,
            modules_by_id,
            external_results,
            visited_module_ids | {module_id},
        )
        if target_gate_result.get("status") != "met":
            return {
                "status": target_gate_result.get("status") or "unknown",
                "reasons": _dedupe_preserve(
                    [f"Module {module_id} result is not eligible because its own activation gate is not met."]
                    + _as_str_list(target_gate_result.get("reasons"))
                ),
            }
        return {"status": "met", "reasons": []}

    external_gate = str(normalized.get("external_gate") or "").strip().lower()
    if external_gate:
        actual = _gate_external_result(external_results, external_gate)
        accepted = _normalize_gate_result_values(normalized.get("accepted_results"))
        if actual in accepted:
            return {"status": "met", "reasons": []}
        return {
            "status": "unmet",
            "reasons": [
                f"External gate {external_gate} result_status={actual}; accepted_results={','.join(accepted)}."
            ],
        }

    operator = str(normalized.get("operator") or "all_of")
    clause_results = [
        _evaluate_activation_gate(
            clause,
            modules_by_id,
            external_results,
            visited_module_ids,
        )
        for clause in normalized.get("clauses") or []
    ]
    if not clause_results:
        return {"status": "met" if operator == "all_of" else "unmet", "reasons": []}

    statuses = [result.get("status") for result in clause_results]
    if operator == "any_of":
        if "met" in statuses:
            return {"status": "met", "reasons": []}
        # Missing optional OR branches are not mandatory when at least one
        # defined branch exists and is simply awaiting a result.
        if "unmet" in statuses:
            status = "unmet"
            reason_results = [result for result in clause_results if result.get("status") == "unmet"]
        else:
            status = "unknown"
            reason_results = clause_results
    else:
        if all(status == "met" for status in statuses):
            return {"status": "met", "reasons": []}
        status = "unmet" if "unmet" in statuses else "unknown"
        reason_results = clause_results
    reasons = _dedupe_preserve(
        reason
        for result in reason_results
        for reason in _as_str_list(result.get("reasons"))
    )
    return {"status": status, "reasons": reasons}


def _parameter_provenance_from_reported_conditions(item: dict) -> List[dict]:
    provenance: List[dict] = []
    for condition in _as_str_list(item.get("reported_conditions"))[:12]:
        lowered = condition.lower()
        status = "unresolved"
        if any(term in lowered for term in ("pilot", "proposed", "optimization", "optimisation")):
            status = "proposed_pilot"
        elif re.search(r"\b(?:PMID|PMCID|DOI)\b", condition, flags=re.IGNORECASE):
            status = "reported"
        source_match = re.search(r"\b(?:PMID|PMCID|DOI)\s*:?[ ]*([^;|,]+)", condition, flags=re.IGNORECASE)
        provenance.append(
            {
                "parameter": "reported condition",
                "value": condition,
                "status": status,
                "source": source_match.group(0).strip() if source_match else "NONE",
                "source_context": "",
                "transfer_rationale": "",
                "pilot_check": condition if status in {"proposed_pilot", "unresolved"} else "",
            }
        )
    return provenance


def _default_experimental_unit(role: str, is_in_vivo: bool = False, is_human: bool = False) -> str:
    if is_human:
        return "One consented participant; site and repeated observations are modeled explicitly when applicable."
    if is_in_vivo:
        return "One animal, with cage or isolator represented as a design and analysis factor when applicable."
    if role in {"direct_production_gate", "indirect_ecology_gate"}:
        return "One independently inoculated culture vessel."
    if role == "host_response_gate":
        return "One independently treated host-model unit nested within an independent experiment, donor, or organoid line."
    return "One independently assigned biological unit."


def _experiment_text_bundle(item: dict, is_in_vivo: bool = False) -> str:
    fields = [
        item.get("aim"), item.get("biological_question"), item.get("model" if is_in_vivo else "model_system"),
        item.get("experimental_material"), item.get("design"), item.get("group_logic"),
        item.get("intervention"), item.get("dose_timing_logic"), item.get("timeline"),
        item.get("primary_endpoint"), item.get("success_threshold"), item.get("data_analysis"),
        item.get("replication_and_sampling"), item.get("sample_size_basis"), item.get("randomization_and_blinding"),
    ]
    list_fields = [
        "hypothesis_tested", "procedure_steps", "controls", "groups", "readouts", "primary_endpoints",
        "secondary_endpoints", "mechanistic_endpoints", "key_confounders", "reported_conditions",
    ]
    fields.extend(" ".join(_as_str_list(item.get(field))) for field in list_fields)
    return " ".join(str(value or "") for value in fields).lower()


def _human_module_is_observational(item: dict) -> bool:
    group_text = " ".join(
        " ".join(
            str(group.get(field) or "")
            for field in ("group_name", "exposure_or_condition", "control_purpose")
        )
        for group in (item.get("groups") or [])
        if isinstance(group, dict)
    )
    text = " ".join(
        [
            _experiment_text_bundle(item, is_in_vivo=True),
            str(item.get("scientific_question") or ""),
            str(item.get("why_needed") or ""),
            str(item.get("study_object") or ""),
            str(item.get("samples_and_timing") or ""),
            str(item.get("claim_boundary") or ""),
            group_text,
        ]
    ).lower()
    explicitly_observational = any(
        term in text for term in ("observational", "cohort", "longitudinal", "cross-sectional")
    )
    interventional_assignment = any(
        term in text
        for term in (
            "randomized trial",
            "randomised trial",
            "assigned to treatment",
            "intervention arm",
            "administered to participants",
            "dose participants",
        )
    )
    claim_boundary = str(item.get("claim_boundary") or "").lower()
    association_only = "association" in claim_boundary and any(
        term in claim_boundary
        for term in ("cannot prove", "does not prove", "no causal", "not causal")
    )
    return not interventional_assignment and (explicitly_observational or association_only)


def _candidate_experiment_completion_issues(item: dict, is_in_vivo: bool = False, is_human: bool = False) -> List[str]:
    issues: List[str] = []
    role = _canonicalize_experiment_role(item.get("experiment_role"))
    model_field = "model" if is_in_vivo else "model_system"
    required_text = {
        "aim": "Missing experiment aim.",
        model_field: "Missing one prespecified primary model.",
        "design": "Missing executable design description.",
        "primary_endpoint": "Missing one prespecified primary endpoint.",
        "success_threshold": "Missing quantitative or source-resolving success threshold.",
        "failure_action": "Missing failure action.",
        "data_analysis": "Missing statistical analysis plan.",
    }
    for field, message in required_text.items():
        if not str(item.get(field) or "").strip():
            issues.append(message)
    if len(_as_str_list(item.get("procedure_steps"))) < 2:
        issues.append("Missing numbered operational procedure.")
    if not str(item.get("experimental_unit") or "").strip():
        issues.append("Missing experimental unit.")
    text = _experiment_text_bundle(item, is_in_vivo=is_in_vivo)
    if not str(item.get("replication_and_sampling") or "").strip() and not any(
        term in text for term in ("independent biological", "biological replicate", "independent donor", "multiple cages")
    ):
        issues.append("Missing biological-versus-technical replication and sampling structure.")
    if not item.get("parameter_provenance"):
        issues.append("Missing line-item parameter provenance or an explicit pilot plan.")
    if not is_in_vivo:
        model_text = str(item.get("model_system") or "").lower()
        if " or " in model_text and sum(term in model_text for term in ("caco", "ht-29", "macroph", "organoid", "epithelial")) > 1:
            issues.append("Primary host model is presented as a menu rather than one executable model.")
    if role == "direct_production_gate":
        checks = {
            "Missing time-zero or baseline medium measurement.": ("time-zero", "time zero", "baseline"),
            "Missing washed-inoculum carryover control.": ("wash", "carryover"),
            "Missing uninoculated or analytical blank control.": ("uninoculated", "blank"),
            "Missing validated LOD/LOQ logic.": ("loq", "limit of quantification", "limit of detection"),
            "Missing source-resolution or isotope linkage logic.": ("isotope", "source-resolv", "source resolv"),
        }
        for message, terms in checks.items():
            if terms == ("wash", "carryover"):
                missing = not all(term in text for term in terms)
            else:
                missing = not any(term in text for term in terms)
            if missing:
                issues.append(message)
    elif role == "indirect_ecology_gate":
        for message, terms in (
            ("Missing a directly validated producer.", ("validated producer", "producer in monoculture")),
            ("Missing complete candidate-by-producer factorial design.", ("2 x 2", "2 × 2", "factorial")),
            ("Missing conditioned-medium residual-cell testing.", ("residual-cell", "residual cell", "sterility")),
            ("Missing pH/osmolarity/nutrient matching.", ("osmolar", "nutrient-matched", "nutrient matched")),
            ("Missing source-resolution logic for ecological production.", ("source-resolv", "isotope", "release-and-uptake")),
        ):
            if not any(term in text for term in terms):
                issues.append(message)
    elif role == "host_response_gate":
        if not any(term in text for term in ("factorial", "interaction")):
            issues.append("Host-response design does not test the microbe-by-metabolite interaction.")
        if not any(term in text for term in ("viability", "toxicity", "non-toxic", "nontoxic")):
            issues.append("Missing host-model toxicity or viability gate.")
    elif role == "effect_interaction_gate":
        if len(_as_str_list(item.get("groups"))) < 4 or not any(term in text for term in ("2 x 2", "2 × 2", "two-factor", "interaction")):
            issues.append("Missing complete microbe-by-metabolite animal factorial core.")
        for message, terms in (
            ("Missing interaction-focused power or precision basis.", ("power", "sample-size", "sample size", "precision")),
            ("Missing cage or cluster allocation.", ("cage", "cluster", "isolator")),
            ("Missing randomization and blinded assessment.", ("random", "blind")),
            ("Missing sex-handling strategy.", ("sex", "male", "female")),
        ):
            if not any(term in text for term in terms):
                issues.append(message)
    elif role == "direct_causal_rescue":
        for message, terms in (
            ("Missing function-loss arm.", ("function-loss", "function loss", "deficient", "knockout")),
            ("Missing functional or genetic complementation.", ("complement" ,)),
            ("Missing metabolite add-back rescue.", ("add-back", "add back", "rescue")),
            ("Missing comparable microbial exposure or colonization check.", ("comparable exposure", "comparable colonization")),
        ):
            if not any(term in text for term in terms):
                issues.append(message)
    elif role == "indirect_ecology_causal":
        for message, terms in (
            ("Missing controlled defined-community model.", ("defined community", "gnotobiotic")),
            ("Missing candidate-by-producer factorial manipulation.", ("2 x 2", "2 × 2", "two-factor", "factorial")),
            ("Missing directly validated producer.", ("validated producer",)),
            ("Missing source-resolved metabolite output.", ("source-resolved", "source resolved", "isotope")),
        ):
            if not any(term in text for term in terms):
                issues.append(message)
    if is_human and role == "human_translation" and not _human_module_is_observational(item):
        issues.append("Initial human stage is not observational or association-only.")
    return _dedupe_preserve(issues)


def _candidate_experiment_decision_issues(
    item: dict,
    is_in_vivo: bool = False,
    is_human: bool = False,
) -> List[str]:
    issues: List[str] = []
    required_text = {
        "scientific_question": "Missing one scientific question.",
        "why_needed": "Missing the module rationale.",
        "study_object": "Missing one prespecified study object.",
        "samples_and_timing": "Missing collected samples or timing logic.",
        "primary_indicator": "Missing one prespecified primary indicator.",
        "positive_gate": "Missing the positive-result gate.",
        "unlock_rule": "Missing the module unlock rule.",
        "claim_boundary": "Missing the permitted claim boundary.",
    }
    for field, message in required_text.items():
        if not str(item.get(field) or "").strip():
            issues.append(message)
    if not isinstance(item.get("groups"), list) or not item.get("groups"):
        issues.append("Missing experimental groups or comparison strata.")
    if not str(item.get("branch_if_positive") or "").strip():
        issues.append("Missing positive-result branch.")
    if not str(item.get("branch_if_negative") or "").strip():
        issues.append("Missing negative-result branch.")
    if is_human:
        if not _human_module_is_observational(item):
            issues.append("Initial human stage is not observational or association-only.")
    return _dedupe_preserve(issues)


def _candidate_role_strategy_issues(item: dict) -> List[str]:
    """Audit high-level causal design without requiring an SOP or fixed parameters."""
    role = _canonicalize_experiment_role(item.get("experiment_role"))
    groups = item.get("groups") if isinstance(item.get("groups"), list) else []
    group_text = " ".join(
        " ".join(
            str(group.get(field) or "")
            for field in ("group_name", "exposure_or_condition", "control_purpose")
        )
        for group in groups
        if isinstance(group, dict)
    ).lower()
    study_object = str(item.get("study_object") or "").lower()
    samples = str(item.get("samples_and_timing") or "").lower()
    primary = str(item.get("primary_indicator") or "").lower()
    secondary = " ".join(_as_str_list(item.get("secondary_indicators"))).lower()
    controls = " ".join(_as_str_list(item.get("key_controls"))).lower()
    positive = str(item.get("positive_gate") or "").lower()
    route_text = " ".join(
        str(item.get(field) or "")
        for field in ("unlock_rule", "branch_if_positive", "branch_if_negative")
    )
    text = " ".join((study_object, group_text, samples, primary, secondary, controls, positive))
    issues: List[str] = []

    if re.search(r"\bM\d+(?:\.\d+)?\b", route_text, flags=re.IGNORECASE):
        issues.append("A branch or unlock rule references a noncanonical M-series module ID.")

    if role == "direct_production_gate":
        if len(groups) < 3:
            issues.append("V1 must include live candidate culture, uninoculated medium, and a justified substrate or precursor control.")
        if not any(term in text for term in ("uninoculated", "no inoculum", "medium blank", "media blank")):
            issues.append("V1 is missing an uninoculated medium control.")
        if not any(term in text for term in ("substrate", "precursor")):
            issues.append("V1 is missing a biologically justified substrate or precursor control.")
        if not any(term in samples for term in ("baseline", "time zero", "time-zero", "0 h", "0-hour")):
            issues.append("V1 samples do not include a baseline or time-zero measurement.")
        if not any(term in secondary for term in ("growth", "abundance", "biomass", "cell count")):
            issues.append("V1 does not pair target-metabolite output with microbial growth or abundance.")
        if not any(term in positive for term in ("newly formed", "net formation", "source-resolved", "source resolved")):
            issues.append("V1 positive gate does not distinguish newly formed metabolite from background or carryover.")
    elif role == "indirect_ecology_gate":
        systems = sum(
            bool(term in study_object)
            for term in ("conditioned medium", "co-culture", "coculture", "defined community")
        )
        if systems > 1 or " or " in study_object:
            issues.append("V2 lists multiple ecological systems as a menu; select one matched system for the primary test.")
        if not any(term in text for term in ("validated producer", "independently confirmed producer", "producer confirmed")):
            issues.append("V2 does not require independent confirmation of the proposed producer.")
        has_candidate_present = any(
            term in group_text
            for term in ("candidate-present", "candidate present", "with candidate", "plus candidate")
        )
        has_candidate_absent = any(
            term in group_text
            for term in ("candidate-absent", "candidate absent", "without candidate", "excluding candidate", "candidate-free")
        )
        if not (has_candidate_present and has_candidate_absent):
            issues.append("V2 groups are not a matched candidate-present versus candidate-absent comparison within one ecological system.")
        if not any(term in secondary for term in ("producer abundance", "producer activity", "substrate transfer", "cross-feeding")):
            issues.append("V2 does not measure the validated producer or substrate-transfer process alongside metabolite output.")
    elif role == "host_response_gate":
        host_terms = ("epithelial", "organoid", "macrophage", "monocyte", "host-cell", "host cell")
        if " or " in study_object and sum(term in study_object for term in host_terms) > 1:
            issues.append("V3 presents host models as a menu instead of selecting one primary model.")
        if not any(term in group_text for term in ("challenge", "challenged", "stimulus", "stimulated", "disease condition")):
            issues.append("V3 lacks an explicit disease-relevant challenge control.")
        if len(groups) < 5 or not any(term in group_text for term in ("interaction", "combined", "both")):
            issues.append("V3 lacks an unchallenged reference plus the complete challenge-condition microbe-by-metabolite factorial core.")
        if " or " in primary or ";" in primary:
            issues.append("V3 primary indicator is a menu; prespecify one directional primary host indicator.")
        if not any(term in positive for term in ("increase", "decrease", "reduce", "restore", "improve", "attenuate")):
            issues.append("V3 positive gate does not state the expected direction relative to its challenge control.")
        if not any(term in text for term in ("viability", "toxicity", "non-toxic", "nontoxic")) or not any(
            term in positive for term in ("without toxicity", "without viability loss", "non-toxic", "nontoxic")
        ):
            issues.append("V3 lacks a no-toxicity condition in its positive gate.")
    elif role == "effect_interaction_gate":
        required_group_signals = (
            ("healthy", "nondisease", "non-disease"),
            ("disease control", "disease baseline"),
            ("microbial main effect", "microbe effect", "disease + microbe"),
            ("metabolite main effect", "metabolite effect", "disease + metabolite"),
            ("interaction", "combined", "both"),
        )
        if len(groups) != 5 or not all(any(term in group_text for term in terms) for terms in required_group_signals):
            issues.append("A1 must contain the healthy reference and all four disease-condition microbe-by-metabolite groups.")
        if not any(term in secondary for term in ("colonization", "microbial exposure", "microbe abundance", "absolute abundance")):
            issues.append("A1 does not confirm candidate-microbe exposure or colonization.")
        if not any(term in secondary for term in ("target metabolite", "metabolite exposure", "metabolite profile", "scfa", "bcfa")):
            issues.append("A1 does not confirm target-metabolite exposure.")
        if "interaction" not in positive or not any(term in positive for term in ("main effect", "microbial effect", "microbe effect")):
            issues.append("A1 positive gate does not separately define microbial, metabolite, and interaction contrasts.")
    elif role == "direct_causal_rescue":
        if not any(term in text for term in ("function-loss", "function loss", "production-deficient", "source-function")):
            issues.append("A2 lacks a candidate-source function-loss or equivalent necessity comparison.")
        if not any(term in text for term in ("complement", "add-back", "add back", "rescue")):
            issues.append("A2 lacks complementation or target-metabolite add-back rescue.")
        if not any(term in secondary for term in ("colonization", "microbial exposure", "comparable abundance")):
            issues.append("A2 does not verify comparable microbial exposure across causal groups.")
    elif role == "indirect_ecology_causal":
        if not any(term in text for term in ("defined community", "gnotobiotic")):
            issues.append("A3 lacks one controlled defined-community study object.")
        if "validated producer" not in text:
            issues.append("A3 does not condition the design on a V2-validated producer.")
        if not any(term in text for term in ("producer-absent", "producer absent", "without producer", "producer-depleted")):
            issues.append("A3 lacks a producer-necessity comparison.")
        if not any(term in text for term in ("add-back", "add back", "rescue")):
            issues.append("A3 lacks a target-metabolite rescue condition.")
    return _dedupe_preserve(issues)


def _prepare_experiment_modules(items: Any, is_in_vivo: bool = False, is_human: bool = False) -> List[dict]:
    normalized = _normalize_experiment_items(items, is_in_vivo=is_in_vivo, is_human=is_human)
    counters: Dict[str, int] = {}
    module_id_map: Dict[str, str] = {}
    valid_roles = (
        {"human_translation"}
        if is_human
        else (
            {"effect_interaction_gate", "direct_causal_rescue", "indirect_ecology_causal"}
            if is_in_vivo
            else {"direct_production_gate", "indirect_ecology_gate", "host_response_gate"}
        )
    )
    for item in normalized:
        role = _canonicalize_experiment_role(item.get("experiment_role"))
        if is_human:
            role = "human_translation"
            item["experiment_role"] = role
        elif role not in valid_roles:
            role = _infer_candidate_experiment_role(item, is_in_vivo=is_in_vivo)
            item["experiment_role"] = role
        item["route_status"] = _normalize_route_status(item.get("route_status"), role)
        base = EXPERIMENT_ROLE_MODULE_BASE.get(role, "A" if is_in_vivo else "V")
        counters[base] = counters.get(base, 0) + 1
        supplied_module_id = str(item.get("module_id") or "").strip()
        canonical_module_id = base if counters[base] == 1 else f"{base}.{counters[base]}"
        item["module_id"] = canonical_module_id
        if supplied_module_id:
            module_id_map.setdefault(supplied_module_id, canonical_module_id)
        activation_gate = _normalize_activation_gate(item.get("activation_gate"))
        if not activation_gate:
            legacy_prerequisites = _as_str_list(item.get("prerequisite_module_ids"))
            activation_gate = (
                _activation_gate_from_legacy_prerequisites(legacy_prerequisites, role=role)
                if legacy_prerequisites
                else _default_activation_gate(role)
            )
        gate_contract = _canonical_gate_contract(role, activation_gate)
        item.update(gate_contract)
        if not str(item.get("claim_scope") or "").strip():
            item["claim_scope"] = DEFAULT_CLAIM_SCOPES.get(role, "")
        if not item.get("experimental_unit"):
            item["experimental_unit"] = _default_experimental_unit(role, is_in_vivo=is_in_vivo, is_human=is_human)
        sampling = item.get("sampling") if isinstance(item.get("sampling"), dict) else {}
        sampling["experimental_unit"] = item.get("experimental_unit") or ""
        item["sampling"] = sampling
        if not item.get("parameter_provenance"):
            item["parameter_provenance"] = _parameter_provenance_from_reported_conditions(item)
        item["parameters"] = _normalize_parameter_provenance(item.get("parameter_provenance"))
        item["result_status"] = _normalize_result_status(item.get("result_status"))
        status = str(item.get("execution_status") or "").strip().lower()
        if status not in VALID_EXECUTION_STATUSES:
            status = "conditional_future" if is_in_vivo or is_human or role == "indirect_ecology_gate" else "ready_now"
        operational_issues = _candidate_experiment_completion_issues(
            item,
            is_in_vivo=is_in_vivo,
            is_human=is_human,
        )
        decision_issues = _candidate_experiment_decision_issues(
            item,
            is_in_vivo=is_in_vivo,
            is_human=is_human,
        )
        if not operational_issues:
            design_status = "execution_complete"
        elif not decision_issues:
            design_status = "decision_complete"
        else:
            design_status = "needs_resolution"
        if status == "ready_now" and design_status == "needs_resolution":
            status = "conditional_future"
        item["execution_status"] = status
        item["design_status"] = design_status
        item["design_issues"] = decision_issues
        item["operationalization_issues"] = operational_issues
        item["completion_issues"] = decision_issues
        item["generation_status"] = "complete" if design_status != "needs_resolution" else "degraded"
    for item in normalized:
        role = _canonicalize_experiment_role(item.get("experiment_role"))
        remapped_gate = _remap_activation_gate_module_ids(item.get("activation_gate"), module_id_map)
        item.update(_canonical_gate_contract(role, remapped_gate))
        for field in (
            "scientific_question",
            "why_needed",
            "samples_and_timing",
            "positive_gate",
            "branch_if_positive",
            "branch_if_negative",
            "prerequisite_result",
            "stage_gate",
            "claim_boundary",
        ):
            item[field] = _canonicalize_module_reference_text(item.get(field))
    return normalized


def _plan_is_candidate_workflow(plan: dict) -> bool:
    candidate = plan.get("candidate") if isinstance(plan.get("candidate"), dict) else {}
    if str(candidate.get("bacteria") or "").strip() and str(candidate.get("metabolite") or "").strip():
        return True
    brief = plan.get("user_brief") if isinstance(plan.get("user_brief"), dict) else {}
    return _question_is_microbe_metabolite_study(
        str(brief.get("research_question") or ""),
        str(brief.get("prompt_constraints") or ""),
    )


def _dedupe_plan_records(items: List[dict], fields: tuple[str, ...]) -> List[dict]:
    deduped: List[dict] = []
    seen: set[tuple[str, ...]] = set()
    for item in items:
        key = tuple(str(item.get(field) or "").strip() for field in fields)
        if key in seen:
            continue
        seen.add(key)
        deduped.append(item)
    return deduped


def _human_study_excluded_by_user(plan: dict) -> bool:
    brief = plan.get("user_brief") if isinstance(plan.get("user_brief"), dict) else {}
    brief_text = f"{brief.get('research_question') or ''} {brief.get('prompt_constraints') or ''}".lower()
    return any(
        term in brief_text
        for term in (
            "no clinical",
            "avoid clinical",
            "preclinical only",
            "only in vitro",
            "in vitro only",
            "cell only",
            "organoid only",
        )
    )


def _human_gate_for_plan(plan: dict, candidate_workflow: bool) -> dict:
    human_items = [item for item in (plan.get("human_plan") or []) if isinstance(item, dict)]
    if _human_study_excluded_by_user(plan):
        return {
            "status": "excluded_by_user",
            "reason": "The user constrained the workflow to preclinical work, so no human experiment is active or scheduled.",
            "future_requirements": [],
        }
    if any(item.get("execution_status") == "ready_now" for item in human_items):
        return {
            "status": "ready_now",
            "reason": "At least one explicit human module passed its policy, prerequisite, and completeness checks.",
            "future_requirements": [],
        }
    if human_items:
        future_requirements = _dedupe_preserve(
            requirement
            for item in human_items
            for requirement in (
                _as_str_list(item.get("activation_gate_reasons"))
                + _as_str_list(item.get("design_issues"))
                + ([str(item.get("prerequisite_result")).strip()] if str(item.get("prerequisite_result") or "").strip() else [])
            )
            if str(requirement or "").strip()
        )[:12]
        return {
            "status": "blocked_preclinical" if candidate_workflow else "conditional_future",
            "reason": (
                "A human observational blueprint is retained for completeness, but recruitment, sampling, or analysis cannot start "
                "until its preclinical, replication, ethics, and design gates are met."
            ),
            "future_requirements": future_requirements,
        }
    if candidate_workflow:
        return {
            "status": "blocked_preclinical",
            "reason": "The human stage is not yet specified and no independently replicated animal result is available.",
            "future_requirements": [
                "Retain a branch-conditioned observational human blueprint as conditional_future work.",
                "Complete a relevant animal effect or causal module and record result_status=positive.",
                "Record independent_replication and ethics_and_analysis_plan as positive external gate results.",
            ],
        }
    profile = _normalize_question_profile(plan.get("question_profile"))
    policy = profile.get("human_study_policy") or "not_recommended"
    return {
        "status": "not_planned",
        "reason": f"No explicit human module was selected under human_study_policy={policy}.",
        "future_requirements": [],
    }


def _parameter_provenance_audit(plan: dict) -> dict:
    module_results: List[dict] = []
    unresolved_items: List[dict] = []
    for plan_key, modality in (
        ("in_vitro_plan", "in_vitro"),
        ("in_vivo_plan", "in_vivo"),
        ("human_plan", "human"),
    ):
        for item in plan.get(plan_key) or []:
            if not isinstance(item, dict) or item.get("execution_status") == "excluded":
                continue
            module_id = str(item.get("module_id") or "unassigned").strip()
            provenance = _normalize_parameter_provenance(item.get("parameter_provenance"))
            module_issues: List[dict] = []
            if not provenance:
                module_issues.append(
                    {
                        "module_id": module_id,
                        "parameter": "all critical parameters",
                        "value": "",
                        "status": "unresolved",
                        "issue": "No line-item parameter provenance or concrete pilot-verification record was supplied.",
                    }
                )
            for entry in provenance:
                status = str(entry.get("status") or "unresolved").strip().lower()
                source = str(entry.get("source") or "NONE").strip()
                issue = ""
                if not str(entry.get("value") or "").strip():
                    issue = "Parameter value or proposed range is missing."
                elif status == "unresolved":
                    issue = "Parameter remains unresolved."
                elif status in {"reported", "adapted"} and source.upper() == "NONE":
                    issue = "A reported or adapted parameter lacks a traceable source."
                elif status == "adapted" and not str(entry.get("transfer_rationale") or "").strip():
                    issue = "An adapted parameter lacks a transfer rationale."
                elif status == "proposed_pilot" and not str(entry.get("pilot_check") or "").strip():
                    issue = "A proposed pilot parameter lacks a verification check."
                if issue:
                    module_issues.append(
                        {
                            "module_id": module_id,
                            "parameter": str(entry.get("parameter") or "unspecified parameter").strip(),
                            "value": str(entry.get("value") or "").strip(),
                            "status": status,
                            "issue": issue,
                        }
                    )
            unresolved_items.extend(module_issues)
            module_results.append(
                {
                    "module_id": module_id,
                    "modality": modality,
                    "status": "pass" if not module_issues else "needs_resolution",
                    "unresolved_count": len(module_issues),
                }
            )
    return {
        "status": "not_applicable" if not module_results else ("needs_resolution" if unresolved_items else "pass"),
        "module_results": module_results,
        "unresolved_items": unresolved_items,
    }


def _ensure_unique_module_ids(plan: dict) -> None:
    seen: Dict[str, int] = {}
    for plan_key, default_prefix in (
        ("in_vitro_plan", "V"),
        ("in_vivo_plan", "A"),
        ("human_plan", "H"),
    ):
        for item in plan.get(plan_key) or []:
            if not isinstance(item, dict):
                continue
            role = _canonicalize_experiment_role(item.get("experiment_role"))
            item["experiment_role"] = role
            base = str(item.get("module_id") or EXPERIMENT_ROLE_MODULE_BASE.get(role) or default_prefix).strip()
            seen[base] = seen.get(base, 0) + 1
            if seen[base] == 1:
                item["module_id"] = base
                continue
            replacement = f"{base}.{seen[base]}"
            item["module_id"] = replacement
            item["completion_issues"] = _dedupe_preserve(
                _as_str_list(item.get("completion_issues"))
                + [f"Duplicate module ID {base} was remapped to {replacement}; verify all intended dependencies."]
            )


def _apply_candidate_dependency_gates(plan: dict) -> None:
    modules = [
        item
        for plan_key in ("in_vitro_plan", "in_vivo_plan", "human_plan")
        for item in (plan.get(plan_key) or [])
        if isinstance(item, dict)
    ]
    by_id = {
        str(item.get("module_id") or "").strip(): item
        for item in modules
        if str(item.get("module_id") or "").strip()
    }
    external_results = plan.get("external_gate_results")
    if not isinstance(external_results, dict):
        external_results = {}
    for item in modules:
        role = _canonicalize_experiment_role(item.get("experiment_role"))
        item["experiment_role"] = role
        item["result_status"] = _normalize_result_status(item.get("result_status"))
        item.update(_canonical_gate_contract(role, item.get("activation_gate")))
        gate = item["activation_gate"]
        evaluation = _evaluate_activation_gate(gate, by_id, external_results)
        item["activation_gate_status"] = evaluation.get("status") or "unknown"
        item["activation_gate_reasons"] = _as_str_list(evaluation.get("reasons"))
        design_status = _normalize_design_status(item.get("design_status"))
        if not design_status:
            design_status = (
                "execution_complete"
                if not _as_str_list(item.get("completion_issues"))
                else "needs_resolution"
            )
            item["design_status"] = design_status
        if item.get("activation_gate_status") == "unknown":
            item["design_issues"] = _dedupe_preserve(
                _as_str_list(item.get("design_issues"))
                + _as_str_list(item.get("activation_gate_reasons"))
            )
            item["design_status"] = "needs_resolution"
            item["generation_status"] = "degraded"
        if item.get("execution_status") != "excluded":
            item["execution_status"] = (
                "ready_now"
                if item.get("activation_gate_status") == "met"
                and item.get("design_status") != "needs_resolution"
                else "conditional_future"
            )


def _design_completeness_audit(plan: dict, strict_candidate_design: bool = False) -> dict:
    module_results: List[dict] = []
    all_issues: List[dict] = []
    for plan_key, modality, is_in_vivo, is_human in (
        ("in_vitro_plan", "in_vitro", False, False),
        ("in_vivo_plan", "in_vivo", True, False),
        ("human_plan", "human", True, True),
    ):
        for item in plan.get(plan_key) or []:
            if not isinstance(item, dict) or item.get("execution_status") == "excluded":
                continue
            operational_issues = _dedupe_preserve(
                _as_str_list(item.get("operationalization_issues"))
                + _candidate_experiment_completion_issues(item, is_in_vivo=is_in_vivo, is_human=is_human)
            )
            decision_issues = _dedupe_preserve(
                _candidate_experiment_decision_issues(
                    item,
                    is_in_vivo=is_in_vivo,
                    is_human=is_human,
                )
                + (
                    _candidate_role_strategy_issues(item)
                    if strict_candidate_design and not is_human
                    else []
                )
            )
            if not operational_issues:
                design_status = "execution_complete"
            elif not decision_issues:
                design_status = "decision_complete"
            else:
                design_status = "needs_resolution"
            item["completion_issues"] = decision_issues
            item["operationalization_issues"] = operational_issues
            item["design_issues"] = decision_issues
            item["design_status"] = design_status
            item["generation_status"] = "complete" if design_status != "needs_resolution" else "degraded"
            if item.get("execution_status") == "ready_now" and design_status == "needs_resolution":
                item["execution_status"] = "conditional_future"
            module_id = str(item.get("module_id") or "unassigned").strip()
            audit_issues = decision_issues
            module_results.append(
                {
                    "module_id": module_id,
                    "modality": modality,
                    "execution_status": str(item.get("execution_status") or "conditional_future"),
                    "design_status": design_status,
                    "status": "pass" if not audit_issues else "needs_resolution",
                    "issues": audit_issues,
                    "deferred_operationalization_issues": (
                        operational_issues if item.get("execution_status") != "ready_now" else []
                    ),
                }
            )
            all_issues.extend(
                {"module_id": module_id, "issue": issue}
                for issue in audit_issues
            )
    return {
        "status": "not_applicable" if not module_results else ("needs_resolution" if all_issues else "pass"),
        "module_results": module_results,
        "issues": all_issues,
    }


def _refresh_plan_level_structure(plan: dict, candidate_workflow: Optional[bool] = None) -> dict:
    """Rebuild execution views and audits after any module-level mutation."""
    is_candidate = _plan_is_candidate_workflow(plan) if candidate_workflow is None else bool(candidate_workflow)
    for key in ("in_vitro_plan", "in_vivo_plan", "human_plan"):
        if not isinstance(plan.get(key), list):
            plan[key] = []
    if _human_study_excluded_by_user(plan):
        for item in plan.get("human_plan") or []:
            if isinstance(item, dict):
                item["execution_status"] = "excluded"

    _ensure_unique_module_ids(plan)
    plan["design_completeness_audit"] = _design_completeness_audit(
        plan, strict_candidate_design=is_candidate
    )
    _apply_candidate_dependency_gates(plan)
    plan["design_completeness_audit"] = _design_completeness_audit(
        plan, strict_candidate_design=is_candidate
    )
    current = {"in_vitro": [], "in_vivo": [], "human": []}
    future = {"in_vitro": [], "in_vivo": [], "human": []}
    selected_modules: List[dict] = []
    omitted_roles: List[dict] = []
    decision_graph: List[dict] = []
    selected_roles: set[str] = set()

    for plan_key, modality in (
        ("in_vitro_plan", "in_vitro"),
        ("in_vivo_plan", "in_vivo"),
        ("human_plan", "human"),
    ):
        for item in plan.get(plan_key) or []:
            if not isinstance(item, dict):
                continue
            status = str(item.get("execution_status") or "conditional_future").strip().lower()
            if status not in VALID_EXECUTION_STATUSES:
                status = "conditional_future"
                item["execution_status"] = status
            role = _canonicalize_experiment_role(item.get("experiment_role"))
            item["experiment_role"] = role
            module_id = str(item.get("module_id") or "unassigned").strip()
            if status == "excluded":
                omitted_roles.append(
                    {
                        "experiment_role": role,
                        "reason": item.get("prerequisite_result") or "The returned module was explicitly classified as excluded.",
                    }
                )
                continue
            selected_roles.add(role)
            target = current if status == "ready_now" else future
            target[modality].append(item)
            selection_reason = (
                "Prerequisites are currently satisfiable and deterministic completeness checks passed."
                if status == "ready_now"
                else "; ".join(_as_str_list(item.get("activation_gate_reasons")))
                or item.get("prerequisite_result")
                or "Retained for a later decision branch; its prerequisites or design details are not yet satisfied."
            )
            selected_modules.append(
                {
                    "module_id": module_id,
                    "experiment_role": role,
                    "execution_status": status,
                    "result_status": _normalize_result_status(item.get("result_status")),
                    "design_status": _normalize_design_status(item.get("design_status")) or "needs_resolution",
                    "activation_gate_status": str(item.get("activation_gate_status") or "not_evaluated"),
                    "selection_reason": str(selection_reason).strip(),
                    "prerequisite_module_ids": _as_str_list(item.get("prerequisite_module_ids")),
                    "activation_gate": _normalize_activation_gate(item.get("activation_gate")),
                }
            )
            decision_graph.append(
                {
                    "module_id": module_id,
                    "execution_status": status,
                    "result_status": _normalize_result_status(item.get("result_status")),
                    "design_status": _normalize_design_status(item.get("design_status")) or "needs_resolution",
                    "activation_gate_status": str(item.get("activation_gate_status") or "not_evaluated"),
                    "prerequisite_module_ids": _as_str_list(item.get("prerequisite_module_ids")),
                    "activation_gate": _normalize_activation_gate(item.get("activation_gate")),
                    "branch_if_positive": str(item.get("branch_if_positive") or "").strip(),
                    "branch_if_negative": str(item.get("branch_if_negative") or "").strip(),
                }
            )

    human_gate = _human_gate_for_plan(plan, is_candidate)
    if is_candidate:
        expected_roles = (
            "direct_production_gate",
            "indirect_ecology_gate",
            "host_response_gate",
            "effect_interaction_gate",
            "direct_causal_rescue",
            "indirect_ecology_causal",
        )
        for role in expected_roles:
            if role not in selected_roles and not any(item.get("experiment_role") == role for item in omitted_roles):
                omitted_roles.append(
                    {
                        "experiment_role": role,
                        "reason": "Not selected in the smallest sufficient current workflow; add it only if the upstream decision branch requires it.",
                    }
                )
        if "human_translation" not in selected_roles and not any(
            item.get("experiment_role") == "human_translation" for item in omitted_roles
        ):
            omitted_roles.append(
                {
                    "experiment_role": "human_translation",
                    "reason": human_gate.get("reason") or "Human translation remains behind the preclinical gate.",
                }
            )
        human_activation_gate = _default_activation_gate("human_translation")
        decision_graph.append(
            {
                "module_id": "H_GATE",
                "execution_status": "excluded" if human_gate.get("status") == "excluded_by_user" else "conditional_future",
                "result_status": "not_run",
                "activation_gate": human_activation_gate,
                "activation_gate_status": (
                    "met" if human_gate.get("status") == "ready_now" else "unmet"
                ),
                "prerequisite_module_ids": _activation_gate_module_references(human_activation_gate),
                "branch_if_positive": "Activate one branch-specific observational human protocol after all gates are documented as positive.",
                "branch_if_negative": "Retain the human blueprint as conditional_future and do not make a causal human interpretation.",
            }
        )

    assessment = _normalize_production_evidence_assessment(plan.get("direct_production_evidence_assessment"))
    limitations = _as_str_list(plan.get("evidence_limitations")) + _as_str_list(assessment.get("evidence_limitations"))
    audit = plan.get("protocol_audit") if isinstance(plan.get("protocol_audit"), dict) else {}
    for bucket, prefix in (
        ("inference_only_claims", "Inference only"),
        ("unsupported_or_overstated_claims", "Unsupported or overstated"),
    ):
        for item in audit.get(bucket) or []:
            if isinstance(item, dict) and str(item.get("claim") or "").strip():
                limitations.append(f"{prefix}: {item.get('claim')}")
    limitations.extend(_as_str_list(plan.get("overall_risk_flags")))

    uncertainties = _as_str_list(plan.get("remaining_uncertainties"))
    for item in plan.get("self_reflection") or []:
        if isinstance(item, dict) and str(item.get("remaining_uncertainty") or "").strip():
            uncertainties.append(str(item.get("remaining_uncertainty")).strip())
    for item in _normalize_hypothesis_branches(plan.get("hypothesis_branches")):
        status = str(item.get("current_evidence_status") or "unresolved").lower()
        if status not in {"supported", "direct_supported", "supported_by_current_culture_evidence"} and item.get("statement"):
            uncertainties.append(str(item.get("statement")).strip())

    plan["module_selection"] = {
        "workflow_type": "microbe_metabolite_disease" if is_candidate else str(plan.get("mode") or "general_validation"),
        "selected_modules": selected_modules,
        "omitted_roles": _dedupe_plan_records(omitted_roles, ("experiment_role", "reason")),
    }
    plan["current_executable_plan"] = current
    plan["conditional_future_plan"] = future
    plan["modules"] = [
        item
        for plan_key in ("in_vitro_plan", "in_vivo_plan", "human_plan")
        for item in (plan.get(plan_key) or [])
        if isinstance(item, dict)
    ]
    plan["human_gate"] = human_gate
    plan["decision_graph"] = decision_graph
    plan["parameter_provenance_audit"] = _parameter_provenance_audit(plan)
    plan["evidence_limitations"] = _dedupe_preserve(limitations)[:24]
    plan["remaining_uncertainties"] = _dedupe_preserve(uncertainties)[:24]
    return plan


def _project_public_experiment_module(item: Any) -> dict:
    source = item if isinstance(item, dict) else {}
    role = _canonicalize_experiment_role(source.get("experiment_role"))
    groups = _normalize_group_definitions(source.get("groups"))
    gate_contract = _canonical_gate_contract(role, source.get("activation_gate"))
    projected = {
        "module_id": str(source.get("module_id") or "").strip(),
        "experiment_role": role,
        "hypothesis_ids": _as_str_list(source.get("hypothesis_ids")),
        "route_status": _normalize_route_status(source.get("route_status"), role),
        "execution_status": str(source.get("execution_status") or "conditional_future").strip().lower(),
        "result_status": _normalize_result_status(source.get("result_status")),
        "design_status": _normalize_design_status(source.get("design_status")) or "needs_resolution",
        "activation_gate": gate_contract["activation_gate"],
        "scientific_question": str(source.get("scientific_question") or "").strip(),
        "why_needed": str(source.get("why_needed") or "").strip(),
        "study_object": str(source.get("study_object") or "").strip(),
        "groups": groups,
        "group_count": len(groups),
        "samples_and_timing": str(source.get("samples_and_timing") or "").strip(),
        "primary_indicator": str(source.get("primary_indicator") or "").strip(),
        "secondary_indicators": _as_str_list(source.get("secondary_indicators")),
        "key_controls": _as_str_list(source.get("key_controls")),
        "positive_gate": str(source.get("positive_gate") or "").strip(),
        "branch_if_positive": str(source.get("branch_if_positive") or "").strip(),
        "branch_if_negative": str(source.get("branch_if_negative") or "").strip(),
        "unlock_rule": gate_contract["unlock_rule"],
        "claim_boundary": str(source.get("claim_boundary") or "").strip(),
        "completion_issues": _as_str_list(source.get("completion_issues")),
    }
    return {field: projected[field] for field in PUBLIC_EXPERIMENT_MODULE_FIELDS}


def _project_public_plan(plan: dict) -> dict:
    projected = dict(plan)
    for plan_key in ("in_vitro_plan", "in_vivo_plan", "human_plan"):
        projected[plan_key] = [
            _project_public_experiment_module(item)
            for item in (plan.get(plan_key) or [])
            if isinstance(item, dict)
        ]
    for view_key in ("current_executable_plan", "conditional_future_plan"):
        source_view = plan.get(view_key) if isinstance(plan.get(view_key), dict) else {}
        projected[view_key] = {
            modality: [
                _project_public_experiment_module(item)
                for item in (source_view.get(modality) or [])
                if isinstance(item, dict)
            ]
            for modality in ("in_vitro", "in_vivo", "human")
        }
    projected["modules"] = (
        projected["in_vitro_plan"]
        + projected["in_vivo_plan"]
        + projected["human_plan"]
    )
    projected.pop("parameter_provenance_audit", None)
    projected.pop("design_completeness_audit", None)
    return projected


def _normalize_reflection_items(items: Any) -> List[dict]:
    normalized = []
    if not isinstance(items, list):
        return normalized
    for item in items[:10]:
        if not isinstance(item, dict):
            continue
        normalized.append(
            {
                "initial_claim": str(item.get("initial_claim") or "").strip(),
                "self_critique": str(item.get("self_critique") or "").strip(),
                "evidence_checked": _as_str_list(item.get("evidence_checked")),
                "revision": str(item.get("revision") or "").strip(),
                "remaining_uncertainty": str(item.get("remaining_uncertainty") or "").strip(),
            }
        )
    return normalized


def _normalize_protocol_audit(audit: Any) -> dict:
    if not isinstance(audit, dict):
        return {
            "overall_assessment": "",
            "supported_claims": [],
            "inference_only_claims": [],
            "unsupported_or_overstated_claims": [],
            "priority_gaps": [],
        }
    return {
        "overall_assessment": str(audit.get("overall_assessment") or "").strip(),
        "supported_claims": _normalize_audit_claims(audit.get("supported_claims")),
        "inference_only_claims": _normalize_audit_claims(audit.get("inference_only_claims"), include_reason=True),
        "unsupported_or_overstated_claims": _normalize_audit_claims(
            audit.get("unsupported_or_overstated_claims"),
            include_reason=True,
        ),
        "priority_gaps": _as_str_list(audit.get("priority_gaps")),
    }


def _normalize_query_log(items: Any) -> List[dict]:
    normalized = []
    if not isinstance(items, list):
        return normalized
    for item in items[:12]:
        if not isinstance(item, dict):
            continue
        normalized.append(
            {
                "round": int(item.get("round") or 0),
                "label": str(item.get("label") or "").strip(),
                "query": str(item.get("query") or "").strip(),
                "rationale": str(item.get("rationale") or "").strip(),
                "pmid_count": int(item.get("pmid_count") or 0),
                "article_count": int(item.get("article_count") or 0),
            }
        )
    return normalized


def _normalize_evidence_strength_map(items: Any) -> List[dict]:
    normalized = []
    if not isinstance(items, list):
        return normalized
    for item in items[:10]:
        if not isinstance(item, dict):
            continue
        normalized.append(
            {
                "topic": str(item.get("topic") or "").strip(),
                "support_level": str(item.get("support_level") or "speculative").strip().lower(),
                "evidence_type": str(item.get("evidence_type") or "unspecified").strip(),
                "rationale": str(item.get("rationale") or "").strip(),
                "representative_pmids": _as_str_list(item.get("representative_pmids")),
            }
        )
    return normalized


def _export_articles(articles: List[dict]) -> List[dict]:
    exported = articles_to_export(articles)
    for idx, article in enumerate(articles[: len(exported)]):
        if idx >= len(exported):
            break
        exported[idx]["query_round"] = article.get("query_round")
        exported[idx]["query_label"] = article.get("query_label")
    return exported


def _empty_plan(
    bacteria: str,
    metabolite: str,
    disease: str,
    candidate_metrics: Dict[str, Any],
    literature: List[dict],
    mode: str,
    research_question: str,
    prompt_constraints: str,
) -> dict:
    return {
        "candidate": {
            "bacteria": bacteria,
            "metabolite": metabolite,
            "disease": disease,
        },
        "candidate_metrics": candidate_metrics,
        "mode": mode,
        "user_brief": {
            "research_question": (research_question or "").strip(),
            "prompt_constraints": (prompt_constraints or "").strip(),
        },
        "question_profile": {},
        "protocol_audit": {
            "overall_assessment": "",
            "supported_claims": [],
            "inference_only_claims": [],
            "unsupported_or_overstated_claims": [],
            "priority_gaps": [],
        },
        "working_hypothesis": {
            "statement": "",
            "direct_support": [],
            "inference_only": [],
        },
        "direct_production_evidence_assessment": {
            "status": "not_assessed",
            "conclusion": "",
            "direct_evidence": [],
            "paper_findings": [],
            "evidence_limitations": [],
        },
        "hypothesis_branches": [],
        "evidence_basis": [],
        "iterative_query_log": [],
        "evidence_strength_map": [],
        "fulltext_method_evidence": {
            "in_vitro": [],
            "in_vivo": [],
        },
        "protocols_io_evidence": {
            "enabled": False,
            "status": "disabled",
            "message": "",
            "queries": [],
            "in_vitro": [],
            "in_vivo": [],
        },
        "validation_summary": {
            "known_evidence": [],
            "direct_production_evidence": [],
            "hypothesis_branches": [],
            "missing_evidence_gap": [],
            "why_this_experiment": [],
            "priority_follow_up_experiments": [],
            "decision_rule": [],
            "future_complete_study": [],
        },
        "validation_protocol_text": "",
        "in_vitro_plan": [],
        "in_vivo_plan": [],
        "human_plan": [],
        "modules": [],
        "external_gate_results": {},
        "module_selection": {
            "workflow_type": "",
            "selected_modules": [],
            "omitted_roles": [],
        },
        "current_executable_plan": {
            "in_vitro": [],
            "in_vivo": [],
            "human": [],
        },
        "conditional_future_plan": {
            "in_vitro": [],
            "in_vivo": [],
            "human": [],
        },
        "human_gate": {
            "status": "not_assessed",
            "reason": "",
            "future_requirements": [],
        },
        "decision_graph": [],
        "parameter_provenance_audit": {
            "status": "not_applicable",
            "module_results": [],
            "unresolved_items": [],
        },
        "design_completeness_audit": {
            "status": "not_applicable",
            "module_results": [],
            "issues": [],
        },
        "evidence_limitations": [],
        "remaining_uncertainties": [],
        "self_reflection": [],
        "overall_risk_flags": [],
        "retrieved_literature": _export_articles(literature),
    }


def _normalize_plan(
    raw: Dict[str, Any],
    bacteria: str,
    metabolite: str,
    disease: str,
    candidate_metrics: Dict[str, Any],
    literature: List[dict],
    mode: str,
    iterative_query_log: List[dict],
    evidence_strength_map: List[dict],
    fulltext_method_evidence: Dict[str, List[dict]],
    protocols_io_evidence: Dict[str, Any],
    research_question: str,
    prompt_constraints: str,
) -> dict:
    plan = _empty_plan(
        bacteria,
        metabolite,
        disease,
        candidate_metrics,
        literature,
        mode,
        research_question,
        prompt_constraints,
    )
    if isinstance(raw, dict):
        plan["validation_protocol_text"] = str(raw.get("validation_protocol_text") or "").strip()
        plan["protocol_audit"] = _normalize_protocol_audit(raw.get("protocol_audit"))

        hypothesis = raw.get("working_hypothesis") or {}
        if isinstance(hypothesis, dict):
            plan["working_hypothesis"] = {
                "statement": str(hypothesis.get("statement") or "").strip(),
                "direct_support": _as_str_list(hypothesis.get("direct_support")),
                "inference_only": _as_str_list(hypothesis.get("inference_only")),
            }

        plan["direct_production_evidence_assessment"] = _normalize_production_evidence_assessment(
            raw.get("direct_production_evidence_assessment")
        )
        plan["hypothesis_branches"] = _normalize_hypothesis_branches(raw.get("hypothesis_branches"))
        plan["evidence_basis"] = _normalize_evidence_items(raw.get("evidence_basis"))
        plan["in_vitro_plan"] = _prepare_experiment_modules(raw.get("in_vitro_plan"))
        plan["in_vivo_plan"] = _prepare_experiment_modules(raw.get("in_vivo_plan"), is_in_vivo=True)
        plan["human_plan"] = _prepare_experiment_modules(raw.get("human_plan"), is_in_vivo=True, is_human=True)
        plan["external_gate_results"] = (
            dict(raw.get("external_gate_results"))
            if isinstance(raw.get("external_gate_results"), dict)
            else {}
        )
        plan["evidence_limitations"] = _as_str_list(raw.get("evidence_limitations"))
        plan["remaining_uncertainties"] = _as_str_list(raw.get("remaining_uncertainties"))
        plan["self_reflection"] = _normalize_reflection_items(raw.get("self_reflection"))
        plan["overall_risk_flags"] = _as_str_list(raw.get("overall_risk_flags"))

    plan["iterative_query_log"] = _normalize_query_log(iterative_query_log)
    plan["evidence_strength_map"] = _normalize_evidence_strength_map(evidence_strength_map)
    plan["fulltext_method_evidence"] = {
        "in_vitro": fulltext_method_evidence.get("in_vitro", [])[:8] if isinstance(fulltext_method_evidence, dict) else [],
        "in_vivo": fulltext_method_evidence.get("in_vivo", [])[:8] if isinstance(fulltext_method_evidence, dict) else [],
    }
    plan["protocols_io_evidence"] = protocols_io_evidence if isinstance(protocols_io_evidence, dict) else {}
    plan = _refresh_plan_level_structure(plan)
    plan["validation_summary"] = _build_validation_summary(plan)
    if not plan["validation_protocol_text"]:
        plan["validation_protocol_text"] = _render_validation_protocol_text(plan)
    return plan


def _strip_reference_noise(text: str) -> str:
    clean = str(text or "")
    clean = re.sub(r"\[PMID:[^\]]+\]", "", clean, flags=re.IGNORECASE)
    clean = re.sub(r"\bPMCID?:\s*\S+", "", clean, flags=re.IGNORECASE)
    clean = re.sub(r"https?://\S+", "", clean)
    clean = re.sub(r"protocols\.io\S*", "", clean, flags=re.IGNORECASE)
    clean = re.sub(r"\s+", " ", clean).strip(" ;|,")
    return clean


def _contains_cjk(text: str) -> bool:
    return bool(re.search(r"[\u4e00-\u9fff]", str(text or "")))


def _question_driven_mode(plan: dict) -> bool:
    return str(plan.get("mode") or "").strip() == "question_driven"


def _normalize_question_profile(raw: Any) -> dict:
    if not isinstance(raw, dict):
        return {
            "question_type": "general_validation",
            "human_feasibility": "uncertain",
            "human_ethics": "uncertain",
            "recommended_focus": "preclinical_first",
            "human_study_policy": "not_recommended",
            "detail_budget_allocation": "balanced_preclinical",
            "reasoning": "",
            "candidate_microbe": "",
            "candidate_metabolite": "",
            "candidate_disease": "",
        }
    return {
        "question_type": str(raw.get("question_type") or "general_validation").strip(),
        "human_feasibility": str(raw.get("human_feasibility") or "uncertain").strip(),
        "human_ethics": str(raw.get("human_ethics") or "uncertain").strip(),
        "recommended_focus": str(raw.get("recommended_focus") or "preclinical_first").strip(),
        "human_study_policy": str(raw.get("human_study_policy") or "not_recommended").strip(),
        "detail_budget_allocation": str(raw.get("detail_budget_allocation") or "balanced_preclinical").strip(),
        "reasoning": str(raw.get("reasoning") or "").strip(),
        "candidate_microbe": str(raw.get("candidate_microbe") or "").strip(),
        "candidate_metabolite": str(raw.get("candidate_metabolite") or "").strip(),
        "candidate_disease": str(raw.get("candidate_disease") or "").strip(),
    }


def _format_question_profile(profile: dict) -> str:
    data = _normalize_question_profile(profile)
    return "\n".join(
        [
            f"- question_type: {data.get('question_type')}",
            f"- human_feasibility: {data.get('human_feasibility')}",
            f"- human_ethics: {data.get('human_ethics')}",
            f"- recommended_focus: {data.get('recommended_focus')}",
            f"- human_study_policy: {data.get('human_study_policy')}",
            f"- detail_budget_allocation: {data.get('detail_budget_allocation')}",
            f"- reasoning: {data.get('reasoning') or 'Not provided'}",
            f"- candidate_microbe: {data.get('candidate_microbe') or 'Not resolved'}",
            f"- candidate_metabolite: {data.get('candidate_metabolite') or 'Not resolved'}",
            f"- candidate_disease: {data.get('candidate_disease') or 'Not resolved'}",
        ]
    )


def classify_question_profile(research_question: str, prompt_constraints: str, disease: str = "") -> dict:
    system_prompt = """You are classifying a standalone experimental planning question so later agents can allocate detail intelligently.

Return JSON with exactly these keys:
- question_type
- human_feasibility
- human_ethics
- recommended_focus
- human_study_policy
- detail_budget_allocation
- reasoning
- candidate_microbe
- candidate_metabolite
- candidate_disease

Allowed values:
- question_type: sensor_discovery, toxic_exposure, beneficial_intervention, active_compound_screening, behavioral_or_observational, general_validation
- human_feasibility: high, medium, low, none
- human_ethics: appropriate, conditional, inappropriate
- recommended_focus: discovery_first, preclinical_first, translational_balanced, human_late_only
- human_study_policy: not_recommended, allow_observational_only, allow_low_risk_behavior_only, allow_late_stage_rct, forbid_direct_exposure
- detail_budget_allocation: discovery_heavy, cell_heavy, animal_heavy, balanced_preclinical, human_late, human_behavior_heavy, observational_human_heavy

Requirements:
- Judge the scientific objective, not just isolated keywords.
- Questions about discovering a sensor protein, receptor, target, or molecular detector should usually be sensor_discovery and discovery_first.
- Questions involving toxicants, hazardous exposures, harmful environmental agents, or unethical direct exposure should usually be toxic_exposure. Direct human exposure should usually be forbidden. Use allow_observational_only only if the human component would be non-interventional exposure assessment or biomonitoring.
- Low-risk behavior or preference questions may use allow_low_risk_behavior_only.
- Beneficial interventions, repurposed drugs, or plausibly safe translational questions may use allow_late_stage_rct, but only when preclinical evidence can reasonably support that escalation.
- Human studies should be recommended only when ethically and practically plausible.
- The detail_budget_allocation should indicate where the answer should spend most of its depth, rather than mechanically adding every stage.
- Keep reasoning short and concrete."""
    system_prompt += """
- If the question explicitly names a microbe, metabolite, or disease, copy each entity into candidate_microbe, candidate_metabolite, and candidate_disease respectively. Use an empty string when an entity is not explicit. Do not substitute a broader class, infer a preferred disease subtype, or return explanatory text in these three fields."""
    user_prompt = f"""Research question:
{research_question or 'Not provided'}

Prompt constraints:
{prompt_constraints or 'Not provided'}

Disease scope:
{disease or 'Not specified'}

Classify the question for downstream planning."""
    try:
        raw = _chat_json(system_prompt, user_prompt, max_tokens=500, temperature=0.1)
        return _normalize_question_profile(raw)
    except Exception:
        text = f"{research_question} {prompt_constraints}".lower()
        if any(term in text for term in ("sensor", "receptor", "odorant-binding protein", "obp", "target identification", "ligand-binding", "binding protein")):
            return _normalize_question_profile(
                {
                    "question_type": "sensor_discovery",
                    "human_feasibility": "none",
                    "human_ethics": "conditional",
                    "recommended_focus": "discovery_first",
                    "human_study_policy": "not_recommended",
                    "detail_budget_allocation": "discovery_heavy",
                    "reasoning": "Target-identification questions should prioritize discovery and functional validation before translation.",
                }
            )
        if any(term in text for term in ("toxic", "toxicant", "poison", "pollutant", "pesticide", "heavy metal", "hazardous", "deltamethrin", "insecticide")):
            return _normalize_question_profile(
                {
                    "question_type": "toxic_exposure",
                    "human_feasibility": "low",
                    "human_ethics": "inappropriate",
                    "recommended_focus": "preclinical_first",
                    "human_study_policy": "forbid_direct_exposure",
                    "detail_budget_allocation": "animal_heavy",
                    "reasoning": "Direct human exposure is ethically problematic for harmful agents, so causal depth should stay preclinical.",
                }
            )
        if any(term in text for term in ("mosquito", "repellent", "clothing color", "colour preference", "color preference", "host-seeking behavior", "landing", "bite avoidance")):
            return _normalize_question_profile(
                {
                    "question_type": "behavioral_or_observational",
                    "human_feasibility": "medium",
                    "human_ethics": "appropriate",
                    "recommended_focus": "translational_balanced",
                    "human_study_policy": "allow_low_risk_behavior_only",
                    "detail_budget_allocation": "human_behavior_heavy",
                    "reasoning": "Low-risk mosquito behavior questions can include controlled human behavioral testing without relying on invasive clinical translation.",
                }
            )
        if any(term in text for term in ("elderly", "adolescent", "blood", "cohort", "screen senescence-associated proteins", "proteomic", "proteomics", "population")):
            return _normalize_question_profile(
                {
                    "question_type": "behavioral_or_observational",
                    "human_feasibility": "high",
                    "human_ethics": "appropriate",
                    "recommended_focus": "translational_balanced",
                    "human_study_policy": "allow_observational_only",
                    "detail_budget_allocation": "observational_human_heavy",
                    "reasoning": "Population biospecimen screening is best handled as an observational human study rather than an intervention trial.",
                }
            )
        if any(term in text for term in ("plant", "extract", "fraction", "active compound", "screen active", "screening active", "natural product")):
            return _normalize_question_profile(
                {
                    "question_type": "active_compound_screening",
                    "human_feasibility": "medium",
                    "human_ethics": "conditional",
                    "recommended_focus": "preclinical_first",
                    "human_study_policy": "allow_low_risk_behavior_only",
                    "detail_budget_allocation": "balanced_preclinical",
                    "reasoning": "Active-compound screening should stay preclinical first, with only low-risk behavioral human confirmation if the intervention becomes sufficiently standardized.",
                }
            )
        return _normalize_question_profile({})


def _extract_question_candidate_entities(
    research_question: str,
    prompt_constraints: str,
    disease: str,
    question_profile: Any,
) -> dict:
    profile = _normalize_question_profile(question_profile)
    raw_text = f"{research_question} {prompt_constraints}".strip()
    lower_text = raw_text.lower()

    def clean(value: Any) -> str:
        candidate = str(value or "").strip()
        return "" if candidate.lower() in {"none", "not specified", "not resolved", "unknown", "n/a"} else candidate

    microbe = clean(profile.get("candidate_microbe"))
    metabolite = clean(profile.get("candidate_metabolite"))
    resolved_disease = clean(disease) or clean(profile.get("candidate_disease"))

    if not microbe:
        if re.search(r"(?<![a-z0-9])akk(?:ermansia(?: muciniphila)?)?(?=$|[^a-z0-9])", lower_text):
            microbe = "Akkermansia muciniphila"
        else:
            binomial_matches = re.findall(
                r"\b([A-Z][a-z]{2,}(?:\s+|_)[a-z][a-z0-9-]{2,})\b",
                raw_text,
            )
            excluded_first_words = {
                "assess", "build", "compare", "determine", "does", "evaluate",
                "how", "investigate", "study", "test", "the", "what", "whether",
            }
            microbe = next(
                (
                    match.replace("_", " ")
                    for match in binomial_matches
                    if match.split()[0].lower() not in excluded_first_words
                ),
                "",
            )

    if not metabolite:
        metabolite_aliases = (
            ("isobutyric acid", ("isobutyric acid", "isobutyrate", "\u5f02\u4e01\u9178")),
            ("butyric acid", ("butyric acid", "butyrate", "\u4e01\u9178")),
            ("propionic acid", ("propionic acid", "propionate", "\u4e19\u9178")),
            ("acetic acid", ("acetic acid", "acetate", "\u4e59\u9178")),
            ("succinic acid", ("succinic acid", "succinate", "\u7425\u73c0\u9178")),
            ("indole", ("indole", "\u5432\u54da")),
            ("bile acid", ("bile acid", "bile acids", "\u80c6\u6c41\u9178")),
        )
        metabolite = next(
            (
                canonical
                for canonical, aliases in metabolite_aliases
                if any(alias in lower_text or alias in raw_text for alias in aliases)
            ),
            "",
        )

    if not resolved_disease:
        mentions_crohn = bool(re.search(r"\bcrohn(?:'s)?(?: disease)?\b|\bcd\b", lower_text))
        mentions_uc = bool("ulcerative colitis" in lower_text or re.search(r"\buc\b", lower_text))
        if mentions_crohn and mentions_uc:
            resolved_disease = "IBD"
        elif mentions_crohn:
            resolved_disease = "Crohn disease"
        elif mentions_uc:
            resolved_disease = "ulcerative colitis"
        elif "ibd" in lower_text or "inflammatory bowel disease" in lower_text:
            resolved_disease = "IBD"

    return {
        "bacteria": microbe,
        "metabolite": metabolite,
        "disease": resolved_disease,
    }


def _list_to_lines(items: List[str], prefix: str = "- ", limit: int = 8) -> List[str]:
    lines = []
    for item in _dedupe_preserve([_strip_reference_noise(v) for v in items])[:limit]:
        if item:
            lines.append(f"{prefix}{item}")
    return lines


def _render_group_lines(item: dict) -> List[str]:
    groups = item.get("groups") if isinstance(item.get("groups"), list) else []
    lines: List[str] = []
    for index, group in enumerate(groups, start=1):
        if isinstance(group, dict):
            name = _strip_reference_noise(group.get("group_name")) or f"Group {index}"
            exposure = _strip_reference_noise(group.get("exposure_or_condition"))
            purpose = _strip_reference_noise(group.get("control_purpose"))
            details = []
            if exposure:
                details.append(exposure)
            if purpose:
                details.append(f"purpose: {purpose}")
            lines.append(f"   - {name}" + (f": {'; '.join(details)}" if details else ""))
            continue
        clean = _strip_reference_noise(group)
        if clean:
            lines.append(f"   - {clean}")
    return lines


def _summarize_groups(item: dict, is_in_vivo: bool = False) -> str:
    group_lines = _render_group_lines(item)
    if group_lines:
        return f"Groups ({len(group_lines)}): " + "; ".join(line[5:] for line in group_lines)
    controls = _as_str_list(item.get("key_controls")) or _as_str_list(item.get("controls"))
    return "Controls: " + "; ".join(controls) if controls else ""


def _format_parameter_provenance_lines(item: dict) -> List[str]:
    provenance = _normalize_parameter_provenance(item.get("parameter_provenance"))
    if not provenance:
        return ["   Parameter provenance: Not provided; exact parameters require source resolution or a prespecified pilot."]

    lines = ["   Parameter provenance:"]
    for entry in provenance[:16]:
        parameter = str(entry.get("parameter") or "unspecified parameter").strip()
        value = str(entry.get("value") or "not specified").strip()
        status = str(entry.get("status") or "unresolved").strip()
        source = str(entry.get("source") or "NONE").strip()
        parts = [f"{parameter}={value}", f"status={status}", f"source={source}"]
        source_context = _strip_reference_noise(entry.get("source_context"))
        transfer_rationale = _strip_reference_noise(entry.get("transfer_rationale"))
        pilot_check = _strip_reference_noise(entry.get("pilot_check"))
        if source_context:
            parts.append(f"source_context={source_context}")
        if transfer_rationale:
            parts.append(f"transfer_rationale={transfer_rationale}")
        if pilot_check:
            parts.append(f"pilot_check={pilot_check}")
        lines.append("   - " + " | ".join(parts))
    return lines


def _format_experiment_block(item: dict, index: int, is_in_vivo: bool = False) -> List[str]:
    module_id = _strip_reference_noise(item.get("module_id")) or "unassigned"
    execution_status = str(item.get("execution_status") or "unclassified").strip().lower()
    prerequisite_module_ids = [
        _strip_reference_noise(value)
        for value in _as_str_list(item.get("prerequisite_module_ids"))
        if _strip_reference_noise(value)
    ]
    aim = _strip_reference_noise(item.get("aim") or ("Animal Experiment" if is_in_vivo else "In Vitro Experiment"))
    lines = [f"{index}. [{module_id}] {aim}"]
    experiment_role = _strip_reference_noise(item.get("experiment_role"))
    hypotheses = "; ".join(
        _strip_reference_noise(v) for v in _as_str_list(item.get("hypothesis_tested"))[:4] if _strip_reference_noise(v)
    )
    prerequisite = _strip_reference_noise(item.get("prerequisite_result"))
    stage_gate = _strip_reference_noise(item.get("stage_gate"))
    branch_if_positive = _strip_reference_noise(item.get("branch_if_positive"))
    branch_if_negative = _strip_reference_noise(item.get("branch_if_negative"))
    question = _strip_reference_noise(item.get("biological_question"))
    rationale = _strip_reference_noise(item.get("priority_rationale"))
    gap_addressed = _strip_reference_noise(item.get("gap_addressed"))
    model = _strip_reference_noise(item.get("model") if is_in_vivo else item.get("model_system"))
    why_this_model = _strip_reference_noise(item.get("why_this_model")) if is_in_vivo else ""
    material = _strip_reference_noise(item.get("experimental_material"))
    key_materials = "; ".join(_strip_reference_noise(v) for v in _as_str_list(item.get("key_materials_equipment")) if _strip_reference_noise(v))
    design = _strip_reference_noise(item.get("design"))
    group_logic = _strip_reference_noise(item.get("group_logic"))
    intervention = _strip_reference_noise(item.get("intervention"))
    route = _strip_reference_noise(item.get("intervention_route")) if is_in_vivo else ""
    dose_timing = _strip_reference_noise(item.get("dose_timing_logic"))
    timeline = _strip_reference_noise(item.get("timeline")) if is_in_vivo else ""
    analysis = _strip_reference_noise(item.get("data_analysis"))
    readouts = "; ".join(_strip_reference_noise(v) for v in _as_str_list(item.get("readouts")) if _strip_reference_noise(v))
    readout_rationale = _strip_reference_noise(item.get("readout_rationale"))
    positive_interpretation = _strip_reference_noise(item.get("positive_result_interpretation"))
    negative_interpretation = _strip_reference_noise(item.get("negative_result_interpretation"))
    confounders = "; ".join(_strip_reference_noise(v) for v in _as_str_list(item.get("key_confounders")) if _strip_reference_noise(v))
    evidence_basis = "; ".join(_strip_reference_noise(v) for v in _as_str_list(item.get("evidence_basis"))[:4] if _strip_reference_noise(v))
    reported_conditions = "; ".join(_strip_reference_noise(v) for v in _as_str_list(item.get("reported_conditions")) if _strip_reference_noise(v))
    source_citations = "; ".join(_strip_reference_noise(v) for v in _as_str_list(item.get("source_citations"))[:4] if _strip_reference_noise(v))
    primary_endpoints = "; ".join(_strip_reference_noise(v) for v in _as_str_list(item.get("primary_endpoints")) if _strip_reference_noise(v)) if is_in_vivo else ""
    secondary_endpoints = "; ".join(_strip_reference_noise(v) for v in _as_str_list(item.get("secondary_endpoints")) if _strip_reference_noise(v)) if is_in_vivo else ""
    mechanistic_endpoints = "; ".join(_strip_reference_noise(v) for v in _as_str_list(item.get("mechanistic_endpoints")) if _strip_reference_noise(v)) if is_in_vivo else ""
    experimental_unit = _strip_reference_noise(item.get("experimental_unit"))
    replication_and_sampling = _strip_reference_noise(item.get("replication_and_sampling"))
    sample_size_basis = _strip_reference_noise(item.get("sample_size_basis"))
    randomization_and_blinding = _strip_reference_noise(item.get("randomization_and_blinding"))
    safety_and_stopping_rules = _strip_reference_noise(item.get("safety_and_stopping_rules"))
    completion_issues = [
        _strip_reference_noise(value)
        for value in _as_str_list(item.get("completion_issues"))
        if _strip_reference_noise(value)
    ]

    lines.append(f"   Module ID: {module_id}")
    lines.append(f"   Execution status: {execution_status}")
    lines.append(
        "   Prerequisite module IDs: "
        + (", ".join(prerequisite_module_ids) if prerequisite_module_ids else "None")
    )
    if question:
        lines.append(f"   Objective: {question}")
    if experiment_role and experiment_role != "unspecified":
        lines.append(f"   Experiment role: {experiment_role}")
    if hypotheses:
        lines.append(f"   Hypothesis tested: {hypotheses}")
    if prerequisite:
        lines.append(f"   Prerequisite result: {prerequisite}")
    if stage_gate:
        lines.append(f"   Stage gate: {stage_gate}")
    if rationale:
        lines.append(f"   Why first: {rationale}")
    if gap_addressed:
        lines.append(f"   Gap addressed: {gap_addressed}")
    if model:
        lines.append(f"   Model: {model}")
    if why_this_model:
        lines.append(f"   Why this model: {why_this_model}")
    if material:
        lines.append(f"   Experimental material: {material}")
    if key_materials:
        lines.append(f"   Key materials / equipment: {key_materials}")
    group_line = _summarize_groups(item, is_in_vivo=is_in_vivo)
    if group_line:
        lines.append(f"   {group_line}")
    if group_logic:
        lines.append(f"   Group logic: {group_logic}")
    if intervention:
        lines.append(f"   Intervention: {intervention}")
    if route:
        lines.append(f"   Route: {route}")
    if dose_timing:
        lines.append(f"   Time / concentration: {dose_timing}")
    if timeline:
        lines.append(f"   Timeline: {timeline}")
    if design:
        lines.append(f"   Design: {design}")
    lines.append(f"   Experimental unit: {experimental_unit or 'Not specified.'}")
    lines.append(f"   Replication and sampling: {replication_and_sampling or 'Not specified.'}")
    lines.append(f"   Sample-size basis: {sample_size_basis or 'Not specified; treat as a pilot until justified.'}")
    lines.append(f"   Randomization and blinding: {randomization_and_blinding or 'Not specified.'}")
    procedure_steps = _list_to_lines(_as_str_list(item.get("procedure_steps")), prefix="   - ", limit=16)
    if procedure_steps:
        lines.append("   Procedure:")
        lines.extend(procedure_steps)
    if readouts:
        lines.append(f"   Readouts: {readouts}")
    if primary_endpoints:
        lines.append(f"   Primary endpoints: {primary_endpoints}")
    if secondary_endpoints:
        lines.append(f"   Secondary endpoints: {secondary_endpoints}")
    if mechanistic_endpoints:
        lines.append(f"   Mechanistic endpoints: {mechanistic_endpoints}")
    if readout_rationale:
        lines.append(f"   Readout rationale: {readout_rationale}")
    if evidence_basis:
        lines.append(f"   Evidence basis: {evidence_basis}")
    if reported_conditions:
        lines.append(f"   Literature-reported conditions: {reported_conditions}")
    if source_citations:
        lines.append(f"   Source citations: {source_citations}")
    lines.extend(_format_parameter_provenance_lines(item))
    lines.append(f"   Safety and stopping rules: {safety_and_stopping_rules or 'Not specified.'}")
    if positive_interpretation:
        lines.append(f"   If positive: {positive_interpretation}")
    if negative_interpretation:
        lines.append(f"   If negative: {negative_interpretation}")
    if branch_if_positive:
        lines.append(f"   Next branch if positive: {branch_if_positive}")
    if branch_if_negative:
        lines.append(f"   Next branch if negative: {branch_if_negative}")
    if confounders:
        lines.append(f"   Alternative explanations / confounders: {confounders}")
    if analysis:
        lines.append(f"   Data analysis: {analysis}")
    if completion_issues:
        lines.append("   Completion issues:")
        lines.extend(f"   - {issue}" for issue in completion_issues[:16])
    else:
        lines.append("   Completion issues: None identified by the current completeness audit.")
    return lines


def _compact_field_text(value: Any) -> str:
    if isinstance(value, list):
        parts = []
        for entry in value:
            if isinstance(entry, dict):
                continue
            clean = _strip_reference_noise(entry)
            if clean:
                parts.append(clean)
        return "; ".join(parts)
    return _strip_reference_noise(value)


def _format_concise_experiment_block(item: dict, index: int, is_in_vivo: bool = False) -> List[str]:
    module_id = _strip_reference_noise(item.get("module_id")) or "unassigned"
    role = _strip_reference_noise(item.get("experiment_role")) or "unspecified"
    hypothesis_ids = _dedupe_preserve(
        hypothesis_id.upper()
        for hypothesis_id in _as_str_list(item.get("hypothesis_ids"))
        if re.fullmatch(r"H\d+", hypothesis_id.strip(), flags=re.IGNORECASE)
    )
    question = _strip_reference_noise(
        item.get("scientific_question") or item.get("biological_question") or item.get("aim")
    ) or "Scientific question not specified."
    route_status = _strip_reference_noise(item.get("route_status")) or "deferred"
    execution_status = _render_execution_status(item)
    result_status = _strip_reference_noise(item.get("result_status")) or "not_run"
    why_needed = _strip_reference_noise(
        item.get("why_needed") or item.get("priority_rationale") or item.get("gap_addressed")
    )
    study_object = _strip_reference_noise(
        item.get("study_object")
        or item.get("model" if is_in_vivo else "model_system")
        or item.get("experimental_material")
    )
    samples_and_timing = _compact_field_text(
        item.get("samples_and_timing") or item.get("replication_and_sampling") or item.get("timeline")
    )
    primary_indicator = _strip_reference_noise(item.get("primary_indicator") or item.get("primary_endpoint"))
    if not primary_indicator:
        primary_indicator = _compact_field_text(item.get("primary_endpoints") or item.get("readouts"))
    secondary_indicators = _compact_field_text(
        item.get("secondary_indicators")
        or item.get("secondary_endpoints")
        or item.get("mechanistic_endpoints")
    )
    key_controls = _compact_field_text(item.get("key_controls") or item.get("controls"))
    positive_gate = _strip_reference_noise(
        item.get("positive_gate") or item.get("success_threshold") or item.get("go_no_go")
    )
    unlock_rule = _strip_reference_noise(
        item.get("unlock_rule") or item.get("prerequisite_result") or item.get("stage_gate")
    )
    branch_if_positive = _strip_reference_noise(
        item.get("branch_if_positive") or item.get("positive_result_interpretation")
    )
    branch_if_negative = _strip_reference_noise(
        item.get("branch_if_negative") or item.get("negative_result_interpretation")
    )
    claim_boundary = _strip_reference_noise(
        item.get("claim_boundary") or item.get("phenomenology_support_criterion")
    )

    lines = [f"{index}. [{module_id}] {question}"]
    lines.append(
        f"   Role: {role} | Route: {route_status} | Status: {execution_status} | Result: {result_status}"
    )
    if hypothesis_ids:
        lines.append(f"   Hypothesis tested: {', '.join(hypothesis_ids)}")
    if why_needed:
        lines.append(f"   Why needed: {why_needed}")
    if study_object:
        lines.append(f"   Study object: {study_object}")
    group_lines = _render_group_lines(item)
    if group_lines:
        lines.append(f"   Groups ({len(group_lines)}):")
        lines.extend(group_lines)
    if samples_and_timing:
        lines.append(f"   Samples and timing: {samples_and_timing}")
    if primary_indicator:
        lines.append(f"   Primary indicator: {primary_indicator}")
    if secondary_indicators:
        lines.append(f"   Secondary indicators: {secondary_indicators}")
    if key_controls:
        lines.append(f"   Key controls: {key_controls}")
    if positive_gate:
        lines.append(f"   Positive gate: {positive_gate}")
    if unlock_rule:
        lines.append(f"   Unlock rule: {unlock_rule}")
    if branch_if_positive:
        lines.append(f"   If positive: {branch_if_positive}")
    if branch_if_negative:
        lines.append(f"   If negative: {branch_if_negative}")
    if claim_boundary:
        lines.append(f"   Claim boundary: {claim_boundary}")
    return lines


def _collect_key_materials(plan: dict) -> List[str]:
    values: List[str] = []
    for item in (plan.get("in_vitro_plan") or []) + (plan.get("in_vivo_plan") or []):
        if not isinstance(item, dict):
            continue
        values.extend(_as_str_list(item.get("key_materials_equipment")))
        if item.get("experimental_material"):
            values.append(str(item.get("experimental_material")))
        if item.get("model"):
            values.append(str(item.get("model")))
        if item.get("model_system"):
            values.append(str(item.get("model_system")))
    return _dedupe_preserve(_strip_reference_noise(v) for v in values if _strip_reference_noise(v))


def _build_background_lines(plan: dict) -> List[str]:
    bacteria = str((plan.get("candidate") or {}).get("bacteria") or "").replace("_", " ")
    metabolite = str((plan.get("candidate") or {}).get("metabolite") or "")
    disease = str((plan.get("candidate") or {}).get("disease") or "")
    question = _strip_reference_noise((plan.get("user_brief") or {}).get("research_question") or "")
    if bacteria and metabolite:
        focus_line = f"This validation protocol focuses on basic verification of the {bacteria}-{metabolite}-{disease} relationship."
    elif disease:
        focus_line = f"This validation protocol focuses on a standalone research question in the context of {disease}."
    else:
        focus_line = "This validation protocol focuses on a standalone user-defined experimental question."

    lines = [
        "Step 3: Validation Protocol",
        "",
        "Background & Objective",
        focus_line,
    ]
    if question:
        if _contains_cjk(question):
            lines.append("Primary question: User-specified validation focus applied.")
        else:
            lines.append(f"Primary question: {question}")
        lines.append(
            "This plan is organized primarily around the user's research question, using literature to support or constrain the proposed experiments."
        )
    lines.append(
        "The plan is limited to fundamental in vitro and in vivo validation tasks, with emphasis on practical grouping, timing, concentration, and readout design."
    )
    return lines


def _format_claim_line(item: dict, include_reason: bool = False, include_support: bool = True) -> str:
    if not isinstance(item, dict):
        return ""
    claim = _strip_reference_noise(item.get("claim"))
    parts = [claim] if claim else []
    tag = _strip_reference_noise(item.get("support_level"))
    if include_support and tag:
        parts.append(f"support={tag}")
    if include_reason:
        reason = _strip_reference_noise(item.get("reason"))
        if reason:
            parts.append(f"reason={reason}")
    pmids = _as_str_list(item.get("pmids"))
    if pmids:
        parts.append("PMIDs: " + ", ".join(pmids[:4]))
    return " | ".join(part for part in parts if part)


def _claim_key(text: str) -> str:
    clean = _strip_reference_noise(text)
    if "|" in clean:
        clean = clean.split("|", 1)[0]
    clean = re.sub(r"[^a-z0-9\s]+", " ", clean.lower())
    return re.sub(r"\s+", " ", clean).strip()


def _has_literature_supported_pmids(item: dict) -> bool:
    if not isinstance(item, dict):
        return False
    pmids = _as_str_list(item.get("pmids"))
    support = _strip_reference_noise(item.get("support_level")).lower()
    if not pmids:
        return False
    blocked_terms = ("speculative", "unsupported", "overstated", "inference", "uncertain", "not supported")
    return not any(term in support for term in blocked_terms)


def _summary_gap_claim_keys(plan: dict) -> set[str]:
    keys: set[str] = set()
    audit = plan.get("protocol_audit") or {}
    for bucket in ("unsupported_or_overstated_claims", "inference_only_claims"):
        for item in audit.get(bucket) or []:
            if isinstance(item, dict):
                key = _claim_key(item.get("claim"))
                if key:
                    keys.add(key)
    for gap in _as_str_list(audit.get("priority_gaps")):
        key = _claim_key(gap)
        if key:
            keys.add(key)
    return keys


def _claim_keys_overlap(left: str, right: str) -> bool:
    left_key = _claim_key(left)
    right_key = _claim_key(right)
    if not left_key or not right_key:
        return False
    return (
        left_key == right_key
        or left_key in right_key
        or right_key in left_key
    )


def _supported_claim_keys(plan: dict) -> set[str]:
    keys: set[str] = set()
    audit = plan.get("protocol_audit") or {}
    for item in audit.get("supported_claims") or []:
        if _has_literature_supported_pmids(item):
            key = _claim_key(item.get("claim"))
            if key:
                keys.add(key)
    for item in plan.get("evidence_basis") or []:
        if isinstance(item, dict) and _has_literature_supported_pmids(item):
            key = _claim_key(item.get("claim"))
            if key:
                keys.add(key)
    return keys


def _is_supported_claim_text(text: str, supported_keys: set[str]) -> bool:
    return any(_claim_keys_overlap(text, supported_key) for supported_key in supported_keys)


def _soften_overclaim_language(text: str) -> str:
    clean = _strip_reference_noise(text)
    modal_replacements = [
        (r"\b(?:can|could|may|might)\s+proves?\b", "may support"),
        (r"\b(?:can|could|may|might)\s+causes?\b", "may be associated with"),
        (r"\b(?:can|could|may|might)\s+drives?\b", "may contribute to"),
        (r"\b(?:can|could|may|might)\s+mediates?\b", "may be involved in"),
    ]
    for pattern, replacement in modal_replacements:
        clean = re.sub(pattern, replacement, clean, flags=re.IGNORECASE)
    replacements = [
        (r"\bproves?\b", "supports"),
        (r"\bcauses?\b", "is associated with"),
        (r"\bdrives?\b", "may contribute to"),
        (r"\bmediates?\b", "may be involved in"),
    ]
    for pattern, replacement in replacements:
        clean = re.sub(pattern, replacement, clean, flags=re.IGNORECASE)
    clean = re.sub(r"\b(?:can|could|may|might)\s+supports\b", "may support", clean, flags=re.IGNORECASE)
    clean = re.sub(r"\b(?:can|could|may|might)\s+is associated with\b", "may be associated with", clean, flags=re.IGNORECASE)
    clean = re.sub(r"\b(can|could|may|might)\s+may\b", "may", clean, flags=re.IGNORECASE)
    clean = re.sub(r"\bmay\s+may\b", "may", clean, flags=re.IGNORECASE)
    return re.sub(r"\s+", " ", clean).strip(" ;|,")


def _soften_process_claim_language(text: str) -> str:
    clean = _strip_reference_noise(text)
    modal_replacements = [
        (r"\b(?:can|could|may|might)\s+(?:directly\s+)?produces?\b", "may be associated with levels of"),
        (r"\b(?:can|could|may|might)\s+(?:directly\s+)?absorbs?\b", "may be associated with responses to"),
    ]
    for pattern, replacement in modal_replacements:
        clean = re.sub(pattern, replacement, clean, flags=re.IGNORECASE)
    replacements = [
        (r"\bis associated with the production of\b", "is associated with levels of"),
        (r"\bis associated with production of\b", "is associated with levels of"),
        (r"\bproduction of\b", "levels of"),
        (r"\b(?:directly\s+)?produces?\b", "is associated with levels of"),
        (r"\bis associated with the absorption of\b", "is associated with responses to"),
        (r"\bis associated with absorption of\b", "is associated with responses to"),
        (r"\babsorption of\b", "responses to"),
        (r"\b(?:directly\s+)?absorbs?\b", "is associated with responses to"),
    ]
    for pattern, replacement in replacements:
        clean = re.sub(pattern, replacement, clean, flags=re.IGNORECASE)
    clean = re.sub(r"\b(?:can|could|may|might)\s+is associated with\b", "may be associated with", clean, flags=re.IGNORECASE)
    clean = re.sub(r"\b(can|could|may|might)\s+may\b", "may", clean, flags=re.IGNORECASE)
    clean = re.sub(r"\bmay\s+may\b", "may", clean, flags=re.IGNORECASE)
    return re.sub(r"\s+", " ", clean).strip(" ;|,")


def _is_generic_gap_text(text: str) -> bool:
    lowered = _strip_reference_noise(text).lower()
    generic_patterns = (
        "association but does not confirm direct causation",
        "suggests an association but lacks direct mechanistic proof of causation",
        "does not support this claim as the primary mechanism",
        "supporting a cure",
        "primary mechanism",
        "direct causation",
        "mechanistic proof of causation",
    )
    return any(pattern in lowered for pattern in generic_patterns)


def _compress_missing_gap_lines(plan: dict, lines: List[str]) -> List[str]:
    cleaned = [line for line in _dedupe_preserve([_strip_reference_noise(line) for line in lines]) if line]
    specific = [line for line in cleaned if not _is_generic_gap_text(line)]
    candidate = plan.get("candidate") or {}
    bacteria = _strip_reference_noise(candidate.get("bacteria"))
    metabolite = _strip_reference_noise(candidate.get("metabolite"))
    disease = _strip_reference_noise(candidate.get("disease"))

    if bacteria and metabolite and disease:
        assessment = _normalize_production_evidence_assessment(
            plan.get("direct_production_evidence_assessment")
        )
        strong_types = {
            str(item.get("evidence_type") or "").strip().lower()
            for item in (assessment.get("direct_evidence") or []) + (assessment.get("paper_findings") or [])
            if isinstance(item, dict) and bool(item.get("citation_eligible"))
        }
        gaps: List[str] = []
        if "direct_monoculture_production" not in strong_types:
            gaps.append(
                f"Whether {bacteria} directly produces newly formed {metabolite} under controlled candidate-specific culture conditions remains unresolved."
            )
        if "indirect_ecological_evidence" not in strong_types:
            gaps.append(
                f"Whether {bacteria} indirectly changes {metabolite} through a validated producer, substrate transfer, or community interaction remains unresolved."
            )
        if "candidate_microbe_disease_intervention" not in strong_types:
            gaps.append(
                f"A direct disease-effect study of {bacteria} in the paper-relevant {disease} model was not found in the current bounded search."
            )
        if "candidate_metabolite_disease_intervention" in strong_types:
            gaps.append(
                f"Whether the paper-reported {metabolite}-to-{disease} effect or mechanism transfers to the {bacteria}-{metabolite} route, and whether {metabolite} is necessary or sufficient for the microbial effect, remains unresolved."
            )
        else:
            gaps.append(
                f"The direct host or disease effect of {metabolite} in a paper-relevant {disease} model was not established in the current bounded search."
            )
            gaps.append(
                f"Whether {metabolite} is necessary or sufficient for any disease effect of {bacteria} remains unresolved."
            )
        for line in specific:
            lowered = line.lower()
            if any(
                phrase in lowered
                for phrase in (
                    "downstream mechanism linking",
                    "mechanism remains to be directly established",
                    "specific mechanism",
                )
            ):
                continue
            gaps.append(line)
        return _dedupe_lines_by_claim(gaps, prefer_longer=True)[:6]

    if specific:
        return specific[: min(3, max(1, len(specific[:3])))]
    return cleaned[:1]


def _soften_known_evidence_claim(claim: str, support_level: str) -> str:
    clean_claim = _strip_reference_noise(claim)
    support = _strip_reference_noise(support_level).lower()
    if not clean_claim:
        return ""
    if "direct" in support:
        return clean_claim
    if re.search(r"\b(produc|absorp)\w*\b", clean_claim, flags=re.IGNORECASE):
        return _soften_process_claim_language(clean_claim)
    return _soften_overclaim_language(clean_claim)


def _dedupe_lines_by_claim(items: List[str], prefer_longer: bool = False) -> List[str]:
    selected: Dict[str, str] = {}
    order: List[str] = []
    for item in items:
        clean = str(item or "").strip()
        if not clean:
            continue
        key = _claim_key(clean)
        if not key:
            key = clean.lower()
        if key not in selected:
            selected[key] = clean
            order.append(key)
            continue
        if prefer_longer and len(clean) > len(selected[key]):
            selected[key] = clean
    return [selected[key] for key in order]


def _summarize_followup_experiment(item: dict, is_in_vivo: bool = False) -> str:
    if not isinstance(item, dict):
        return ""
    aim = _strip_reference_noise(item.get("aim") or ("Animal follow-up" if is_in_vivo else "In vitro follow-up"))
    model = _strip_reference_noise(item.get("model") if is_in_vivo else item.get("model_system"))
    readouts = [v for v in (_strip_reference_noise(x) for x in _as_str_list(item.get("readouts"))[:3]) if v]
    support = _strip_reference_noise(item.get("support_result"))
    parts = [aim]
    if model:
        parts.append(f"model={model}")
    if readouts:
        parts.append("readouts=" + ", ".join(readouts))
    if support:
        parts.append(f"supports_if={support}")
    return " | ".join(parts)


def _build_known_evidence_lines(plan: dict) -> List[str]:
    candidate = plan.get("candidate") or {}
    if candidate.get("bacteria") and candidate.get("metabolite"):
        assessment = _normalize_production_evidence_assessment(
            plan.get("direct_production_evidence_assessment")
        )
        allowed_types = {
            "direct_monoculture_production",
            "direct_monoculture_nonproduction",
            "candidate_microbe_disease_intervention",
            "candidate_metabolite_disease_intervention",
            "indirect_ecological_evidence",
        }
        paper_lines: List[str] = []
        unsupported_boundaries = {
            "direct_monoculture_production": (
                "disease benefit or proof that the metabolite mediates the microbe's disease effect."
            ),
            "direct_monoculture_nonproduction": (
                "absence of production under untested culture conditions or absence of indirect ecological regulation."
            ),
            "candidate_microbe_disease_intervention": (
                "production of the named metabolite or mediation of the microbial effect by that metabolite."
            ),
            "candidate_metabolite_disease_intervention": (
                "the named microbe as the direct or indirect metabolite source, or metabolite mediation of that microbe's effect."
            ),
            "indirect_ecological_evidence": (
                "direct production by the candidate microbe or causal disease mediation."
            ),
        }
        source_items = (assessment.get("direct_evidence") or []) + (assessment.get("paper_findings") or [])
        seen_pmids = set()
        for item in source_items:
            if not isinstance(item, dict):
                continue
            pmid = str(item.get("pmid") or "").strip()
            evidence_type = str(item.get("evidence_type") or "").strip().lower()
            relevance = str(item.get("candidate_relevance") or "").strip().lower()
            if (
                not pmid
                or pmid in seen_pmids
                or evidence_type not in allowed_types
                or not bool(item.get("citation_eligible"))
            ):
                continue
            if relevance != "direct":
                continue
            claim = _strip_reference_noise(item.get("claim"))
            if not claim:
                continue
            scope = _strip_reference_noise(item.get("claim_scope"))
            title = _strip_reference_noise(item.get("title"))
            parts = [
                f"PMID {pmid} ({title or 'title unavailable'})",
                f"Supports: {claim}",
                f"Does not support: {unsupported_boundaries[evidence_type]}",
            ]
            if scope:
                parts.append(f"Evidence scope: {scope}")
            paper_lines.append(" | ".join(parts))
            seen_pmids.add(pmid)
        edge_labels = {
            "direct_monoculture_production": "direct microbe-to-metabolite production",
            "direct_monoculture_nonproduction": "direct microbe-to-metabolite testing with a negative result",
            "indirect_ecological_evidence": "candidate-dependent ecological metabolite regulation",
            "candidate_microbe_disease_intervention": "named-microbe-to-disease intervention",
            "candidate_metabolite_disease_intervention": "named-metabolite-to-disease intervention",
        }
        required_edges = (
            "direct_monoculture_production",
            "indirect_ecological_evidence",
            "candidate_microbe_disease_intervention",
            "candidate_metabolite_disease_intervention",
        )
        found_types = {
            str(item.get("evidence_type") or "").strip().lower()
            for item in source_items
            if isinstance(item, dict) and bool(item.get("citation_eligible"))
        }
        found_labels = [edge_labels[edge] for edge in required_edges if edge in found_types]
        missing_labels = [edge_labels[edge] for edge in required_edges if edge not in found_types]
        coverage = (
            "Search coverage: qualifying evidence was found for "
            + (", ".join(found_labels) if found_labels else "none of the four prespecified candidate edges")
            + "; no qualifying direct paper was found in the current bounded search for "
            + (", ".join(missing_labels) if missing_labels else "none of the prespecified edges")
            + ". Not found in this search is not evidence of absence."
        )
        if paper_lines:
            return paper_lines[:7] + [coverage]
        return [
            "No strongly candidate-relevant paper met the citation threshold; weak co-occurrence, background, and analogous-method papers were intentionally excluded.",
            coverage,
        ]

    lines: List[str] = []
    audit = plan.get("protocol_audit") or {}
    for item in audit.get("supported_claims") or []:
        if not _has_literature_supported_pmids(item):
            continue
        line = _format_claim_line(
            {**item, "claim": _soften_known_evidence_claim(item.get("claim"), item.get("support_level"))},
            include_support=False,
        )
        if line:
            lines.append(line)
    for item in plan.get("evidence_basis") or []:
        if not isinstance(item, dict):
            continue
        if not _has_literature_supported_pmids(item):
            continue
        claim = _soften_known_evidence_claim(item.get("claim"), item.get("support_level"))
        if not claim:
            continue
        summary = _strip_reference_noise(item.get("evidence_summary"))
        pmids = _as_str_list(item.get("pmids"))
        line_parts = [claim]
        if summary:
            line_parts.append(summary)
        if pmids:
            line_parts.append("PMIDs: " + ", ".join(pmids[:4]))
        lines.append(" | ".join(line_parts))
    deduped = _dedupe_lines_by_claim([line for line in lines if line], prefer_longer=True)
    return deduped[:6] or ["No clearly supported literature-backed claim met the inclusion threshold after consistency review."]


def _build_direct_production_evidence_lines(plan: dict) -> List[str]:
    assessment = _normalize_production_evidence_assessment(plan.get("direct_production_evidence_assessment"))
    status = assessment.get("status") or "not_assessed"
    lines = [f"Status: {status}. {assessment.get('conclusion') or 'No conclusion was returned.'}"]
    seen_pmids = set()
    if status == "direct_supported":
        for item in assessment.get("direct_evidence") or []:
            pmid = str(item.get("pmid") or "").strip()
            lines.append(
                f"Direct culture evidence — PMID {item.get('pmid')}: {item.get('claim')} "
                f"| model={item.get('model_system') or 'not stated'} | measured={item.get('measured_output') or 'not stated'}"
            )
            if pmid:
                seen_pmids.add(pmid)
    else:
        lines.append(
            "No retrieved paper met the direct-production threshold unless listed above; this means not found in the current search, not biologically impossible."
        )
    for item in (assessment.get("paper_findings") or [])[:12]:
        if not bool(item.get("citation_eligible")):
            continue
        pmid = str(item.get("pmid") or "").strip()
        if pmid and pmid in seen_pmids:
            continue
        title = _strip_reference_noise(item.get("title"))
        lines.append(
            f"PMID {item.get('pmid')} ({title or 'title unavailable'}) individually supports: {item.get('claim')} "
            f"| evidence_type={item.get('evidence_type')} "
            f"| candidate_relevance={item.get('candidate_relevance')}"
        )
    for limitation in (assessment.get("evidence_limitations") or [])[:4]:
        lines.append(f"Limitation: {limitation}")
    return _dedupe_preserve([line for line in lines if line])[:14]


def _build_hypothesis_branch_lines(plan: dict) -> List[str]:
    lines = []
    for item in _normalize_hypothesis_branches(plan.get("hypothesis_branches")):
        parts = [
            f"{item.get('hypothesis_id')}: {item.get('statement')}",
            f"current_status={item.get('current_evidence_status')}",
        ]
        if item.get("discriminating_prediction"):
            parts.append(f"prediction={item.get('discriminating_prediction')}")
        if item.get("in_vitro_gate"):
            parts.append(f"first_gate={item.get('in_vitro_gate')}")
        if item.get("falsification_or_redirection"):
            parts.append(f"redirect_if_not_supported={item.get('falsification_or_redirection')}")
        lines.append(" | ".join(parts))
    return lines[:6] or ["Competing hypotheses were not returned; do not collapse direct production, indirect production, and independent effects into one claim."]


def _build_missing_gap_lines(plan: dict) -> List[str]:
    lines: List[str] = []
    audit = plan.get("protocol_audit") or {}
    supported_keys = _supported_claim_keys(plan)
    production_assessment = _normalize_production_evidence_assessment(
        plan.get("direct_production_evidence_assessment")
    )
    candidate = plan.get("candidate") or {}
    if (
        candidate.get("bacteria")
        and candidate.get("metabolite")
        and production_assessment.get("status") != "direct_supported"
    ):
        lines.append(
            f"No single retrieved PMID directly established controlled-culture production of {candidate.get('metabolite')} by {candidate.get('bacteria')}; direct production remains an experimental hypothesis."
        )
    for item in audit.get("unsupported_or_overstated_claims") or []:
        reason = _strip_reference_noise(item.get("reason"))
        if not reason or _is_supported_claim_text(reason, supported_keys):
            continue
        lines.append(reason)
    for item in audit.get("inference_only_claims") or []:
        reason = _strip_reference_noise(item.get("reason"))
        if not reason or _is_supported_claim_text(reason, supported_keys):
            continue
        lines.append(reason)
    for gap in _as_str_list(audit.get("priority_gaps")):
        if _is_supported_claim_text(gap, supported_keys):
            continue
        gap_clean = _strip_reference_noise(gap)
        if gap_clean:
            lines.append(gap_clean)
    for item in plan.get("self_reflection") or []:
        if not isinstance(item, dict):
            continue
        remaining = _strip_reference_noise(item.get("remaining_uncertainty"))
        if _is_supported_claim_text(remaining, supported_keys):
            continue
        if remaining:
            lines.append(remaining)
    for flag in _as_str_list(plan.get("overall_risk_flags")):
        flag_clean = _strip_reference_noise(flag)
        if _is_supported_claim_text(flag_clean, supported_keys):
            continue
        if flag_clean:
            lines.append(flag_clean)
    compressed = _compress_missing_gap_lines(plan, [line for line in lines if line])
    return compressed[:6] or ["No explicit missing-evidence gap was returned; treat the mechanism and causal direction as not yet directly established."]


def _build_why_this_experiment_lines(plan: dict) -> List[str]:
    candidate = plan.get("candidate") or {}
    bacteria = _strip_reference_noise(candidate.get("bacteria"))
    metabolite = _strip_reference_noise(candidate.get("metabolite"))
    disease = _strip_reference_noise(candidate.get("disease"))
    question = _strip_reference_noise((plan.get("user_brief") or {}).get("research_question"))
    profile = _normalize_question_profile(plan.get("question_profile"))
    question_type = profile.get("question_type", "general_validation")
    policy = profile.get("human_study_policy", "not_recommended")
    gaps = plan.get("validation_summary", {}).get("missing_evidence_gap") if isinstance(plan.get("validation_summary"), dict) else []
    gap_hint = gaps[0] if gaps else ""

    lines: List[str] = []
    if bacteria and metabolite and disease:
        lines.append(
            f"The proposed experiment set is intended to test whether the {bacteria}-{metabolite} signal can be linked to {disease} with direct experimental support, rather than relying on correlation or literature inference alone."
        )
        lines.append(
            "The sequence keeps direct production, indirect ecological production, and parallel disease effects as separate hypotheses: controlled culture decides the first branch, branch-specific animal work tests causality, and human sampling is designed only after the animal result is interpretable."
        )
    elif question:
        lines.append(
            "The proposed experiment set is intended to answer the user's validation question with direct operational evidence, rather than a broad descriptive protocol."
        )
    else:
        lines.append(
            "The proposed experiment set is intended to close the highest-priority evidence gaps before escalating to a larger study."
        )
    if gap_hint:
        lines.append(f"The immediate reason for these experiments is the current gap: {gap_hint}")
    if question_type == "sensor_discovery":
        lines.append(
            "The sequence is organized as a gatekeeper workflow: identify the candidate sensor or target first, then test necessity, sufficiency, and pathway relevance before any broader translational claim."
        )
    elif question_type == "toxic_exposure":
        lines.append(
            "The sequence is organized as a gatekeeper workflow: establish exposure-linked biological signal and mechanism preclinically first, because direct human causal testing would be ethically weak or inappropriate."
        )
    elif policy == "allow_observational_only":
        lines.append(
            "The sequence is organized as a gatekeeper workflow: establish the preclinical logic first, then use observational human data only to test translational relevance rather than direct efficacy."
        )
    elif not (bacteria and metabolite and disease):
        lines.append(
            "The sequence is organized as a gatekeeper workflow: verify the most fragile link first, then escalate only if the basic biological signal is reproducible."
        )
    return lines[:3]


def _build_priority_followup_lines(plan: dict) -> List[str]:
    lines: List[str] = []
    for item in (plan.get("in_vitro_plan") or [])[:3]:
        line = _summarize_followup_experiment(item, is_in_vivo=False)
        if line:
            lines.append(line)
    for item in (plan.get("in_vivo_plan") or [])[:3]:
        line = _summarize_followup_experiment(item, is_in_vivo=True)
        if line:
            lines.append(line)
    return _dedupe_preserve(lines)[:6] or ["No concrete follow-up experiment block was returned; start with a minimal culture or host-cell verification study."]


def _backfill_conditional_causal_followups(plan: dict) -> dict:
    """Retain branch-specific A2/A3 decision cards without rebuilding a fixed plan."""
    candidate = plan.get("candidate") if isinstance(plan.get("candidate"), dict) else {}
    bacteria = _strip_reference_noise(candidate.get("bacteria")) or "the candidate microbe"
    metabolite = _strip_reference_noise(candidate.get("metabolite")) or "the target metabolite"
    disease = _strip_reference_noise(candidate.get("disease")) or "the target disease"
    vitro_items = [item for item in (plan.get("in_vitro_plan") or []) if isinstance(item, dict)]
    vivo_items = [item for item in (plan.get("in_vivo_plan") or []) if isinstance(item, dict)]
    vitro_roles = {_canonicalize_experiment_role(item.get("experiment_role")) for item in vitro_items}
    vivo_roles = {_canonicalize_experiment_role(item.get("experiment_role")) for item in vivo_items}
    if "effect_interaction_gate" not in vivo_roles:
        return plan

    a1 = next(
        (item for item in vivo_items if _canonicalize_experiment_role(item.get("experiment_role")) == "effect_interaction_gate"),
        {},
    )
    a1_primary = str(item_value or "").strip() if (
        item_value := a1.get("primary_indicator") or a1.get("primary_endpoint")
    ) else "The same prespecified primary disease indicator used in A1"
    additions: List[dict] = []

    if "direct_production_gate" in vitro_roles and "direct_causal_rescue" not in vivo_roles:
        additions.append(
            {
                "module_id": "A2",
                "experiment_role": "direct_causal_rescue",
                "hypothesis_ids": ["H1"],
                "route_status": "next_if_positive",
                "result_status": "not_run",
                "scientific_question": (
                    f"Is {metabolite} production necessary for the disease effect of {bacteria}, and can the effect be rescued?"
                ),
                "why_needed": (
                    "Separate direct-branch mediation from parallel effects after both direct production and the A1 disease effect are positive."
                ),
                "study_object": (
                    f"The A1 {disease} model with a candidate-specific {metabolite}-source function perturbation selected only after V1."
                ),
                "groups": [
                    {
                        "group_name": "Disease control",
                        "exposure_or_condition": "Disease condition without candidate microbe or metabolite rescue",
                        "control_purpose": "Define the untreated disease reference",
                    },
                    {
                        "group_name": "Source-function competent microbe",
                        "exposure_or_condition": f"Disease condition plus source-function competent {bacteria}",
                        "control_purpose": "Reproduce the A1 microbial effect",
                    },
                    {
                        "group_name": "Source-function loss",
                        "exposure_or_condition": f"Disease condition plus {bacteria} with specific loss of the V1-supported source function",
                        "control_purpose": "Test necessity of the direct source function",
                    },
                    {
                        "group_name": "Function complementation",
                        "exposure_or_condition": "The function-loss condition with candidate-specific functional complementation when feasible",
                        "control_purpose": "Exclude unrelated perturbation effects",
                    },
                    {
                        "group_name": "Metabolite add-back rescue",
                        "exposure_or_condition": f"The function-loss condition plus {metabolite} add-back",
                        "control_purpose": "Test metabolite-specific rescue",
                    },
                ],
                "samples_and_timing": (
                    "Match A1 baseline and endpoint sampling; collect material needed to verify microbial exposure, target-metabolite exposure, and the primary disease outcome."
                ),
                "primary_indicator": a1_primary,
                "secondary_indicators": [
                    "Comparable colonization or microbial exposure across microbe-containing groups",
                    f"Target-metabolite exposure and source-function readout for {metabolite}",
                    "A small set of disease, barrier, inflammatory, and safety indicators aligned with A1",
                ],
                "key_controls": [
                    "Source-function competent comparator",
                    "Perturbation-specific complementation when feasible",
                    "Matched metabolite vehicle and add-back control",
                ],
                "positive_gate": (
                    "The A1 microbial benefit is lost after specific source-function loss, microbial exposure remains comparable, and the effect is restored by complementation or target-metabolite add-back."
                ),
                "branch_if_positive": "Support a direct-branch causal contribution within the tested model and conditions.",
                "branch_if_negative": (
                    "Interpret the microbial and metabolite effects as parallel, partly overlapping, or mediated by another function; do not claim metabolite mediation."
                ),
                "claim_boundary": DEFAULT_CLAIM_SCOPES["direct_causal_rescue"],
            }
        )

    if "indirect_ecology_gate" in vitro_roles and "indirect_ecology_causal" not in vivo_roles:
        additions.append(
            {
                "module_id": "A3",
                "experiment_role": "indirect_ecology_causal",
                "hypothesis_ids": ["H2"],
                "route_status": "next_if_positive",
                "result_status": "not_run",
                "scientific_question": (
                    f"Is the V2-validated ecological producer route necessary for the disease effect associated with {bacteria} and {metabolite}?"
                ),
                "why_needed": (
                    "Test ecological mediation only after V2 identifies a candidate-dependent producer route and A1 establishes a disease effect."
                ),
                "study_object": (
                    f"The A1 {disease} model with one controlled defined community containing the V2-validated producer."
                ),
                "groups": [
                    {
                        "group_name": "Disease control",
                        "exposure_or_condition": "Disease condition without the defined community",
                        "control_purpose": "Define the untreated disease reference",
                    },
                    {
                        "group_name": "Candidate-absent community",
                        "exposure_or_condition": "Validated-producer community without candidate microbe",
                        "control_purpose": "Estimate producer-community activity without the candidate",
                    },
                    {
                        "group_name": "Candidate-present community",
                        "exposure_or_condition": f"The matched validated-producer community plus {bacteria}",
                        "control_purpose": "Test the candidate-dependent ecological effect",
                    },
                    {
                        "group_name": "Producer-absent community",
                        "exposure_or_condition": f"The matched community plus {bacteria} but without the validated producer",
                        "control_purpose": "Test producer necessity",
                    },
                    {
                        "group_name": "Producer-absent metabolite rescue",
                        "exposure_or_condition": f"The producer-absent condition plus {metabolite} add-back rescue",
                        "control_purpose": "Test whether the target metabolite rescues loss of the ecological route",
                    },
                ],
                "samples_and_timing": (
                    "Match A1 baseline and endpoint sampling; collect material needed to verify candidate and producer abundance, target-metabolite exposure, and the primary disease outcome."
                ),
                "primary_indicator": a1_primary,
                "secondary_indicators": [
                    "Candidate and validated-producer abundance or activity",
                    f"Source-resolved or otherwise producer-linked {metabolite} exposure",
                    "A small set of disease, barrier, inflammatory, and safety indicators aligned with A1",
                ],
                "key_controls": [
                    "Matched candidate-present and candidate-absent communities",
                    "Producer-present and producer-absent communities",
                    "Matched metabolite vehicle and add-back rescue",
                ],
                "positive_gate": (
                    "The candidate-dependent A1 benefit and metabolite signal require the validated producer and are restored by target-metabolite add-back after producer removal."
                ),
                "branch_if_positive": "Support an ecological causal contribution within the tested defined community and disease model.",
                "branch_if_negative": (
                    "Do not claim ecological mediation; retain parallel-effect or alternative-community explanations."
                ),
                "claim_boundary": DEFAULT_CLAIM_SCOPES["indirect_ecology_causal"],
            }
        )

    if additions:
        plan["in_vivo_plan"] = _prepare_experiment_modules(
            vivo_items + additions,
            is_in_vivo=True,
        )
    return plan


def _build_question_driven_human_experiments(plan: dict) -> List[dict]:
    brief = plan.get("user_brief") or {}
    question = _strip_reference_noise(brief.get("research_question"))
    constraints = _strip_reference_noise(brief.get("prompt_constraints"))
    combined_brief = f"{question} {constraints}".lower()
    candidate = plan.get("candidate") or {}
    disease = _strip_reference_noise(candidate.get("disease")) or "the target condition"
    bacteria = _strip_reference_noise(candidate.get("bacteria"))
    metabolite = _strip_reference_noise(candidate.get("metabolite"))
    mode = str(plan.get("mode") or "").strip()
    profile = _normalize_question_profile(plan.get("question_profile"))
    question_type = profile.get("question_type", "general_validation")
    human_feasibility = profile.get("human_feasibility", "uncertain")
    human_ethics = profile.get("human_ethics", "uncertain")
    recommended_focus = profile.get("recommended_focus", "preclinical_first")
    human_study_policy = profile.get("human_study_policy", "not_recommended")
    detail_budget = profile.get("detail_budget_allocation", "balanced_preclinical")
    toxic_like = question_type == "toxic_exposure"
    sensor_like = question_type == "sensor_discovery"
    beneficial_like = question_type in ("beneficial_intervention", "behavioral_or_observational", "active_compound_screening")
    behavior_only = human_study_policy == "allow_low_risk_behavior_only"
    observational_only = human_study_policy == "allow_observational_only"
    allows_rct = human_study_policy == "allow_late_stage_rct"
    direct_exposure_forbidden = human_study_policy == "forbid_direct_exposure"
    wants_human = any(
        term in combined_brief
        for term in ("human", "patient", "clinical", "trial", "rct", "randomized", "efficacy", "translate", "translational")
    )
    preclinical_only = any(
        term in combined_brief
        for term in ("only in vitro", "in vitro only", "cell only", "organoid only", "avoid clinical", "no clinical", "preclinical only")
    )

    if _plan_is_candidate_workflow(plan):
        if preclinical_only:
            return []
        bacteria_label = bacteria or "the candidate microbe"
        metabolite_label = metabolite or "the target metabolite"
        disease_label = disease or "the target disease"
        disease_context = f"{combined_brief} {disease_label}".lower()
        disease_name_is_crohn = bool(
            re.search(r"\bcrohn(?:'s)?(?: disease)?\b", disease_label.lower())
            or re.fullmatch(r"cd", disease_label.strip(), flags=re.IGNORECASE)
        )
        broad_ibd_brief = bool(
            re.search(
                r"(?:\buc\b|ulcerative colitis).{0,40}(?:\bcd\b|crohn)|"
                r"(?:\bcd\b|crohn).{0,40}(?:\buc\b|ulcerative colitis)",
                combined_brief,
            )
        )
        brief_explicitly_crohn = bool(
            re.search(r"\bcrohn(?:'s)?(?: disease)?\b", combined_brief)
            or re.search(
                r"\b(?:focus(?:ed)? on|restricted to|specifically in|patients? with)\s+cd\b",
                combined_brief,
            )
        )
        crohn_specific = disease_name_is_crohn or (
            brief_explicitly_crohn and not broad_ibd_brief
        )
        ibd_like = crohn_specific or any(
            term in disease_context
            for term in ("\bibd\b", "inflammatory bowel", "ulcerative colitis", "colitis")
        )
        if "ibd" in disease_context:
            ibd_like = True
        stratum_label = "Crohn disease (CD)" if crohn_specific else disease_label
        groups = [
            {
                "group_name": f"Active {stratum_label}",
                "exposure_or_condition": f"Participants with active {stratum_label}",
                "control_purpose": "Estimate association during active disease",
            },
            {
                "group_name": f"Remission or stable {stratum_label}",
                "exposure_or_condition": f"Participants with clinically defined remission or stable {stratum_label}",
                "control_purpose": "Separate disease-state association from diagnosis alone",
            },
            {
                "group_name": "Healthy control",
                "exposure_or_condition": "Participants without the target disease",
                "control_purpose": "Provide the nondisease reference distribution",
            },
        ]
        secondary_indicators = [
            f"Absolute abundance of {bacteria_label} in stool and disease-relevant tissue when clinically obtained",
            f"Quantitative {metabolite_label} and related metabolite profile",
            "Inflammatory and clinical disease-activity indicators",
            "Direction and temporal coherence in a prespecified longitudinal subset",
        ]
        if ibd_like:
            secondary_indicators.append(
                "Fecal calprotectin, CRP, and clinical or endoscopic disease activity when available"
            )
        return [
            {
                "module_id": "H1",
                "experiment_role": "human_translation",
                "hypothesis_ids": ["HUMAN_BRANCH_PENDING"],
                "route_status": "conditional",
                "execution_status": "conditional_future",
                "result_status": "not_run",
                "design_status": "decision_complete",
                "activation_gate": _default_activation_gate("human_translation"),
                "scientific_question": (
                    f"Is the animal-supported relationship among {bacteria_label}, {metabolite_label}, "
                    f"and {disease_label} directionally and temporally compatible with human observations?"
                ),
                "why_needed": (
                    "Test human relevance of one preclinically supported branch without using association "
                    "to repair incomplete causal evidence."
                ),
                "study_object": f"Observational cohort spanning activity states of {stratum_label}",
                "groups": groups,
                "samples_and_timing": (
                    "Collect stool and blood at the prespecified study visit; use disease-relevant tissue only "
                    "when obtained for clinical care, and repeat the same samples in a longitudinal subset."
                ),
                "primary_indicator": (
                    f"Prespecified association of quantitative {metabolite_label} with disease activity, "
                    f"conditioned on the animal-supported {bacteria_label} branch"
                ),
                "secondary_indicators": secondary_indicators,
                "key_controls": [
                    "Dietary protein and fiber intake",
                    "Antibiotics, probiotics, and disease-directed medication",
                    "Smoking, BMI, disease subtype, and anatomical disease location",
                    "Collection timing, sample processing, storage, and analytical batch",
                ],
                "positive_gate": (
                    "The prespecified direction is consistent after major confounder control and is supported "
                    "by the longitudinal subset or an independent cohort."
                ),
                "branch_if_positive": (
                    "Support human relevance of the selected branch while retaining an association-only claim."
                ),
                "branch_if_negative": (
                    "Do not infer human relevance or mediation; reassess transferability, measurement, and branch selection."
                ),
                "unlock_rule": (
                    "Require a positive A1, A2, or A3 result, positive independent replication, and positive "
                    "ethics_and_analysis_plan gates before recruitment or sampling."
                ),
                "claim_boundary": (
                    "Association, direction, and temporal compatibility only; this module cannot prove microbial "
                    "production, causal mediation, or disease causation."
                ),
            }
        ]

    if "mosquito" in combined_brief or "repellent" in combined_brief or "biting" in combined_brief or "landing" in combined_brief:
        primary_endpoint = "Human landing or bite-avoidance endpoint under standardized exposure"
        mediator_endpoint = "Exposure-time mediator and dose-response linkage to the repellency endpoint"
        active_arm = "Active plant extract, fraction, or compound arm"
        control_arm = "Vehicle control and positive repellent control arm"
    elif observational_only:
        primary_endpoint = "Prespecified human observational endpoint linked to the biological question"
        mediator_endpoint = "Exposure-biomarker or phenotype-biomarker relationship under non-interventional sampling"
        active_arm = "Higher-exposure or higher-signal stratum"
        control_arm = "Lower-exposure or matched reference stratum"
    elif question_type == "active_compound_screening":
        primary_endpoint = "Human translational endpoint linked to the prioritized active fraction or compound"
        mediator_endpoint = "Mediator and exposure-response relationship for the prioritized active component"
        active_arm = "Active fraction or prioritized compound arm"
        control_arm = "Matched placebo or comparator arm"
    else:
        primary_endpoint = "Primary human efficacy endpoint"
        mediator_endpoint = "Candidate mediator or biomarker linked to the mechanism"
        active_arm = "Active intervention arm"
        control_arm = "Matched placebo or standard-control arm"

    human_ethically_plausible = False
    if human_ethics == "appropriate" and human_feasibility in ("high", "medium"):
        human_ethically_plausible = True
    elif human_ethics == "conditional" and human_feasibility in ("high", "medium") and human_study_policy in (
        "allow_observational_only",
        "allow_low_risk_behavior_only",
        "allow_late_stage_rct",
    ):
        human_ethically_plausible = True
    elif beneficial_like and not toxic_like and not sensor_like and allows_rct:
        human_ethically_plausible = True

    if mode == "question_driven" and (preclinical_only or direct_exposure_forbidden or human_study_policy == "not_recommended"):
        return []
    if mode == "question_driven" and not human_ethically_plausible and not (bacteria and metabolite and disease):
        return []
    if mode == "question_driven" and detail_budget in ("discovery_heavy", "cell_heavy", "animal_heavy") and not wants_human:
        return []

    if mode == "question_driven":
        if observational_only:
            primary_objective = (
                f"Can a non-interventional human cohort or case-control design test whether the main human endpoint covaries with the exposure or biomarker relevant to {disease or 'the target condition'}?"
                if not question
                else f"Can the user question be narrowed in humans through observational sampling rather than intervention? {question}"
            )
            mechanistic_objective = f"Do serial or stratified biospecimens show a human mediator pattern aligned with the observed phenotype in {disease or 'the target condition'}?"
            efficacy_material = "Human biospecimens, exposure metadata, matched phenotype groups, and validated biomarker assays"
            mechanistic_material = "Longitudinal or stratified biospecimen collection, phenotype metadata, and mediator assays"
            efficacy_intervention = "No active intervention; perform stratified observational sampling and prespecified covariate adjustment"
        elif behavior_only:
            primary_objective = (
                f"Can a low-risk human behavioral study test the operational endpoint relevant to {disease or 'the user question'}?"
                if not question
                else f"Can the user question be answered in a low-risk controlled human behavioral setting? {question}"
            )
            mechanistic_objective = "Do repeated behavioral measurements and exposure verification support a mediator-linked interpretation rather than a one-off preference effect?"
            efficacy_material = "Standardized exposure setup, matched controls, exposure verification, and safety-compatible operational measurements"
            mechanistic_material = "Repeated-measures behavioral setup, serial sampling when relevant, and mediator or exposure assays"
            efficacy_intervention = "Apply the low-risk behavioral or exposure-control protocol under tightly standardized conditions"
        else:
            primary_objective = (
                f"Can the intervention change the main human endpoint relevant to {disease} in a randomized setting?"
                if not question
                else f"Can the proposed intervention answer the user question in humans under randomized conditions? {question}"
            )
            mechanistic_objective = (
                f"Does change in the proposed mediator temporally track change in the human endpoint relevant to {disease}?"
            )
            efficacy_material = "Standardized intervention/exposure, matched placebo or control, adherence and exposure assays"
            mechanistic_material = "Intervention, matched control condition, serial biospecimen collection, and mediator assays"
            efficacy_intervention = "Administer the candidate intervention or exposure-control strategy under a prespecified protocol"
    elif bacteria and metabolite and disease:
        primary_objective = (
            f"Does the {bacteria}-{metabolite} axis improve the main human endpoint relevant to {disease} in a randomized translational study?"
        )
        mechanistic_objective = (
            f"Do changes in {metabolite} temporally track human endpoint changes after {bacteria} or metabolite-targeted intervention in {disease}?"
        )
        efficacy_material = f"{bacteria}-targeted intervention or standardized {metabolite}-modulating exposure, matched placebo/control, adherence and exposure assays"
        mechanistic_material = f"Intervention, matched control condition, serial biospecimen collection, {metabolite} quantification, and mediator assays"
        efficacy_intervention = f"Administer a prespecified {bacteria}- or {metabolite}-targeted intervention under a randomized translational protocol"
    else:
        return []

    if observational_only:
        return [
            {
                "aim": "To test the human relevance signal through an observational study rather than direct intervention.",
                "biological_question": primary_objective,
                "model": f"Cross-sectional or longitudinal human observational study relevant to {disease or 'the target condition'}",
                "experimental_material": efficacy_material,
                "groups": [
                    active_arm,
                    control_arm,
                    "Matching or covariate adjustment for age, sex, baseline phenotype, and major background exposures",
                ],
                "intervention": efficacy_intervention,
                "intervention_route": "No intervention route; sampling and metadata collection only",
                "dose_timing_logic": "Use predefined sampling windows and exposure metadata capture so temporal ambiguity is minimized where possible",
                "timeline": "Single-visit case-control design or repeated biospecimen sampling if longitudinal follow-up is feasible",
                "design": "Observational human design focused on relevance, stratification, confounder control, and mediator-phenotype coupling rather than efficacy testing",
                "readouts": [
                    primary_endpoint,
                    mediator_endpoint,
                    "Key covariates, medication use, and background exposure metrics",
                ],
                "priority_rationale": "Use a human block here only because non-interventional sampling can improve relevance without creating unethical exposure.",
                "gap_addressed": "Whether the proposed signal appears in humans at all, and whether it aligns with the mediator pattern predicted from preclinical work.",
                "group_logic": "Matched or adjusted strata reduce obvious bias while preserving the main question of whether phenotype and biomarker co-vary in humans.",
                "readout_rationale": "Pair the phenotype with a plausible mediator so the study can separate simple association from a more coherent translational signal.",
                "primary_endpoint": primary_endpoint,
                "success_threshold": "A prespecified and directionally consistent association that survives major confounder adjustment and is reproducible across strata or time points.",
                "failure_action": "Do not claim human relevance; refine biomarker selection, improve exposure measurement, or return to stronger mechanistic preclinical work.",
                "positive_result_interpretation": "A positive result supports human relevance and prioritizes a cleaner mechanistic or interventional next step if ethics later allow it.",
                "negative_result_interpretation": "A null result suggests weak translational relevance, poor exposure discrimination, or the wrong mediator.",
                "mechanism_support_criterion": "Mechanistic support requires that the mediator tracks the phenotype in the predicted direction after prespecified confounder adjustment.",
                "phenomenology_support_criterion": "If only the phenotype association appears without mediator alignment, treat the result as observational relevance only.",
                "key_confounders": [
                    "Reverse causation and temporal ambiguity",
                    "Medication, diet, or behavioral co-exposures",
                    "Selection bias and residual confounding",
                ],
                "go_no_go": "Go only if the association is directionally consistent, robust to adjustment, and supported by sampling quality rather than one unstable subgroup.",
                "decision_impact": "A passing result justifies denser longitudinal sampling or a lower-risk causal test where ethically feasible.",
                "data_analysis": "Multivariable regression, stratified sensitivity analysis, repeated-measures modeling when available, and robustness checks for confounding and missingness",
            }
        ]

    first_block = {
        "aim": "To test the human efficacy signal in a pilot randomized controlled trial." if allows_rct else "To test the human operational signal in a low-risk controlled study.",
        "biological_question": primary_objective,
        "model": (
            f"Pilot randomized, placebo-controlled human study in participants relevant to {disease}"
            if allows_rct
            else f"Controlled low-risk human behavioral study relevant to {disease or 'the target question'}"
        ),
        "experimental_material": efficacy_material,
        "groups": [
            active_arm,
            control_arm,
            "Stratification by baseline severity or exposure level when appropriate",
        ],
        "intervention": efficacy_intervention,
        "intervention_route": "Route matched to the intended human use case" if allows_rct else "Exposure format matched to the low-risk operational setting",
        "dose_timing_logic": (
            "Use a literature-supported starting regimen with predefined adherence checks and escalation only if safety permits"
            if allows_rct
            else "Use short, standardized exposure windows with prespecified repetition and safety-compatible operational timing"
        ),
        "timeline": (
            "Parallel-group 4-12 week study or shorter proof-of-concept window matched to endpoint kinetics"
            if allows_rct
            else "Repeated-session crossover or balanced within-subject comparison matched to short-term behavioral endpoint kinetics"
        ),
        "design": (
            "Randomized, assessor-blinded design with predefined primary endpoint, safety monitoring, and exposure verification"
            if allows_rct
            else "Randomized or counterbalanced low-risk behavioral design with predefined operational endpoint and exposure verification"
        ),
        "readouts": [
            primary_endpoint,
            mediator_endpoint,
            "Safety, tolerability, and adherence" if allows_rct else "Operational consistency, carryover checks, and exposure fidelity",
        ],
        "priority_rationale": "Use this only after preclinical or operational evidence is strong enough that a human proof-of-concept decision is justified.",
        "gap_addressed": "Whether the intervention produces a reproducible human efficacy or operational signal and whether the proposed mediator moves in the same direction.",
        "group_logic": "The comparator arm estimates nonspecific change, while stratification protects against imbalance in baseline severity or exposure.",
        "readout_rationale": "Pair one prespecified human endpoint with one mediator measurement so the study can distinguish efficacy or operational signal from simple exposure without mechanism.",
        "primary_endpoint": primary_endpoint,
        "success_threshold": "A prespecified clinically or operationally meaningful improvement versus control with consistent exposure verification.",
        "failure_action": "Do not escalate to a larger trial; either refine the intervention, shorten the mechanistic claim, or return to preclinical optimization.",
        "positive_result_interpretation": "A positive signal supports human-level relevance and justifies a tighter mechanism-oriented follow-up, but does not by itself prove mediation.",
        "negative_result_interpretation": "A null result argues against immediate translation and suggests either insufficient exposure, wrong population, or a weak underlying biological effect.",
        "mechanism_support_criterion": "The human endpoint and mediator should shift in the predicted direction with temporal alignment and exposure verification.",
        "phenomenology_support_criterion": "If the endpoint improves but the mediator does not move as predicted, count this as efficacy or operational support without mechanism confirmation.",
        "key_confounders": [
            "Baseline severity imbalance",
            "Adherence failure or exposure misclassification",
            "Diet, co-medication, or background behavior changes",
        ],
        "go_no_go": "Go only if the prespecified primary endpoint changes in the expected direction with acceptable safety and exposure confirmation.",
        "decision_impact": "A passing signal justifies a larger human study or a mechanistic crossover design with denser mediator sampling.",
        "data_analysis": "Intention-to-treat comparison with mixed-effects models, prespecified subgroup analysis, mediator-outcome coupling checks, and effect-size estimation",
    }

    second_block = {
        "aim": "To test the causal mechanism in a deeper human follow-up study.",
        "biological_question": mechanistic_objective,
        "model": (
            f"Crossover, randomized-withdrawal, or washout-rechallenge human mechanistic study in {disease}"
            if allows_rct
            else f"Repeated-measures low-risk human mechanism study in {disease or 'the target question'}"
        ),
        "experimental_material": mechanistic_material,
        "groups": [
            "Intervention-first sequence",
            "Control-first sequence",
            "Washout or withdrawal phase when ethically feasible",
        ],
        "intervention": "Apply timed intervention periods with repeated mediator and endpoint sampling to test temporal precedence",
        "intervention_route": "Same route as the intended clinical use case with exposure confirmation" if allows_rct else "Same operational exposure format with repeated verification",
        "dose_timing_logic": "Include repeated-measures sampling, dose-response contrasts when feasible, and predefined washout timing",
        "timeline": "Multi-period design with baseline, on-treatment, washout, and rechallenge or withdrawal assessments",
        "design": "Mechanistic human study focused on temporal ordering, dose-response, and mediation rather than efficacy alone",
        "readouts": [
            f"{primary_endpoint} over time",
            f"{mediator_endpoint} over time",
            "Within-subject coupling between mediator change and endpoint change",
        ],
        "priority_rationale": "Escalate to this only when a simpler efficacy or operational signal already exists and the main open question is causal ordering rather than basic activity.",
        "gap_addressed": "Whether the proposed mediator truly sits on the causal path between intervention and human outcome.",
        "group_logic": "Sequence switching or rechallenge helps separate stable subject effects from treatment-linked shifts in mediator and endpoint.",
        "readout_rationale": "Dense repeated measures are required to test temporal precedence, mediator movement, and within-subject reversibility.",
        "primary_endpoint": f"{primary_endpoint} over time",
        "success_threshold": "A reproducible within-subject or sequence-specific endpoint shift accompanied by aligned mediator change during intervention and reversal during washout or withdrawal.",
        "failure_action": "Downgrade the mechanism claim and treat any efficacy or operational signal as phenomenological unless another mediator or exposure model is better supported.",
        "positive_result_interpretation": "A positive result strengthens the causal interpretation because timing, reversibility, and mediator coupling all point in the same direction.",
        "negative_result_interpretation": "A null or uncoupled result argues that the proposed mediator is incomplete, secondary, or not causal in humans.",
        "mechanism_support_criterion": "Mechanism support requires temporal precedence, intervention-linked mediator change, and mediator-outcome coupling under repeated measures.",
        "phenomenology_support_criterion": "If endpoints move without mediator coupling or reversibility, count the result as phenomenological rather than causal.",
        "key_confounders": [
            "Carryover effects between periods",
            "Time-varying co-exposures or adherence drift",
            "Regression to the mean and spontaneous fluctuation",
        ],
        "go_no_go": "Go only if the repeated-measures design shows the predicted temporal ordering and mediator-endpoint coupling, not just endpoint change alone.",
        "decision_impact": "A passing result supports a stronger causal claim and a more definitive confirmatory human study.",
        "data_analysis": "Mixed-effects longitudinal modeling, mediation analysis, carryover sensitivity analysis, and within-subject response classification",
    }
    return [first_block, second_block]


def _render_execution_status(item: dict) -> str:
    status = str(item.get("execution_status") or "").strip().lower()
    if status in VALID_EXECUTION_STATUSES:
        return status
    # Unknown-status blocks must not be promoted into the executable plan.
    return "conditional_future"


def _iter_render_experiment_modules(plan: dict) -> List[tuple[str, dict, bool]]:
    modules: List[tuple[str, dict, bool]] = []
    seen = set()
    sources: List[tuple[str, Any, bool]] = [
        ("In Vitro / Cell Experiments", plan.get("in_vitro_plan") or [], False),
        ("In Vivo / Animal Experiments", plan.get("in_vivo_plan") or [], True),
        ("Human / Translational Experiments", plan.get("human_plan") or [], True),
    ]
    partition_fields = (
        ("current_executable_plan", "ready_now"),
        ("conditional_future_plan", "conditional_future"),
    )
    category_fields = (
        ("In Vitro / Cell Experiments", "in_vitro", False),
        ("In Vivo / Animal Experiments", "in_vivo", True),
        ("Human / Translational Experiments", "human", True),
    )
    for field, _expected_status in partition_fields:
        bundle = plan.get(field) or {}
        if not isinstance(bundle, dict):
            continue
        for category, key, is_in_vivo in category_fields:
            sources.append((category, bundle.get(key) or [], is_in_vivo))

    for category, items, is_in_vivo in sources:
        if not isinstance(items, list):
            continue
        for item_index, item in enumerate(items):
            if not isinstance(item, dict):
                continue
            module_id = str(item.get("module_id") or "").strip()
            fallback_identity = str(item.get("aim") or "").strip() or f"anonymous-{item_index}-{len(modules)}"
            identity = (category, module_id or fallback_identity)
            if identity in seen:
                continue
            seen.add(identity)
            modules.append((category, item, is_in_vivo))
    return modules


def _render_mapping_summary(item: dict, preferred_fields: Optional[List[str]] = None) -> str:
    if not isinstance(item, dict):
        return _strip_reference_noise(item)
    fields = preferred_fields or list(item.keys())
    parts: List[str] = []
    for field in fields:
        value = item.get(field)
        if isinstance(value, list):
            rendered = ", ".join(
                _strip_reference_noise(entry)
                for entry in value
                if not isinstance(entry, (dict, list)) and _strip_reference_noise(entry)
            )
        elif isinstance(value, dict):
            rendered = ""
        else:
            rendered = _strip_reference_noise(value)
        if rendered:
            parts.append(f"{field}={rendered}")
    return " | ".join(parts)


def _render_loose_value(value: Any, label: str = "", limit: int = 16) -> List[str]:
    lines: List[str] = []
    if isinstance(value, str):
        clean = _strip_reference_noise(value)
        if clean:
            lines.append(f"{label}: {clean}" if label else clean)
    elif isinstance(value, list):
        for entry in value[:limit]:
            if isinstance(entry, dict):
                rendered = _render_mapping_summary(entry)
            else:
                rendered = _strip_reference_noise(entry)
            if rendered:
                lines.append(f"{label}: {rendered}" if label else rendered)
    elif isinstance(value, dict):
        for key, entry in list(value.items())[:limit]:
            if isinstance(entry, list):
                lines.extend(_render_loose_value(entry, label=str(key), limit=limit))
            elif isinstance(entry, dict):
                rendered = _render_mapping_summary(entry)
                if rendered:
                    lines.append(f"{key}: {rendered}")
            else:
                clean = _strip_reference_noise(entry)
                if clean:
                    lines.append(f"{key}: {clean}")
    return lines


def _render_priority_followup_section(plan: dict, execution_status: Optional[str] = None) -> List[str]:
    grouped: Dict[str, List[tuple[dict, bool]]] = {}
    for category, item, is_in_vivo in _iter_render_experiment_modules(plan):
        status = _render_execution_status(item)
        if execution_status and status != execution_status:
            continue
        grouped.setdefault(category, []).append((item, is_in_vivo))

    lines: List[str] = []
    for category in (
        "In Vitro / Cell Experiments",
        "In Vivo / Animal Experiments",
        "Human / Translational Experiments",
    ):
        items = grouped.get(category) or []
        if not items:
            continue
        if lines:
            lines.append("")
        lines.append(category)
        for idx, (item, is_in_vivo) in enumerate(items, start=1):
            lines.extend(_format_concise_experiment_block(item, idx, is_in_vivo=is_in_vivo))

    if lines:
        return lines
    if execution_status == "ready_now":
        return ["No module currently meets the ready_now completeness and prerequisite criteria."]
    if execution_status == "conditional_future":
        return ["No conditional_future module is currently specified."]
    if execution_status == "excluded":
        return ["No experiment module is explicitly excluded."]
    return ["No experiment module was returned."]


def _build_candidate_decision_question_lines(plan: dict) -> List[str]:
    candidate = plan.get("candidate") or {}
    bacteria = _strip_reference_noise(candidate.get("bacteria"))
    metabolite = _strip_reference_noise(candidate.get("metabolite"))
    disease = _strip_reference_noise(candidate.get("disease"))
    question = _strip_reference_noise((plan.get("user_brief") or {}).get("research_question"))
    constraints = _strip_reference_noise((plan.get("user_brief") or {}).get("prompt_constraints"))
    lines: List[str] = []
    candidate_parts = []
    if bacteria:
        candidate_parts.append(f"microbe={bacteria}")
    if metabolite:
        candidate_parts.append(f"metabolite={metabolite}")
    if disease:
        candidate_parts.append(f"disease={disease}")
    if candidate_parts:
        lines.append("Candidate: " + " | ".join(candidate_parts))
    if bacteria and metabolite:
        disease_clause = f" in the context of {disease}" if disease else ""
        lines.append(
            f"Decision question: Does {bacteria} directly produce {metabolite}, regulate it indirectly through an ecological route, or act independently or interactively with it{disease_clause}?"
        )
        lines.append(
            f"Core uncertainty: which of these source and host-effect explanations is supported, without presuming that {metabolite} is the only or principal mediator of {bacteria}."
        )
    elif question and not re.search(r"[\u3400-\u9fff]", question):
        lines.append(f"Decision question: {question}")
    else:
        hypothesis = _strip_reference_noise(plan.get("working_hypothesis"))
        if hypothesis and not re.search(r"[\u3400-\u9fff]", hypothesis):
            lines.append(f"Decision question: {hypothesis}")
        else:
            lines.append("Decision question: Resolve the most consequential uncertainty in the user-defined biological relationship.")
    if constraints and not re.search(r"[\u3400-\u9fff]", constraints):
        lines.append(f"User constraints: {constraints}")
    return lines


def _build_evidence_limitation_render_lines(plan: dict, direct_summary: List[str]) -> List[str]:
    lines = [
        _strip_reference_noise(item[len("Limitation:") :])
        for item in direct_summary
        if str(item).strip().lower().startswith("limitation:")
        and _strip_reference_noise(item[len("Limitation:") :])
    ]
    lines.extend(_as_str_list(plan.get("evidence_limitations")))
    assessment = _normalize_production_evidence_assessment(plan.get("direct_production_evidence_assessment"))
    lines.extend(_as_str_list(assessment.get("evidence_limitations")))
    audit = plan.get("protocol_audit") or {}
    for key in ("unsupported_or_overstated_claims", "inference_only_claims"):
        for item in audit.get(key) or []:
            if isinstance(item, dict):
                claim = _strip_reference_noise(item.get("claim"))
                reason = _strip_reference_noise(item.get("reason"))
                rendered = " | ".join(value for value in (claim, reason) if value)
            else:
                rendered = _strip_reference_noise(item)
            if rendered:
                lines.append(rendered)
    return _dedupe_preserve(_strip_reference_noise(item) for item in lines if _strip_reference_noise(item))[:16]


def _build_module_selection_render_lines(plan: dict) -> List[str]:
    selection = plan.get("module_selection") or {}
    lines: List[str] = []
    if isinstance(selection, dict):
        workflow_type = _strip_reference_noise(selection.get("workflow_type"))
        if workflow_type:
            lines.append(f"Workflow type: {workflow_type}")
        for item in selection.get("selected_modules") or []:
            if not isinstance(item, dict):
                continue
            rendered = _render_mapping_summary(
                item,
                [
                    "module_id",
                    "experiment_role",
                    "execution_status",
                    "selection_reason",
                    "prerequisite_module_ids",
                ],
            )
            if rendered:
                lines.append(rendered)
        for item in selection.get("omitted_roles") or []:
            if isinstance(item, dict):
                rendered = _render_mapping_summary(item, ["experiment_role", "reason"])
            else:
                rendered = _strip_reference_noise(item)
            if rendered:
                lines.append(f"Omitted role: {rendered}")
    elif selection:
        lines.extend(_render_loose_value(selection))

    modules = _iter_render_experiment_modules(plan)
    status_counts = {status: 0 for status in ("ready_now", "conditional_future", "excluded")}
    derived_lines: List[str] = []
    for category, item, _ in modules:
        status = _render_execution_status(item)
        status_counts[status] += 1
        module_id = _strip_reference_noise(item.get("module_id")) or "unassigned"
        role = _strip_reference_noise(item.get("experiment_role")) or "unspecified"
        prerequisites = ", ".join(_as_str_list(item.get("prerequisite_module_ids"))) or "none"
        derived_lines.append(
            f"{module_id} | layer={category} | role={role} | status={status} | prerequisites={prerequisites}"
        )
    if modules:
        lines.insert(
            0,
            "Status totals: "
            + ", ".join(f"{status}={status_counts[status]}" for status in status_counts),
        )
    lines.extend(derived_lines)
    return _dedupe_preserve(lines) or ["No module was selected."]


def _build_human_translation_render_lines(plan: dict) -> List[str]:
    gate = plan.get("human_gate") or {}
    lines: List[str] = []
    if isinstance(gate, dict):
        status = _strip_reference_noise(gate.get("status"))
        reason = _strip_reference_noise(gate.get("reason"))
        if status:
            lines.append(f"Human gate status: {status}")
        if reason:
            lines.append(f"Reason: {reason}")
        requirements = gate.get("future_requirements")
        if requirements is None:
            requirements = gate.get("future_requirement")
        requirement_items = _as_str_list(requirements)
        if isinstance(requirements, str) and requirements.strip():
            requirement_items = [requirements.strip()]
        for requirement in requirement_items:
            clean = _strip_reference_noise(requirement)
            if clean:
                lines.append(f"Future requirement: {clean}")
        known_fields = {"status", "reason", "future_requirements", "future_requirement"}
        remaining_gate = {key: value for key, value in gate.items() if key not in known_fields}
        lines.extend(_render_loose_value(remaining_gate))
    elif gate:
        lines.extend(_render_loose_value(gate, label="Human gate"))

    human_items = [
        item
        for category, item, _is_in_vivo in _iter_render_experiment_modules(plan)
        if category == "Human / Translational Experiments"
    ]
    if human_items:
        status_counts = {status: 0 for status in VALID_EXECUTION_STATUSES}
        for item in human_items:
            status_counts[_render_execution_status(item)] += 1
        lines.append(
            "Human modules: "
            + ", ".join(f"{status}={status_counts[status]}" for status in ("ready_now", "conditional_future", "excluded"))
        )
    elif not lines:
        lines.append("No active human module is present; human translation remains gated or outside the current scope.")
    return _dedupe_preserve(lines)


def _build_decision_graph_render_lines(plan: dict, decision_rules: List[str]) -> List[str]:
    lines = list(decision_rules)
    graph = plan.get("decision_graph") or []
    if isinstance(graph, list):
        for item in graph:
            if not isinstance(item, dict):
                clean = _strip_reference_noise(item)
                if clean:
                    lines.append(clean)
                continue
            rendered = _render_mapping_summary(
                item,
                [
                    "module_id",
                    "execution_status",
                    "prerequisite_module_ids",
                    "branch_if_positive",
                    "branch_if_negative",
                ],
            )
            if rendered:
                lines.append(rendered)
    elif graph:
        lines.extend(_render_loose_value(graph))

    if not graph:
        for _, item, _ in _iter_render_experiment_modules(plan):
            module_id = _strip_reference_noise(item.get("module_id")) or "unassigned"
            prerequisites = ", ".join(_as_str_list(item.get("prerequisite_module_ids"))) or "none"
            positive = _strip_reference_noise(item.get("branch_if_positive")) or "not specified"
            negative = _strip_reference_noise(item.get("branch_if_negative")) or "not specified"
            lines.append(
                f"{module_id} | status={_render_execution_status(item)} | prerequisites={prerequisites} "
                f"| if_positive={positive} | if_negative={negative}"
            )
    return _dedupe_preserve(lines)[:24]


def _build_parameter_provenance_audit_render_lines(plan: dict) -> List[str]:
    audit = plan.get("parameter_provenance_audit") or {}
    lines: List[str] = []
    if isinstance(audit, dict):
        status = _strip_reference_noise(audit.get("status"))
        if status:
            lines.append(f"Audit status: {status}")
        for key in ("module_results", "unresolved_items"):
            lines.extend(_render_loose_value(audit.get(key), label=key))
        remaining = {key: value for key, value in audit.items() if key not in {"status", "module_results", "unresolved_items"}}
        lines.extend(_render_loose_value(remaining))
    elif audit:
        lines.extend(_render_loose_value(audit))

    for _, item, _ in _iter_render_experiment_modules(plan):
        module_id = _strip_reference_noise(item.get("module_id")) or "unassigned"
        provenance = _normalize_parameter_provenance(item.get("parameter_provenance"))
        if not provenance:
            lines.append(f"{module_id}: no line-item parameter provenance was supplied.")
            continue
        counts = {status: 0 for status in ("reported", "adapted", "proposed_pilot", "unresolved")}
        for entry in provenance:
            counts[entry.get("status") or "unresolved"] += 1
        lines.append(
            f"{module_id}: " + ", ".join(f"{status}={counts[status]}" for status in counts)
        )
        for entry in provenance:
            if entry.get("status") not in {"adapted", "proposed_pilot", "unresolved"}:
                continue
            parameter = str(entry.get("parameter") or "unspecified parameter").strip()
            value = str(entry.get("value") or "not specified").strip()
            source = str(entry.get("source") or "NONE").strip()
            rationale = _strip_reference_noise(entry.get("transfer_rationale"))
            pilot = _strip_reference_noise(entry.get("pilot_check"))
            detail = [f"{module_id}.{parameter}={value}", f"status={entry.get('status')}", f"source={source}"]
            if rationale:
                detail.append(f"transfer_rationale={rationale}")
            if pilot:
                detail.append(f"pilot_check={pilot}")
            lines.append(" | ".join(detail))
    return _dedupe_preserve(lines)[:32] or ["No parameter-provenance audit was returned."]


def _build_design_completeness_audit_render_lines(plan: dict) -> List[str]:
    audit = plan.get("design_completeness_audit") or {}
    lines: List[str] = []
    if isinstance(audit, dict):
        status = _strip_reference_noise(audit.get("status"))
        if status:
            lines.append(f"Audit status: {status}")
        for key in ("module_results", "issues"):
            lines.extend(_render_loose_value(audit.get(key), label=key))
        remaining = {key: value for key, value in audit.items() if key not in {"status", "module_results", "issues"}}
        lines.extend(_render_loose_value(remaining))
    elif audit:
        lines.extend(_render_loose_value(audit))

    for _, item, _ in _iter_render_experiment_modules(plan):
        module_id = _strip_reference_noise(item.get("module_id")) or "unassigned"
        issues = [
            _strip_reference_noise(value)
            for value in _as_str_list(item.get("completion_issues"))
            if _strip_reference_noise(value)
        ]
        if issues:
            lines.append(
                f"{module_id}: status={_render_execution_status(item)} | incomplete=" + "; ".join(issues)
            )
        else:
            lines.append(f"{module_id}: status={_render_execution_status(item)} | no completion issue recorded.")
    return _dedupe_preserve(lines)[:32] or ["No design-completeness audit was returned."]


def _build_remaining_uncertainty_render_lines(plan: dict, deferred_roadmap: List[str]) -> List[str]:
    lines = _as_str_list(plan.get("remaining_uncertainties"))
    for item in plan.get("self_reflection") or []:
        if isinstance(item, dict):
            uncertainty = _strip_reference_noise(item.get("remaining_uncertainty"))
            if uncertainty:
                lines.append(uncertainty)
    lines.extend(deferred_roadmap)
    return _dedupe_preserve(_strip_reference_noise(item) for item in lines if _strip_reference_noise(item))[:16]


def _build_decision_rule_lines(plan: dict) -> List[str]:
    lines: List[str] = []
    for item in (plan.get("in_vitro_plan") or [])[:3]:
        if not isinstance(item, dict):
            continue
        aim = _strip_reference_noise(item.get("aim") or "In vitro gate")
        endpoint = _strip_reference_noise(item.get("primary_endpoint")) or (
            _strip_reference_noise((_as_str_list(item.get("readouts")) or [""])[0])
        )
        threshold = _strip_reference_noise(item.get("success_threshold"))
        go_no_go = _strip_reference_noise(item.get("go_no_go"))
        failure_action = _strip_reference_noise(item.get("failure_action"))
        mechanism = _strip_reference_noise(item.get("mechanism_support_criterion"))
        phenomenology = _strip_reference_noise(item.get("phenomenology_support_criterion"))
        impact = _strip_reference_noise(item.get("decision_impact"))
        parts = [aim]
        if endpoint:
            parts.append(f"endpoint={endpoint}")
        if threshold:
            parts.append(f"threshold={threshold}")
        if go_no_go:
            parts.append(f"rule={go_no_go}")
        if failure_action:
            parts.append(f"if_fail={failure_action}")
        if mechanism:
            parts.append(f"mechanism_support={mechanism}")
        if phenomenology:
            parts.append(f"phenomenology_only_if={phenomenology}")
        if impact:
            parts.append(f"if_pass={impact}")
        if len(parts) > 1:
            lines.append(" | ".join(parts))
    for item in (plan.get("in_vivo_plan") or [])[:3]:
        if not isinstance(item, dict):
            continue
        aim = _strip_reference_noise(item.get("aim") or "Animal escalation gate")
        endpoint = _strip_reference_noise(item.get("primary_endpoint")) or (
            _strip_reference_noise((_as_str_list(item.get("primary_endpoints")) or [""])[0])
        )
        threshold = _strip_reference_noise(item.get("success_threshold")) or _strip_reference_noise(item.get("go_no_go_threshold"))
        go_no_go = _strip_reference_noise(item.get("go_no_go"))
        failure_action = _strip_reference_noise(item.get("failure_action"))
        mechanism = _strip_reference_noise(item.get("mechanism_support_criterion"))
        phenomenology = _strip_reference_noise(item.get("phenomenology_support_criterion"))
        impact = _strip_reference_noise(item.get("decision_impact"))
        parts = [aim]
        if endpoint:
            parts.append(f"endpoint={endpoint}")
        if threshold:
            parts.append(f"threshold={threshold}")
        if go_no_go:
            parts.append(f"rule={go_no_go}")
        if failure_action:
            parts.append(f"if_fail={failure_action}")
        if mechanism:
            parts.append(f"mechanism_support={mechanism}")
        if phenomenology:
            parts.append(f"phenomenology_only_if={phenomenology}")
        if impact:
            parts.append(f"if_pass={impact}")
        if len(parts) > 1:
            lines.append(" | ".join(parts))
    for item in (plan.get("human_plan") or [])[:2]:
        if not isinstance(item, dict):
            continue
        aim = _strip_reference_noise(item.get("aim") or "Human escalation gate")
        endpoint = _strip_reference_noise(item.get("primary_endpoint")) or (
            _strip_reference_noise((_as_str_list(item.get("readouts")) or [""])[0])
        )
        threshold = _strip_reference_noise(item.get("success_threshold"))
        go_no_go = _strip_reference_noise(item.get("go_no_go"))
        failure_action = _strip_reference_noise(item.get("failure_action"))
        mechanism = _strip_reference_noise(item.get("mechanism_support_criterion"))
        phenomenology = _strip_reference_noise(item.get("phenomenology_support_criterion"))
        impact = _strip_reference_noise(item.get("decision_impact"))
        parts = [aim]
        if endpoint:
            parts.append(f"endpoint={endpoint}")
        if threshold:
            parts.append(f"threshold={threshold}")
        if go_no_go:
            parts.append(f"rule={go_no_go}")
        if failure_action:
            parts.append(f"if_fail={failure_action}")
        if mechanism:
            parts.append(f"mechanism_support={mechanism}")
        if phenomenology:
            parts.append(f"phenomenology_only_if={phenomenology}")
        if impact:
            parts.append(f"if_pass={impact}")
        if len(parts) > 1:
            lines.append(" | ".join(parts))
    return _dedupe_preserve(lines)[:8] or ["Advance only if the first gatekeeper experiment shows a reproducible effect in the predicted direction; otherwise stop escalation and revise the hypothesis."]


def _build_future_complete_study_lines(plan: dict) -> List[str]:
    lines: List[str] = []
    profile = _normalize_question_profile(plan.get("question_profile"))
    policy = profile.get("human_study_policy", "not_recommended")
    detail_budget = profile.get("detail_budget_allocation", "balanced_preclinical")
    question_type = profile.get("question_type", "general_validation")
    for item in (plan.get("in_vivo_plan") or [])[:3]:
        if not isinstance(item, dict):
            continue
        endpoints = [v for v in (_strip_reference_noise(x) for x in _as_str_list(item.get("mechanistic_endpoints"))[:3]) if v]
        if endpoints:
            lines.append("If the gatekeeper study passes, expand to mechanism-oriented endpoints: " + ", ".join(endpoints))
    remaining_uncertainties = []
    for item in (plan.get("self_reflection") or [])[:3]:
        if isinstance(item, dict):
            remaining = _strip_reference_noise(item.get("remaining_uncertainty"))
            if remaining:
                remaining_uncertainties.append(remaining)
    if remaining_uncertainties:
        lines.append("The complete follow-up study should resolve the remaining uncertainties: " + "; ".join(remaining_uncertainties[:3]))
    lines.append("A complete study should add dose-response, time-course, replication in an independent system, and deeper mechanism mapping only after the gatekeeper experiments succeed.")
    if policy == "allow_observational_only":
        lines.append(
            "The human extension should remain observational: expand to an independent cohort, denser longitudinal biospecimen sampling, stronger covariate control, and external replication before making any causal human claim."
        )
    elif policy == "allow_low_risk_behavior_only":
        lines.append(
            "The late-stage human package should stay within low-risk operational testing: add standardized crossover sessions, independent replication sites, exposure verification, and repeated mediator or behavioral sampling before making a stronger translational claim."
        )
    elif policy == "allow_late_stage_rct":
        lines.append(
            "If the preclinical signal is reproducible and the safety profile is acceptable, escalate to a randomized, placebo-controlled human study with predefined primary endpoints, exposure verification, and mediator measurement."
        )
        lines.append(
            "The late-stage causal package should add crossover or washout-rechallenge testing, longitudinal mediator sampling, and formal mediation analysis to test whether the human endpoint tracks the proposed mechanism."
        )
    elif question_type == "sensor_discovery" or detail_budget == "discovery_heavy":
        lines.append(
            "The complete study should deepen target-confidence rather than human translation: extend to orthogonal target-identification assays, rescue logic, structure-function perturbation, and in vivo necessity/sufficiency testing in an independent system."
        )
    elif question_type == "toxic_exposure":
        lines.append(
            "The late-stage package should avoid direct human exposure and instead add orthogonal exposure assessment, environmental measurement, biomonitoring, and stronger causal triangulation across cell, animal, and observational datasets."
        )
    return _dedupe_preserve([line for line in lines if line])[:4]


def _build_validation_summary(plan: dict) -> dict:
    summary = {
        "known_evidence": _build_known_evidence_lines(plan),
        "direct_production_evidence": _build_direct_production_evidence_lines(plan),
        "hypothesis_branches": _build_hypothesis_branch_lines(plan),
        "missing_evidence_gap": _build_missing_gap_lines(plan),
        "why_this_experiment": [],
        "priority_follow_up_experiments": _build_priority_followup_lines(plan),
        "decision_rule": _build_decision_rule_lines(plan),
        "future_complete_study": _build_future_complete_study_lines(plan),
    }
    summary["why_this_experiment"] = _build_why_this_experiment_lines({**plan, "validation_summary": summary})
    return _run_final_summary_consistency_review(summary)


def _run_final_summary_consistency_review(summary: dict) -> dict:
    reviewed: dict[str, List[str]] = {}
    for key, value in (summary or {}).items():
        if isinstance(value, list):
            cleaned_items: List[str] = []
            for item in value:
                cleaned = (
                    _strip_reference_noise(item)
                    if key in {"known_evidence", "direct_production_evidence", "hypothesis_branches"}
                    else _soften_overclaim_language(item)
                )
                if cleaned:
                    cleaned_items.append(cleaned)
            reviewed[key] = _dedupe_preserve(cleaned_items)
        else:
            reviewed[key] = value

    known_keys = {_claim_key(item) for item in reviewed.get("known_evidence", []) if _claim_key(item)}
    reviewed["missing_evidence_gap"] = [
        item for item in reviewed.get("missing_evidence_gap", []) if _claim_key(item) not in known_keys
    ]
    if not reviewed.get("known_evidence"):
        reviewed["known_evidence"] = [
            "No clearly supported literature-backed claim met the inclusion threshold after consistency review."
        ]

    for section in (
        "direct_production_evidence",
        "hypothesis_branches",
        "missing_evidence_gap",
        "why_this_experiment",
        "priority_follow_up_experiments",
        "decision_rule",
        "future_complete_study",
    ):
        limit = 14 if section == "direct_production_evidence" else 8
        reviewed[section] = reviewed.get(section, [])[:limit]
    known_evidence = reviewed.get("known_evidence", [])
    coverage_lines = [line for line in known_evidence if line.startswith("Search coverage:")]
    if coverage_lines and coverage_lines[0] not in known_evidence[:6]:
        reviewed["known_evidence"] = known_evidence[:5] + [coverage_lines[0]]
    else:
        reviewed["known_evidence"] = known_evidence[:6]
    return reviewed


def _activation_gate_summary(gate: Any) -> str:
    if not isinstance(gate, dict) or not gate:
        return ""
    clauses = gate.get("clauses")
    operator = str(gate.get("operator") or "all_of").strip().lower()
    if isinstance(clauses, list):
        rendered = [value for value in (_activation_gate_summary(clause) for clause in clauses) if value]
        if not rendered:
            return ""
        joiner = " AND " if operator == "all_of" else " OR "
        return "(" + joiner.join(rendered) + ")"
    accepted = _as_str_list(gate.get("accepted_results"))
    expected = "/".join(accepted) or "a required result"
    module_id = _strip_reference_noise(gate.get("module_id"))
    external_gate = _strip_reference_noise(gate.get("external_gate"))
    if module_id:
        return f"{module_id}={expected}"
    if external_gate:
        return f"{external_gate}={expected}"
    return ""


def _render_stage_modules(plan: dict, plan_key: str, is_in_vivo: bool = False) -> List[str]:
    items = [
        item
        for item in (plan.get(plan_key) or [])
        if isinstance(item, dict) and _render_execution_status(item) != "excluded"
    ]
    if not items:
        return ["No module was selected for this stage under the current scientific question."]
    lines: List[str] = []
    for index, item in enumerate(items, start=1):
        if lines:
            lines.append("")
        lines.extend(_format_concise_experiment_block(item, index, is_in_vivo=is_in_vivo))
    return lines


def _build_route_priority_render_lines(plan: dict) -> List[str]:
    lines: List[str] = []
    for plan_key in ("in_vitro_plan", "in_vivo_plan", "human_plan"):
        for item in plan.get(plan_key) or []:
            if not isinstance(item, dict) or _render_execution_status(item) == "excluded":
                continue
            module_id = _strip_reference_noise(item.get("module_id")) or "unassigned"
            route_status = _strip_reference_noise(item.get("route_status")) or "deferred"
            question = _strip_reference_noise(
                item.get("scientific_question") or item.get("biological_question") or item.get("aim")
            )
            hypothesis_ids = _dedupe_preserve(
                hypothesis_id.upper()
                for hypothesis_id in _as_str_list(item.get("hypothesis_ids"))
                if re.fullmatch(r"H\d+", hypothesis_id.strip(), flags=re.IGNORECASE)
            )
            unlock = _canonical_gate_contract(
                item.get("experiment_role"), item.get("activation_gate")
            )["unlock_rule"]
            line = f"{module_id} | {route_status} | {question or 'question not specified'}"
            if hypothesis_ids:
                line += f" | hypothesis: {', '.join(hypothesis_ids)}"
            if unlock:
                line += f" | unlock: {unlock}"
            lines.append(line)
    return _dedupe_preserve(lines)[:12] or ["No experimental route was selected."]


def _build_hypothesis_experiment_mapping_lines(plan: dict) -> List[str]:
    modules_by_hypothesis: Dict[str, List[str]] = {}
    for plan_key in ("in_vitro_plan", "in_vivo_plan"):
        for item in plan.get(plan_key) or []:
            if not isinstance(item, dict) or _render_execution_status(item) == "excluded":
                continue
            module_id = _strip_reference_noise(item.get("module_id")) or "unassigned"
            role = _strip_reference_noise(item.get("experiment_role")) or "unspecified"
            status = _render_execution_status(item)
            label = f"{module_id} ({role}; {status})"
            for hypothesis_id in _as_str_list(item.get("hypothesis_ids")):
                canonical_id = hypothesis_id.strip().upper()
                if re.fullmatch(r"H\d+", canonical_id):
                    modules_by_hypothesis.setdefault(canonical_id, []).append(label)

    lines: List[str] = []
    for branch in _normalize_hypothesis_branches(plan.get("hypothesis_branches")):
        hypothesis_id = str(branch.get("hypothesis_id") or "").strip().upper()
        if not re.fullmatch(r"H\d+", hypothesis_id):
            continue
        modules = _dedupe_preserve(modules_by_hypothesis.get(hypothesis_id, []))
        if modules:
            lines.append(
                f"Hypothesis-to-experiment mapping: {hypothesis_id} -> {', '.join(modules)}."
            )
        else:
            lines.append(
                f"Hypothesis-to-experiment mapping: {hypothesis_id} -> no experimental module was returned; keep this branch unresolved rather than inferring support."
            )
    return lines


def _build_evidence_gap_and_hypothesis_lines(summary: dict, plan: Optional[dict] = None) -> List[str]:
    lines = [f"Evidence gap: {item}" for item in _as_str_list(summary.get("missing_evidence_gap"))[:8]]
    lines.extend(
        f"Competing hypothesis: {item}" for item in _as_str_list(summary.get("hypothesis_branches"))[:6]
    )
    if isinstance(plan, dict):
        lines.extend(_build_hypothesis_experiment_mapping_lines(plan))
    return _dedupe_preserve(lines) or ["The evidence gaps and competing hypotheses were not resolved."]


def _build_decision_alternative_render_lines(plan: dict) -> List[str]:
    lines: List[str] = []
    for plan_key in ("in_vitro_plan", "in_vivo_plan", "human_plan"):
        for item in plan.get(plan_key) or []:
            if not isinstance(item, dict) or _render_execution_status(item) == "excluded":
                continue
            module_id = _strip_reference_noise(item.get("module_id")) or "unassigned"
            positive_gate = _strip_reference_noise(
                item.get("positive_gate") or item.get("success_threshold") or item.get("go_no_go")
            )
            positive = _strip_reference_noise(
                item.get("branch_if_positive") or item.get("positive_result_interpretation")
            )
            negative = _strip_reference_noise(
                item.get("branch_if_negative") or item.get("negative_result_interpretation")
            )
            parts = [module_id]
            if positive_gate:
                parts.append(f"positive gate: {positive_gate}")
            if positive:
                parts.append(f"if positive: {positive}")
            if negative:
                parts.append(f"if negative: {negative}")
            lines.append(" | ".join(parts))
    return _dedupe_preserve(lines)[:12] or ["No branch-specific decision rule was returned."]


def _build_core_uncertainty_render_lines(plan: dict) -> List[str]:
    lines = _as_str_list(plan.get("remaining_uncertainties"))
    for plan_key in ("in_vitro_plan", "in_vivo_plan", "human_plan"):
        for item in plan.get(plan_key) or []:
            if not isinstance(item, dict):
                continue
            module_id = _strip_reference_noise(item.get("module_id")) or "unassigned"
            route_status = str(item.get("route_status") or "").strip().lower()
            if route_status not in {"alternative", "deferred"}:
                continue
            question = _strip_reference_noise(
                item.get("scientific_question") or item.get("biological_question") or item.get("aim")
            )
            lines.append(f"Deferred branch {module_id}: {question or 'retain only if an upstream result requires it.'}")
    if not lines:
        lines.append("The main remaining uncertainty is which source or host-effect branch will pass its first decisive gate.")
    return _dedupe_preserve(lines)[:10]


def _render_validation_protocol_text(plan: dict) -> str:
    summary = plan.get("validation_summary") or _build_validation_summary(plan)
    if not isinstance(summary, dict):
        summary = _build_validation_summary(plan)

    human_lines = _build_human_translation_render_lines(plan)
    human_modules = _render_stage_modules(plan, "human_plan", is_in_vivo=True)
    if plan.get("human_plan"):
        human_lines.extend([""] + human_modules)

    sections: List[tuple[str, List[str], bool]] = [
        ("Core scientific question", _build_candidate_decision_question_lines(plan), True),
        (
            "Known evidence and supporting literature",
            _as_str_list(summary.get("known_evidence")),
            True,
        ),
        (
            "Evidence gaps and competing hypotheses",
            _build_evidence_gap_and_hypothesis_lines(summary, plan),
            True,
        ),
        ("Overall experimental route and priorities", _build_route_priority_render_lines(plan), True),
        ("In vitro and cell studies", _render_stage_modules(plan, "in_vitro_plan"), False),
        ("Animal studies", _render_stage_modules(plan, "in_vivo_plan", is_in_vivo=True), False),
        ("Human study, conditional", human_lines, False),
        (
            "Decision rules, alternatives, and deferred studies",
            _build_decision_alternative_render_lines(plan),
            True,
        ),
        ("Core uncertainty and future roadmap", _build_core_uncertainty_render_lines(plan), True),
    ]

    lines: List[str] = []
    for heading, items, use_bullets in sections:
        if lines:
            lines.append("")
        lines.append(heading)
        if use_bullets:
            lines.extend(_list_to_lines(items, prefix="- ", limit=16) or ["- Not available."])
        else:
            lines.extend(items or ["Not available."])
    return "\n".join(line.rstrip() for line in lines if line is not None).strip()


def _quoted_or_clause(values: List[str]) -> str:
    terms = [f'"{item}"' for item in _dedupe_preserve(values) if item]
    return " OR ".join(terms[:6]) or '""'


def _build_round_one_query_specs(bacteria: str, metabolite: str, disease: str) -> List[dict]:
    b_expr = build_pubmed_query([bacteria])
    m_expr = build_pubmed_query([metabolite])
    d_expr = build_pubmed_query([disease])
    return [
        {
            "round": 1,
            "label": "Direct monoculture production evidence",
            "rationale": "Find a single study that directly measures the metabolite after culturing the named microbe under controlled culture conditions.",
            "query": f"{b_expr} AND {m_expr} AND (monoculture OR \"pure culture\" OR \"defined medium\" OR \"culture supernatant\" OR fermentation) AND (production OR produced OR quantified OR metabolomics OR \"mass spectrometry\")",
        },
        {
            "round": 1,
            "label": "Indirect ecological production",
            "rationale": "Find evidence for co-culture, cross-feeding, substrate release, conditioned-medium effects, or community remodeling that could change metabolite production indirectly.",
            "query": f"{b_expr} AND {m_expr} AND (\"co-culture\" OR coculture OR \"cross-feeding\" OR \"conditioned medium\" OR \"defined community\" OR consortium OR substrate OR mucin)",
        },
        {
            "round": 1,
            "label": "Metabolite and disease intervention evidence",
            "rationale": "Find direct named-metabolite exposure in a disease-relevant cell, tissue, or animal model.",
            "query": f"{m_expr} AND {d_expr} AND (administered OR treated OR supplementation OR exposure OR addition) AND (\"in vitro\" OR epithelial OR macrophage OR organoid OR mouse OR mice OR animal)",
        },
        {
            "round": 1,
            "label": "Microbe and disease intervention evidence",
            "rationale": "Identify whether the microbe itself changes disease outcomes independently of the candidate metabolite.",
            "query": f"{b_expr} AND {d_expr} AND (intervention OR administration OR supplementation OR mouse OR mice OR colitis)",
        },
    ]


def _build_candidate_missing_edge_query_specs(
    bacteria: str,
    metabolite: str,
    disease: str,
    evidence_strength_map: List[dict],
) -> List[dict]:
    """Run one targeted second-round query for each unresolved candidate edge."""
    support = {
        str(item.get("topic") or ""): str(item.get("support_level") or "speculative")
        for item in evidence_strength_map
        if isinstance(item, dict)
    }
    base_specs = _build_round_one_query_specs(bacteria, metabolite, disease)
    topic_by_label = {
        "Direct monoculture production evidence": "Direct Monoculture Production Evidence",
        "Indirect ecological production": "Indirect Ecological Production Evidence",
        "Metabolite and disease intervention evidence": "Metabolite -> Disease Host Response",
        "Microbe and disease intervention evidence": "Microbe -> Disease Intervention",
    }
    discriminator_by_label = {
        "Direct monoculture production evidence": "(time course OR baseline OR uninoculated OR isotope OR precursor)",
        "Indirect ecological production": "(producer OR depletion OR removal OR tracing OR substrate transfer)",
        "Metabolite and disease intervention evidence": "(administered OR treated OR dosed OR exposed)",
        "Microbe and disease intervention evidence": "(administered OR gavage OR colonized OR supplemented)",
    }
    specs: List[dict] = []
    for base in base_specs:
        label = str(base.get("label") or "")
        topic = topic_by_label.get(label)
        if not topic or support.get(topic) == "direct":
            continue
        specs.append(
            {
                "round": 2,
                "label": f"Unresolved edge follow-up: {label}",
                "rationale": f"No qualifying direct evidence was identified for {topic}; repeat with stricter manipulation and control terms.",
                "query": f"({base.get('query')}) AND {discriminator_by_label[label]}",
            }
        )
    return specs


def _fallback_question_query(research_question: str, disease: str) -> str:
    clean_question = _strip_reference_noise(research_question).replace('"', "")
    keyword_tokens = _extract_question_keywords(clean_question)
    method_block = "(\"in vitro\" OR cell OR animal OR mouse OR methods OR protocol OR cohort OR patient)"
    validation_block = "(validation OR experiment OR protocol OR methods)"
    keyword_query = build_pubmed_query(keyword_tokens[:5])
    disease_block = build_pubmed_query([disease]) if str(disease or "").strip() else ""

    if keyword_query and disease_block:
        return f"({keyword_query}) AND {disease_block} AND {method_block}"
    if keyword_query:
        return f"({keyword_query}) AND {method_block}"
    if disease_block:
        return f"{disease_block} AND {validation_block} AND {method_block}"
    return f"{validation_block} AND {method_block}"


def build_question_only_query_specs(
    research_question: str,
    prompt_constraints: str,
    disease: str = "",
) -> List[dict]:
    fallback_query = _fallback_question_query(research_question, disease)
    keyword_specs = _build_question_keyword_query_specs(
        research_question=research_question,
        prompt_constraints=prompt_constraints,
        disease=disease,
        max_specs=8,
    )
    specs: List[dict] = []
    seen_queries: set[str] = set()
    _append_query_spec(
        specs,
        seen_queries,
        fallback_query,
        "Question-driven evidence search",
        "Search for literature directly addressing keywords from the user's experimental question.",
    )
    for spec in keyword_specs:
        _append_query_spec(
            specs,
            seen_queries,
            spec.get("query", ""),
            str(spec.get("label") or "Question keyword search"),
            str(spec.get("rationale") or "Keyword-derived literature search."),
            round_id=int(spec.get("round") or 1),
        )
        if len(specs) >= 8:
            break
    return specs[:8]


def search_literature_online(query_specs: List[dict], max_results: int = MAX_QUERY_RESULTS) -> dict:
    articles_by_key: Dict[str, dict] = {}
    query_log: List[dict] = []
    for spec in query_specs[:8]:
        query = str(spec.get("query") or "").strip()
        if not query:
            continue
        try:
            pmids = search_pubmed(query, max_results=max_results)
            abstracts = fetch_abstracts(pmids[:max_results])
        except Exception:
            pmids = []
            abstracts = []
        query_log.append(
            {
                "round": int(spec.get("round") or 0),
                "label": str(spec.get("label") or "").strip(),
                "query": query,
                "rationale": str(spec.get("rationale") or "").strip(),
                "pmid_count": len(pmids),
                "article_count": len(abstracts),
            }
        )
        for article in abstracts:
            key = str(article.get("pmid") or article.get("title") or "").strip()
            if not key:
                continue
            enriched = dict(article)
            enriched["topic"] = spec.get("label") or article.get("topic") or "online_search"
            enriched["query_round"] = int(spec.get("round") or 0)
            enriched["query_label"] = str(spec.get("label") or "").strip()
            enriched["query"] = query
            articles_by_key.setdefault(key, enriched)
    return {
        "articles": list(articles_by_key.values()),
        "query_log": query_log,
    }


def _merge_articles(*article_lists: List[dict]) -> List[dict]:
    merged = []
    seen = set()
    for article_list in article_lists:
        for article in article_list or []:
            key = str(article.get("pmid") or article.get("title") or "").strip().lower()
            if not key or key in seen:
                continue
            seen.add(key)
            merged.append(article)
    return merged


def rank_evidence_strength(
    articles: List[dict],
    bacteria: str,
    metabolite: str,
    disease: str,
) -> List[dict]:
    bacteria_aliases = _citation_entity_aliases(bacteria, is_microbe=True)
    metabolite_aliases = _citation_entity_aliases(metabolite)
    disease_aliases = _citation_disease_aliases(disease)
    facets = {
        "Direct Monoculture Production Evidence": {"pmids": [], "direct": 0, "indirect": 0, "speculative": 0, "note": ""},
        "Indirect Ecological Production Evidence": {"pmids": [], "direct": 0, "indirect": 0, "speculative": 0, "note": ""},
        "Microbe -> Metabolite": {"pmids": [], "direct": 0, "indirect": 0, "speculative": 0, "note": ""},
        "Microbe -> Disease Intervention": {"pmids": [], "direct": 0, "indirect": 0, "speculative": 0, "note": ""},
        "Metabolite -> Disease Host Response": {"pmids": [], "direct": 0, "indirect": 0, "speculative": 0, "note": ""},
        "Cell / In Vitro Evidence": {"pmids": [], "direct": 0, "indirect": 0, "speculative": 0, "note": ""},
        "Animal / In Vivo Evidence": {"pmids": [], "direct": 0, "indirect": 0, "speculative": 0, "note": ""},
        "Contradictory / Null Evidence": {"pmids": [], "direct": 0, "indirect": 0, "speculative": 0, "note": ""},
    }

    for article in articles[:24]:
        text = f"{article.get('title', '')} {article.get('abstract', '')}"
        has_bacteria = _contains_entity(text, bacteria_aliases)
        has_metabolite = _contains_entity(text, metabolite_aliases)
        has_disease = _contains_entity(text, disease_aliases)
        has_cell = _contains_any(text, CELL_TERMS)
        has_animal = _contains_any(text, ANIMAL_TERMS)
        is_contradictory = _contains_any(text, CONTRADICT_TERMS)
        pmid = str(article.get("pmid") or "").strip()
        evidence_type, _, citation_eligible = _classify_candidate_article_evidence(
            text,
            bacteria_aliases,
            metabolite_aliases,
            disease_aliases,
        )

        if evidence_type == "direct_monoculture_production":
            facets["Microbe -> Metabolite"]["direct"] += 1
            facets["Direct Monoculture Production Evidence"]["direct"] += 1
            if pmid:
                facets["Microbe -> Metabolite"]["pmids"].append(pmid)
                facets["Direct Monoculture Production Evidence"]["pmids"].append(pmid)
        elif evidence_type == "direct_monoculture_nonproduction":
            facets["Contradictory / Null Evidence"]["direct"] += 1
            if pmid:
                facets["Contradictory / Null Evidence"]["pmids"].append(pmid)
        elif evidence_type == "indirect_ecological_evidence":
            facets["Microbe -> Metabolite"]["indirect"] += 1
            facets["Indirect Ecological Production Evidence"]["direct"] += 1
            if pmid:
                facets["Microbe -> Metabolite"]["pmids"].append(pmid)
                facets["Indirect Ecological Production Evidence"]["pmids"].append(pmid)
        elif has_bacteria and has_metabolite:
            facets["Microbe -> Metabolite"]["speculative"] += 1

        if evidence_type == "candidate_microbe_disease_intervention":
            facets["Microbe -> Disease Intervention"]["direct"] += 1
            if pmid:
                facets["Microbe -> Disease Intervention"]["pmids"].append(pmid)
        elif has_bacteria and has_disease:
            facets["Microbe -> Disease Intervention"]["speculative"] += 1

        if evidence_type == "candidate_metabolite_disease_intervention":
            facets["Metabolite -> Disease Host Response"]["direct"] += 1
            if pmid:
                facets["Metabolite -> Disease Host Response"]["pmids"].append(pmid)
        elif has_metabolite and has_disease:
            facets["Metabolite -> Disease Host Response"]["speculative"] += 1
        if has_cell and citation_eligible:
            level = "direct" if has_disease else "indirect"
            facets["Cell / In Vitro Evidence"][level] += 1
            if pmid:
                facets["Cell / In Vitro Evidence"]["pmids"].append(pmid)
        if has_animal and citation_eligible:
            level = "direct" if has_disease else "indirect"
            facets["Animal / In Vivo Evidence"][level] += 1
            if pmid:
                facets["Animal / In Vivo Evidence"]["pmids"].append(pmid)
        if is_contradictory and evidence_type != "direct_monoculture_nonproduction":
            facets["Contradictory / Null Evidence"]["direct"] += 1
            if pmid:
                facets["Contradictory / Null Evidence"]["pmids"].append(pmid)

    evidence_map = []
    for topic, bucket in facets.items():
        if bucket["direct"] > 0:
            support_level = "direct"
        elif bucket["indirect"] > 0:
            support_level = "indirect"
        else:
            support_level = "speculative"
        rationale = f"Direct={bucket['direct']}, indirect={bucket['indirect']}, speculative={bucket['speculative']} across retrieved literature."
        if bucket["note"]:
            rationale = f"{rationale} {bucket['note']}"
        evidence_map.append(
            {
                "topic": topic,
                "support_level": support_level,
                "evidence_type": topic.lower().replace(" / ", " ").replace("->", "to"),
                "rationale": rationale,
                "representative_pmids": _dedupe_preserve(bucket["pmids"])[:6],
            }
        )
    return evidence_map


def _format_article_context(articles: List[dict], limit: int = 10) -> str:
    if not articles:
        return "No literature abstracts available."
    chunks = []
    for article in articles[:limit]:
        label = article.get("query_label") or article.get("topic") or "literature"
        round_label = f"round {article.get('query_round')}" if article.get("query_round") else "seed"
        chunks.append(
            f"[{label}; {round_label}] PMID {article.get('pmid', '')}: {article.get('title', '')} "
            f"({article.get('journal', '')}, {article.get('year', '')})\n"
            f"{_truncate(article.get('abstract', ''), 360)}"
        )
    return "\n\n".join(chunks)


def _compact_protocol_text(protocol_text: str, limit: int = 1600) -> str:
    return _truncate(_strip_reference_noise(protocol_text), limit)


def _compact_experiment_item(item: dict, is_in_vivo: bool = False) -> dict:
    if not isinstance(item, dict):
        return {}
    keep_fields = [
        "module_id",
        "experiment_role",
        "hypothesis_ids",
        "route_status",
        "execution_status",
        "result_status",
        "design_status",
        "activation_gate",
        "scientific_question",
        "why_needed",
        "study_object",
        "groups",
        "samples_and_timing",
        "primary_indicator",
        "secondary_indicators",
        "key_controls",
        "positive_gate",
        "branch_if_positive",
        "branch_if_negative",
        "unlock_rule",
        "claim_boundary",
        "source_citations",
    ]
    compact = {}
    for field in keep_fields:
        value = item.get(field)
        if isinstance(value, dict):
            if value:
                compact[field] = value
            continue
        if isinstance(value, list):
            if field == "groups":
                compact[field] = value[:12]
            else:
                compact[field] = _dedupe_preserve(
                    [_truncate(str(v), 160) for v in value if str(v).strip()]
                )[:8]
        elif str(value or "").strip():
            compact[field] = _truncate(str(value), 320)
    return compact


def _compact_experiment_bundle(items: List[dict], is_in_vivo: bool = False, limit: int = 2) -> str:
    compact_items = [_compact_experiment_item(item, is_in_vivo=is_in_vivo) for item in items[:limit] if isinstance(item, dict)]
    if not compact_items:
        return "[]"
    return json.dumps(compact_items, ensure_ascii=False, indent=2)


def _compact_protocols_io_summary(protocols_io_evidence: Dict[str, Any], limit: int = 2) -> str:
    if not isinstance(protocols_io_evidence, dict):
        return "{}"
    compact = {
        "status": protocols_io_evidence.get("status"),
        "message": _truncate(protocols_io_evidence.get("message") or "", 240),
        "queries": _dedupe_preserve([str(v) for v in (protocols_io_evidence.get("queries") or [])])[:6],
        "all": [],
    }
    for item in (protocols_io_evidence.get("all") or [])[:limit]:
        if not isinstance(item, dict):
            continue
        compact["all"].append(
            {
                "citation": _truncate(item.get("citation") or "", 120),
                "title": _truncate(item.get("title") or "", 140),
                "url": _truncate(item.get("url") or "", 180),
                "materials": _dedupe_preserve([str(v) for v in (item.get("materials") or [])])[:5],
                "steps": _dedupe_preserve([_truncate(str(v), 140) for v in (item.get("steps") or [])])[:5],
            }
        )
    return json.dumps(compact, ensure_ascii=False, indent=2)


def _format_query_log(query_log: List[dict]) -> str:
    if not query_log:
        return "No iterative queries were executed."
    return "\n".join(
        f"Round {item.get('round')}: {item.get('label')} | PMIDs={item.get('pmid_count')} | Articles={item.get('article_count')} | Query={item.get('query')}"
        for item in query_log[:10]
    )


def _format_evidence_strength_map(evidence_strength_map: List[dict]) -> str:
    if not evidence_strength_map:
        return "No structured evidence strength map available."
    return "\n".join(
        f"- {item.get('topic')}: {item.get('support_level')} | PMIDs: {', '.join(item.get('representative_pmids') or ['none'])} | {item.get('rationale')}"
        for item in evidence_strength_map[:8]
    )


def _format_production_evidence_assessment(value: Any) -> str:
    assessment = _normalize_production_evidence_assessment(value)
    lines = [
        f"Status: {assessment.get('status')}",
        f"Conclusion: {assessment.get('conclusion') or 'Not available'}",
    ]
    if assessment.get("direct_evidence"):
        lines.append("Direct culture evidence:")
        for item in assessment["direct_evidence"][:6]:
            lines.append(
                f"- PMID {item.get('pmid')}: {item.get('claim')} | model={item.get('model_system') or 'not stated'} | "
                f"measurement={item.get('measured_output') or 'not stated'} | reason={item.get('why_this_support_level') or 'not stated'}"
            )
    if assessment.get("paper_findings"):
        lines.append("Paper-by-paper findings (one paper per item):")
        for item in assessment["paper_findings"][:10]:
            lines.append(
                f"- PMID {item.get('pmid')}: {item.get('claim')} | evidence_type={item.get('evidence_type')} | "
                f"reason={item.get('why_this_support_level') or 'not stated'}"
            )
    if assessment.get("evidence_limitations"):
        lines.append("Limitations:")
        lines.extend(f"- {item}" for item in assessment["evidence_limitations"][:6])
    return "\n".join(lines)


def _format_hypothesis_branches(items: Any) -> str:
    branches = _normalize_hypothesis_branches(items)
    if not branches:
        return "No structured competing hypotheses are available."
    chunks = []
    for item in branches[:6]:
        chunks.append(
            f"- {item.get('hypothesis_id')}: {item.get('statement')}\n"
            f"  evidence_status={item.get('current_evidence_status')}\n"
            f"  evidence_basis={'; '.join(item.get('evidence_basis') or ['none'])}\n"
            f"  discriminating_prediction={item.get('discriminating_prediction') or 'not specified'}\n"
            f"  in_vitro_gate={item.get('in_vitro_gate') or 'not specified'}\n"
            f"  animal_gate={item.get('animal_gate') or 'not specified'}\n"
            f"  human_gate={item.get('human_gate') or 'not specified'}\n"
            f"  redirect_if_not_supported={item.get('falsification_or_redirection') or 'not specified'}"
        )
    return "\n".join(chunks)


def _has_user_research_question(research_question: str) -> bool:
    return bool((research_question or "").strip())


def _format_user_brief(research_question: str, prompt_constraints: str) -> str:
    question = (research_question or "").strip()[:4000]
    constraints = (prompt_constraints or "").strip()[:6000]
    if _has_user_research_question(research_question):
        priority = (
            "Instruction priority: The user's research question is the primary organizing principle for this validation plan. "
            "Structure audit priorities, in vitro/in vivo experiments, and conclusions mainly to address that question within "
            "the candidate microbe-metabolite-disease context. Use retrieved literature to support, refine, or constrain the "
            "design, but do not replace the user's focus with a generic validation template when literature is sparse. "
            "Never invent citations or unsupported experimental details. If evidence is weak, state that explicitly while "
            "still framing experiments around the user's question."
        )
    else:
        priority = (
            "Instruction priority: retrieved literature and scientific accuracy override the user brief. "
            "Use the brief to set focus, scope, preferred models, exclusions, and output emphasis, but never "
            "invent experimental details or citations to satisfy it."
        )
    return (
        f"Research question:\n{question or 'No additional research question provided.'}\n\n"
        f"Prompt constraints:\n{constraints or 'No additional prompt constraints provided.'}\n\n"
        f"{priority}"
    )


def _user_question_plan_directive(research_question: str) -> str:
    if not _has_user_research_question(research_question):
        return ""
    return (
        "\n\nUser-question-first mode:\n"
        "- Organize all outputs primarily to answer the research question above.\n"
        "- Choose models, readouts, grouping, and go/no-go criteria based on what the question asks.\n"
        "- Use literature to operationalize the question (doses, models, timing), not to redirect away from it.\n"
        "- Keep the candidate context, but do not default to generic validation when the user asked something specific."
    )


def _vitro_priority_block(research_question: str) -> str:
    if _has_user_research_question(research_question):
        return (
            "1. experiments that directly answer the user's research question in vitro or in cell models\n"
            "2. supporting microbe-to-metabolite or metabolite-to-host checks only when needed for that question\n"
            "3. controls required to interpret the user-focused readouts, not a generic validation checklist"
        )
    return (
        "1. direct microbe-to-metabolite production gate in controlled culture\n"
        "2. an indirect-production branch using conditioned medium, co-culture, cross-feeding, or a defined community if the direct gate fails\n"
        "3. a separate metabolite-to-host response experiment that does not assume mediation\n"
        "4. only the controls needed to distinguish the competing hypotheses and exclude media carryover or nonspecific effects"
    )


def _vivo_priority_block(research_question: str) -> str:
    if _has_user_research_question(research_question):
        return (
            "- explain whether an animal model is justified for answering the user's research question\n"
            "- define the exact no-go threshold tied to the user's question, not generic disease outcomes alone\n"
            "- prioritize models and readouts from retrieved literature that best operationalize the user's question"
        )
    return (
        "- first test effects and interaction without claiming that the metabolite mediates the microbe's action\n"
        "- only after a branch-specific in vitro gate passes, test causality using function loss plus rescue or a defined-community ecological design\n"
        "- require comparable microbial exposure or colonization when comparing strains or communities\n"
        "- define the exact no-go threshold that would stop animal escalation or redirect to another hypothesis\n"
        "- prioritize models and operational conditions used in the retrieved literature when possible"
    )


def build_user_focus_query_spec(
    bacteria: str,
    metabolite: str,
    disease: str,
    research_question: str,
    prompt_constraints: str,
) -> Optional[dict]:
    if not (research_question or "").strip():
        return None

    entity_query = build_pubmed_query([bacteria, metabolite, disease])
    source_text = " ".join(
        part for part in [research_question, prompt_constraints] if str(part or "").strip()
    )
    keyword_query = build_pubmed_query(_extract_question_keywords(source_text, limit=4))
    query = " AND ".join(part for part in [entity_query, keyword_query] if part)
    if not query:
        return None
    return {
        "round": 1,
        "label": "User-defined research focus",
        "rationale": "Search for evidence directly addressing the user's Step 3 question.",
        "query": query[:600],
    }


def rank_question_evidence_strength(
    articles: List[dict],
    disease: str,
) -> List[dict]:
    disease_aliases = _disease_aliases(disease) if disease else []
    method_terms = (
        "method",
        "methods",
        "protocol",
        "concentration",
        "dose",
        "dosing",
        "duration",
        "time course",
        "cell line",
        "animal model",
        "treatment",
    )
    facets = {
        "Question-Targeted Evidence": {"pmids": [], "direct": 0, "indirect": 0, "speculative": 0, "note": ""},
        "Cell / In Vitro Evidence": {"pmids": [], "direct": 0, "indirect": 0, "speculative": 0, "note": ""},
        "Animal / In Vivo Evidence": {"pmids": [], "direct": 0, "indirect": 0, "speculative": 0, "note": ""},
        "Methods / Operational Detail": {"pmids": [], "direct": 0, "indirect": 0, "speculative": 0, "note": ""},
        "Contradictory / Null Evidence": {"pmids": [], "direct": 0, "indirect": 0, "speculative": 0, "note": ""},
    }
    for article in articles[:24]:
        text = f"{article.get('title', '')} {article.get('abstract', '')}"
        has_disease = _contains_entity(text, disease_aliases) if disease_aliases else True
        has_cell = _contains_any(text, CELL_TERMS)
        has_animal = _contains_any(text, ANIMAL_TERMS)
        has_method = _contains_any(text, method_terms)
        is_contradictory = _contains_any(text, CONTRADICT_TERMS)
        pmid = str(article.get("pmid") or "").strip()
        if has_disease:
            level = "direct" if (has_cell or has_animal or has_method) else "indirect"
            facets["Question-Targeted Evidence"][level] += 1
            if pmid:
                facets["Question-Targeted Evidence"]["pmids"].append(pmid)
        if has_cell:
            level = "direct" if has_disease else "indirect"
            facets["Cell / In Vitro Evidence"][level] += 1
            if pmid:
                facets["Cell / In Vitro Evidence"]["pmids"].append(pmid)
        if has_animal:
            level = "direct" if has_disease else "indirect"
            facets["Animal / In Vivo Evidence"][level] += 1
            if pmid:
                facets["Animal / In Vivo Evidence"]["pmids"].append(pmid)
        if has_method:
            level = "direct" if (has_cell or has_animal or has_disease) else "indirect"
            facets["Methods / Operational Detail"][level] += 1
            if pmid:
                facets["Methods / Operational Detail"]["pmids"].append(pmid)
        if is_contradictory:
            facets["Contradictory / Null Evidence"]["direct"] += 1
            if pmid:
                facets["Contradictory / Null Evidence"]["pmids"].append(pmid)

    evidence_map = []
    for topic, bucket in facets.items():
        if bucket["direct"] > 0:
            support_level = "direct"
        elif bucket["indirect"] > 0:
            support_level = "indirect"
        else:
            support_level = "speculative"
        rationale = f"Direct={bucket['direct']}, indirect={bucket['indirect']}, speculative={bucket['speculative']} across retrieved literature."
        evidence_map.append(
            {
                "topic": topic,
                "support_level": support_level,
                "evidence_type": topic.lower().replace(" / ", " ").replace("->", "to"),
                "rationale": rationale,
                "representative_pmids": _dedupe_preserve(bucket["pmids"])[:6],
            }
        )
    return evidence_map


def audit_question_scope(
    research_question: str,
    prompt_constraints: str,
    literature_context: str,
    evidence_strength_map: List[dict],
    disease: str = "",
) -> dict:
    system_prompt = """You are auditing a standalone validation-planning request for scientific overreach.

Return JSON with exactly these keys:
- protocol_audit
- claims_to_verify
- high_risk_claims
- priority_questions
- experimental_jumps

protocol_audit must include:
- overall_assessment
- supported_claims
- inference_only_claims
- unsupported_or_overstated_claims
- priority_gaps

Requirements:
- Treat the user question as a draft planning objective, not ground truth.
- Use the actual literature content provided.
- Downgrade claims aggressively if direct cell or animal evidence is missing.
- If direct experimental validation for verbs such as produce, absorb, drive, or cause is missing, rewrite them in cautious language such as improve, decrease, is associated with, or may contribute to.
- Avoid absolute mechanistic wording unless the cited experiments directly support it.
- Prefer priority questions that can be answered by basic cell experiments, microbial culture, or animal studies.
- Do not recommend asking the human user for more inputs."""
    user_prompt = f"""Standalone validation question:
{_format_user_brief(research_question, prompt_constraints)}
{_user_question_plan_directive(research_question)}

Disease scope:
{disease or 'Not specified'}

Evidence strength map:
{_format_evidence_strength_map(evidence_strength_map)}

Literature excerpts:
{literature_context}

Audit the question scope and literature support. The weakest claims should be pushed to the front of the later validation plan.
Priority questions should be phrased as internal evidence questions for the system to resolve in later rounds."""
    try:
        raw = _chat_json(system_prompt, user_prompt, max_tokens=2200, temperature=0.15)
    except Exception:
        raw = {
            "protocol_audit": {
                "overall_assessment": "Automated question-scope audit was unavailable; deterministic evidence and completeness safeguards remain active.",
                "supported_claims": [],
                "inference_only_claims": [],
                "unsupported_or_overstated_claims": [],
                "priority_gaps": ["Recheck every claim against the retrieved primary literature before interpretation."],
            },
            "claims_to_verify": [],
            "high_risk_claims": [],
            "priority_questions": [],
            "experimental_jumps": [],
        }
    return {
        "protocol_audit": _normalize_protocol_audit(raw.get("protocol_audit")),
        "claims_to_verify": _as_str_list(raw.get("claims_to_verify")),
        "high_risk_claims": _as_str_list(raw.get("high_risk_claims")),
        "priority_questions": _as_str_list(raw.get("priority_questions")),
        "experimental_jumps": _as_str_list(raw.get("experimental_jumps")),
    }


def audit_protocol_claims(
    bacteria: str,
    metabolite: str,
    disease: str,
    protocol_text: str,
    candidate_context: str,
    literature_context: str,
    evidence_strength_map: List[dict],
    research_question: str,
    prompt_constraints: str,
) -> dict:
    system_prompt = """You are auditing a post-protocol validation draft for scientific overreach.

Return JSON with exactly these keys:
- protocol_audit
- claims_to_verify
- high_risk_claims
- priority_questions
- experimental_jumps

protocol_audit must include:
- overall_assessment
- supported_claims
- inference_only_claims
- unsupported_or_overstated_claims
- priority_gaps

Each supported_claims item must include: claim, support_level, pmids
Each inference_only_claims item must include: claim, reason, pmids
Each unsupported_or_overstated_claims item must include: claim, reason, pmids

Requirements:
- Treat the protocol as a draft, not ground truth.
- Use the actual literature content provided.
- Evaluate every PMID independently. Do not combine a culture method from one paper, an association from another, and a disease effect from a third into a single direct mechanistic claim.
- A direct microbe-to-metabolite production claim requires that the same cited paper reports controlled culture of the named microbe and direct measurement of the named metabolite as an output.
- If no such paper is present, list the narrower claim that each paper actually supports and state that direct production was not found in the current search.
- Downgrade claims aggressively if direct cell or animal evidence is missing.
- If direct experimental validation for verbs such as produce, absorb, drive, or cause is missing, rewrite them in cautious language such as improve, decrease, is associated with, or may contribute to.
- Avoid absolute mechanistic wording unless the cited experiments directly support it.
- Prefer questions that can be answered by cell experiments, microbial culture, or animal studies.
- Do not recommend asking the human user for more inputs."""
    if _has_user_research_question(research_question):
        system_prompt += (
            "\n- When a user research question is provided, prioritize whether the protocol and evidence gaps "
            "serve that question. Flag content that is generic validation but does not address the user's focus."
        )

    user_prompt = f"""Candidate:
Microbe: {bacteria}
Metabolite: {metabolite}
Disease: {disease}

Candidate context:
{candidate_context}

User-defined Step 3 brief:
{_format_user_brief(research_question, prompt_constraints)}
{_user_question_plan_directive(research_question)}

Evidence strength map:
{_format_evidence_strength_map(evidence_strength_map)}

Literature excerpts:
{literature_context}

Protocol under audit:
{protocol_text or 'No protocol text provided.'}

Audit the protocol. The weakest claims should be pushed to the front of the later validation plan.
Priority questions should be phrased as internal evidence questions for the system to resolve in later rounds."""
    if _has_user_research_question(research_question):
        user_prompt += (
            "\nWhen auditing, treat the user's research question as the main thread. "
            "Prioritize gaps that would block answering that question."
        )

    try:
        raw = _chat_json(system_prompt, user_prompt, max_tokens=2200, temperature=0.15)
    except Exception:
        raw = {
            "protocol_audit": {
                "overall_assessment": "Automated protocol audit was unavailable; deterministic paper grounding and experiment-completeness safeguards remain active.",
                "supported_claims": [],
                "inference_only_claims": [],
                "unsupported_or_overstated_claims": [],
                "priority_gaps": ["Do not treat the draft protocol as evidence; adjudicate candidate links paper by paper."],
            },
            "claims_to_verify": [],
            "high_risk_claims": [],
            "priority_questions": [],
            "experimental_jumps": [],
        }
    return {
        "protocol_audit": _normalize_protocol_audit(raw.get("protocol_audit")),
        "claims_to_verify": _as_str_list(raw.get("claims_to_verify")),
        "high_risk_claims": _as_str_list(raw.get("high_risk_claims")),
        "priority_questions": _as_str_list(raw.get("priority_questions")),
        "experimental_jumps": _as_str_list(raw.get("experimental_jumps")),
    }


def assess_direct_production_evidence_and_hypotheses(
    bacteria: str,
    metabolite: str,
    disease: str,
    literature: List[dict],
    literature_context: str,
    fulltext_method_evidence: Dict[str, List[dict]],
    evidence_strength_map: List[dict],
) -> dict:
    system_prompt = """You are the evidence adjudication and hypothesis-branching agent for a microbe-metabolite-disease validation workflow.

Return JSON with exactly these keys:
- direct_production_evidence_assessment
- hypothesis_branches

direct_production_evidence_assessment must include:
- status: direct_supported, indirect_only, not_found, or conflicting
- conclusion
- direct_evidence
- paper_findings
- evidence_limitations

Each direct_evidence item must include:
- pmid
- title
- claim
- evidence_type
- candidate_relevance
- claim_scope
- model_system
- measured_output
- why_this_support_level

Each paper_findings item must include the same fields and must describe exactly one paper identified by one PMID.

Allowed evidence_type values only:
- direct_monoculture_production
- direct_monoculture_nonproduction
- candidate_microbe_disease_intervention
- candidate_metabolite_disease_intervention
- candidate_pair_association
- indirect_ecological_evidence
- analogous_method_only
- background_only
- unrelated

candidate_relevance must be one of: direct, partial, analogous, unrelated.

hypothesis_branches must contain 2 to 4 competing hypotheses. Each item must include:
- hypothesis_id
- statement
- current_evidence_status
- evidence_basis
- discriminating_prediction
- in_vitro_gate
- animal_gate
- human_gate
- falsification_or_redirection

Evidence adjudication rules:
- Direct production means that a single cited study itself reports a controlled monoculture, pure-culture, axenic-culture, or otherwise defined culture experiment in which the named microbe was cultured and the named metabolite was directly measured as an output, preferably against an uninoculated control, baseline, time course, substrate control, or isotope-tracing design.
- A controlled candidate-specific culture study that directly measures the metabolite but reports no production is strong negative evidence. Classify it as direct_monoculture_nonproduction, preserve its PMID and negative result, and never use it to support direct production.
- A paper that only reports abundance-metabolite correlation, mixed-community fermentation, animal fecal concentrations, disease association, genomic potential, pathway annotation, or a host response is not direct production evidence.
- Never combine separate facts from different papers into one stronger claim. One paper may establish culture conditions, another an association, and another a disease effect; together they still do not become direct production evidence.
- direct_evidence items must each be defensible from that same PMID alone. If this threshold is not met, direct_evidence must be empty and status cannot be direct_supported.
- paper_findings must state only what each paper itself supports and must retain the PMID.
- Include a separate paper_findings item for every retrieved PMID that materially informs the candidate relationship; omit clearly irrelevant search hits rather than forcing a claim.
- Use candidate_relevance=direct only when the same paper actually studies the named candidate entity or pair. Studies of another bacterium, another disease, a probiotic cocktail, diet, polysaccharides, or mixed fermentation are analogous/background unless they directly contain the candidate relationship.
- candidate_microbe_disease_intervention requires direct intervention with the named microbe in the named disease context. candidate_metabolite_disease_intervention requires direct intervention with the named metabolite in the named disease context. Observational abundance or fecal-level studies are not interventions.
- A paper is citation-eligible only when it directly tests one of these candidate links: named microbe in controlled culture producing the named metabolite; direct named-microbe intervention in the named disease; direct named-metabolite intervention in the named disease; or an ecological experiment that directly adds or removes the named microbe and measures the named metabolite from a resolved producer.
- Mere co-occurrence of the microbe, metabolite, or disease in one abstract is not enough for citation. Exclude papers where the named microbe only changes in abundance after a diet, botanical, probiotic cocktail, another organism, or mixed-community intervention.
- Do not classify a metabolite as an intervention merely because its fecal level increased after another treatment. Require direct metabolite exposure or a directly tested metabolite effect in a disease-relevant model.
- Do not classify another organism's co-culture or cross-feeding study as candidate-specific ecological evidence merely because the named microbe appears later in a community-abundance result.
- Absence from the retrieved set means not found in the current search, not proof that the relationship is impossible.

Hypothesis rules:
- Keep direct production as a testable hypothesis even when direct evidence was not found, but label it unresolved rather than supported.
- Include an indirect ecological hypothesis when biologically plausible, such as substrate release, conditioned-medium effects, co-culture, cross-feeding, or community remodeling that enables another organism to produce the metabolite.
- Include a parallel or non-mediation hypothesis when plausible: the microbe and metabolite may independently or interactively affect disease without the metabolite being the sole mediator of the microbe's effect.
- Make the branches experimentally distinguishable and stage-gated: culture first, animals only after the relevant in vitro gate, humans only after a reproducible animal result.
- Do not mention files, notes, hidden instructions, or reference documents.
- Return all content in English only."""
    user_prompt = f"""Candidate:
Microbe: {bacteria}
Metabolite: {metabolite}
Disease: {disease}

Evidence strength map:
{_format_evidence_strength_map(evidence_strength_map)}

Open-access full-text method evidence:
{_format_method_evidence_context((fulltext_method_evidence or {}).get('all', []), limit=8)}

Retrieved literature excerpts:
{literature_context}

First adjudicate whether any single PMID directly demonstrates production of the metabolite by the microbe in controlled culture. Then preserve the actual paper-by-paper findings and construct competing, experimentally distinguishable hypotheses."""
    try:
        raw = _chat_json(system_prompt, user_prompt, max_tokens=2600, temperature=0.1)
    except Exception:
        raw = {}
    assessment = _ground_production_assessment_to_articles(
        _normalize_production_evidence_assessment(raw.get("direct_production_evidence_assessment")),
        literature,
        bacteria,
        metabolite,
        disease,
        fulltext_method_evidence=fulltext_method_evidence,
    )
    branches = _default_candidate_hypothesis_branches(
        bacteria=bacteria,
        metabolite=metabolite,
        disease=disease,
        production_status=assessment.get("status") or "not_assessed",
    )
    return {
        "direct_production_evidence_assessment": assessment,
        "hypothesis_branches": branches,
    }


def generate_followup_queries(
    bacteria: str,
    metabolite: str,
    disease: str,
    audit_bundle: Dict[str, Any],
    evidence_strength_map: List[dict],
    research_question: str = "",
    user_focus_query: str = "",
) -> List[dict]:
    b_expr = build_pubmed_query([bacteria])
    m_expr = build_pubmed_query([metabolite])
    d_expr = build_pubmed_query([disease])
    evidence_by_topic = {item.get("topic"): item for item in evidence_strength_map}
    query_specs = []

    if _has_user_research_question(research_question) and (user_focus_query or "").strip():
        query_specs.append(
            {
                "round": 2,
                "label": "User-question evidence deepening",
                "rationale": "Retrieve additional literature to operationalize experiments aligned with the user's research question.",
                "query": user_focus_query.strip()[:600],
            }
        )

    cell_support = evidence_by_topic.get("Cell / In Vitro Evidence", {})
    animal_support = evidence_by_topic.get("Animal / In Vivo Evidence", {})
    contradiction_support = evidence_by_topic.get("Contradictory / Null Evidence", {})
    direct_production_support = evidence_by_topic.get("Direct Monoculture Production Evidence", {})
    high_risk_text = " ".join(audit_bundle.get("high_risk_claims") or []).lower()
    gap_text = " ".join(audit_bundle.get("protocol_audit", {}).get("priority_gaps") or []).lower()

    if direct_production_support.get("support_level") != "direct":
        query_specs.append(
            {
                "round": 2,
                "label": "Direct production evidence deepening",
                "rationale": "Search specifically for a single controlled monoculture or defined-culture study that directly measured metabolite production by the named microbe.",
                "query": f"{b_expr} AND {m_expr} AND (monoculture OR \"pure culture\" OR \"defined medium\" OR \"culture supernatant\") AND (quantified OR \"mass spectrometry\" OR metabolomics OR production OR produced)",
            }
        )

    if cell_support.get("support_level") != "direct" or "cell" in high_risk_text or "culture" in gap_text:
        query_specs.append(
            {
                "round": 2,
                "label": "Targeted cell evidence check",
                "rationale": "Search for organoid, epithelial, macrophage, or co-culture evidence before escalating to animals.",
                "query": f"{b_expr} AND {m_expr} AND {d_expr} AND (organoid OR epithelial OR macrophage OR co-culture OR barrier OR inflammation)",
            }
        )
    if animal_support.get("support_level") == "speculative" or "animal" in high_risk_text or "in vivo" in gap_text:
        query_specs.append(
            {
                "round": 2,
                "label": "Targeted animal validation",
                "rationale": "Find whether there is any real animal model support and what model/readouts are used.",
                "query": f"({b_expr} OR {m_expr}) AND {d_expr} AND (mouse OR mice OR murine OR animal OR colitis OR gnotobiotic)",
            }
        )
    if contradiction_support.get("support_level") == "speculative" or "contradict" in gap_text or "overstate" in high_risk_text:
        query_specs.append(
            {
                "round": 2,
                "label": "Contradictory or null evidence",
                "rationale": "Actively search for negative or non-replicating findings for later self-reflection and claim revision.",
                "query": f"({b_expr} OR {m_expr}) AND {d_expr} AND (\"no significant\" OR contradictory OR negative OR null OR inconsistent)",
            }
        )
    if not query_specs:
        query_specs.append(
            {
                "round": 2,
                "label": "Mechanistic pathway follow-up",
                "rationale": "Use a fallback second-round mechanistic query when no other follow-up trigger is strong enough.",
                "query": f"{m_expr} AND {d_expr} AND (mechanism OR receptor OR signaling OR pathway OR inflammation)",
            }
        )
    elif len(query_specs) < MAX_QUERY_ROUNDS:
        query_specs.append(
            {
                "round": 3,
                "label": "Mechanistic pathway follow-up",
                "rationale": "Backfill mechanistic experiments only after second-round culture/cell/animal gaps are checked.",
                "query": f"{m_expr} AND {d_expr} AND (mechanism OR receptor OR signaling OR pathway OR inflammation)",
            }
        )
    return query_specs[:MAX_QUERY_ROUNDS]


def generate_question_followup_queries(
    research_question: str,
    prompt_constraints: str,
    seed_query: str,
    evidence_strength_map: List[dict],
) -> List[dict]:
    del prompt_constraints
    evidence_by_topic = {item.get("topic"): item for item in evidence_strength_map}
    query_specs = [
        {
            "round": 2,
            "label": "Methods and reported conditions",
            "rationale": "Find papers that report experimental conditions, models, doses, concentrations, and timing relevant to the question.",
            "query": f"({seed_query}) AND (methods OR protocol OR concentration OR dose OR dosing OR time OR duration OR cell line OR animal model)",
        }
    ]
    wants_cell = _contains_any(research_question, ("cell", "cells", "in vitro", "organoid", "culture", "epithelial", "macrophage"))
    wants_animal = _contains_any(research_question, ("animal", "in vivo", "mouse", "mice", "murine", "rat", "colitis"))
    cell_support = evidence_by_topic.get("Cell / In Vitro Evidence", {})
    animal_support = evidence_by_topic.get("Animal / In Vivo Evidence", {})
    contradiction_support = evidence_by_topic.get("Contradictory / Null Evidence", {})

    if wants_cell or cell_support.get("support_level") != "direct":
        query_specs.append(
            {
                "round": 2,
                "label": "Targeted cell evidence check",
                "rationale": "Deepen the cell or in vitro evidence base for the standalone question.",
                "query": f"({seed_query}) AND (\"in vitro\" OR cell OR epithelial OR macrophage OR organoid OR co-culture)",
            }
        )
    if wants_animal or animal_support.get("support_level") != "direct":
        query_specs.append(
            {
                "round": 2,
                "label": "Targeted animal validation",
                "rationale": "Find whether there is real animal-model support and what models, routes, and timelines are used.",
                "query": f"({seed_query}) AND (mouse OR mice OR murine OR animal OR colitis OR gnotobiotic)",
            }
        )
    if contradiction_support.get("support_level") != "direct":
        query_specs.append(
            {
                "round": 2,
                "label": "Contradictory or null evidence",
                "rationale": "Actively search for negative, null, or inconsistent findings for later self-reflection.",
                "query": f"({seed_query}) AND (\"no significant\" OR contradictory OR negative OR null OR inconsistent)",
            }
        )
    return query_specs[:MAX_QUERY_ROUNDS]


def expand_in_vitro_plan(
    bacteria: str,
    metabolite: str,
    disease: str,
    protocol_text: str,
    audit_bundle: Dict[str, Any],
    literature_context: str,
    evidence_strength_map: List[dict],
    query_log: List[dict],
    fulltext_method_evidence: List[dict],
    protocols_io_evidence: List[dict],
    research_question: str,
    prompt_constraints: str,
    direct_production_evidence_assessment: Optional[dict] = None,
    hypothesis_branches: Optional[List[dict]] = None,
    question_profile: Optional[dict] = None,
) -> dict:
    system_prompt = f"""You are designing the in vitro / cell experiment package for a validation plan.

Return JSON with exactly one key:
- in_vitro_plan

Each in_vitro_plan item must include:
- module_id
- experiment_role
- hypothesis_ids
- route_status
- execution_status
- result_status
- design_status
- activation_gate
- scientific_question
- why_needed
- study_object
- groups
- samples_and_timing
- primary_indicator
- secondary_indicators
- key_controls
- positive_gate
- branch_if_positive
- branch_if_negative
- unlock_rule
- claim_boundary

Requirements:
- Return all content in English only.
- Do not output any Chinese characters.
- Return the smallest sufficient set of modules. Do not force unrelated roles or use a numerical quota, but retain V2 as an explicit conditional branch when indirect ecology remains a live competing explanation.
- For groups, return objects with group_name, exposure_or_condition, and control_purpose. Never return a separate group count.
- Use result_status=not_run for a prospective plan. Do not infer that an upstream experiment is positive from its execution status.
- Put direct source testing first when unresolved, keep host-cell activity available in parallel, and make indirect ecology conditional on a negative or inconclusive direct-source result unless ecological amplification is an explicit objective.
- Describe only the study object, groups, samples and timing, indicators, controls, and branch logic needed to answer the scientific question.
- Do not return equipment, numbered procedures, detailed concentrations, full analysis models, or parameter-provenance records.
{UNIVERSAL_EXPERIMENT_DESIGN_PROMPT_RULES}
{MICROBE_METABOLITE_IN_VITRO_PROMPT_RULES}
- If the question is not a microbe-metabolite-disease workflow, adapt the same concise decision-first structure to the actual question."""

    if bacteria and metabolite:
        system_prompt += """

Candidate hypothesis-coverage requirement:
- The supplied H1, H2, and H3 branches are active competing explanations, not optional prose. Return one complete V1 card with hypothesis_ids=["H1"], one complete conditional V2 card with hypothesis_ids=["H2"], and one complete parallel V3 card with hypothesis_ids=["H3"].
- A conditional card must still contain the study object, explicit groups and control purposes, samples and timing, one primary indicator, secondary indicators, key controls, a directional positive gate, both result branches, and a claim boundary.
- Do not merge H1/H2/H3 into one omnibus experiment, do not duplicate a canonical role, and do not create decimal or M-series module aliases."""

    user_prompt = f"""Candidate:
Microbe: {bacteria}
Metabolite: {metabolite}
Disease: {disease}

Protocol under review:
{protocol_text}

User-defined Step 3 brief:
{_format_user_brief(research_question, prompt_constraints)}
{_user_question_plan_directive(research_question)}

Question profile:
{_format_question_profile(question_profile or {})}

Direct production evidence assessment:
{_format_production_evidence_assessment(direct_production_evidence_assessment or {})}

Competing hypothesis branches:
{_format_hypothesis_branches(hypothesis_branches or [])}

Protocol audit:
{audit_bundle.get('protocol_audit')}

Evidence strength map:
{_format_evidence_strength_map(evidence_strength_map)}

Query log:
{_format_query_log(query_log)}

Open-access full-text methods evidence:
{_format_method_evidence_context(fulltext_method_evidence)}

protocols.io operational references:
{_format_protocol_evidence_context(protocols_io_evidence)}

Literature excerpts:
{literature_context}

Return the smallest sufficient set of in vitro / cell experiment blocks, prioritizing:
{_vitro_priority_block(research_question)}

The blocks must form an explicit decision tree. State which hypothesis each block tests, what result unlocks the next block, and where a negative result redirects the program.

If the user question can be answered without a generic multi-step package, return the smallest literature-supported in vitro set that directly answers it.

Evidence-use requirement:
- Use candidate-relevant methods evidence only to choose a defensible study object, sample type, timing level, or indicator. Do not copy operational details from unrelated candidates and do not expand the module into an SOP."""

    if not (bacteria or metabolite):
        user_prompt += """

Question-driven safeguards:
- Do not allocate most of the answer to human studies by default.
- If the question involves toxicants, hazardous exposures, or other ethically problematic interventions, do not propose direct human exposure trials; keep the plan preclinical or observational.
- Only propose human trials when the intervention is plausibly ethical and feasible, such as beneficial interventions, repurposed drugs, low-risk behavioral designs, or clearly observational human work.
- If the question is mainly about finding a sensor protein or molecular target, spend most of the detail on discovery and functional validation rather than downstream human translation."""

    try:
        raw = _chat_json(system_prompt, user_prompt, max_tokens=3200, temperature=0.2)
    except Exception:
        raw = {}
    items = raw.get("in_vitro_plan")
    items = _backfill_method_fields(items if isinstance(items, list) else [], fulltext_method_evidence)
    items = _backfill_protocol_fields(items, protocols_io_evidence)
    if bacteria and metabolite:
        items = _assign_candidate_experiment_roles(items, is_in_vivo=False)
    return {"in_vitro_plan": _prepare_experiment_modules(items)}


def expand_in_vivo_plan(
    bacteria: str,
    metabolite: str,
    disease: str,
    protocol_text: str,
    audit_bundle: Dict[str, Any],
    literature_context: str,
    evidence_strength_map: List[dict],
    query_log: List[dict],
    fulltext_method_evidence: List[dict],
    protocols_io_evidence: List[dict],
    research_question: str,
    prompt_constraints: str,
    direct_production_evidence_assessment: Optional[dict] = None,
    hypothesis_branches: Optional[List[dict]] = None,
    question_profile: Optional[dict] = None,
) -> dict:
    system_prompt = f"""You are designing the in vivo / animal experiment package for a validation plan.

Return JSON with exactly one key:
- in_vivo_plan

Each in_vivo_plan item must include:
- module_id
- experiment_role
- hypothesis_ids
- route_status
- execution_status
- result_status
- design_status
- activation_gate
- scientific_question
- why_needed
- study_object
- groups
- samples_and_timing
- primary_indicator
- secondary_indicators
- key_controls
- positive_gate
- branch_if_positive
- branch_if_negative
- unlock_rule
- claim_boundary

Requirements:
- Return all content in English only.
- Do not output any Chinese characters.
- Return the smallest sufficient animal set. For a microbe-metabolite-disease workflow, return A1 plus brief conditional A2 and/or A3 blueprints when the corresponding V1 and/or V2 source branch is present. Their inclusion records the decision route and does not make them executable before their code-owned gates pass.
- For groups, return objects with group_name, exposure_or_condition, and control_purpose. Never return a separate group count.
- Use result_status=not_run for a prospective plan. execution_status must not stand in for an experimental result.
- Make the first animal effect module a healthy reference plus the disease-condition microbe-by-metabolite 2 x 2 core when that design answers the question; keep branch-specific mediation modules conditional.
- State the study object, explicit groups, baseline and endpoint samples, exposure indicators, one primary disease indicator, a short secondary-indicator list, and positive/negative branches.
- Do not return equipment, numbered procedures, detailed dose schedules, full analysis models, or parameter-provenance records.
{UNIVERSAL_EXPERIMENT_DESIGN_PROMPT_RULES}
{MICROBE_METABOLITE_IN_VIVO_PROMPT_RULES}
- Treat a positive effect/interaction study as efficacy or interaction evidence only; mediation requires the appropriate conditional causal module."""
    if bacteria and metabolite:
        system_prompt += "\n" + MICROBE_METABOLITE_HUMAN_PROMPT_RULES + """

Candidate hypothesis-coverage requirement:
- Return one complete A1 effect-and-interaction card with hypothesis_ids=["H3"].
- When V1 and A1 are present, retain a complete conditional A2 decision card with hypothesis_ids=["H1"]; when V2 and A1 are present, retain a complete conditional A3 decision card with hypothesis_ids=["H2"].
- Conditional A2/A3 cards require explicit necessity and rescue comparisons, samples and timing, one primary indicator, exposure/process indicators, key controls, a positive gate, both result branches, and a claim boundary. Do not invent candidate-specific strains, functions, producers, doses, or disease models."""

    user_prompt = f"""Candidate:
Microbe: {bacteria}
Metabolite: {metabolite}
Disease: {disease}

Protocol under review:
{protocol_text}

User-defined Step 3 brief:
{_format_user_brief(research_question, prompt_constraints)}
{_user_question_plan_directive(research_question)}

Question profile:
{_format_question_profile(question_profile or {})}

Direct production evidence assessment:
{_format_production_evidence_assessment(direct_production_evidence_assessment or {})}

Competing hypothesis branches:
{_format_hypothesis_branches(hypothesis_branches or [])}

Protocol audit:
{audit_bundle.get('protocol_audit')}

Evidence strength map:
{_format_evidence_strength_map(evidence_strength_map)}

Query log:
{_format_query_log(query_log)}

Open-access full-text methods evidence:
{_format_method_evidence_context(fulltext_method_evidence)}

protocols.io operational references:
{_format_protocol_evidence_context(protocols_io_evidence)}

Literature excerpts:
{literature_context}

Return the smallest sufficient set of animal experiment blocks. The first returned block should act as a gatekeeper:
{_vivo_priority_block(research_question)}

The animal blocks must state the prerequisite in vitro result and the branch-specific interpretation of positive and negative outcomes. Human work must remain downstream of a reproducible, interpretable animal result.

If the user question does not truly require an animal step at this stage, return no animal block rather than forcing a generic one.

Evidence-use requirement:
- Use candidate-relevant methods evidence only to choose a defensible animal object, sample type, timing level, or indicator. Do not copy operational details from unrelated candidates and do not expand the module into an SOP."""

    if not (bacteria or metabolite):
        user_prompt += """

Question-driven safeguards:
- Do not allocate most of the answer to human studies by default.
- If the question involves toxicants, hazardous exposures, or other ethically problematic interventions, do not propose direct human exposure trials; keep the plan preclinical or observational.
- Only propose human trials when the intervention is plausibly ethical and feasible, such as beneficial interventions, repurposed drugs, low-risk behavioral designs, or clearly observational human work.
- If the question is mainly about finding a sensor protein or molecular target, spend most of the detail on in vivo necessity/sufficiency logic only after the target-identification stage is well defined."""

    try:
        raw = _chat_json(system_prompt, user_prompt, max_tokens=3200, temperature=0.2)
    except Exception:
        raw = {}
    items = raw.get("in_vivo_plan")
    items = _backfill_method_fields(items if isinstance(items, list) else [], fulltext_method_evidence)
    items = _backfill_protocol_fields(items, protocols_io_evidence)
    if bacteria and metabolite:
        items = _assign_candidate_experiment_roles(items, is_in_vivo=True)
    return {"in_vivo_plan": _prepare_experiment_modules(items, is_in_vivo=True)}


def run_self_reflection(
    bacteria: str,
    metabolite: str,
    disease: str,
    protocol_text: str,
    audit_bundle: Dict[str, Any],
    literature_context: str,
    evidence_strength_map: List[dict],
    query_log: List[dict],
    protocols_io_evidence: Dict[str, Any],
    in_vitro_plan: List[dict],
    in_vivo_plan: List[dict],
    research_question: str,
    prompt_constraints: str,
    direct_production_evidence_assessment: Optional[dict] = None,
    hypothesis_branches: Optional[List[dict]] = None,
) -> dict:
    system_prompt = """You are the self-reflection and critical revision agent for a post-protocol validation plan.

Return JSON with exactly these keys:
- working_hypothesis
- evidence_basis
- self_reflection
- revised_in_vitro_plan
- revised_in_vivo_plan
- overall_risk_flags

working_hypothesis must include:
- statement
- direct_support
- inference_only

Each evidence_basis item must include:
- claim_type
- claim
- support_level
- pmids
- model_type
- evidence_summary

Each self_reflection item must include:
- initial_claim
- self_critique
- evidence_checked
- revision
- remaining_uncertainty

Requirements:
- Return all content in English only.
- Do not output any Chinese characters.
- This is self-reflection, not a request to invent falsification experiments or falsification conditions.
- If microbe and metabolite are not specified, treat this as a standalone question-driven planning request and critique whether the draft actually answers the user brief.
- Critique the generated plan's own claims, citations, method transfers, and confidence.
- Check whether correlation was overstated as causation.
- Recheck the direct-production decision paper by paper. Never upgrade indirect evidence to direct production by combining multiple PMIDs.
- If the assessment says direct production was not found, preserve that conclusion unless one specific cited study independently meets the controlled-culture and direct-measurement threshold.
- Check whether direct-production, indirect-production, and parallel-effect hypotheses remain separate and experimentally distinguishable throughout the plan.
- Do not preserve a fixed numerical quota of experiment roles. For a microbe-metabolite-disease workflow, retain brief A2/A3 conditional decision cards whenever A1 and the corresponding V1/V2 source branch are present; do not expand them into executable protocols before their gates pass.
- A downstream module with an unmet prerequisite must remain in branch logic or be marked conditional_future; it must not be presented as ready_now.
- Audit every returned module for one clear scientific question, a justified study object, explicit structured groups, samples and timing, one primary indicator, key controls, a positive gate, branch logic, and an explicit claim boundary.
- Keep the revision at research-strategy level. Remove equipment lists, numbered procedures, exhaustive dose schedules, full statistical models, and line-item parameter-provenance records.
- Confirm that concentration change is not called direct production, efficacy is not called mediation, and human observation is not used to repair incomplete preclinical causality.
- Remove citations based only on co-occurrence or keyword overlap. Retain only candidate-specific controlled culture, direct named intervention, or directly manipulated ecological evidence.
- Check that every animal experiment names the in vitro prerequisite that unlocks it, and that every human experiment is conditional on an interpretable animal result.
- Check that a metabolite effect, a microbial effect, or a combination effect is not mislabeled as mediation without function-loss and rescue or an equivalent causal design.
- Rewrite unsupported absolute verbs such as produce, absorb, drive, or cause into cautious language such as improve, decrease, is associated with, or may contribute to unless the cited experiments directly justify the stronger wording.
- If an exact concentration, model, route, duration, or material is retained, check that it is directly relevant and supported; otherwise replace it with a high-level optimization note.
- Do not treat protocols.io operational references as candidate-specific scientific evidence.
- Check whether the plan complied with the user brief by overstating evidence or inventing details.
- Use contradictory or null evidence to revise or downgrade claims.
- If direct evidence is absent, say so explicitly and revise the claim.
- Resolve contradictions across sections so that supported claims, missing gaps, and follow-up experiments do not conflict with one another.
- Record what changed after reflection and what remains uncertain."""
    if bacteria and metabolite:
        system_prompt += "\n" + MICROBE_METABOLITE_HUMAN_PROMPT_RULES
    if _has_user_research_question(research_question):
        system_prompt += (
            "\n- When a user research question is provided, check whether the plan addresses it as the main thread. "
            "Revise generic validation content that drifts away from the user's question."
        )

    user_prompt = f"""Candidate:
Microbe: {bacteria}
Metabolite: {metabolite}
Disease: {disease}

Protocol under review:
{_compact_protocol_text(protocol_text, 1200)}

User-defined Step 3 brief:
{_format_user_brief(research_question, prompt_constraints)}
{_user_question_plan_directive(research_question)}

Direct production evidence assessment:
{_format_production_evidence_assessment(direct_production_evidence_assessment or {})}

Competing hypothesis branches:
{_format_hypothesis_branches(hypothesis_branches or [])}

Protocol audit:
{audit_bundle.get('protocol_audit')}

High-risk claims:
{audit_bundle.get('high_risk_claims')}

Priority questions:
{audit_bundle.get('priority_questions')}

Evidence strength map:
{_format_evidence_strength_map(evidence_strength_map)}

Query log:
{_format_query_log(query_log)}

Draft in vitro plan:
{_compact_experiment_bundle(in_vitro_plan, is_in_vivo=False, limit=3)}

Draft in vivo plan:
{_compact_experiment_bundle(in_vivo_plan, is_in_vivo=True, limit=3)}

protocols.io retrieval summary:
{_compact_protocols_io_summary(protocols_io_evidence, limit=2)}

Literature excerpts:
{literature_context}

Reflect on the draft and revise overconfident conclusions without expanding its technical detail.
Return complete revised experiment arrays using the same compact fields supplied in the draft. Preserve module IDs, roles, hypothesis IDs, route/result/design status, activation gates, structured groups, samples, indicators, branch rules, and claim boundaries."""

    try:
        raw = _chat_json(system_prompt, user_prompt, max_tokens=3500, temperature=0.15)
    except Exception:
        raw = {
            "working_hypothesis": {},
            "evidence_basis": [],
            "self_reflection": [],
            "revised_in_vitro_plan": in_vitro_plan,
            "revised_in_vivo_plan": in_vivo_plan,
            "overall_risk_flags": ["Automated reflection was unavailable; deterministic completeness and evidence filters were retained."],
        }
    return {
        "working_hypothesis": raw.get("working_hypothesis") or {},
        "evidence_basis": _normalize_evidence_items(raw.get("evidence_basis")),
        "self_reflection": _normalize_reflection_items(raw.get("self_reflection")),
        "revised_in_vitro_plan": _prepare_experiment_modules(
            _preserve_experiment_source_fields(raw.get("revised_in_vitro_plan"), in_vitro_plan)
        ),
        "revised_in_vivo_plan": _prepare_experiment_modules(
            _preserve_experiment_source_fields(raw.get("revised_in_vivo_plan"), in_vivo_plan),
            is_in_vivo=True,
        ),
        "overall_risk_flags": _as_str_list(raw.get("overall_risk_flags")),
    }


def assemble_validation_plan(
    bacteria: str,
    metabolite: str,
    disease: str,
    candidate_metrics: Dict[str, Any],
    literature: List[dict],
    mode: str,
    iterative_query_log: List[dict],
    evidence_strength_map: List[dict],
    fulltext_method_evidence: Dict[str, List[dict]],
    protocols_io_evidence: Dict[str, Any],
    audit_bundle: Dict[str, Any],
    vitro_bundle: Dict[str, Any],
    vivo_bundle: Dict[str, Any],
    reflection_bundle: Dict[str, Any],
    hypothesis_bundle: Dict[str, Any],
    research_question: str,
    prompt_constraints: str,
) -> dict:
    raw = {
        "protocol_audit": audit_bundle.get("protocol_audit"),
        "working_hypothesis": reflection_bundle.get("working_hypothesis"),
        "direct_production_evidence_assessment": hypothesis_bundle.get("direct_production_evidence_assessment"),
        "hypothesis_branches": hypothesis_bundle.get("hypothesis_branches"),
        "evidence_basis": reflection_bundle.get("evidence_basis"),
        "in_vitro_plan": reflection_bundle.get("revised_in_vitro_plan") or vitro_bundle.get("in_vitro_plan"),
        "in_vivo_plan": reflection_bundle.get("revised_in_vivo_plan") or vivo_bundle.get("in_vivo_plan"),
        "self_reflection": reflection_bundle.get("self_reflection"),
        "overall_risk_flags": reflection_bundle.get("overall_risk_flags"),
    }
    plan = _normalize_plan(
        raw=raw,
        bacteria=bacteria,
        metabolite=metabolite,
        disease=disease,
        candidate_metrics=candidate_metrics,
        literature=literature,
        mode=mode,
        iterative_query_log=iterative_query_log,
        evidence_strength_map=evidence_strength_map,
        fulltext_method_evidence=fulltext_method_evidence,
        protocols_io_evidence=protocols_io_evidence,
        research_question=research_question,
        prompt_constraints=prompt_constraints,
    )
    if bacteria and metabolite:
        assessment = _normalize_production_evidence_assessment(
            hypothesis_bundle.get("direct_production_evidence_assessment")
        )
        plan["direct_production_evidence_assessment"] = assessment
        plan["hypothesis_branches"] = _default_candidate_hypothesis_branches(
            bacteria=bacteria,
            metabolite=metabolite,
            disease=disease,
            production_status=assessment.get("status") or "not_assessed",
        )
        plan["in_vitro_plan"] = _assign_candidate_experiment_roles(
            plan.get("in_vitro_plan") or [],
            is_in_vivo=False,
        )
        plan["in_vivo_plan"] = _assign_candidate_experiment_roles(
            plan.get("in_vivo_plan") or [],
            is_in_vivo=True,
        )
        plan["in_vitro_plan"] = _prepare_experiment_modules(plan.get("in_vitro_plan") or [])
        plan["in_vivo_plan"] = _prepare_experiment_modules(plan.get("in_vivo_plan") or [], is_in_vivo=True)
        plan = _backfill_conditional_causal_followups(plan)
        plan = _sanitize_candidate_plan_references(plan)
        plan = _attach_strong_evidence_to_plan(plan)
        plan["human_plan"] = _prepare_experiment_modules(
            plan.get("human_plan") or _build_question_driven_human_experiments(plan),
            is_in_vivo=True,
            is_human=True,
        )
        plan = _refresh_plan_level_structure(plan, candidate_workflow=True)
        plan["validation_summary"] = _build_validation_summary(plan)
        plan["validation_protocol_text"] = _render_validation_protocol_text(plan)
    return plan


def generate_validation_plan(
    bacteria: str,
    metabolite: str,
    disease: str = "IBD",
    run_id: Optional[str] = None,
    mechanism_summary: str = "",
    mode: str = "evidence_self_reflection",
    protocol_text: str = "",
    research_question: str = "",
    prompt_constraints: str = "",
    progress_callback: Optional[Callable[[str, str, int, Optional[Dict[str, Any]]], None]] = None,
) -> dict:
    if mode == "question_driven" and not (bacteria and metabolite):
        return generate_question_driven_validation_plan(
            research_question=research_question,
            prompt_constraints=prompt_constraints,
            disease=disease,
            progress_callback=progress_callback,
        )

    _report_progress(progress_callback, "prepare_support", "Loading candidate context and seed literature.", 5)
    support = _prepare_protocol_support(
        bacteria=bacteria,
        metabolite=metabolite,
        disease=disease,
        run_id=run_id,
        mechanism_summary=mechanism_summary,
    )
    candidate_context = support["candidate_context"]
    literature = support["literature"]
    candidate_metrics = support["candidate_metrics"]

    _report_progress(progress_callback, "round1_search", "Running first-round PubMed queries.", 22)
    round_one_specs = _build_round_one_query_specs(bacteria, metabolite, disease)
    user_focus_spec = build_user_focus_query_spec(
        bacteria=bacteria,
        metabolite=metabolite,
        disease=disease,
        research_question=research_question,
        prompt_constraints=prompt_constraints,
    )
    if user_focus_spec:
        round_one_specs.insert(0, user_focus_spec)
    round_one = search_literature_online(round_one_specs)
    initial_literature = _merge_articles(literature, round_one["articles"])
    initial_strength_map = rank_evidence_strength(initial_literature, bacteria, metabolite, disease)
    _report_progress(
        progress_callback,
        "audit_protocol",
        "Auditing protocol claims against first-round evidence.",
        35,
        {
            "initial_articles": len(initial_literature),
            "round1_queries": len(round_one["query_log"]),
        },
    )

    audit_bundle = audit_protocol_claims(
        bacteria=bacteria,
        metabolite=metabolite,
        disease=disease,
        protocol_text=protocol_text,
        candidate_context=candidate_context,
        literature_context=_format_article_context(initial_literature, limit=6),
        evidence_strength_map=initial_strength_map,
        research_question=research_question,
        prompt_constraints=prompt_constraints,
    )

    _report_progress(progress_callback, "followup_queries", "Generating follow-up literature questions for weak claims.", 48)
    followup_specs = generate_followup_queries(
        bacteria=bacteria,
        metabolite=metabolite,
        disease=disease,
        audit_bundle=audit_bundle,
        evidence_strength_map=initial_strength_map,
        research_question=research_question,
        user_focus_query=(user_focus_spec or {}).get("query", ""),
    )
    _report_progress(progress_callback, "round2_search", "Running follow-up PubMed queries and merging evidence.", 58)
    followup_round = search_literature_online(followup_specs)
    merged_literature = _merge_articles(initial_literature, followup_round["articles"])
    query_log = round_one["query_log"] + followup_round["query_log"]
    evidence_strength_map = rank_evidence_strength(merged_literature, bacteria, metabolite, disease)
    literature_context = _format_article_context(merged_literature, limit=8)
    _report_progress(progress_callback, "fulltext_methods", "Retrieving open-access full-text methods for concentrations and model details.", 66)
    fulltext_method_evidence = collect_fulltext_method_evidence(merged_literature, bacteria, metabolite, disease)
    _report_progress(
        progress_callback,
        "protocols_io_search",
        "Searching protocols.io for operational procedures and materials.",
        69,
    )
    protocols_io_evidence = collect_protocol_evidence(
        bacteria=bacteria,
        metabolite=metabolite,
        disease=disease,
        research_question=research_question,
    )
    _report_progress(
        progress_callback,
        "protocols_io_search",
        protocols_io_evidence.get("message") or "protocols.io search finished.",
        71,
        {
            "status": protocols_io_evidence.get("status"),
            "protocol_count": len(protocols_io_evidence.get("all", [])),
        },
    )

    _report_progress(
        progress_callback,
        "evidence_adjudication",
        "Adjudicating direct culture evidence paper by paper and building competing hypotheses.",
        72,
    )
    hypothesis_bundle = assess_direct_production_evidence_and_hypotheses(
        bacteria=bacteria,
        metabolite=metabolite,
        disease=disease,
        literature=merged_literature,
        literature_context=_format_article_context(merged_literature, limit=12),
        fulltext_method_evidence=fulltext_method_evidence,
        evidence_strength_map=evidence_strength_map,
    )

    _report_progress(
        progress_callback,
        "in_vitro_design",
        "Expanding in vitro and cell experiment package.",
        73,
        {
            "total_articles": len(merged_literature),
            "total_queries": len(query_log),
            "fulltext_method_hits": len(fulltext_method_evidence.get("all", [])),
            "protocols_io_hits": len(protocols_io_evidence.get("all", [])),
            "protocols_io_status": protocols_io_evidence.get("status"),
        },
    )
    vitro_bundle = expand_in_vitro_plan(
        bacteria=bacteria,
        metabolite=metabolite,
        disease=disease,
        protocol_text=protocol_text,
        audit_bundle=audit_bundle,
        literature_context=literature_context,
        evidence_strength_map=evidence_strength_map,
        query_log=query_log,
        fulltext_method_evidence=fulltext_method_evidence.get("in_vitro", []),
        protocols_io_evidence=protocols_io_evidence.get("in_vitro", []),
        research_question=research_question,
        prompt_constraints=prompt_constraints,
        direct_production_evidence_assessment=hypothesis_bundle.get("direct_production_evidence_assessment"),
        hypothesis_branches=hypothesis_bundle.get("hypothesis_branches"),
    )
    _report_progress(progress_callback, "in_vivo_design", "Expanding animal experiment package and escalation gates.", 83)
    vivo_bundle = expand_in_vivo_plan(
        bacteria=bacteria,
        metabolite=metabolite,
        disease=disease,
        protocol_text=protocol_text,
        audit_bundle=audit_bundle,
        literature_context=literature_context,
        evidence_strength_map=evidence_strength_map,
        query_log=query_log,
        fulltext_method_evidence=fulltext_method_evidence.get("in_vivo", []),
        protocols_io_evidence=protocols_io_evidence.get("in_vivo", []),
        research_question=research_question,
        prompt_constraints=prompt_constraints,
        direct_production_evidence_assessment=hypothesis_bundle.get("direct_production_evidence_assessment"),
        hypothesis_branches=hypothesis_bundle.get("hypothesis_branches"),
    )
    _report_progress(progress_callback, "self_reflection", "Critiquing claims, citations, method transfers, and confidence.", 92)
    reflection_bundle = run_self_reflection(
        bacteria=bacteria,
        metabolite=metabolite,
        disease=disease,
        protocol_text=protocol_text,
        audit_bundle=audit_bundle,
        literature_context=literature_context,
        evidence_strength_map=evidence_strength_map,
        query_log=query_log,
        protocols_io_evidence=protocols_io_evidence,
        in_vitro_plan=vitro_bundle.get("in_vitro_plan") or [],
        in_vivo_plan=vivo_bundle.get("in_vivo_plan") or [],
        research_question=research_question,
        prompt_constraints=prompt_constraints,
        direct_production_evidence_assessment=hypothesis_bundle.get("direct_production_evidence_assessment"),
        hypothesis_branches=hypothesis_bundle.get("hypothesis_branches"),
    )

    _report_progress(progress_callback, "assemble", "Assembling final validation plan payload.", 97)
    plan = assemble_validation_plan(
        bacteria=bacteria,
        metabolite=metabolite,
        disease=disease,
        candidate_metrics=candidate_metrics,
        literature=merged_literature,
        mode=mode,
        iterative_query_log=query_log,
        evidence_strength_map=evidence_strength_map,
        fulltext_method_evidence=fulltext_method_evidence,
        protocols_io_evidence=protocols_io_evidence,
        audit_bundle=audit_bundle,
        vitro_bundle=vitro_bundle,
        vivo_bundle=vivo_bundle,
        reflection_bundle=reflection_bundle,
        hypothesis_bundle=hypothesis_bundle,
        research_question=research_question,
        prompt_constraints=prompt_constraints,
    )
    plan = _project_public_plan(plan)
    _report_progress(progress_callback, "completed", "Validation plan is ready.", 100)
    return plan


def generate_question_driven_validation_plan(
    research_question: str,
    prompt_constraints: str = "",
    disease: str = "",
    progress_callback: Optional[Callable[[str, str, int, Optional[Dict[str, Any]]], None]] = None,
) -> dict:
    protocol_text = "Standalone question-driven request. There is no prior protocol draft. Build the protocol directly from the user brief and retrieved evidence."
    _report_progress(progress_callback, "prepare_question", "Parsing standalone validation question and building the search focus.", 5)
    question_profile = classify_question_profile(research_question=research_question, prompt_constraints=prompt_constraints, disease=disease)
    candidate_like_question = _question_is_microbe_metabolite_study(
        research_question, prompt_constraints
    )
    extracted_candidate = _extract_question_candidate_entities(
        research_question,
        prompt_constraints,
        disease,
        question_profile,
    )
    candidate_entities_resolved = bool(
        candidate_like_question
        and extracted_candidate.get("bacteria")
        and extracted_candidate.get("metabolite")
    )
    design_bacteria = (
        extracted_candidate.get("bacteria")
        or ("the candidate microbe specified in the research question" if candidate_like_question else "")
    )
    design_metabolite = (
        extracted_candidate.get("metabolite")
        or ("the candidate metabolite specified in the research question" if candidate_like_question else "")
    )
    design_disease = extracted_candidate.get("disease") or disease
    question_specs = build_question_only_query_specs(
        research_question=research_question,
        prompt_constraints=prompt_constraints,
        disease=disease,
    )
    round_one_specs: List[dict] = []
    seen_queries: set[str] = set()
    for spec in question_specs:
        query = str(spec.get("query") or "").strip()
        if not query or query.lower() in seen_queries:
            continue
        seen_queries.add(query.lower())
        round_one_specs.append(spec)

    _report_progress(progress_callback, "round1_search", "Running first-round PubMed queries for the standalone question.", 22)
    round_one = search_literature_online(round_one_specs, max_results=QUESTION_MAX_QUERY_RESULTS)
    initial_literature = _filter_question_relevant_articles(
        round_one["articles"], research_question, design_disease
    )
    initial_strength_map = (
        rank_evidence_strength(
            initial_literature,
            extracted_candidate.get("bacteria") or "",
            extracted_candidate.get("metabolite") or "",
            design_disease,
        )
        if candidate_entities_resolved
        else rank_question_evidence_strength(initial_literature, design_disease)
    )

    _report_progress(
        progress_callback,
        "audit_scope",
        "Auditing the standalone question scope against first-round evidence.",
        35,
        {
            "initial_articles": len(initial_literature),
            "round1_queries": len(round_one["query_log"]),
        },
    )
    audit_bundle = audit_question_scope(
        research_question=research_question,
        prompt_constraints=prompt_constraints,
        literature_context=_format_article_context(initial_literature, limit=6),
        evidence_strength_map=initial_strength_map,
        disease=design_disease,
    )

    _report_progress(progress_callback, "followup_queries", "Generating follow-up literature questions for the standalone request.", 48)
    followup_specs = generate_question_followup_queries(
        research_question=research_question,
        prompt_constraints=prompt_constraints,
        seed_query=(round_one_specs[0].get("query") if round_one_specs else "") or _fallback_question_query(research_question, disease),
        evidence_strength_map=initial_strength_map,
    )
    _report_progress(progress_callback, "round2_search", "Running follow-up PubMed queries and merging evidence.", 58)
    followup_round = search_literature_online(followup_specs, max_results=QUESTION_MAX_QUERY_RESULTS)
    merged_literature = _filter_question_relevant_articles(
        _merge_articles(initial_literature, followup_round["articles"]),
        research_question,
        design_disease,
    )
    query_log = round_one["query_log"] + followup_round["query_log"]
    evidence_strength_map = (
        rank_evidence_strength(
            merged_literature,
            extracted_candidate.get("bacteria") or "",
            extracted_candidate.get("metabolite") or "",
            design_disease,
        )
        if candidate_entities_resolved
        else rank_question_evidence_strength(merged_literature, design_disease)
    )
    literature_context = _format_article_context(merged_literature, limit=8)

    _report_progress(progress_callback, "fulltext_methods", "Retrieving open-access full-text methods for concentrations and model details.", 66)
    fulltext_method_evidence = collect_fulltext_method_evidence(
        merged_literature,
        extracted_candidate.get("bacteria") or "",
        extracted_candidate.get("metabolite") or "",
        design_disease,
    )
    _report_progress(progress_callback, "protocols_io_search", "Searching protocols.io for operational procedures and materials.", 69)
    protocols_io_evidence = collect_protocol_evidence(
        bacteria=extracted_candidate.get("bacteria") or "",
        metabolite=extracted_candidate.get("metabolite") or "",
        disease=design_disease,
        research_question=research_question,
    )
    hypothesis_bundle = (
        assess_direct_production_evidence_and_hypotheses(
            bacteria=extracted_candidate.get("bacteria") or "",
            metabolite=extracted_candidate.get("metabolite") or "",
            disease=design_disease,
            literature=merged_literature,
            literature_context=_format_article_context(merged_literature, limit=12),
            fulltext_method_evidence=fulltext_method_evidence,
            evidence_strength_map=evidence_strength_map,
        )
        if candidate_entities_resolved
        else {
            "direct_production_evidence_assessment": _normalize_production_evidence_assessment({}),
            "hypothesis_branches": [],
        }
    )
    _report_progress(
        progress_callback,
        "protocols_io_search",
        protocols_io_evidence.get("message") or "protocols.io search finished.",
        71,
        {
            "status": protocols_io_evidence.get("status"),
            "protocol_count": len(protocols_io_evidence.get("all", [])),
        },
    )
    _report_progress(
        progress_callback,
        "in_vitro_design",
        "Expanding in vitro and cell experiment package.",
        73,
        {
            "total_articles": len(merged_literature),
            "total_queries": len(query_log),
            "fulltext_method_hits": len(fulltext_method_evidence.get("all", [])),
            "protocols_io_hits": len(protocols_io_evidence.get("all", [])),
            "protocols_io_status": protocols_io_evidence.get("status"),
        },
    )
    vitro_bundle = expand_in_vitro_plan(
        bacteria=design_bacteria,
        metabolite=design_metabolite,
        disease=design_disease,
        protocol_text=protocol_text,
        audit_bundle=audit_bundle,
        literature_context=literature_context,
        evidence_strength_map=evidence_strength_map,
        query_log=query_log,
        fulltext_method_evidence=fulltext_method_evidence.get("in_vitro", []),
        protocols_io_evidence=protocols_io_evidence.get("in_vitro", []),
        research_question=research_question,
        prompt_constraints=prompt_constraints,
        direct_production_evidence_assessment=hypothesis_bundle.get("direct_production_evidence_assessment"),
        hypothesis_branches=hypothesis_bundle.get("hypothesis_branches"),
        question_profile=question_profile,
    )
    _report_progress(progress_callback, "in_vivo_design", "Expanding in vivo and animal experiment package.", 82)
    vivo_bundle = expand_in_vivo_plan(
        bacteria=design_bacteria,
        metabolite=design_metabolite,
        disease=design_disease,
        protocol_text=protocol_text,
        audit_bundle=audit_bundle,
        literature_context=literature_context,
        evidence_strength_map=evidence_strength_map,
        query_log=query_log,
        fulltext_method_evidence=fulltext_method_evidence.get("in_vivo", []),
        protocols_io_evidence=protocols_io_evidence.get("in_vivo", []),
        research_question=research_question,
        prompt_constraints=prompt_constraints,
        direct_production_evidence_assessment=hypothesis_bundle.get("direct_production_evidence_assessment"),
        hypothesis_branches=hypothesis_bundle.get("hypothesis_branches"),
        question_profile=question_profile,
    )
    if candidate_like_question:
        vitro_bundle["in_vitro_plan"] = _assign_candidate_experiment_roles(
            vitro_bundle.get("in_vitro_plan") or [],
            is_in_vivo=False,
        )
        vivo_bundle["in_vivo_plan"] = _assign_candidate_experiment_roles(
            vivo_bundle.get("in_vivo_plan") or [],
            is_in_vivo=True,
        )
    _report_progress(progress_callback, "self_reflection", "Critiquing claims, citations, method transfers, and confidence.", 92)
    reflection_bundle = run_self_reflection(
        bacteria=design_bacteria,
        metabolite=design_metabolite,
        disease=design_disease,
        protocol_text=protocol_text,
        audit_bundle=audit_bundle,
        literature_context=literature_context,
        evidence_strength_map=evidence_strength_map,
        query_log=query_log,
        protocols_io_evidence=protocols_io_evidence,
        in_vitro_plan=vitro_bundle.get("in_vitro_plan") or [],
        in_vivo_plan=vivo_bundle.get("in_vivo_plan") or [],
        research_question=research_question,
        prompt_constraints=prompt_constraints,
        direct_production_evidence_assessment=hypothesis_bundle.get("direct_production_evidence_assessment"),
        hypothesis_branches=hypothesis_bundle.get("hypothesis_branches"),
    )
    _report_progress(progress_callback, "assemble", "Assembling final validation protocol.", 97)
    plan = assemble_validation_plan(
        bacteria=extracted_candidate.get("bacteria") if candidate_entities_resolved else "",
        metabolite=extracted_candidate.get("metabolite") if candidate_entities_resolved else "",
        disease=design_disease,
        candidate_metrics={},
        literature=merged_literature,
        mode="question_driven",
        iterative_query_log=query_log,
        evidence_strength_map=evidence_strength_map,
        fulltext_method_evidence=fulltext_method_evidence,
        protocols_io_evidence=protocols_io_evidence,
        audit_bundle=audit_bundle,
        vitro_bundle=vitro_bundle,
        vivo_bundle=vivo_bundle,
        reflection_bundle=reflection_bundle,
        hypothesis_bundle=hypothesis_bundle,
        research_question=research_question,
        prompt_constraints=prompt_constraints,
    )
    plan["question_profile"] = question_profile
    plan["human_plan"] = _prepare_experiment_modules(
        plan.get("human_plan") or _build_question_driven_human_experiments(plan),
        is_in_vivo=True,
        is_human=True,
    )
    plan = _sanitize_question_plan_references(plan)
    plan = _refresh_plan_level_structure(
        plan,
        candidate_workflow=candidate_like_question,
    )
    plan["validation_summary"] = _build_validation_summary(plan)
    plan["validation_protocol_text"] = _render_validation_protocol_text(plan)
    plan = _project_public_plan(plan)
    _report_progress(progress_callback, "completed", "Standalone validation protocol is ready.", 100)
    return plan
