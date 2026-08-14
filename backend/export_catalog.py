"""
Archive pipeline JSON outputs into data/results_catalog/ for browse & reuse.
"""

from __future__ import annotations

import json
import re
import shutil
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, List, Optional

BACKEND_DIR = Path(__file__).parent
PROJECT_DIR = BACKEND_DIR.parent
CATALOG_DIR = PROJECT_DIR / "data" / "results_catalog"
RUNS_DIR = CATALOG_DIR / "runs"
INDEX_PATH = CATALOG_DIR / "index.json"

QUADRANT_ORDER = {"I": 0, "II": 1, "III": 2, "IV": 3}

ARCHIVE_FILES = [
    "step1_candidates.json",
    "step2_chain_novelty.json",
    "step2b_agent_evidence.json",
    "step3_quadrant.json",
    "pipeline_status.json",
    "knowledge_graph.gml",
]


def _slug(text: str, max_len: int = 24) -> str:
    s = re.sub(r"[^\w]+", "_", (text or "run").strip()).strip("_")
    return (s[:max_len] if s else "run").lower()


def _portable_source_path(path: str | Path) -> str:
    source = Path(path)
    try:
        return source.resolve().relative_to(PROJECT_DIR.resolve()).as_posix()
    except ValueError:
        return source.as_posix()


def _load_json(path: Path) -> Any:
    if not path.exists():
        return None
    with open(path, "r", encoding="utf-8") as f:
        return json.load(f)


def _save_json(path: Path, data: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8") as f:
        json.dump(data, f, indent=2, ensure_ascii=False)


def _infer_bacteria_hint(output_dir: Path, bacteria_filter: Optional[str]) -> str:
    if bacteria_filter:
        return bacteria_filter
    step3 = _load_json(output_dir / "step3_quadrant.json")
    if isinstance(step3, list) and step3:
        return step3[0].get("bacteria", "mixed")
    step1 = _load_json(output_dir / "step1_candidates.json")
    if isinstance(step1, list) and step1:
        return step1[0].get("bacteria", "mixed")
    return "mixed"


def make_run_id(disease: str, bacteria_hint: str) -> str:
    ts = datetime.now().strftime("%Y%m%d_%H%M%S")
    return f"run_{ts}_{_slug(bacteria_hint, 20)}_{_slug(disease, 12)}"


def _index_record(run_id: str, item: dict, meta: dict) -> dict:
    cand = item.get("candidate") or {}
    return {
        "run_id": run_id,
        "bacteria": item.get("bacteria") or cand.get("bacteria", ""),
        "metabolite": item.get("metabolite") or cand.get("metabolite", ""),
        "disease": item.get("disease") or meta.get("disease", ""),
        "quadrant": item.get("quadrant", ""),
        "quadrant_label": item.get("quadrant_label", ""),
        "composite_score": item.get("composite_score"),
        "mmsage_norm": item.get("mmsage_norm"),
        "chain_novelty": item.get("chain_novelty"),
        "chain_count": item.get("chain_count"),
        "is_dark_matter": item.get("is_dark_matter", False),
        "pair_bm_exp": item.get("pair_bm_exp"),
        "pair_md_exp": item.get("pair_md_exp"),
        "evidence_foundation": item.get("evidence_foundation"),
        "has_path": item.get("has_path", False),
        "computed_at": meta.get("created_at", ""),
    }


def rebuild_index() -> dict:
    """Rebuild catalog index.json from all archived runs."""
    records: List[dict] = []
    runs_meta: List[dict] = []

    if RUNS_DIR.exists():
        for run_dir in sorted(RUNS_DIR.iterdir()):
            if not run_dir.is_dir():
                continue
            meta = _load_json(run_dir / "meta.json") or {}
            meta["run_id"] = meta.get("run_id") or run_dir.name
            runs_meta.append(meta)

            candidates = _load_json(run_dir / "candidates.json")
            if not isinstance(candidates, list):
                continue
            for item in candidates:
                if isinstance(item, dict):
                    records.append(_index_record(run_dir.name, item, meta))

    records.sort(
        key=lambda r: (r.get("computed_at") or "", r.get("composite_score") or 0),
        reverse=True,
    )

    index = {
        "version": 1,
        "updated_at": datetime.now(timezone.utc).isoformat(),
        "total_records": len(records),
        "run_count": len(runs_meta),
        "runs": runs_meta,
        "records": records,
    }
    _save_json(INDEX_PATH, index)
    return index


def archive_run(
    output_dir: str | Path,
    disease: str = "IBD",
    coordinates_file: Optional[str] = None,
    bacteria_filter: Optional[str] = None,
    pipeline_status: Optional[dict] = None,
    run_id: Optional[str] = None,
) -> str:
    """
    Copy pipeline outputs into data/results_catalog/runs/<run_id>/ and refresh index.
    Returns the run_id used.
    """
    out = Path(output_dir)
    if not (out / "step3_quadrant.json").exists():
        raise FileNotFoundError(f"step3_quadrant.json not found in {out}")

    bacteria_hint = _infer_bacteria_hint(out, bacteria_filter)
    rid = run_id or make_run_id(disease, bacteria_hint)

    run_dir = RUNS_DIR / rid
    if run_dir.exists():
        shutil.rmtree(run_dir)
    run_dir.mkdir(parents=True, exist_ok=True)

    copied = []
    for fname in ARCHIVE_FILES:
        src = out / fname
        if src.exists():
            shutil.copy2(src, run_dir / fname)
            copied.append(fname)

    step3_path = run_dir / "step3_quadrant.json"
    if step3_path.exists():
        shutil.copy2(step3_path, run_dir / "candidates.json")

    candidates = _load_json(run_dir / "candidates.json") or []
    status = pipeline_status or _load_json(out / "pipeline_status.json") or {}

    coord_name = ""
    if coordinates_file:
        coord_name = _portable_source_path(coordinates_file)

    meta = {
        "run_id": rid,
        "disease": disease,
        "bacteria_hint": bacteria_hint,
        "bacteria_filter": bacteria_filter,
        "coordinates_file": coord_name,
        "candidate_count": len(candidates) if isinstance(candidates, list) else 0,
        "archived_files": copied,
        "pipeline_status": status.get("status"),
        "total_duration_s": status.get("total_duration_s"),
        "quadrant_counts": (status.get("steps") or {}).get("step3", {}).get("quadrant_counts"),
        "created_at": datetime.now(timezone.utc).isoformat(),
    }
    _save_json(run_dir / "meta.json", meta)

    rebuild_index()
    print(f"[catalog] Archived run {rid} ({meta['candidate_count']} candidates)")
    return rid


def list_runs() -> List[dict]:
    index = _load_json(INDEX_PATH)
    if isinstance(index, dict) and index.get("runs"):
        return index["runs"]
    runs = []
    if RUNS_DIR.exists():
        for run_dir in sorted(RUNS_DIR.iterdir(), reverse=True):
            if run_dir.is_dir():
                meta = _load_json(run_dir / "meta.json")
                if meta:
                    runs.append(meta)
    return runs


def _match_text(hay: str, needle: str) -> bool:
    return needle.lower() in (hay or "").lower()


def search_records(
    q: str = "",
    bacteria: str = "",
    metabolite: str = "",
    quadrant: str = "",
    disease: str = "",
    run_id: str = "",
    limit: int = 100,
    offset: int = 0,
) -> dict:
    index = _load_json(INDEX_PATH) or rebuild_index()
    records = index.get("records") or []

    filtered = []
    for r in records:
        if run_id and r.get("run_id") != run_id:
            continue
        if quadrant and r.get("quadrant") != quadrant:
            continue
        if disease and not _match_text(r.get("disease", ""), disease):
            continue
        if bacteria and not _match_text(r.get("bacteria", ""), bacteria):
            continue
        if metabolite and not _match_text(r.get("metabolite", ""), metabolite):
            continue
        if q:
            blob = f"{r.get('bacteria','')} {r.get('metabolite','')} {r.get('disease','')}"
            if not _match_text(blob, q):
                continue
        filtered.append(r)

    filtered.sort(
        key=lambda r: (
            QUADRANT_ORDER.get(r.get("quadrant"), 9),
            -(r.get("mmsage_norm") or 0),
        )
    )

    total = len(filtered)
    page = filtered[offset : offset + limit]
    return {
        "total": total,
        "offset": offset,
        "limit": limit,
        "records": page,
    }


def get_run_candidates(run_id: str) -> List[dict]:
    path = RUNS_DIR / run_id / "candidates.json"
    data = _load_json(path)
    return data if isinstance(data, list) else []


def _pair_keys(row: dict) -> tuple:
    """Resolve bacteria / metabolite from a pipeline row or nested candidate."""
    b = row.get("bacteria")
    m = row.get("metabolite")
    cand = row.get("candidate")
    if isinstance(cand, dict):
        b = b or cand.get("bacteria")
        m = m or cand.get("metabolite")
    sr = row.get("scoresp_result")
    if isinstance(sr, dict):
        c2 = sr.get("candidate")
        if isinstance(c2, dict):
            b = b or c2.get("bacteria")
            m = m or c2.get("metabolite")
    return b, m


def _find_in_list(rows: Any, bacteria: str, metabolite: str) -> Optional[dict]:
    if not isinstance(rows, list):
        return None
    for row in rows:
        if not isinstance(row, dict):
            continue
        b, m = _pair_keys(row)
        if b == bacteria and m == metabolite:
            return row
    return None


def _index_record_for_pair(run_id: str, bacteria: str, metabolite: str) -> Optional[dict]:
    index = _load_json(INDEX_PATH)
    if not isinstance(index, dict):
        return None
    for rec in index.get("records") or []:
        if (
            rec.get("run_id") == run_id
            and rec.get("bacteria") == bacteria
            and rec.get("metabolite") == metabolite
        ):
            return rec
    return None


def build_record_export(run_id: str, bacteria: str, metabolite: str) -> Optional[dict]:
    """
    Extract this microbe–metabolite pair from archived run JSON files and merge
    into one downloadable document.
    """
    run_dir = RUNS_DIR / run_id
    if not run_dir.is_dir():
        return None

    step3_row = _find_in_list(_load_json(run_dir / "candidates.json"), bacteria, metabolite)
    if not step3_row:
        step3_path = run_dir / "step3_quadrant.json"
        if step3_path.exists():
            step3_row = _find_in_list(_load_json(step3_path), bacteria, metabolite)
    if not step3_row and not _find_in_list(_load_json(run_dir / "step1_candidates.json"), bacteria, metabolite):
        return None

    meta = _load_json(run_dir / "meta.json") or {}
    pipeline_status = _load_json(run_dir / "pipeline_status.json")

    sources: Dict[str, Any] = {
        "step1_candidates": _find_in_list(_load_json(run_dir / "step1_candidates.json"), bacteria, metabolite),
        "step2_chain_novelty": _find_in_list(_load_json(run_dir / "step2_chain_novelty.json"), bacteria, metabolite),
        "step2b_agent_evidence": _find_in_list(_load_json(run_dir / "step2b_agent_evidence.json"), bacteria, metabolite),
        "step3_quadrant": step3_row or _find_in_list(_load_json(run_dir / "candidates.json"), bacteria, metabolite),
    }

    merged = get_record(run_id, bacteria, metabolite)
    if merged and "step3_detail" in merged:
        merged = {k: v for k, v in merged.items() if k != "step3_detail"}

    return {
        "export_version": 1,
        "exported_at": datetime.now(timezone.utc).isoformat(),
        "run_id": run_id,
        "bacteria": bacteria,
        "metabolite": metabolite,
        "disease": (merged or {}).get("disease") or meta.get("disease", ""),
        "run_meta": meta,
        "pipeline_status": pipeline_status,
        "catalog_index_record": _index_record_for_pair(run_id, bacteria, metabolite),
        "source_files": {
            **{
                name: f"runs/{run_id}/{name}.json" if (run_dir / f"{name}.json").exists() else None
                for name in [
                    "step1_candidates",
                    "step2_chain_novelty",
                    "step2b_agent_evidence",
                    "candidates",
                    "pipeline_status",
                    "meta",
                ]
            },
            "knowledge_graph": f"runs/{run_id}/knowledge_graph.gml" if (run_dir / "knowledge_graph.gml").exists() else None,
        },
        "extracted": {k: v for k, v in sources.items() if v is not None},
        "merged_record": merged,
    }


def export_download_filename(run_id: str, bacteria: str, metabolite: str) -> str:
    base = _slug(f"{bacteria}_{metabolite}_{run_id}", 80)
    return f"catalog_{base}.json"


def _merge_step2_fields(item: dict, step2_row: dict) -> None:
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
        if step2_row.get(key) is not None:
            item[key] = step2_row[key]


def _build_step3_detail(item: dict) -> dict:
    """Structured evidence payload for catalog detail UI (full archived content)."""
    ae = item.get("agent_evidence") or {}
    if not isinstance(ae, dict):
        ae = {}

    return {
        "recommendation": ae.get("recommendation") or item.get("agent_recommendation"),
        "scoring": {
            "chain_count": ae.get("chain_count", item.get("chain_count")),
            "chain_count_raw": ae.get("chain_count_raw"),
            "chain_novelty": ae.get("chain_novelty", item.get("chain_novelty")),
            "bottleneck": ae.get("bottleneck"),
            "db_bonus": ae.get("db_bonus"),
            "hop_counts": ae.get("hop_counts") or {},
            "db_hits": ae.get("db_hits") or {},
        },
        "sources": ae.get("sources") or [],
        "hop_evidence": ae.get("hop_evidence") or [],
        "db_details": ae.get("db_details") or {},
        "pairwise_counts": item.get("pairwise_counts") or {},
        "edge_cooccurrences": item.get("edge_cooccurrences") or [],
        "bottleneck_edge": item.get("bottleneck_edge"),
        "chain_query": item.get("chain_query"),
        "chain_path_str": item.get("chain_path_str"),
        "chain_path": item.get("chain_path") or [],
        "query_terms": item.get("query_terms") or [],
        "all_paths_info": item.get("all_paths_info") or [],
    }


def get_record(run_id: str, bacteria: str, metabolite: str) -> Optional[dict]:
    candidates = get_run_candidates(run_id)
    item = None
    for c in candidates:
        if c.get("bacteria") == bacteria and c.get("metabolite") == metabolite:
            item = dict(c)
            break
    if not item:
        return None

    meta = _load_json(RUNS_DIR / run_id / "meta.json") or {}
    if not item.get("disease"):
        item["disease"] = meta.get("disease", "")

    step2 = _load_json(RUNS_DIR / run_id / "step2_chain_novelty.json")
    if isinstance(step2, list):
        for s in step2:
            if s.get("bacteria") == bacteria and s.get("metabolite") == metabolite:
                _merge_step2_fields(item, s)
                break

    if not item.get("agent_evidence"):
        step2b = _load_json(RUNS_DIR / run_id / "step2b_agent_evidence.json")
        if isinstance(step2b, list):
            for s in step2b:
                if s.get("bacteria") == bacteria and s.get("metabolite") == metabolite:
                    item["agent_evidence"] = s
                    if s.get("recommendation"):
                        item.setdefault("agent_recommendation", s["recommendation"])
                    break

    ae = item.get("agent_evidence") or {}
    if isinstance(ae, dict) and ae.get("chain_novelty") is not None:
        item["agent_chain_novelty"] = ae["chain_novelty"]
    item["step3_detail"] = _build_step3_detail(item)
    return item


def run_dir_path(run_id: str) -> Path:
    return RUNS_DIR / run_id


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Archive pipeline outputs to results catalog")
    parser.add_argument("--output-dir", default=str(PROJECT_DIR / "outputs"))
    parser.add_argument("--disease", default="IBD")
    parser.add_argument("--coordinates-file", default=None)
    parser.add_argument("--bacteria-filter", default=None)
    parser.add_argument("--rebuild-only", action="store_true")
    args = parser.parse_args()

    if args.rebuild_only:
        idx = rebuild_index()
        print(f"Rebuilt index: {idx['total_records']} records, {idx['run_count']} runs")
    else:
        rid = archive_run(
            args.output_dir,
            disease=args.disease,
            coordinates_file=args.coordinates_file,
            bacteria_filter=args.bacteria_filter,
        )
        print(f"Done: {rid}")
