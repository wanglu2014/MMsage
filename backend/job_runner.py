from __future__ import annotations

import argparse
import threading
import time
from pathlib import Path

from runtime_state import job_output_dir, job_status_path, load_json, save_json


def update_status(job_id: str, **fields) -> dict:
    path = job_status_path(job_id)
    current = load_json(path)
    if not isinstance(current, dict):
        current = {"job_id": job_id}
    current.update(fields)
    save_json(path, current)
    return current


def infer_step(output_dir: Path) -> str:
    if (output_dir / "step3_quadrant.json").exists():
        return "step3 (quadrant assignment)"
    if (output_dir / "step2b_agent_evidence.json").exists():
        return "step2b (multi-agent evidence)"
    if (output_dir / "step2_chain_novelty.json").exists():
        return "step2b (multi-agent evidence)"
    if (output_dir / "step1_candidates.json").exists():
        return "step2 (PubMed queries, may take 3-5 min)"
    return "step1 (MMSage signal)"


def run_job(
    job_id: str,
    csv_path: str,
    disease: str,
    bacteria_filter: str | None = None,
    max_queries: int = 15,
    max_articles_per_query: int = 20,
) -> dict:
    from run_pipeline import run_pipeline

    output_dir = job_output_dir(job_id)
    output_dir.mkdir(parents=True, exist_ok=True)

    update_status(
        job_id,
        status="running",
        step="step1 (MMSage signal)",
        disease=disease,
        run_id="",
        output_dir=str(output_dir),
        error=None,
    )

    result_holder: dict = {"result": None, "error": None}

    def _inner() -> None:
        try:
            result_holder["result"] = run_pipeline(
                coordinates_file=csv_path,
                output_dir=str(output_dir),
                top_n=200,
                disease=disease,
                max_depth=3,
                bacteria_filter=bacteria_filter,
                max_queries=max_queries,
                max_articles_per_query=max_articles_per_query,
            )
        except Exception as exc:  # pragma: no cover - subprocess runtime guard
            result_holder["error"] = exc

    worker = threading.Thread(target=_inner, daemon=True)
    worker.start()

    while worker.is_alive():
        time.sleep(2)
        update_status(job_id, status="running", step=infer_step(output_dir))

    worker.join()

    if result_holder["error"] is not None:
        raise result_holder["error"]

    status_blob = load_json(output_dir / "pipeline_status.json") or {}
    run_id = status_blob.get("catalog_run_id", "")
    candidates = load_json(output_dir / "step3_quadrant.json") or []

    final = update_status(
        job_id,
        status="done",
        step="complete",
        disease=disease,
        run_id=run_id,
        candidates=len(candidates) if isinstance(candidates, list) else 0,
        total_duration_s=status_blob.get("total_duration_s"),
        error=None,
    )
    return final


def main() -> int:
    parser = argparse.ArgumentParser(description="Run an isolated MMSage pipeline job")
    parser.add_argument("--job-id", required=True)
    parser.add_argument("--csv-path", required=True)
    parser.add_argument("--disease", default="IBD")
    parser.add_argument("--bacteria-filter", default=None)
    parser.add_argument("--max-queries", type=int, default=15)
    parser.add_argument("--max-articles-per-query", type=int, default=20)
    args = parser.parse_args()

    try:
        run_job(
            job_id=args.job_id,
            csv_path=args.csv_path,
            disease=args.disease,
            bacteria_filter=args.bacteria_filter,
            max_queries=args.max_queries,
            max_articles_per_query=args.max_articles_per_query,
        )
        return 0
    except Exception as exc:  # pragma: no cover - subprocess runtime guard
        update_status(
            args.job_id,
            status="error",
            step="",
            disease=args.disease,
            run_id="",
            error=str(exc),
        )
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
