"""
Pipeline Orchestrator
======================
Run step1 -> step2 -> step2b -> step3 sequentially.
"""

from __future__ import annotations

import io
import json
import socket
import subprocess
import sys
import time
from pathlib import Path
from typing import Optional

socket.setdefaulttimeout(30.0)


if sys.stdout.encoding and sys.stdout.encoding.lower() != "utf-8":
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")
    sys.stderr = io.TextIOWrapper(sys.stderr.buffer, encoding="utf-8", errors="replace")


from agents.master_agent import MasterAgent
from step1_mmsage_signal import run_step1
from step2_chain_novelty import run_step2
from step3_quadrant import run_step3


DEFAULT_COORDINATES_DIR = str(Path(__file__).parent.parent / "data" / "sample_coordinates")
DEFAULT_COORDINATES_FILE = str(
    Path(__file__).parent.parent
    / "data"
    / "sample_coordinates"
    / "pluscombno1V0317_min0_Akkermansia_muciniphila_seed_1_dim_10_neighbor_2_dist_0.4_metric_euclidean_rank_rootknow_cor0303_1_top_50_Pthre_0.1_pair.csv_clu_1_coordinates_tunek.csv"
)
DEFAULT_OUTPUT_DIR = "outputs"
DEFAULT_DISEASE = "IBD"
AUTO_GML_REL = Path("data") / "knowledge_graph" / "auto_built_kg.gml"
TEST_GML_REL = Path("data") / "knowledge_graph" / "ibd_test_kg.gml"


def resolve_default_gml(project_dir: Path) -> str:
    auto_gml = project_dir / AUTO_GML_REL
    test_gml = project_dir / TEST_GML_REL
    return str(auto_gml if auto_gml.exists() else test_gml)


def run_pipeline(
    coordinates_dir: Optional[str] = None,
    gml_path: Optional[str] = None,
    output_dir: Optional[str] = None,
    top_n: int = 200,
    disease: str = DEFAULT_DISEASE,
    max_depth: int = 3,
    bacteria_filter: Optional[str] = None,
    coordinates_file: Optional[str] = None,
    max_queries: int = 15,
    max_articles_per_query: int = 20,
) -> dict:
    backend_dir = Path(__file__).parent
    project_dir = backend_dir.parent

    coords = coordinates_dir or DEFAULT_COORDINATES_DIR
    coord_file = coordinates_file or DEFAULT_COORDINATES_FILE
    gml = gml_path or resolve_default_gml(project_dir)
    out_dir = Path(output_dir or (project_dir / DEFAULT_OUTPUT_DIR))
    out_dir.mkdir(parents=True, exist_ok=True)

    step1_out = str(out_dir / "step1_candidates.json")
    step2_out = str(out_dir / "step2_chain_novelty.json")
    step2b_out = str(out_dir / "step2b_agent_evidence.json")
    step3_out = str(out_dir / "step3_quadrant.json")

    status = {
        "status": "running",
        "current_step": 0,
        "steps": {},
        "start_time": time.time(),
        "disease": disease,
    }

    print("\n" + "=" * 70)
    print("  MMSage x Chain Novelty DUAL-AXIS PIPELINE")
    print("=" * 70)
    print(f"\nCoordinates: {coord_file or coords}")
    print(f"Knowledge Graph: {gml}")
    print(f"Disease context: {disease}")
    print(f"Output: {out_dir}")
    print(f"Top-N: {top_n}, Max-Depth: {max_depth}")

    try:
        status["current_step"] = 1
        t0 = time.time()
        candidates = run_step1(
            coords,
            step1_out,
            top_n,
            bacteria_filter,
            coordinates_file=coord_file,
        )
        status["steps"]["step1"] = {
            "status": "completed",
            "candidates": len(candidates),
            "duration_s": round(time.time() - t0, 1),
        }

        print(f"\n[Step 1.5] Building a run-specific knowledge graph for {disease}...")
        target_gml_path = str(out_dir / "knowledge_graph.gml")
        try:
            subprocess.run(
                [
                    sys.executable,
                    str(backend_dir / "build_kg.py"),
                    "candidates",
                    "--candidates",
                    step1_out,
                    "--disease",
                    disease,
                    "--output",
                    target_gml_path,
                    "--max-queries",
                    str(max_queries),
                    "--max-articles",
                    str(max_articles_per_query),
                ],
                check=True,
            )
            print("[OK] Run-specific knowledge graph build completed.")
            gml = target_gml_path
        except subprocess.CalledProcessError as e:
            print(f"[WARN] Run-specific knowledge graph build failed: {e}")
            print("[INFO] Falling back to the existing graph for the remaining steps.")

        status["current_step"] = 2
        t0 = time.time()
        novelty_results = run_step2(
            step1_out,
            gml,
            step2_out,
            disease=disease,
            max_depth=max_depth,
        )
        has_path_n = sum(1 for row in novelty_results if row.get("has_path"))
        dark_n = sum(1 for row in novelty_results if row.get("chain_count", 0) == 0)
        status["steps"]["step2"] = {
            "status": "completed",
            "scored": len(novelty_results),
            "has_path": has_path_n,
            "dark_matter_chains": dark_n,
            "duration_s": round(time.time() - t0, 1),
        }

        status["current_step"] = "2b"
        t0 = time.time()
        master = MasterAgent(c_max=500, parallel=True)
        agent_results = []
        for idx, row in enumerate(novelty_results, start=1):
            bacteria = row.get("bacteria", row.get("candidate", {}).get("bacteria", ""))
            metabolite = row.get("metabolite", row.get("candidate", {}).get("metabolite", ""))
            if not bacteria or not metabolite:
                continue

            print(f"  [Step 2b] Calling AI for candidate {idx}/{len(novelty_results)}: {bacteria} x {metabolite} ...")
            agent_res = master.run(bacteria, metabolite, disease)
            row["agent_evidence"] = agent_res
            row["agent_chain_novelty"] = agent_res["chain_novelty"]
            row["agent_chain_count"] = agent_res["chain_count"]
            row["agent_bottleneck"] = agent_res["bottleneck"]
            row["agent_hop_counts"] = agent_res["hop_counts"]
            row["agent_db_bonus"] = agent_res["db_bonus"]
            row["agent_sources"] = agent_res["sources"]
            row["agent_recommendation"] = agent_res["recommendation"]
            agent_results.append(agent_res)

        with open(step2b_out, "w", encoding="utf-8") as f:
            json.dump(agent_results, f, indent=2, ensure_ascii=False)

        with open(step2_out, "w", encoding="utf-8") as f:
            json.dump(novelty_results, f, indent=2, ensure_ascii=False)

        status["steps"]["step2b"] = {
            "status": "completed",
            "agents_run": len(agent_results),
            "duration_s": round(time.time() - t0, 1),
        }

        status["current_step"] = 3
        t0 = time.time()
        quadrants = run_step3(step2_out, step3_out, disease=disease)
        quadrant_counts = {}
        for row in quadrants:
            qid = row.get("quadrant", "?")
            quadrant_counts[qid] = quadrant_counts.get(qid, 0) + 1

        status["steps"]["step3"] = {
            "status": "completed",
            "quadrant_counts": quadrant_counts,
            "dark_matter": sum(1 for row in quadrants if row.get("is_dark_matter")),
            "duration_s": round(time.time() - t0, 1),
        }

        status["status"] = "completed"
        status["current_step"] = 3
        status["total_duration_s"] = round(time.time() - status["start_time"], 1)

        print("\n" + "=" * 70)
        print("  PIPELINE COMPLETE")
        print("=" * 70)
        print(f"\nTotal duration: {status['total_duration_s']}s")
        for step_name, step_info in status["steps"].items():
            print(f"  {step_name}: {step_info['duration_s']}s")
        print(f"\nQuadrant distribution: {quadrant_counts}")

        try:
            from export_catalog import archive_run

            catalog_run_id = archive_run(
                output_dir=out_dir,
                disease=disease,
                coordinates_file=coord_file,
                bacteria_filter=bacteria_filter,
                pipeline_status=status,
            )
            status["catalog_run_id"] = catalog_run_id
        except Exception as arch_e:
            print(f"[WARN] Results catalog archive failed: {arch_e}")

    except Exception as exc:
        status["status"] = "error"
        status["error"] = str(exc)
        print(f"\n[ERROR] Pipeline failed at step {status['current_step']}: {exc}")
        import traceback

        traceback.print_exc()

    status_path = out_dir / "pipeline_status.json"
    with open(status_path, "w", encoding="utf-8") as f:
        json.dump(status, f, indent=2)

    return status


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="MMSage x Chain Novelty Pipeline")
    parser.add_argument("--coordinates-dir", default=None)
    parser.add_argument("--gml", default=None)
    parser.add_argument("--output-dir", default=None)
    parser.add_argument("--top-n", type=int, default=200)
    parser.add_argument("--disease", default="IBD", help="Disease context for PubMed queries")
    parser.add_argument("--max-depth", type=int, default=3, help="Max KG path depth")
    parser.add_argument("--bacteria-filter", default=None, help="Only process files containing this bacteria name")
    parser.add_argument("--coordinates-file", default=None, help="Process a single coordinates CSV file instead of directory")
    args = parser.parse_args()

    run_pipeline(
        coordinates_dir=args.coordinates_dir,
        gml_path=args.gml,
        output_dir=args.output_dir,
        top_n=args.top_n,
        disease=args.disease,
        max_depth=args.max_depth,
        bacteria_filter=args.bacteria_filter,
        coordinates_file=args.coordinates_file,
    )
