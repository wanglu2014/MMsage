from __future__ import annotations
import csv
import os
import shutil
import subprocess
from pathlib import Path
from runtime_state import job_output_dir, job_status_path, save_json
import json
import zipfile
import pandas as pd
from scipy.stats import spearmanr
from reference_inputs import stage_inputs


def save_status(path, data):
    save_json(path, data)


def run_correlation_job(job_id: str, microbe_path: str, metabolite_path: str,
                        taxonomy_path: str | None = None, correlation_path: str | None = None,
                        pvalues_path: str | None = None, microbe_filter: str | None = None,
                        p_threshold: float = 0.05, top_n: int = 50, configuration: dict | None = None) -> None:
    status_path = job_status_path(job_id)
    out_dir = job_output_dir(job_id).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)
    scripts = Path(__file__).resolve().parents[1] / "r_pipeline"
    def status(state, progress, **fields):
        save_status(status_path, dict(job_id=job_id, status=state, progress=progress, **fields))
    try:
        status("running", "Validating input tables")
        work_dir = out_dir.parent / "analysis"
        bundled_taxonomy = Path(__file__).resolve().parents[1] / "data" / "microbe_taxonomy.tsv"
        stage_inputs(work_dir, microbe_path, metabolite_path, taxonomy_path, correlation_path, pvalues_path,
                     taxonomy_source=bundled_taxonomy)
        env = os.environ.copy()
        env["MMSAGE_JOB_DIR"] = str(work_dir)
        env["MMSAGE_OUTPUT_DIR"] = str(out_dir)
        env["MMSAGE_ROOT_P_THRESHOLD"] = str(p_threshold)
        env["MMSAGE_MIN_METABOLITE_MEAN"] = str(os.environ.get("MMSAGE_MIN_METABOLITE_MEAN", "0.0026"))
        if microbe_filter:
            env["MMSAGE_MICROBE_FILTER"] = microbe_filter
        for key, env_name in (("num_dim", "MMSAGE_NUM_DIM"), ("neighbors", "MMSAGE_UMAP_NEIGHBORS"),
                              ("min_dist", "MMSAGE_UMAP_MIN_DIST"), ("metric", "MMSAGE_UMAP_METRIC"),
                              ("cluster_res", "MMSAGE_CLUSTER_RESOLUTION")):
            if configuration and configuration.get(key) is not None:
                env[env_name] = str(configuration[key])
        env.setdefault("MMSAGE_WORKERS", "1")
        rscript = os.environ.get("MMSAGE_RSCRIPT") or shutil.which("Rscript")
        if not rscript:
            raise RuntimeError("Rscript not found. Install R or set MMSAGE_RSCRIPT.")
        status("running", "Checking R dependencies")
        env["MMSAGE_CORE_PATH"] = str(scripts / "trajectory_core" / "runtime.R")
        env["TREEMM_R_DIR"] = str(Path(__file__).resolve().parents[1] / "TreeMM" / "R")
        check = subprocess.run([rscript, "-e", 'source(Sys.getenv("MMSAGE_CORE_PATH"))'], capture_output=True, text=True, errors="replace", env=env)
        if check.returncode:
            raise RuntimeError(check.stderr.strip())
        stages = [("M1_targets_clr.R", "Resolving target CAGs and applying CLR"),
                  ("M2_root_select.R", "Selecting roots at P < 0.05 (NArank 1, 2, 3, 5)"),
                  ("M3_pluscombno1.R", "Building interaction matrices"),
                  ("M4_trajectory.R", "Computing configured trajectory configurations per CAG")]
        for script, label in stages:
            status("running", label)
            log = out_dir / (script + ".log")
            with log.open("w", encoding="utf-8") as stream:
                process = subprocess.run([rscript, str(scripts / script)], stdout=stream, stderr=subprocess.STDOUT, env=env)
            if process.returncode:
                raise RuntimeError(script + ": " + log.read_text(encoding="utf-8", errors="replace")[-3000:])
        # Always provide a compact, directly usable ranking artifact.  This is
        # independent of the optional trajectory calculation and works for any
        # feature identifiers supplied by the user.
        mic = pd.read_csv(microbe_path, sep=None, engine="python", index_col=0)
        met = pd.read_csv(metabolite_path, sep=None, engine="python", index_col=0)
        shared = [s for s in mic.columns if s in met.columns]
        ranking_rows = []
        selected_microbes = list(mic.index)
        if microbe_filter:
            query = str(microbe_filter).strip().casefold()
            selected_microbes = [name for name in mic.index if query in str(name).casefold()]
            if not selected_microbes:
                raise ValueError(f"Requested microbe or taxon is absent from the input: {microbe_filter}")
        for microbe in selected_microbes:
            for metabolite in met.index:
                rho, p = spearmanr(mic.loc[microbe, shared], met.loc[metabolite, shared], nan_policy="omit")
                if pd.notna(rho):
                    ranking_rows.append({"Microbe": microbe, "Metabolite": metabolite,
                                 "Correlation": float(rho), "P_value": float(p) if pd.notna(p) else None})
        ranking = pd.DataFrame(ranking_rows)
        if not ranking.empty:
            ranking["Significant"] = ranking["P_value"] < float(p_threshold)
            ranking["Rank"] = ranking.groupby("Microbe")["Correlation"].rank(method="first", ascending=False).astype(int)
            ranking = ranking.sort_values(["Microbe", "Correlation"], ascending=[True, False]).groupby("Microbe", group_keys=False).head(int(top_n))
        ranking.to_csv(out_dir / "correlation_ranking.csv", index=False)
        coordinates = sorted((out_dir / "coords").glob("coords_*.csv"))
        legacy_coordinates = sorted((out_dir / "coords").glob("*_newseed_*_coordinates.csv"))
        if not coordinates:
            raise RuntimeError("No trajectory coordinate files were produced")
        # Materialize a stable ranking artifact from the TreeMM trajectory
        # itself.  Lower pseudotime is the first ranked metabolite.
        trajectory_rows = []
        for path in coordinates:
            frame = pd.read_csv(path)
            if {"metabolite", "Pseudotime"}.issubset(frame.columns):
                frame = frame.sort_values("Pseudotime", kind="stable").copy()
                frame["Rank"] = range(1, len(frame) + 1)
                frame["Configuration"] = path.stem
                trajectory_rows.append(frame[["metabolite", "Pseudotime", "Rank", "Configuration"]])
        if trajectory_rows:
            pd.concat(trajectory_rows, ignore_index=True).to_csv(out_dir / "trajectory_ranking.csv", index=False)
        rows = 0
        for path in coordinates:
            with path.open(encoding="utf-8", newline="") as stream:
                rows += sum(1 for _ in csv.DictReader(stream))
        with zipfile.ZipFile(out_dir / "trajectory_results.zip", "w", zipfile.ZIP_DEFLATED) as bundle:
            for path in [out_dir / "correlation_ranking.csv", out_dir / "trajectory_ranking.csv"] + coordinates + legacy_coordinates + sorted((out_dir / "roots").glob("*.csv")) + sorted((out_dir / "logs").glob("*")) + sorted(out_dir.glob("*.R.log")):
                if path.is_file():
                    bundle.write(path, path.relative_to(out_dir).as_posix())
        status("completed", "Complete", rows=rows, coordinate_files=len(coordinates), ranking_file="correlation_ranking.csv", download=f"/api/correlation/jobs/{job_id}/download")
    except Exception as exc:
        status("error", "Failed", error=str(exc))
