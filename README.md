# MMSage Platform

MMSage Platform is a runnable release of the MMSage mechanism-discovery system. It connects rank-based microbe-metabolite signals with literature-derived knowledge graphs, PubMed co-occurrence search, multi-agent evidence aggregation, and dual-axis candidate prioritization.

The goal is to help identify microbe-metabolite-disease mechanisms that are biologically plausible, literature-grounded, and experimentally actionable.

This directory is not a placeholder. It contains a FastAPI backend, static web frontend, knowledge-graph construction scripts, pipeline modules, test code, local example data, and generated outputs.

## What the platform does

Given a MMSage-ranked microbe-metabolite pair, MMSage asks whether the pair may participate in a disease-relevant mechanism chain, and whether that chain is already known, weakly supported, or potentially novel dark matter.

The platform produces:

- MMSage candidate lists extracted from coordinate CSV files.
- Chain Novelty scores based on knowledge-graph paths and PubMed co-occurrence counts.
- Multi-agent evidence summaries for microbe-disease, microbe-metabolite, and metabolite-disease links.
- Dual-axis quadrant assignments using MMSage signal and Chain Novelty.
- A local web dashboard for candidate review, evidence inspection, and knowledge-graph exploration.

## Repository layout

```text
0810_MMSage_Platform/
├── backend/
│   ├── api_server.py                  # FastAPI server and static frontend entry point
│   ├── correlation_job_runner.py      # Uploaded abundance matrices -> R M1–M4 -> ranking CSV/ZIP
│   ├── reference_inputs.py            # Validate/align matrices and derive Spearman reference tables
│   ├── run_pipeline.py                # Full Step 1 -> KG build -> Step 2 -> Step 2b -> Step 3 pipeline
│   ├── step1_mmsage_signal.py         # Extract candidate pairs and Rank values from MMSage coordinates
│   ├── step2_chain_novelty.py         # Chain Novelty scoring using KG paths and PubMed co-occurrence
│   ├── step2_scoresp_energy.py        # Optional ScoreSP energy scoring module
│   ├── step3_quadrant.py              # Dual-axis quadrant assignment
│   ├── step4_causal_chain.py          # Legacy optional module; see compatibility note below
│   ├── validation_planner.py          # Evidence-driven validation planning (new)
│   ├── protocol_refiner.py            # Multi-role protocol review and synthesis (new)
│   ├── protocols_io_tool.py           # Read-only protocols.io retrieval (new)
│   ├── export_catalog.py              # Result catalogue and export (new)
│   ├── job_runner.py                  # Background job execution, file-based state (new)
│   ├── runtime_state.py               # Runtime state helpers (new)
│   ├── build_kg.py                    # PubMed + LLM relation extraction into GML knowledge graphs
│   ├── knowledge_pump.py              # KG loading, path search, and context construction
│   ├── knowledge_formatter.py         # Evidence formatting utilities
│   ├── llm_reasoning.py               # LLM reasoning interface
│   ├── agents/                        # Multi-agent literature evidence modules
│   ├── db_checkers/                   # Disbiome, CTD, ChEBI, BacDive, SapBERT, and related checkers
│   └── requirements.txt               # One line: `-r ../requirements.txt`
├── frontend/
│   ├── MMSage_Dashboard.html          # Main dashboard
│   ├── dashboard.compiled.js          # Dashboard logic, A4/150 dpi render sizes (new)
│   ├── index.html                     # Entry page
│   ├── pipeline.html                  # Abundance upload, M1–M4 progress and result downloads
│   ├── evidence.html                  # Agent evidence browser
│   ├── graph.html                     # Knowledge-graph explorer
│   ├── candidate.html                 # Single-candidate detail page
│   ├── browse.html                    # Result browser (new)
│   ├── browse_record.html             # Single-record view (new)
│   └── kg_evidence.js                 # Knowledge-graph evidence view (new)
├── r_pipeline/
│   ├── M1_targets_clr.R               # Resolve target CAGs/metabolites and create the CLR cache
│   ├── M2_root_select.R               # Select Noback metabolite roots from Spearman matrices
│   ├── M3_pluscombno1.R               # Calculate sample-level CAG-metabolite interaction scores
│   ├── M4_trajectory.R                # Run trajectory and pseudotime analysis
│   └── trajectory_core/               # Bundled matrix trajectory runtime and its dependencies
├── data/
│   ├── sample_coordinates/            # Example coordinate input
│   ├── knowledge_graph/               # GML knowledge graphs
│   └── databases/                     # Empty; populate it yourself (see Known limitations)
├── outputs/                           # Generated pipeline outputs
├── tests/test_pipeline.py             # Backend unit and regression tests
└── requirements.txt                   # Python dependencies
```

`backend/db_checkers/sapbert_index/` is not shipped (769 MB); rebuild it with `python backend/db_checkers/build_sapbert_index.py` when SapBERT-based post-retrieval alignment is needed.

## R preprocessing and trajectory pipeline

The `/pipeline` page accepts two abundance matrices and starts an independent background job through `backend/correlation_job_runner.py`. Python validates and aligns the tables, computes Spearman reference matrices, and invokes the four R scripts sequentially using `subprocess` and `Rscript`:

```text
Two CSV/TSV abundance tables -> reference_inputs.py -> M1 -> M2 -> M3 -> M4
                                                           -> ranking CSVs and ZIP
```

This workflow is separate from the downstream literature/KG workflow in `backend/run_pipeline.py`.

| Stage | Script | Purpose | Main outputs under the job output directory |
|---|---|---|---|
| M1 | `M1_targets_clr.R` | Resolve matching input microbe IDs; normalize within samples and apply RCLR; filter metabolites by mean normalized abundance (default `0.0026`). | `M1_cache.rds` |
| M2 | `M2_root_select.R` | Select positive Spearman edges below the requested P threshold, rank eligible metabolites by mean CLR abundance, and form root sets for NArank 1, 2, 3 and 5. If no edge qualifies, try finite positive edges, then finite edges. | `roots/all_roots_noback.csv` |
| M3 | `M3_pluscombno1.R` | Calculate sample-level `sCor = CLR(metabolite) * CLR(microbe)` interaction matrices. | `pluscomb/pluscombno1_<CAG>.csv` |
| M4 | `M4_trajectory.R` | Run PCA/UMAP, clustering, graph learning and root-based pseudotime ordering using the bundled matrix runtime. | `coords/*.csv`, `logs/*` |

### Upload format and controls

- Use UTF-8 CSV or TSV tables with a header row, feature IDs in the first column and sample IDs in the remaining columns. Rows represent microbes or metabolites; values must be finite numeric abundances.
- Feature and sample IDs must be nonempty and unique. At least three sample IDs must be shared; the staged tables follow the metabolite table's sample order.
- `Target Microbe (Optional)` filters input row names. Leave blank for all microbes, enter a complete ID for a specific microbe, or enter `Cyanobacteria` to match IDs containing that text (for example `rplo 1 (Cyanobacteria)`). Matching is based on the input names, not an external taxonomy lookup.
- `Top-N Results per Microbe` (page default 50) limits the correlation ranking for each microbe; it does not truncate the trajectory coordinates.
- `Significance Threshold (P-value)` (page default 0.05) controls root eligibility and the correlation table's `Significant` flag. The correlation CSV retains nonsignificant rows within Top-N. The API default, when omitted, is 1.0.
- `Configuration` uses `num_dim,neighbors,min_dist,metric,cluster_res`; the page defaults to `20,3,0.5,cosine,0.1`.

### Job outputs and interpretation

Each task is stored under `data/jobs/<job_id>/`: `analysis/` holds staged input/reference tables, `outputs/` holds results, and `status.json` holds progress.

| Output under `outputs/` | Contents |
|---|---|
| `correlation_ranking.csv` | `Microbe`, `Metabolite`, Spearman `Correlation`, `P_value`, `Significant`, and descending-correlation `Rank`; calculated from uploaded abundances. |
| `trajectory_ranking.csv` | `metabolite`, `Pseudotime`, `Rank`, `Configuration`; sorted by pseudotime separately for each coordinate configuration. |
| `coords/coords_<CAG>_seed1_<config>_NArank<N>.csv` | `metabolite`, `rowcolumn`, `UMAP1`, `UMAP2`, `Pseudotime`. |
| `coords/<CAG>_newseed_1_<config>_NArank<N>_coordinates.csv` | TreeMM-style export with `Row.names`, `CellID`, `Pseudotime`, and PC columns. The exporter pads to at least 18 PC columns; unavailable components are `NA`. |
| `roots/all_roots_noback.csv` | Selected roots and their NArank, correlation, P-value and CLR mean. |
| `logs/M4_summary.log`, `*.R.log` | Per-configuration outcomes and R diagnostics. |
| `trajectory_results.zip` | Ranking tables, coordinate exports, root tables and logs. |

Correlation rank and pseudotime rank answer different questions. TreeMM-style exports describe the current run; compatible column names do not guarantee reproduction of the bundled TreeMM example's values or complete PC values. Very small inputs can use a deterministic coordinate fallback rather than the full trajectory model. If some configurations fail but coordinate files exist, M4 currently allows the job to complete; inspect `logs/M4_summary.log` for coverage before interpreting the result as complete for every requested CAG. Compatibility exports are generated on the regular trajectory path.

When a task completes, the page displays **Download results ZIP**, **Download ranking CSV**, and **Output location** with the job ID. Existing downloaded ZIPs do not change after a deployment; submit a new task to obtain the current export format.

### R environment

Use Python 3.10+ (the deployment uses 3.11) and an R installation compatible with the following packages (the deployment uses R 4.5.3). Python discovers `Rscript` on `PATH`, or uses the explicit `MMSAGE_RSCRIPT` environment variable. R packages are installed separately from `requirements.txt`.

Run these commands in R:

```r
install.packages(c(
  "BiocManager", "Matrix", "tidyr", "data.table", "assertthat", "digest",
  "dplyr", "igraph", "irlba", "leidenbase", "pbapply", "proxy", "RANN",
  "RcppAnnoy", "RcppHNSW", "RhpcBLASctl", "uwot", "openssl", "plyr", "future"
))
BiocManager::install("SingleCellExperiment", ask = FALSE, update = FALSE)
source("r_pipeline/trajectory_core/runtime.R")
```

The runtime checks its required packages on each job. Linux source installations may additionally require C/C++/Fortran compilers and system development libraries; follow the failing package's installation diagnostics. The web M1–M4 path uses the bundled `trajectory_core` and a base-R CLR implementation. The separate original `TreeMM/R/run_pipeline.R` still calls `SpiecEasi`; its dependencies should be assessed separately if running that entry point directly.

Python sets `MMSAGE_JOB_DIR`, `MMSAGE_OUTPUT_DIR`, root threshold, target filter and trajectory configuration for each R child process. For manual execution, explicitly set the job/output directories and stage the reference tables first; the scripts no longer use a hard-coded data root. Optional deployment controls include `MMSAGE_MIN_METABOLITE_MEAN` (default 0.0026) and `MMSAGE_WORKERS` (web default 1).

### Linux deployment and upload smoke test

From `/srv/mmsage`, install the Python requirements in the service's virtual environment and validate the R runtime:

```bash
cd /srv/mmsage
venv/bin/python -m pip install -r requirements.txt
Rscript -e 'source("r_pipeline/trajectory_core/runtime.R")'
# Foreground launch for a new installation; avoid starting a second server on the same port.
venv/bin/python -m uvicorn backend.api_server:app --host 127.0.0.1 --port 8000
```

On the existing systemd deployment, use `systemctl restart mmsage.service` after replacing Python code, and verify with `systemctl is-active mmsage.service`. Keep the working directory at `/srv/mmsage`, preserve write access to `data/jobs/`, and retain the existing reverse proxy for public access. Do not restart while jobs are running: background jobs are threads inside the service process.

```bash
curl -sS -F 'microbes=@TreeMM/Example/input/microbes.tsv' \
  -F 'metabolites=@TreeMM/Example/input/metabolites.tsv' \
  'http://127.0.0.1:8000/api/correlation/jobs?microbe_filter=Cyanobacteria&top_n=50&p_threshold=0.05&num_dim=20&neighbors=3&min_dist=0.5&metric=cosine&cluster_res=0.1'
# Replace <job_id> with the returned identifier and poll until terminal status.
curl 'http://127.0.0.1:8000/api/correlation/jobs/<job_id>'
# After completion:
curl -f 'http://127.0.0.1:8000/api/correlation/jobs/<job_id>/download' -o trajectory_results.zip
curl -f 'http://127.0.0.1:8000/api/correlation/jobs/<job_id>/ranking' -o correlation_ranking.csv
```

Validate ZIP contents and `M4_summary.log` as well as HTTP status. The abundance/trajectory workflow does not require an LLM API key.

## Packaged results

The package includes the MMSage example input, generated Steps 1–3 outputs, knowledge graphs, and one archived catalogue copy with the same 146 candidates and Step 1–3 core fields. The archived Step 2 copy additionally contains the separately recorded Step 2b evidence. Existing scientific values are retained as the packaged results.

| File | Current content |
|---|---:|
| `outputs/step1_candidates.json` | 146 candidate pairs |
| `outputs/step2_chain_novelty.json` | 146 Chain Novelty results |
| `outputs/step2b_agent_evidence.json` | 146 multi-agent evidence results |
| `outputs/step3_quadrant.json` | 146 quadrant results |
| `data/results_catalog/runs/run_20260602_164050_akkermansia_muciniph_ibd/candidates.json` | 146 archived Step 1–3 catalogue records |
| `data/knowledge_graph/auto_built_kg.gml` | 403 nodes, 447 edges |
| `data/knowledge_graph/direct_kg.gml` | 768 nodes, 1136 edges |
| `data/knowledge_graph/ibd_test_kg.gml` | 30 nodes, 46 edges |

The current `outputs/pipeline_status.json` reports a completed run with this quadrant distribution:

| Quadrant | Interpretation | Count |
|---|---|---:|
| I | High MMSage signal + high Chain Novelty; priority dark-matter candidates | 55 |
| II | Low MMSage signal + high Chain Novelty; novel but weaker trajectory signal | 54 |
| III | High MMSage signal + low Chain Novelty; more literature-supported relationships | 18 |
| IV | Low MMSage signal + low Chain Novelty; lower-priority candidates | 19 |

## Installation

Run commands from the project root so that relative input and output paths resolve correctly.

```powershell
cd 0810_MMSage_Platform
python -m venv .venv
.\.venv\Scripts\Activate.ps1
python -m pip install --upgrade pip
python -m pip install -r requirements.txt
```

On macOS or Linux, activate the environment with `source .venv/bin/activate`. If you already have a conda or micromamba environment, install `requirements.txt` in that environment instead.

## Start the web app

```bash
cd 0810_MMSage_Platform
python backend/api_server.py --port 8000
```

Open these local pages after the server starts:

- `http://localhost:8000/` for the main dashboard.
- `http://localhost:8000/pipeline` for abundance CSV/TSV upload and M1–M4 execution.
- `http://localhost:8000/evidence` for multi-agent evidence review.
- `http://localhost:8000/graph` for knowledge-graph exploration.
- `http://localhost:8000/docs` for FastAPI-generated API documentation.

The frontend API-key box syncs the key only to the local backend runtime. Do not commit real API keys into repository files.

## Run the full pipeline

```bash
cd 0810_MMSage_Platform
python backend/run_pipeline.py --disease IBD --top-n 200 --max-depth 3
```

Default inputs and outputs:

- Sample coordinate directory: `data/sample_coordinates/`
- Default run input: `data/sample_coordinates/pluscombno1V0317_min0_Akkermansia_muciniphila_seed_1_dim_10_neighbor_2_dist_0.4_metric_euclidean_rank_rootknow_cor0303_1_top_50_Pthre_0.1_pair.csv_clu_1_coordinates_tunek.csv`
- Default disease context: `IBD`
- Fallback KG: `data/knowledge_graph/auto_built_kg.gml`, then `data/knowledge_graph/ibd_test_kg.gml`; a successful run-specific build is written to `outputs/knowledge_graph.gml` and used for that run
- Default output directory: `outputs/`

Example with explicit options:

```bash
python backend/run_pipeline.py \
  --coordinates-file "data/sample_coordinates/pluscombno1V0317_min0_Akkermansia_muciniphila_seed_1_dim_10_neighbor_2_dist_0.4_metric_euclidean_rank_rootknow_cor0303_1_top_50_Pthre_0.1_pair.csv_clu_1_coordinates_tunek.csv" \
  --disease IBD \
  --top-n 50 \
  --max-depth 3 \
  --output-dir outputs
```

`run_pipeline.py` calls `backend/build_kg.py candidates` after Step 1 to build a targeted graph from the current candidate list. This step depends on PubMed access and LLM API access, so it can be slow or fail if the network or API is unavailable.

## Run pipeline stages manually

Use this mode when you want to reuse the existing KG and existing candidate files.

```bash
# Step 1: coordinate files -> candidate pairs
python backend/step1_mmsage_signal.py \
  --coordinates-file "data/sample_coordinates/pluscombno1V0317_min0_Akkermansia_muciniphila_seed_1_dim_10_neighbor_2_dist_0.4_metric_euclidean_rank_rootknow_cor0303_1_top_50_Pthre_0.1_pair.csv_clu_1_coordinates_tunek.csv" \
  --output outputs/step1_candidates.json \
  --top-n 200

# Step 2: Chain Novelty scoring
python backend/step2_chain_novelty.py \
  --candidates outputs/step1_candidates.json \
  --gml data/knowledge_graph/auto_built_kg.gml \
  --output outputs/step2_chain_novelty.json \
  --disease IBD \
  --max-depth 3

# Step 3: Dual-axis quadrant assignment
python backend/step3_quadrant.py \
  --input outputs/step2_chain_novelty.json \
  --output outputs/step3_quadrant.json

```

### Legacy optional Step 4 source

`backend/step4_causal_chain.py` is retained for historical reference, but it expects the legacy nested `scoresp_result.candidate` schema and is not part of the current Step 1–3 pipeline. It must not be run directly against the shipped `outputs/step3_quadrant.json` without a schema adapter. No Step 4 output is included in the package or archived by `export_catalog.py`.

## Build knowledge graphs separately

`backend/build_kg.py` supports two modes.

### Direct mode: species plus disease

```bash
python backend/build_kg.py direct \
  --species "Akkermansia muciniphila" \
  --disease IBD \
  --output data/knowledge_graph/direct_kg.gml \
  --max-articles 40 \
  --max-queries 10
```

### Candidates mode: Step 1 candidates to KG

```bash
python backend/build_kg.py candidates \
  --candidates outputs/step1_candidates.json \
  --disease IBD \
  --output data/knowledge_graph/auto_built_kg.gml \
  --max-articles 40 \
  --max-queries 100
```

Important node and edge fields:

| Object | Field | Purpose |
|---|---|---|
| node | `label` | Node identifier used by NetworkX after GML loading |
| node | `type` / `node_type` | Entity type, such as microbe, metabolite, pathway, or disease |
| node | `description` | LLM-extracted entity description |
| edge | `relation` | Relation type, such as produces, metabolizes, inhibits, or promotes |
| edge | `pmid` | PubMed IDs supporting the edge |
| edge | `edge_impact_factor` | Journal influence field used by ScoreSP |
| edge | `edge_citation_count` | Citation-count field used by ScoreSP |
| edge | `publication_year` | Publication year |

## API endpoints

Common backend endpoints:

| Method | Endpoint | Description |
|---|---|---|
| `POST` | `/api/correlation/jobs` | Upload multipart `microbes` and `metabolites` files; start M1–M4 with query parameters |
| `GET` | `/api/correlation/jobs/{job_id}` | Read job status, progress, output counts and download path |
| `GET` | `/api/correlation/jobs/{job_id}/ranking` | Download the correlation ranking CSV |
| `GET` | `/api/correlation/jobs/{job_id}/download` | Download the result ZIP after completion |
| `POST` | `/api/pipeline/run` | Start the pipeline in a background thread |
| `GET` | `/api/pipeline/status` | Get current pipeline status |
| `GET` | `/api/step1/candidates` | Get Step 1 candidates |
| `GET` | `/api/step2/chain_novelty` | Get Chain Novelty results |
| `GET` | `/api/step2b/agent_evidence` | Get multi-agent evidence results |
| `GET` | `/api/step3/quadrant` | Get quadrant results |
| `GET` | `/api/graph/stats` | Get KG size and node-type statistics |
| `GET` | `/api/graph/full` | Get the full KG in Cytoscape.js format |
| `GET` | `/api/graph/subgraph` | Extract a bacteria/metabolite-filtered subgraph |
| `GET` | `/api/graph/chain` | Extract a focused bacteria-metabolite-disease chain subgraph |
| `POST` | `/api/upload-coordinates` | Upload a coordinate CSV and start the pipeline |
| `GET` | `/api/pipeline-progress` | Poll progress after CSV upload |
| `POST` | `/api/generate-protocol` | Generate an experimental validation protocol for one candidate |

## Output files

| File | Produced by | Key fields |
|---|---|---|
| `outputs/pipeline_status.json` | `run_pipeline.py` | `status`, `steps`, `total_duration_s` |
| `outputs/step1_candidates.json` | Step 1 | `bacteria`, `metabolite`, `mmsage_norm`, `rank_in_microbe`, `source_file` |
| `outputs/step2_chain_novelty.json` | Step 2 | `chain_count`, `chain_novelty`, `has_path`, `pairwise_counts` |
| `outputs/step2b_agent_evidence.json` | Step 2b | `hop_counts`, `bottleneck`, `db_bonus`, `sources`, `recommendation` |
| `outputs/step3_quadrant.json` | Step 3 | `quadrant`, `quadrant_label`, `is_dark_matter`, `evidence_foundation` |
| `data/results_catalog/runs/*/candidates.json` | Archived Step 1–3 catalogue snapshot | Candidate, Chain Novelty, quadrant, and evidence fields |

Current candidate and stage outputs use `mmsage_norm`. Candidate-facing views display Rank, and downstream summaries use Prioritization.

## Tests and validation

Run the backend tests with:

```bash
cd 0810_MMSage_Platform
python -m pytest tests/test_pipeline.py -q
```

The suite validates the 146-candidate sample fixture, current `mmsage_norm` schema,
knowledge-graph structure, stage-to-stage candidate consistency, and Chain Novelty formula.
Required fixtures fail explicitly when missing instead of being silently skipped.

PubMed entity terms are queried exactly as provided, with underscores converted to spaces.
SapBERT is not used to expand PubMed queries; it is reserved for post-retrieval entity
alignment against controlled databases and downloaded abstracts.

For lightweight service validation:

```bash
python backend/api_server.py --port 8000
# In another terminal:
curl http://localhost:8000/api/pipeline/status
curl http://localhost:8000/api/graph/stats
curl http://localhost:8000/api/step3/quadrant
```

## Known limitations

- The full pipeline accesses PubMed and LLM APIs. Dynamic KG construction, agent evidence, and protocol generation may fail or run slowly when network or API access is unavailable.
- `run_pipeline.py` launches dynamic graph construction with a relative `backend/build_kg.py` path, so it should be run from the project root.
- PubMed entity terms are not expanded. SapBERT aliases must never be written back into search queries.
- The frontend is static HTML plus JavaScript, not an npm application. CSS and plotting libraries are loaded from CDNs.
- `backend/db_checkers/sapbert_index/` is not shipped (769 MB). Entity normalisation and SapBERT retrieval are unavailable until it is rebuilt with `python backend/db_checkers/build_sapbert_index.py`, which itself needs the reference tables below.
- `data/databases/` is present but empty. Populating it requires the five reference sources used by `backend/db_checkers/`: CTD chemical-disease associations, BacDive edges, BindeD annotations, ChEBI relationships, and Disbiome edges. Each checker module documents the columns it expects.
- No API credentials are included. `llm_reasoning.py` reads `OPENAI_API_KEY`, `OPENAI_API_BASE` and `OPENAI_MODEL` from the environment; `build_kg.py` reads a rotating DeepSeek key CSV supplied via `--keys-csv`.
