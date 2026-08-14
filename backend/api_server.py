"""
FastAPI Backend Server
======================
REST API for the MMSage x Chain Novelty Dual-Axis Decision System.
"""

import json
import asyncio
import threading
import subprocess
import re
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, List, Optional
from urllib.parse import urlencode
from uuid import uuid4

from fastapi import FastAPI, BackgroundTasks, Query, UploadFile, File, HTTPException
from fastapi.staticfiles import StaticFiles
from fastapi.responses import HTMLResponse, FileResponse, PlainTextResponse, JSONResponse
from fastapi.middleware.cors import CORSMiddleware
from pydantic import BaseModel, Field
import networkx as nx
import uvicorn

import sys

# Paths
BACKEND_DIR = Path(__file__).parent
sys.path.append(str(BACKEND_DIR))
PROJECT_DIR = BACKEND_DIR.parent
OUTPUTS_DIR = PROJECT_DIR / "outputs"
DATA_DIR = PROJECT_DIR / "data"
CATALOG_DIR = DATA_DIR / "results_catalog"
FRONTEND_DIR = PROJECT_DIR / "frontend"
# Prefer auto-built KG, fallback to test KG
AUTO_GML = DATA_DIR / "knowledge_graph" / "auto_built_kg.gml"
TEST_GML = DATA_DIR / "knowledge_graph" / "ibd_test_kg.gml"
DEFAULT_GML = AUTO_GML if AUTO_GML.exists() else TEST_GML
DEFAULT_COORDINATES_DIR = str(PROJECT_DIR / "data" / "sample_coordinates")
DEFAULT_COORDINATES_FILE = str(PROJECT_DIR / "data" / "sample_coordinates" / "pluscombno1V0317_min0_Akkermansia_muciniphila_seed_1_dim_10_neighbor_2_dist_0.4_metric_euclidean_rank_rootknow_cor0303_1_top_50_Pthre_0.1_pair.csv_clu_1_coordinates_tunek.csv")

from runtime_state import (
    APP_STATE_PATH,
    JOBS_DIR,
    job_dir,
    job_inputs_dir,
    job_output_dir,
    job_status_path,
    load_app_state,
    load_json as load_runtime_json,
    save_app_state,
    save_json as save_runtime_json,
)

app = FastAPI(title="MMSage x Chain Novelty Dual-Axis System", version="2.0.0")

# CORS for local development
app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_methods=["*"],
    allow_headers=["*"],
)

# Pipeline state
pipeline_state = {
    "status": "idle",
    "current_step": 0,
    "progress": "",
    "error": None,
}


# --- Models ---

class PipelineConfig(BaseModel):
    coordinates_dir: Optional[str] = None
    coordinates_file: Optional[str] = None
    gml_path: Optional[str] = None
    output_dir: Optional[str] = None
    top_n: int = 200
    disease: str = "IBD"
    max_depth: int = 3
    bacteria_filter: Optional[str] = None


class KGRequest(BaseModel):
    bacteria: str = ""
    metabolite: str = ""
    include_fella_network: bool = False


class ProtocolRequest(BaseModel):
    bacteria: str = ""
    metabolite: str = ""
    mechanism_summary: str = ""
    disease: str = "IBD"
    run_id: Optional[str] = None
    research_question: str = ""
    prompt_constraints: str = ""


class ValidationPlanRequest(BaseModel):
    bacteria: str = ""
    metabolite: str = ""
    mechanism_summary: str = ""
    disease: str = "IBD"
    run_id: Optional[str] = None
    mode: str = "evidence_self_reflection"
    protocol_text: str = ""
    research_question: str = ""
    prompt_constraints: str = ""


class ProtocolRefineStartRequest(BaseModel):
    bacteria: str = ""
    metabolite: str = ""
    mechanism_summary: str = ""
    disease: str = "IBD"
    run_id: Optional[str] = None
    protocol_text: str = ""


class ProtocolRefineConfirmRequest(BaseModel):
    session_id: str
    modules: List[Dict[str, Any]] = Field(default_factory=list)
    extra_requirements: str = ""


class ProtocolRefineRunRequest(BaseModel):
    session_id: str


# --- Helpers ---

def read_json(filename: str) -> list:
    """Read JSON output file."""
    path = OUTPUTS_DIR / filename
    if not path.exists():
        return []
    with open(path, 'r', encoding='utf-8') as f:
        return json.load(f)


def read_json_for_run(filename: str, run_id: str = "") -> list:
    resolved_run_id = resolve_requested_run_id(run_id)
    if resolved_run_id:
        try:
            from export_catalog import run_dir_path

            path = run_dir_path(resolved_run_id) / filename
            if path.exists():
                with open(path, "r", encoding="utf-8") as f:
                    return json.load(f)
        except Exception:
            pass
    return read_json(filename)


_OBSERVATIONAL_PROMOTE_MARKERS = (
    "upregulated",
    "increased",
    "increase in",
    "enriched",
    "abundance",
    "correlated",
    "correlation",
    "associated",
    "dysbiosis",
    "dss mice",
    "in dss",
    "tnbs",
    "mouse model",
    "mice",
    "patients",
    "patient",
)


def _normalize_edge_relation(
    G: nx.DiGraph, src: str, tgt: str, edata: dict
) -> tuple[str, bool]:
    """Downgrade suspicious observational microbe->IBD promotes edges."""
    relation = str(edata.get("relation", "") or "")
    if relation.lower() != "promotes":
        return relation, False

    src_type = str(G.nodes[src].get("type", "") or "").lower()
    tgt_type = str(G.nodes[tgt].get("type", "") or "").lower()
    if src_type != "microbe" or tgt_type != "disease":
        return relation, False

    target_text = " ".join([
        str(G.nodes[tgt].get("label", "") or "").lower().replace("_", " "),
        str(G.nodes[tgt].get("description", "") or "").lower(),
    ])
    if not any(term in target_text for term in (
        "ibd", "inflammatory bowel", "colitis", "enteritis", "intestinal inflammation"
    )):
        return relation, False

    desc = str(edata.get("description", "") or "").lower()
    if any(marker in desc for marker in _OBSERVATIONAL_PROMOTE_MARKERS):
        return "associated_with", True

    return relation, False


def _load_pipeline_status_from_outputs() -> dict:
    path = OUTPUTS_DIR / "pipeline_status.json"
    if not path.exists():
        return {}
    with open(path, "r", encoding="utf-8") as f:
        data = json.load(f)
    return data if isinstance(data, dict) else {}


def ensure_default_dashboard_run_id() -> str:
    state = load_app_state()
    run_id = state.get("default_dashboard_run_id", "")
    if run_id:
        try:
            from export_catalog import run_dir_path

            if run_dir_path(run_id).is_dir():
                return run_id
        except Exception:
            pass

    try:
        from export_catalog import archive_run, list_runs

        if (OUTPUTS_DIR / "step3_quadrant.json").exists():
            status = _load_pipeline_status_from_outputs()
            disease = status.get("disease") or "IBD"
            run_id = archive_run(output_dir=str(OUTPUTS_DIR), disease=disease)
            state["default_dashboard_run_id"] = run_id
            save_app_state(state)
            return run_id

        runs = list_runs()
        if runs:
            run_id = runs[0].get("run_id", "")
            if run_id:
                state["default_dashboard_run_id"] = run_id
                save_app_state(state)
                return run_id
    except Exception:
        pass

    return ""


def resolve_requested_run_id(run_id: str = "") -> str:
    return run_id or ensure_default_dashboard_run_id()


def get_run_meta(run_id: str = "") -> dict:
    resolved_run_id = resolve_requested_run_id(run_id)
    if resolved_run_id:
        try:
            from export_catalog import run_dir_path

            path = run_dir_path(resolved_run_id) / "meta.json"
            if path.exists():
                with open(path, "r", encoding="utf-8") as f:
                    data = json.load(f)
                if isinstance(data, dict):
                    return data
        except Exception:
            pass
    return {}


def resolve_gml_path(run_id: str = "") -> Path:
    resolved_run_id = resolve_requested_run_id(run_id)
    if resolved_run_id:
        try:
            from export_catalog import run_dir_path

            run_dir = run_dir_path(resolved_run_id)
            for name in ("knowledge_graph.gml", "auto_built_kg.gml"):
                path = run_dir / name
                if path.exists():
                    return path
        except Exception:
            pass
    return DEFAULT_GML


def load_job_status(job_id: str) -> dict:
    data = load_runtime_json(job_status_path(job_id))
    return data if isinstance(data, dict) else {"job_id": job_id, "status": "not_found"}


def _now_iso() -> str:
    return datetime.now(timezone.utc).isoformat()


def load_validation_plan_job_status(job_id: str) -> dict:
    data = load_runtime_json(job_status_path(job_id))
    return data if isinstance(data, dict) else {"job_id": job_id, "status": "not_found"}


def save_validation_plan_job_status(job_id: str, data: dict) -> None:
    save_runtime_json(job_status_path(job_id), data)


def create_validation_plan_job_status(job_id: str, req: "ValidationPlanRequest") -> dict:
    mode = str(req.mode or "evidence_self_reflection").strip() or "evidence_self_reflection"
    queued_message = (
        "Standalone question-driven validation protocol queued."
        if mode == "question_driven"
        else "Step 3 validation protocol queued."
    )
    payload = {
        "job_id": job_id,
        "job_type": "validation_plan",
        "mode": mode,
        "status": "queued",
        "progress_percent": 0,
        "current_stage": "queued",
        "current_message": queued_message,
        "created_at": _now_iso(),
        "updated_at": _now_iso(),
        "candidate": {
            "bacteria": req.bacteria,
            "metabolite": req.metabolite,
            "disease": req.disease,
            "run_id": resolve_requested_run_id(req.run_id),
        },
        "user_brief": {
            "research_question": req.research_question,
            "prompt_constraints": req.prompt_constraints,
        },
        "logs": [],
        "result": None,
        "error": "",
    }
    save_validation_plan_job_status(job_id, payload)
    return payload


def append_validation_plan_job_log(
    job_id: str,
    stage: str,
    message: str,
    progress_percent: int,
    extra: Optional[dict] = None,
) -> None:
    print(f"[validation-plan:{job_id}] {progress_percent}% {stage}: {message}")
    status = load_validation_plan_job_status(job_id)
    logs = status.get("logs") if isinstance(status.get("logs"), list) else []
    logs.append(
        {
            "timestamp": _now_iso(),
            "stage": stage,
            "message": message,
            "progress_percent": progress_percent,
            "extra": extra or {},
        }
    )
    status["logs"] = logs[-80:]
    status["status"] = "running" if stage != "completed" else "completed"
    status["progress_percent"] = progress_percent
    status["current_stage"] = stage
    status["current_message"] = message
    status["updated_at"] = _now_iso()
    save_validation_plan_job_status(job_id, status)


def complete_validation_plan_job(job_id: str, result: dict) -> None:
    print(f"[validation-plan:{job_id}] 100% completed")
    status = load_validation_plan_job_status(job_id)
    status["status"] = "completed"
    status["progress_percent"] = 100
    status["current_stage"] = "completed"
    status["current_message"] = "Validation plan completed."
    status["result"] = result
    status["updated_at"] = _now_iso()
    save_validation_plan_job_status(job_id, status)


def fail_validation_plan_job(job_id: str, error: str) -> None:
    print(f"[validation-plan:{job_id}] error: {error}")
    status = load_validation_plan_job_status(job_id)
    logs = status.get("logs") if isinstance(status.get("logs"), list) else []
    logs.append(
        {
            "timestamp": _now_iso(),
            "stage": "error",
            "message": error,
            "progress_percent": status.get("progress_percent", 0),
            "extra": {},
        }
    )
    status["logs"] = logs[-80:]
    status["status"] = "error"
    status["error"] = error
    status["current_stage"] = "error"
    status["current_message"] = error
    status["updated_at"] = _now_iso()
    save_validation_plan_job_status(job_id, status)


def run_validation_plan_job(job_id: str, req: "ValidationPlanRequest") -> None:
    try:
        from validation_planner import generate_validation_plan as build_validation_plan

        start_message = (
            "Standalone question-driven validation worker started."
            if str(req.mode or "").strip() == "question_driven"
            else "Step 3 validation protocol worker started."
        )
        append_validation_plan_job_log(job_id, "started", start_message, 2)

        def progress_callback(stage: str, message: str, percent: int, extra: Optional[dict] = None) -> None:
            append_validation_plan_job_log(job_id, stage, message, percent, extra)

        result = build_validation_plan(
            bacteria=req.bacteria,
            metabolite=req.metabolite,
            disease=req.disease,
            run_id=resolve_requested_run_id(req.run_id),
            mechanism_summary=req.mechanism_summary,
            mode=req.mode,
            protocol_text=req.protocol_text,
            research_question=req.research_question,
            prompt_constraints=req.prompt_constraints,
            progress_callback=progress_callback,
        )
        complete_validation_plan_job(job_id, result)
    except Exception as e:
        current_stage = load_validation_plan_job_status(job_id).get("current_stage") or "unknown"
        fail_validation_plan_job(
            job_id,
            f"Validation plan generation failed during '{current_stage}': {e}",
        )


def run_pipeline_background(config: PipelineConfig):
    """Run pipeline in background thread."""
    global pipeline_state
    pipeline_state = {"status": "running", "current_step": 1, "progress": "Starting...", "error": None}

    try:
        from run_pipeline import run_pipeline
        result = run_pipeline(
            coordinates_dir=config.coordinates_dir,
            gml_path=config.gml_path,
            output_dir=config.output_dir,
            top_n=config.top_n,
            disease=config.disease,
            max_depth=config.max_depth,
            bacteria_filter=config.bacteria_filter,
            coordinates_file=config.coordinates_file,
        )
        pipeline_state["status"] = result.get("status", "completed")
        pipeline_state["current_step"] = 3
        pipeline_state["progress"] = "Pipeline completed"
    except Exception as e:
        pipeline_state["status"] = "error"
        pipeline_state["error"] = str(e)


# --- Pipeline Endpoints ---

@app.post("/api/pipeline/run")
async def trigger_pipeline(config: PipelineConfig = PipelineConfig()):
    """Trigger pipeline execution in background."""
    if pipeline_state["status"] == "running":
        return {"status": "already_running", "current_step": pipeline_state["current_step"]}

    thread = threading.Thread(target=run_pipeline_background, args=(config,), daemon=True)
    thread.start()
    return {"status": "started"}


@app.get("/api/pipeline/status")
async def get_pipeline_status():
    """Get current pipeline status."""
    # Also check if output files exist
    steps_completed = []
    for i, fname in enumerate(["step1_candidates.json", "step2_chain_novelty.json",
                                "step3_quadrant.json"], 1):
        if (OUTPUTS_DIR / fname).exists():
            steps_completed.append(i)
    if (OUTPUTS_DIR / "step2b_agent_evidence.json").exists():
        steps_completed.append("2b")

    return {
        **pipeline_state,
        "steps_completed": steps_completed,
        "output_dir": str(OUTPUTS_DIR),
    }


# --- Results Endpoints ---

@app.get("/api/step1/candidates")
async def get_step1_candidates(run_id: str = Query(default="")):
    """Get step 1 candidates."""
    return read_json_for_run("step1_candidates.json", run_id)


@app.get("/api/step2/chain_novelty")
async def get_step2_chain_novelty(run_id: str = Query(default="")):
    """Get step 2 Chain Novelty results."""
    return read_json_for_run("step2_chain_novelty.json", run_id)


@app.get("/api/step3/quadrant")
async def get_step3_quadrant(run_id: str = Query(default="")):
    """Get step 3 quadrant results."""
    return read_json_for_run("step3_quadrant.json", run_id)


# --- Multi-Agent Evidence Endpoints ---

@app.get("/api/step2b/agent_evidence")
async def get_agent_evidence(run_id: str = Query(default="")):
    """Get step 2b multi-agent evidence results."""
    return read_json_for_run("step2b_agent_evidence.json", run_id)


@app.get("/api/agent/query")
async def query_agent_evidence(
    bacteria: str = Query(default=""),
    metabolite: str = Query(default=""),
    disease: str = Query(default="IBD"),
):
    """Run multi-agent evidence query for a single candidate (on-demand)."""
    if not bacteria or not metabolite:
        return {"error": "bacteria and metabolite are required"}

    from agents.master_agent import MasterAgent
    master = MasterAgent(c_max=500, parallel=True)
    result = master.run(bacteria, metabolite, disease)
    return result


# --- Knowledge Graph Endpoints ---

@app.get("/api/graph/stats")
async def get_graph_stats(run_id: str = Query(default="")):
    """Get knowledge graph statistics."""
    gml_path = resolve_gml_path(run_id)
    if not gml_path.exists():
        return {"error": "GML file not found"}

    G = nx.read_gml(str(gml_path))
    node_types = {}
    for _, d in G.nodes(data=True):
        t = d.get('type', 'unknown')
        node_types[t] = node_types.get(t, 0) + 1

    return {
        "nodes": G.number_of_nodes(),
        "edges": G.number_of_edges(),
        "node_types": node_types,
        "gml_path": str(gml_path),
    }


@app.get("/api/graph/subgraph")
async def get_graph_subgraph(
    bacteria: str = Query(default=""),
    metabolite: str = Query(default=""),
    max_hop: int = Query(default=2),
    run_id: str = Query(default=""),
):
    """Get KG subgraph in Cytoscape.js JSON format."""
    gml_path = resolve_gml_path(run_id)
    if not gml_path.exists():
        return {"error": "GML not found"}

    G = nx.read_gml(str(gml_path))

    # If no filter, return full graph
    if not bacteria and not metabolite:
        subgraph = G
    else:
        # Find matching nodes by label
        seeds = set()
        for nid, ndata in G.nodes(data=True):
            label = ndata.get('label', nid).lower()
            nid_lower = nid.lower()
            if bacteria and (bacteria.lower() in label or bacteria.lower() in nid_lower):
                seeds.add(nid)
            if metabolite and (metabolite.lower() in label or metabolite.lower() in nid_lower):
                seeds.add(nid)

        if seeds:
            # Extract k-hop subgraph
            neighbors = set(seeds)
            for _ in range(max_hop):
                new_neighbors = set()
                for n in neighbors:
                    new_neighbors.update(G.predecessors(n))
                    new_neighbors.update(G.successors(n))
                neighbors.update(new_neighbors)
            subgraph = G.subgraph(neighbors)
        else:
            subgraph = G

    # Convert to Cytoscape.js JSON
    cy_nodes = []
    for nid, ndata in subgraph.nodes(data=True):
        cy_nodes.append({
            "data": {
                "id": nid,
                "label": ndata.get('label', nid),
                "type": ndata.get('type', 'unknown'),
                **{k: v for k, v in ndata.items() if k not in ('label', 'type')},
            }
        })

    cy_edges = []
    for src, tgt, edata in subgraph.edges(data=True):
        relation, downgraded = _normalize_edge_relation(subgraph, src, tgt, edata)
        cy_edges.append({
            "data": {
                "source": src,
                "target": tgt,
                "relation": relation,
                "relation_original": edata.get('relation', ''),
                "relation_downgraded": downgraded,
                "description": edata.get('description', ''),
                "impact_factor": edata.get('edge_impact_factor', 0),
                "citation_count": edata.get('edge_citation_count', 0),
                "pmid": edata.get('pmid', ''),
            }
        })

    return {"nodes": cy_nodes, "edges": cy_edges}


@app.get("/api/graph/full")
async def get_full_graph(run_id: str = Query(default="")):
    """Get full KG as Cytoscape.js JSON."""
    return await get_graph_subgraph(bacteria="", metabolite="", max_hop=2, run_id=run_id)


@app.get("/api/graph/chain")
async def get_chain_subgraph(
    bacteria: str = Query(...),
    metabolite: str = Query(...),
    disease: str = Query(default="IBD"),
    max_neighbors: int = Query(default=6),
    focused: bool = Query(default=False),
    run_id: str = Query(default=""),
):
    """Get a focused single-chain subgraph: bacteria -> metabolite -> disease.

    Only extracts 1-hop neighbors of the three core nodes, limited to
    `max_neighbors` per core node to keep the graph clean (10-20 nodes total).
    Core nodes and chain edges are marked for frontend highlighting.

    If focused=true, only core nodes and direct edges between them are returned
    (avoids repeated neighbor clutter for the same microbe across catalog records).
    """
    gml_path = resolve_gml_path(run_id)
    if not gml_path.exists():
        return {"error": "GML not found"}

    G = nx.read_gml(str(gml_path))

    # --- Find core nodes by fuzzy label match ---
    def find_node(keyword):
        keyword_lower = keyword.lower().replace("_", " ")
        best, best_score = None, 0
        for nid, ndata in G.nodes(data=True):
            label = ndata.get('label', nid).lower().replace("_", " ")
            nid_lower = nid.lower().replace("_", " ")
            # Exact match
            if keyword_lower == label or keyword_lower == nid_lower:
                return nid
            # Substring match — prefer shorter labels (more specific)
            if keyword_lower in label or keyword_lower in nid_lower:
                score = len(keyword_lower) / max(len(label), 1)
                if score > best_score:
                    best, best_score = nid, score
        return best

    core_bacteria = find_node(bacteria)
    core_metabolite = find_node(metabolite)
    core_disease = find_node(disease)

    core_ids = set()
    core_labels = {}
    for role, nid in [("bacteria", core_bacteria), ("metabolite", core_metabolite), ("disease", core_disease)]:
        if nid:
            core_ids.add(nid)
            core_labels[nid] = role

    if not core_ids:
        return {"error": f"No matching nodes found for: {bacteria}, {metabolite}, {disease}",
                "nodes": [], "edges": []}

    core_list = [core_bacteria, core_metabolite, core_disease]
    selected_nodes = set(core_ids)

    if focused:
        for i in range(len(core_list)):
            for j in range(len(core_list)):
                if i != j and core_list[i] and core_list[j] and G.has_edge(core_list[i], core_list[j]):
                    selected_nodes.add(core_list[i])
                    selected_nodes.add(core_list[j])
    else:
        pass  # neighbor expansion below

    # --- Collect 1-hop neighbors per core node (limited) ---
    if not focused:
        for core_nid in core_ids:
            neighbors = set()
            neighbors.update(G.predecessors(core_nid))
            neighbors.update(G.successors(core_nid))
            neighbors -= core_ids

            filtered_neighbors = []
            for n in neighbors:
                node_type = G.nodes[n].get('type', 'unknown')
                if node_type == 'disease' and n != core_disease:
                    continue
                filtered_neighbors.append(n)

            ranked = sorted(filtered_neighbors, key=lambda n: G.degree(n), reverse=True)
            selected_nodes.update(ranked[:max_neighbors])

    subgraph = G.subgraph(selected_nodes)

    # --- Identify chain edges (direct links between core nodes) ---
    chain_edges = set()
    for i in range(len(core_list)):
        for j in range(len(core_list)):
            if i != j and core_list[i] and core_list[j]:
                if subgraph.has_edge(core_list[i], core_list[j]):
                    chain_edges.add((core_list[i], core_list[j]))

    # --- Build Cytoscape.js JSON with core/chain markers ---
    cy_nodes = []
    for nid, ndata in subgraph.nodes(data=True):
        is_core = nid in core_ids
        cy_nodes.append({
            "data": {
                "id": nid,
                "label": ndata.get('label', nid),
                "type": ndata.get('type', 'unknown'),
                "is_core": is_core,
                "core_role": core_labels.get(nid, ""),
                **{k: v for k, v in ndata.items() if k not in ('label', 'type')},
            }
        })

    cy_edges = []
    for src, tgt, edata in subgraph.edges(data=True):
        is_chain = (src, tgt) in chain_edges
        relation, downgraded = _normalize_edge_relation(subgraph, src, tgt, edata)
        cy_edges.append({
            "data": {
                "source": src,
                "target": tgt,
                "relation": relation,
                "relation_original": edata.get('relation', ''),
                "relation_downgraded": downgraded,
                "description": edata.get('description', ''),
                "impact_factor": edata.get('edge_impact_factor', 0),
                "citation_count": edata.get('edge_citation_count', 0),
                "pmid": edata.get('pmid', ''),
                "is_chain": is_chain,
            }
        })

    return {
        "nodes": cy_nodes,
        "edges": cy_edges,
        "core": {
            "bacteria": core_bacteria,
            "metabolite": core_metabolite,
            "disease": core_disease,
        },
    }


HOP_TYPE_LABELS = {
    "microbe_metabolite": "Microbe–Metabolite",
    "metabolite_disease": "Metabolite–Disease",
    "microbe_disease": "Microbe–Disease",
}


def _parse_pmids(pmid_field) -> list:
    if not pmid_field:
        return []
    raw = str(pmid_field).strip()
    parts = raw.replace(";", ",").split(",")
    return [p.strip() for p in parts if p.strip()]


def _load_hop_evidence_for_pair(bacteria: str, metabolite: str, run_id: str = "") -> list:
    if run_id:
        try:
            from export_catalog import get_record
            rec = get_record(run_id, bacteria, metabolite)
            if rec:
                detail = rec.get("step3_detail") or {}
                hops = detail.get("hop_evidence")
                if hops:
                    return hops
                ae = rec.get("agent_evidence") or {}
                if isinstance(ae, dict) and ae.get("hop_evidence"):
                    return ae["hop_evidence"]
        except Exception:
            pass

    step2b = read_json_for_run("step2b_agent_evidence.json")
    if isinstance(step2b, list):
        for item in step2b:
            if item.get("bacteria") == bacteria and item.get("metabolite") == metabolite:
                return item.get("hop_evidence") or []
    return []


def _hop_types_for_node(node_type: str, core_role: str, label: str, bacteria: str, metabolite: str, disease: str) -> list:
    types = []
    nt = (node_type or "").lower()
    role = (core_role or "").lower()
    label_l = (label or "").lower().replace("_", " ")
    bac = (bacteria or "").lower().replace("_", " ")
    met = (metabolite or "").lower()
    dis = (disease or "").lower()

    if role == "bacteria" or nt == "microbe" or (bac and bac in label_l):
        types.extend(["microbe_metabolite", "microbe_disease"])
    if role == "metabolite" or nt == "metabolite" or (met and met.lower() in label_l):
        types.extend(["microbe_metabolite", "metabolite_disease"])
    if role == "disease" or nt == "disease" or (dis and dis in label_l):
        types.append("metabolite_disease")

    seen = set()
    out = []
    for t in types:
        if t not in seen:
            seen.add(t)
            out.append(t)
    return out


@app.get("/api/graph/node-evidence")
async def get_node_evidence(
    node_id: str = Query(...),
    bacteria: str = Query(default=""),
    metabolite: str = Query(default=""),
    disease: str = Query(default="IBD"),
    run_id: str = Query(default=""),
    node_type: str = Query(default=""),
    core_role: str = Query(default=""),
    include_abstracts: bool = Query(default=True),
):
    """Literature evidence for a KG node: incident edges (PMIDs) + hop-level PubMed co-occurrence."""
    gml_path = resolve_gml_path(run_id)
    if not gml_path.exists():
        return {"error": "GML not found"}

    G = nx.read_gml(str(gml_path))
    if node_id not in G:
        return {"error": f"Node not found: {node_id}"}

    ndata = dict(G.nodes[node_id])
    label = ndata.get("label", node_id)
    ntype = node_type or ndata.get("type", "unknown")

    edge_evidence = []
    pmid_list = []
    for src, tgt, edata in G.edges(data=True):
        if src != node_id and tgt != node_id:
            continue
        pmid_raw = edata.get("pmid", "")
        pmids = _parse_pmids(pmid_raw)
        pmid_list.extend(pmids)
        other = tgt if src == node_id else src
        relation, downgraded = _normalize_edge_relation(G, src, tgt, edata)
        edge_evidence.append({
            "source": src,
            "target": tgt,
            "source_label": G.nodes[src].get("label", src),
            "target_label": G.nodes[tgt].get("label", tgt),
            "neighbor_id": other,
            "neighbor_label": G.nodes[other].get("label", other),
            "direction": "outgoing" if src == node_id else "incoming",
            "relation": relation,
            "relation_original": edata.get("relation", ""),
            "relation_downgraded": downgraded,
            "description": edata.get("description", ""),
            "pmid": pmid_raw,
            "pmids": pmids,
            "impact_factor": edata.get("edge_impact_factor", 0),
            "citation_count": edata.get("edge_citation_count", 0),
        })

    edge_evidence.sort(
        key=lambda e: (len(e.get("pmids") or []), e.get("citation_count") or 0),
        reverse=True,
    )

    hop_types = _hop_types_for_node(ntype, core_role, label, bacteria, metabolite, disease)
    all_hops = _load_hop_evidence_for_pair(bacteria, metabolite, run_id) if bacteria and metabolite else []
    hop_evidence = [h for h in all_hops if h.get("hop_type") in hop_types] if hop_types else []

    articles = []
    unique_pmids = list(dict.fromkeys(pmid_list))[:10]
    if include_abstracts and unique_pmids:
        try:
            from build_kg import fetch_abstracts
            articles = fetch_abstracts(unique_pmids)
        except Exception:
            articles = [{"pmid": p} for p in unique_pmids]

    return {
        "node": {
            "id": node_id,
            "label": label,
            "type": ntype,
            "core_role": core_role,
            "description": ndata.get("description", ""),
            "kegg_id": ndata.get("kegg_id", ""),
        },
        "edge_evidence": edge_evidence,
        "hop_evidence": hop_evidence,
        "hop_types_considered": hop_types,
        "articles": articles,
        "pmid_count": len(unique_pmids),
    }


# --- Config Endpoint ---

@app.get("/api/config")
async def get_config():
    """Get current pipeline configuration."""
    default_run_id = ensure_default_dashboard_run_id()
    return {
        "coordinates_dir": DEFAULT_COORDINATES_DIR,
        "gml_path": str(resolve_gml_path(default_run_id)),
        "output_dir": str(OUTPUTS_DIR),
        "top_n": 200,
        "max_depth": 3,
        "disease": "IBD",
        "default_dashboard_run_id": default_run_id,
    }


# --- Compatibility endpoints for MMSage_Dashboard.html ---

@app.get("/api/dual-axis-candidates")
async def get_dual_axis_candidates(run_id: str = Query(default="")):
    """Return step3 quadrant data in the format expected by MMSage_Dashboard.html."""
    import math

    resolved_run_id = resolve_requested_run_id(run_id)
    quadrants = read_json_for_run("step3_quadrant.json", resolved_run_id)
    if not quadrants:
        return {
            "candidates": [],
            "thresholds": {"mmsage": 0.5, "evidence": 0},
            "run_id": resolved_run_id,
            "default_dashboard_run_id": ensure_default_dashboard_run_id(),
            "is_default_dashboard_run": resolved_run_id == ensure_default_dashboard_run_id(),
        }

    candidates = []
    for q in quadrants:
        bm = q.get("pair_bm_exp", 0)
        md = q.get("pair_md_exp", 0)
        mms = q.get("mmsage_norm", 0)
        # Composite: mmsage * (1 + log2(1+bm_exp)) * (1 + log2(1+md_exp))
        composite = mms * (1 + math.log2(1 + bm)) * (1 + math.log2(1 + md))

        candidates.append({
            "bacteria": q.get("bacteria", ""),
            "metabolite": q.get("metabolite", ""),
            "mmsage_norm": mms,
            "chain_count": q.get("chain_count", 0),
            "chain_novelty": q.get("chain_novelty", 1.0),
            "novelty_score": q.get("chain_novelty", 1.0),
            "evidence_foundation": q.get("evidence_foundation", 0),
            "pair_bm_exp": bm,
            "pair_md_exp": md,
            "composite_score": round(composite, 4),
            "quadrant": q.get("quadrant", ""),
            "quadrant_label": q.get("quadrant_label", ""),
            "is_dark_matter": q.get("is_dark_matter", False),
            "pairwise_counts": q.get("candidate", {}).get("pairwise_counts", {}),
        })

    # Group by quadrant (I > II > III > IV), then sort by mmsage_norm descending within each group
    _quadrant_priority = {"I": 0, "II": 1, "III": 2, "IV": 3}
    candidates.sort(
        key=lambda c: (
            _quadrant_priority.get(c.get("quadrant"), 9),
            -c["mmsage_norm"],
        )
    )

    # Compute actual thresholds (median)
    import statistics
    mms_vals = [c["mmsage_norm"] for c in candidates]
    ef_vals = [c["pair_bm_exp"] + c["pair_md_exp"] for c in candidates]
    thresholds = {
        "mmsage": statistics.median(mms_vals) if mms_vals else 0.5,
        "evidence": statistics.median(ef_vals) if ef_vals else 0,
    }

    meta = get_run_meta(resolved_run_id)
    default_run_id = ensure_default_dashboard_run_id()
    return {
        "candidates": candidates,
        "thresholds": thresholds,
        "run_id": resolved_run_id,
        "default_dashboard_run_id": default_run_id,
        "is_default_dashboard_run": resolved_run_id == default_run_id,
        "disease": meta.get("disease", "IBD"),
    }


@app.get("/api/heatmap")
async def get_heatmap():
    """Return heatmap data from coordinates CSV."""
    import pandas as pd

    coord_file = Path(DEFAULT_COORDINATES_FILE)
    if not coord_file.exists():
        return {"matrix": [], "metabolites": [], "microbes": [], "top_pairs": []}

    df = pd.read_csv(coord_file)
    rowname_col = "Row.names" if "Row.names" in df.columns else df.columns[1]

    # Row.names format: "MetaboliteName-Bacteria_name"
    # Bacteria name contains underscores, metabolite may contain hyphens.
    # Split by finding the bacteria suffix (last occurrence of "-" followed by a word with "_")
    rows = []
    for _, row in df.iterrows():
        name = str(row[rowname_col])
        # Try to find bacteria suffix: look for "-SomeWord_SomeWord" pattern at end
        # The bacteria name always contains "_" (genus_species)
        metabolite, bacteria = name, "Unknown"
        # Find the last "-" where the right part contains "_" (bacteria genus_species)
        for i in range(len(name) - 1, 0, -1):
            if name[i] == '-' and '_' in name[i+1:]:
                metabolite = name[:i]
                bacteria = name[i+1:]
                break
        rows.append({
            "Metabolite": metabolite,
            "Microbe": bacteria,
            "Pseudotime": float(row["Pseudotime"]) if "Pseudotime" in df.columns else 0,
        })

    top_pairs = sorted(rows, key=lambda r: r["Pseudotime"])
    microbes = list(set(r["Microbe"] for r in rows))
    metabolites = [r["Metabolite"] for r in sorted(rows, key=lambda r: r["Pseudotime"])]
    matrix = [[r["Pseudotime"] for r in sorted(rows, key=lambda r: r["Pseudotime"])]]

    return {
        "matrix": matrix,
        "metabolites": metabolites,
        "microbes": microbes,
        "top_pairs": [{"Microbe": r["Microbe"], "Metabolite": r["Metabolite"],
                        "Pseudotime": r["Pseudotime"]} for r in top_pairs[:50]],
    }


@app.post("/api/knowledge-graph")
async def post_knowledge_graph(req: KGRequest):
    """Return KG subgraph in Cytoscape.js format (POST, for old dashboard)."""
    return await get_graph_subgraph(
        bacteria=req.bacteria, metabolite=req.metabolite, max_hop=2)


@app.post("/api/generate-protocol")
async def generate_protocol(req: ProtocolRequest):
    """Generate through the same evidence-grounded validation planner as the Step 3 endpoint."""
    try:
        from validation_planner import generate_validation_plan as build_validation_plan

        plan = build_validation_plan(
            bacteria=req.bacteria,
            metabolite=req.metabolite,
            disease=req.disease,
            run_id=resolve_requested_run_id(req.run_id),
            mechanism_summary=req.mechanism_summary,
            mode="evidence_self_reflection",
            protocol_text="",
            research_question=req.research_question,
            prompt_constraints=req.prompt_constraints,
        )
        protocol_text = str(plan.get("validation_protocol_text") or "").strip()
        if not protocol_text:
            raise RuntimeError("Validation planner returned an empty protocol.")
        return {
            "protocol": protocol_text,
            "validation_protocol_text": protocol_text,
            "generated_by": "validation_planner",
        }
    except Exception as e:
        raise HTTPException(
            status_code=503,
            detail=f"Validation protocol generation failed: {e}",
        )


@app.post("/api/validation-plan/start")
async def start_validation_plan(req: ValidationPlanRequest):
    """Start validation-plan generation as a background job with progress logs."""
    job_id = f"vp-{uuid4().hex[:12]}"
    create_validation_plan_job_status(job_id, req)
    thread = threading.Thread(target=run_validation_plan_job, args=(job_id, req), daemon=True)
    thread.start()
    return {
        "job_id": job_id,
        "status": "started",
        "message": "Validation plan started. Poll /api/validation-plan/status/{job_id}.",
    }


@app.get("/api/validation-plan/status/{job_id}")
async def get_validation_plan_status(job_id: str):
    """Get validation-plan progress, logs, and final result."""
    status = load_validation_plan_job_status(job_id)
    if status.get("status") == "not_found":
        raise HTTPException(status_code=404, detail=f"Validation plan job '{job_id}' not found.")
    return status


@app.post("/api/validation-plan")
async def generate_validation_plan(req: ValidationPlanRequest):
    """Compatibility endpoint: generate validation plan synchronously."""
    try:
        from validation_planner import generate_validation_plan as build_validation_plan

        return build_validation_plan(
            bacteria=req.bacteria,
            metabolite=req.metabolite,
            disease=req.disease,
            run_id=resolve_requested_run_id(req.run_id),
                mechanism_summary=req.mechanism_summary,
                mode=req.mode,
                protocol_text=req.protocol_text,
                research_question=req.research_question,
                prompt_constraints=req.prompt_constraints,
            )
    except Exception as e:
        raise HTTPException(
            status_code=503,
            detail=f"Validation plan generation failed: {e}"
        )


@app.post("/api/protocol/refine/start")
async def protocol_refine_start(req: ProtocolRefineStartRequest):
    """Create a file-based Step 3 refinement session from the current protocol draft."""
    try:
        from protocol_refiner import start_refinement_session

        return start_refinement_session(
            bacteria=req.bacteria,
            metabolite=req.metabolite,
            disease=req.disease,
            run_id=resolve_requested_run_id(req.run_id),
            mechanism_summary=req.mechanism_summary,
            protocol_text=req.protocol_text,
        )
    except FileNotFoundError as e:
        raise HTTPException(status_code=404, detail=str(e))
    except Exception as e:
        raise HTTPException(status_code=503, detail=f"Failed to start protocol refinement: {e}")


@app.post("/api/protocol/refine/confirm")
async def protocol_refine_confirm(req: ProtocolRefineConfirmRequest):
    """Save user-confirmed modules and one round of extra requirements."""
    try:
        from protocol_refiner import confirm_refinement_session

        return confirm_refinement_session(
            session_id=req.session_id,
            modules=req.modules,
            extra_requirements=req.extra_requirements,
        )
    except FileNotFoundError as e:
        raise HTTPException(status_code=404, detail=str(e))
    except Exception as e:
        raise HTTPException(status_code=400, detail=f"Failed to confirm refinement session: {e}")


@app.post("/api/protocol/refine/run")
async def protocol_refine_run(req: ProtocolRefineRunRequest):
    """Run the simplified multi-agent Step 3 refinement loop."""
    try:
        from protocol_refiner import run_refinement_session

        return run_refinement_session(req.session_id)
    except FileNotFoundError as e:
        raise HTTPException(status_code=404, detail=str(e))
    except Exception as e:
        raise HTTPException(status_code=503, detail=f"Protocol refinement failed: {e}")


@app.get("/api/protocol/refine/result")
async def protocol_refine_result(session_id: str = Query(...)):
    """Fetch the current or completed Step 3 refinement session payload."""
    from protocol_refiner import get_refinement_session_result

    result = get_refinement_session_result(session_id)
    if not result:
        raise HTTPException(status_code=404, detail="Refinement session not found")
    return result


# --- Coordinates File Upload (async background pipeline) ---

import threading

_pipeline_progress = {"status": "idle", "step": "", "disease": "", "error": None}
_last_job_id = ""

def _run_pipeline_bg(csv_path: str, disease: str, output_dir: str):
    """Run pipeline in background thread, updating progress."""
    global _pipeline_progress
    _pipeline_progress = {"status": "running", "step": "step1 (MMSage signal)", "disease": disease, "error": None}
    try:
        import time as _time
        from pathlib import Path as _P

        # Monitor file creation to update progress
        step1_path = _P(output_dir) / "step1_candidates.json"
        step2_path = _P(output_dir) / "step2_chain_novelty.json"
        step2b_path = _P(output_dir) / "step2b_agent_evidence.json"
        step3_path = _P(output_dir) / "step3_quadrant.json"

        # Remove stale stage outputs so the new disease run starts from scratch.
        for p in [step1_path, step2_path, step2b_path, step3_path]:
            if p.exists():
                p.unlink()

        from run_pipeline import run_pipeline

        # Run in a sub-thread so we can monitor progress
        result_holder = [None, None]  # [result, error]
        def _inner():
            try:
                result_holder[0] = run_pipeline(
                    coordinates_file=csv_path,
                    output_dir=output_dir,
                    top_n=200,
                    disease=disease,
                    max_depth=3,
                )
            except Exception as e:
                result_holder[1] = e

        t = threading.Thread(target=_inner, daemon=True)
        t.start()

        # Monitor progress while pipeline runs
        while t.is_alive():
            _time.sleep(2)
            if step3_path.exists():
                _pipeline_progress["step"] = "step3 (quadrant assignment)"
            elif step2_path.exists():
                _pipeline_progress["step"] = "step2b (multi-agent evidence)"
            elif (_P(output_dir) / "step1_candidates.json").exists():
                _pipeline_progress["step"] = "step2 (PubMed queries, may take 3-5 min)"

        t.join()

        if result_holder[1]:
            raise result_holder[1]

        step3 = read_json("step3_quadrant.json")
        n = len(step3) if step3 else 0
        status_path = Path(output_dir) / "pipeline_status.json"
        catalog_run_id = None
        if status_path.exists():
            try:
                with open(status_path, encoding="utf-8") as f:
                    catalog_run_id = json.load(f).get("catalog_run_id")
            except Exception:
                pass
        _pipeline_progress = {
            "status": "done",
            "step": "complete",
            "disease": disease,
            "candidates": n,
            "catalog_run_id": catalog_run_id,
            "error": None,
        }
    except Exception as e:
        import traceback
        traceback.print_exc()
        _pipeline_progress = {"status": "error", "step": "", "disease": disease, "error": str(e)}


@app.post("/api/upload-coordinates")
async def upload_coordinates(
    file: UploadFile = File(...),
    disease: str = Query(default="IBD"),
):
    """Upload CSV and start an isolated pipeline job."""
    global _last_job_id

    import re

    safe_name = re.sub(r"[^\w\-_. ]", "_", Path(file.filename).name)
    if not safe_name:
        return JSONResponse(status_code=400, content={"error": "Invalid filename"})

    job_id = f"job_{uuid4().hex[:12]}"
    inputs_dir = job_inputs_dir(job_id)
    inputs_dir.mkdir(parents=True, exist_ok=True)
    dest = inputs_dir / safe_name

    content = await file.read()
    with open(dest, "wb") as f:
        f.write(content)

    save_runtime_json(
        job_status_path(job_id),
        {
            "job_id": job_id,
            "status": "queued",
            "step": "queued",
            "disease": disease,
            "filename": safe_name,
            "run_id": "",
            "error": None,
        },
    )

    subprocess.Popen(
        [
            sys.executable,
            str(BACKEND_DIR / "job_runner.py"),
            "--job-id",
            job_id,
            "--csv-path",
            str(dest),
            "--disease",
            disease,
        ],
        cwd=str(PROJECT_DIR),
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
    )
    _last_job_id = job_id

    return {
        "status": "started",
        "job_id": job_id,
        "filename": safe_name,
        "disease": disease,
        "message": "Pipeline started in isolated background job. Poll /api/jobs/{job_id}/status.",
    }


@app.get("/api/jobs/{job_id}/status")
async def job_progress(job_id: str):
    """Poll status for a specific isolated pipeline job."""
    status = load_job_status(job_id)
    if status.get("status") == "not_found":
        raise HTTPException(status_code=404, detail=f"Job not found: {job_id}")
    return status


@app.get("/api/pipeline-progress")
async def pipeline_progress(job_id: str = Query(default="")):
    """Compatibility polling endpoint; prefer /api/jobs/{job_id}/status."""
    target_job_id = job_id or _last_job_id
    if target_job_id:
        return load_job_status(target_job_id)
    return _pipeline_progress


# --- Results Catalog (JSON archive) ---

@app.get("/api/catalog/runs")
async def catalog_list_runs():
    """List archived pipeline runs."""
    from export_catalog import list_runs
    return list_runs()


@app.get("/api/catalog/search")
async def catalog_search(
    q: str = Query(default=""),
    bacteria: str = Query(default=""),
    metabolite: str = Query(default=""),
    quadrant: str = Query(default=""),
    disease: str = Query(default=""),
    run_id: str = Query(default=""),
    limit: int = Query(default=100, ge=1, le=500),
    offset: int = Query(default=0, ge=0),
):
    """Search flattened index across all archived runs."""
    from export_catalog import search_records
    return search_records(
        q=q, bacteria=bacteria, metabolite=metabolite,
        quadrant=quadrant, disease=disease, run_id=run_id,
        limit=limit, offset=offset,
    )


@app.get("/api/catalog/candidates")
async def catalog_candidates(run_id: str = Query(...)):
    """Full candidate list for one archived run."""
    from export_catalog import get_run_candidates
    items = get_run_candidates(run_id)
    if not items:
        raise HTTPException(status_code=404, detail=f"Run not found: {run_id}")
    return items


@app.get("/api/catalog/record")
async def catalog_record(
    run_id: str = Query(...),
    bacteria: str = Query(...),
    metabolite: str = Query(...),
):
    """Single candidate detail from an archived run."""
    from export_catalog import get_record
    item = get_record(run_id, bacteria, metabolite)
    if not item:
        raise HTTPException(status_code=404, detail="Record not found in catalog")
    return item


@app.get("/api/catalog/download")
async def catalog_download(
    run_id: str = Query(...),
    bacteria: str = Query(...),
    metabolite: str = Query(...),
):
    """Download merged JSON export for one catalog record."""
    from export_catalog import build_record_export, export_download_filename

    payload = build_record_export(run_id, bacteria, metabolite)
    if not payload:
        raise HTTPException(status_code=404, detail="Record not found in catalog")

    filename = export_download_filename(run_id, bacteria, metabolite)
    return JSONResponse(
        content=payload,
        media_type="application/json",
        headers={
            "Content-Disposition": f'attachment; filename="{filename}"',
        },
    )


@app.post("/api/catalog/archive")
async def catalog_archive_now(disease: str = Query(default="IBD")):
    """Manually archive current outputs/ into the results catalog."""
    from export_catalog import archive_run
    if not (OUTPUTS_DIR / "step3_quadrant.json").exists():
        raise HTTPException(status_code=400, detail="No step3_quadrant.json in outputs/")
    run_id = archive_run(output_dir=str(OUTPUTS_DIR), disease=disease)
    return {"status": "archived", "run_id": run_id}


@app.get("/api/dashboard/default-run")
async def dashboard_default_run():
    """Return the pinned default dashboard run."""
    run_id = ensure_default_dashboard_run_id()
    meta = get_run_meta(run_id)
    return {
        "run_id": run_id,
        "disease": meta.get("disease", "IBD"),
    }


# --- Serve MMSage_Dashboard.html ---

@app.get("/dashboard", response_class=HTMLResponse)
async def serve_dashboard():
    """Serve the MMSage_Dashboard.html."""
    path = FRONTEND_DIR / "MMSage_Dashboard.html"
    if path.exists():
        return HTMLResponse(path.read_text(encoding='utf-8'))
    return HTMLResponse("<h1>MMSage Dashboard not found</h1>")


# --- Frontend Serving ---

# Mount static files
if FRONTEND_DIR.exists():
    app.mount("/static", StaticFiles(directory=str(FRONTEND_DIR)), name="static")


@app.get("/", response_class=HTMLResponse)
async def serve_index():
    """Serve dashboard page."""
    index_path = FRONTEND_DIR / "index.html"
    if index_path.exists():
        return HTMLResponse(index_path.read_text(encoding='utf-8'))
    return HTMLResponse("<h1>MMSage x Chain Novelty API</h1><p>Frontend not found. API available at /docs</p>")


@app.get("/pipeline", response_class=HTMLResponse)
async def serve_pipeline():
    path = FRONTEND_DIR / "pipeline.html"
    if path.exists():
        return HTMLResponse(path.read_text(encoding='utf-8'))
    return HTMLResponse("<h1>Pipeline page not found</h1>")


@app.get("/candidate", response_class=HTMLResponse)
async def serve_candidate():
    path = FRONTEND_DIR / "candidate.html"
    if path.exists():
        return HTMLResponse(path.read_text(encoding='utf-8'))
    return HTMLResponse("<h1>Candidate page not found</h1>")


@app.get("/graph", response_class=HTMLResponse)
async def serve_graph():
    path = FRONTEND_DIR / "graph.html"
    if path.exists():
        return HTMLResponse(path.read_text(encoding='utf-8'))
    return HTMLResponse("<h1>Graph page not found</h1>")


@app.get("/evidence", response_class=HTMLResponse)
async def serve_evidence():
    path = FRONTEND_DIR / "evidence.html"
    if path.exists():
        return HTMLResponse(path.read_text(encoding='utf-8'))
    return HTMLResponse("<h1>Evidence page not found</h1>")


@app.get("/browse", response_class=HTMLResponse)
async def serve_browse():
    path = FRONTEND_DIR / "browse.html"
    if path.exists():
        return HTMLResponse(path.read_text(encoding='utf-8'))
    return HTMLResponse("<h1>Browse page not found</h1>")


@app.get("/browse/record", response_class=HTMLResponse)
async def serve_browse_record():
    path = FRONTEND_DIR / "browse_record.html"
    if path.exists():
        return HTMLResponse(path.read_text(encoding='utf-8'))
    return HTMLResponse("<h1>Record page not found</h1>")


@app.get("/protocol-refiner", response_class=HTMLResponse)
async def serve_protocol_refiner():
    path = FRONTEND_DIR / "protocol_refiner.html"
    if path.exists():
        return HTMLResponse(path.read_text(encoding='utf-8'))
    return HTMLResponse("<h1>Protocol refiner page not found</h1>")


if __name__ == "__main__":
    CATALOG_DIR.mkdir(parents=True, exist_ok=True)
    OUTPUTS_DIR.mkdir(parents=True, exist_ok=True)
    JOBS_DIR.mkdir(parents=True, exist_ok=True)
    print(f"Starting server...")
    print(f"Project: {PROJECT_DIR}")
    print(f"Frontend: {FRONTEND_DIR}")
    print(f"API docs: http://localhost:8000/docs")
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--port', type=int, default=8000)
    args = parser.parse_args()
    uvicorn.run(app, host="0.0.0.0", port=args.port)
