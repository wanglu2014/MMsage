from __future__ import annotations

import json
from pathlib import Path
from typing import Any


BACKEND_DIR = Path(__file__).parent
PROJECT_DIR = BACKEND_DIR.parent
DATA_DIR = PROJECT_DIR / "data"
OUTPUTS_DIR = PROJECT_DIR / "outputs"
JOBS_DIR = DATA_DIR / "jobs"
APP_STATE_PATH = DATA_DIR / "app_state.json"


def load_json(path: Path) -> Any:
    if not path.exists():
        return None
    with open(path, "r", encoding="utf-8") as f:
        return json.load(f)


def save_json(path: Path, data: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8") as f:
        json.dump(data, f, indent=2, ensure_ascii=False)


def load_app_state() -> dict:
    data = load_json(APP_STATE_PATH)
    return data if isinstance(data, dict) else {}


def save_app_state(state: dict) -> None:
    save_json(APP_STATE_PATH, state)


def job_dir(job_id: str) -> Path:
    return JOBS_DIR / job_id


def job_output_dir(job_id: str) -> Path:
    return job_dir(job_id) / "outputs"


def job_inputs_dir(job_id: str) -> Path:
    return job_dir(job_id) / "inputs"


def job_status_path(job_id: str) -> Path:
    return job_dir(job_id) / "status.json"
