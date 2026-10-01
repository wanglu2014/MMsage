from __future__ import annotations

import json
import os
import tempfile
import threading
import time
from pathlib import Path
from typing import Any


BACKEND_DIR = Path(__file__).parent
PROJECT_DIR = BACKEND_DIR.parent
DATA_DIR = PROJECT_DIR / "data"
OUTPUTS_DIR = PROJECT_DIR / "outputs"
JOBS_DIR = DATA_DIR / "jobs"
APP_STATE_PATH = DATA_DIR / "app_state.json"
_json_lock = threading.RLock()


def load_json(path: Path) -> Any:
    with _json_lock:
        if not path.exists():
            return None
        with open(path, "r", encoding="utf-8") as f:
            return json.load(f)


def save_json(path: Path, data: Any) -> None:
    path = Path(path)
    payload = json.dumps(data, indent=2, ensure_ascii=False)
    with _json_lock:
        path.parent.mkdir(parents=True, exist_ok=True)
        # Close the temporary file before replacement on Windows. Unique names
        # also prevent separate workers from overwriting each other's temp file.
        fd, name = tempfile.mkstemp(prefix=path.name + ".", suffix=".tmp", dir=path.parent)
        temporary = Path(name)
        try:
            with os.fdopen(fd, "w", encoding="utf-8") as f:
                f.write(payload)
            for attempt in range(20):
                try:
                    os.replace(temporary, path)
                    break
                except PermissionError:
                    if attempt == 19:
                        raise
                    # Other processes/Windows scanners may briefly hold a handle.
                    time.sleep(0.05)
        finally:
            temporary.unlink(missing_ok=True)


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
