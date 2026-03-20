import json
import os
import subprocess
from pathlib import Path
from typing import Dict, List


def get_repo_root() -> Path:
    return Path(__file__).resolve().parents[3]


def load_config(config_path: str) -> Dict:
    cfg_path = Path(config_path)
    if not cfg_path.is_absolute():
        cfg_path = get_repo_root() / cfg_path
    with open(cfg_path, "r") as f:
        raw = f.read()
    try:
        import yaml

        return yaml.safe_load(raw)
    except ImportError:
        pass

    # If YAML support is unavailable in the current interpreter, attempt JSON.
    try:
        return json.loads(raw)
    except json.JSONDecodeError:
        # Last resort: use "python" from PATH (often different from notebook kernel)
        # to parse YAML when PyYAML is installed there.
        cmd = [
            "python",
            "-c",
            (
                "import json, sys, yaml; "
                "print(json.dumps(yaml.safe_load(open(sys.argv[1]).read())))"
            ),
            str(cfg_path),
        ]
        proc = subprocess.run(
            cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            universal_newlines=True,
        )
        if proc.returncode != 0:
            raise RuntimeError(
                "Unable to parse config. Install PyYAML in current environment or "
                f"use JSON config. Parser stderr: {proc.stderr.strip()}"
            )
        return json.loads(proc.stdout)


def resolve_path(path_like: str) -> Path:
    path = Path(path_like)
    if path.is_absolute():
        return path
    return get_repo_root() / path


def read_id_list(path_like: str) -> List[str]:
    path = resolve_path(path_like)
    with open(path, "r") as f:
        rows = [line.strip() for line in f if line.strip()]
    return rows


def ensure_dir(path_like: str) -> Path:
    path = resolve_path(path_like)
    os.makedirs(path, exist_ok=True)
    return path

