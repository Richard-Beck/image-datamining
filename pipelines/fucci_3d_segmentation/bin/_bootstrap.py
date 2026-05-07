from __future__ import annotations

import sys
from pathlib import Path


def add_src_to_path() -> Path:
    repo_root = Path(__file__).resolve().parents[3]
    src = repo_root / "src"
    sys.path.insert(0, str(src))
    return repo_root

