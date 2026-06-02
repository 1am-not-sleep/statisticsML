#!/usr/bin/env python
"""Check whether the required PBMC project dependencies import successfully."""

from __future__ import annotations

import importlib
import os
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
os.environ.setdefault("MPLCONFIGDIR", str(ROOT / ".matplotlib_cache"))
os.environ.setdefault("NUMBA_CACHE_DIR", str(ROOT / ".numba_cache"))

PACKAGES = [
    "numpy",
    "pandas",
    "scipy",
    "sklearn",
    "matplotlib",
    "scanpy",
    "anndata",
    "umap",
    "networkx",
]


def main() -> None:
    missing: list[str] = []
    for package in PACKAGES:
        try:
            module = importlib.import_module(package)
            version = getattr(module, "__version__", "installed")
            print(f"{package}: {version}")
        except Exception as exc:  # pragma: no cover - diagnostic script
            missing.append(package)
            print(f"{package}: MISSING ({exc})")
    if missing:
        raise SystemExit(f"Missing dependencies: {', '.join(missing)}")


if __name__ == "__main__":
    main()
