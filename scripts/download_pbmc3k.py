#!/usr/bin/env python
"""Download the public Scanpy PBMC3k raw h5ad file into data/."""

from __future__ import annotations

import urllib.request
import ssl
from pathlib import Path

import certifi

URL = "https://exampledata.scverse.org/scanpy/pbmc3k_raw.h5ad"
ROOT = Path(__file__).resolve().parents[1]
TARGET = ROOT / "data" / "pbmc3k_raw.h5ad"


def main() -> None:
    TARGET.parent.mkdir(parents=True, exist_ok=True)
    print(f"Downloading {URL}")
    print(f"Target: {TARGET}")
    context = ssl.create_default_context(cafile=certifi.where())
    request = urllib.request.Request(
        URL,
        headers={"User-Agent": "pip/25.3 Python/3.13 PBMC3k-course-project"},
    )
    with urllib.request.urlopen(request, context=context) as response, TARGET.open("wb") as handle:
        handle.write(response.read())
    print("Done.")


if __name__ == "__main__":
    main()
