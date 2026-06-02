"""Project constants."""

from __future__ import annotations

KMEANS_RANGE = tuple(range(4, 13))
HIERARCHICAL_RANGE = tuple(range(4, 13))
LOUVAIN_RESOLUTIONS = (0.4, 0.8, 1.2)

MARKER_SETS: dict[str, list[str]] = {
    "T cells": ["CD3D", "CD3E", "CD2", "IL7R"],
    "CD4 T cells": ["IL7R", "CCR7", "LTB"],
    "CD8 T cells": ["CD8A", "CD8B", "GZMK"],
    "B cells": ["MS4A1", "CD79A", "CD79B"],
    "NK cells": ["GNLY", "NKG7", "KLRD1"],
    "Monocytes": ["LYZ", "S100A8", "S100A9", "FCGR3A", "MS4A7"],
    "Dendritic cells": ["FCER1A", "CST3"],
    "Platelets": ["PPBP", "PF4"],
}


def format_resolution(resolution: float) -> str:
    """Format a Louvain resolution for a stable AnnData obs key."""
    return str(resolution).replace(".", "_")
