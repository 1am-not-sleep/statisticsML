"""Marker-gene based interpretation utilities."""

from __future__ import annotations

import numpy as np
import pandas as pd
from scipy import sparse

from pbmc_project.config import MARKER_SETS


def marker_gene_tables(adata, cluster_key: str, marker_sets: dict[str, list[str]] = MARKER_SETS):
    """Compute marker gene mean expression, fraction expressed, and cell-type scores."""
    genes = sorted({gene for marker_list in marker_sets.values() for gene in marker_list})
    source = adata.raw.to_adata() if adata.raw is not None else adata
    available = [gene for gene in genes if gene in source.var_names]
    if not available:
        raise ValueError("None of the configured marker genes are present in this dataset.")

    means, fractions = _cluster_gene_summary(source, cluster_key, available, adata.obs[cluster_key])

    scores = pd.DataFrame(index=means.index)
    available_by_type: dict[str, list[str]] = {}
    for cell_type, marker_list in marker_sets.items():
        present = [gene for gene in marker_list if gene in means.columns]
        available_by_type[cell_type] = present
        scores[cell_type] = means[present].mean(axis=1) if present else np.nan

    annotations = []
    for cluster, row in scores.iterrows():
        if row.notna().any():
            best_type = str(row.idxmax())
            best_score = float(row.max())
            marker_text = ", ".join(available_by_type[best_type])
        else:
            best_type = "Unknown"
            best_score = float("nan")
            marker_text = ""
        annotations.append(
            {
                "cluster": cluster,
                "suggested_cell_type": best_type,
                "marker_score": best_score,
                "supporting_markers": marker_text,
            }
        )

    return means, fractions, scores, pd.DataFrame(annotations)


def _cluster_gene_summary(source_adata, cluster_key: str, genes: list[str], labels: pd.Series):
    sub = source_adata[:, genes]
    X = sub.X
    cluster_labels = labels.astype(str).to_numpy()
    clusters = sorted(np.unique(cluster_labels), key=_cluster_sort_key)

    mean_rows = []
    fraction_rows = []
    for cluster in clusters:
        mask = cluster_labels == cluster
        cluster_matrix = X[mask]
        mean_rows.append(_mean_vector(cluster_matrix))
        fraction_rows.append(_fraction_positive(cluster_matrix))

    means = pd.DataFrame(mean_rows, index=clusters, columns=genes)
    fractions = pd.DataFrame(fraction_rows, index=clusters, columns=genes)
    means.index.name = cluster_key
    fractions.index.name = cluster_key
    return means, fractions


def _mean_vector(matrix) -> np.ndarray:
    values = matrix.mean(axis=0)
    if sparse.issparse(values):
        values = values.A1
    return np.asarray(values).ravel()


def _fraction_positive(matrix) -> np.ndarray:
    positive = matrix > 0
    values = positive.mean(axis=0)
    if sparse.issparse(values):
        values = values.A1
    return np.asarray(values).ravel()


def _cluster_sort_key(value: str):
    return (0, int(value)) if value.isdigit() else (1, value)
