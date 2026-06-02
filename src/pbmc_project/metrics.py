"""Clustering quality and stability metrics."""

from __future__ import annotations

import math
import re

import numpy as np
import pandas as pd

from pbmc_project.clustering import get_pca_matrix, infer_method_and_param, networkx_louvain_labels


def evaluate_clusterings(adata, keys: list[str], n_pcs: int) -> pd.DataFrame:
    """Compute internal clustering metrics for each label column."""
    from sklearn.metrics import (
        calinski_harabasz_score,
        davies_bouldin_score,
        silhouette_score,
    )

    X = get_pca_matrix(adata, n_pcs)
    rows: list[dict[str, object]] = []
    for key in keys:
        labels = adata.obs[key].astype(str).to_numpy()
        method, param = infer_method_and_param(key)
        n_clusters = len(np.unique(labels))
        row: dict[str, object] = {
            "key": key,
            "method": method,
            "param": param,
            "n_clusters": n_clusters,
            "min_cluster_size": int(pd.Series(labels).value_counts().min()),
            "max_cluster_size": int(pd.Series(labels).value_counts().max()),
        }
        if 1 < n_clusters < len(labels):
            row["silhouette"] = float(silhouette_score(X, labels))
            row["calinski_harabasz"] = float(calinski_harabasz_score(X, labels))
            row["davies_bouldin"] = float(davies_bouldin_score(X, labels))
        else:
            row["silhouette"] = math.nan
            row["calinski_harabasz"] = math.nan
            row["davies_bouldin"] = math.nan
        rows.append(row)
    return pd.DataFrame(rows)


def choose_best_key(metrics: pd.DataFrame, method: str) -> str | None:
    """Choose the best key within a method by silhouette score."""
    subset = metrics[(metrics["method"] == method) & metrics["silhouette"].notna()]
    if subset.empty:
        return None
    return str(subset.sort_values("silhouette", ascending=False).iloc[0]["key"])


def compute_stability(
    adata,
    keys: list[str],
    n_pcs: int,
    repeats: int,
    sample_fraction: float,
    random_state: int,
) -> pd.DataFrame:
    """Estimate clustering stability by refitting on random cell subsamples."""
    if repeats <= 0:
        return pd.DataFrame(columns=["key", "ari_mean", "ari_sd", "nmi_mean", "nmi_sd"])

    from sklearn.cluster import AgglomerativeClustering, KMeans
    from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score

    X = get_pca_matrix(adata, n_pcs)
    rng = np.random.default_rng(random_state)
    rows: list[dict[str, object]] = []

    for key in keys:
        ref_labels = adata.obs[key].astype(str).to_numpy()
        ari_values: list[float] = []
        nmi_values: list[float] = []

        for repeat in range(repeats):
            sample_size = max(10, int(len(ref_labels) * sample_fraction))
            sample_idx = np.sort(rng.choice(len(ref_labels), size=sample_size, replace=False))
            sampled_labels = _refit_on_subset(
                adata=adata,
                X=X,
                key=key,
                sample_idx=sample_idx,
                n_pcs=n_pcs,
                random_state=random_state + repeat + 1,
                kmeans_cls=KMeans,
                hierarchical_cls=AgglomerativeClustering,
            )
            ari_values.append(adjusted_rand_score(ref_labels[sample_idx], sampled_labels))
            nmi_values.append(normalized_mutual_info_score(ref_labels[sample_idx], sampled_labels))

        rows.append(
            {
                "key": key,
                "ari_mean": float(np.mean(ari_values)),
                "ari_sd": float(np.std(ari_values, ddof=1)) if len(ari_values) > 1 else 0.0,
                "nmi_mean": float(np.mean(nmi_values)),
                "nmi_sd": float(np.std(nmi_values, ddof=1)) if len(nmi_values) > 1 else 0.0,
            }
        )

    return pd.DataFrame(rows)


def _refit_on_subset(
    adata,
    X,
    key: str,
    sample_idx,
    n_pcs: int,
    random_state: int,
    kmeans_cls,
    hierarchical_cls,
) -> np.ndarray:
    if key.startswith("kmeans_k"):
        k = int(key.removeprefix("kmeans_k"))
        return kmeans_cls(n_clusters=k, n_init=20, random_state=random_state).fit_predict(X[sample_idx])

    if key.startswith("hierarchical_k"):
        k = int(key.removeprefix("hierarchical_k"))
        return hierarchical_cls(n_clusters=k, linkage="ward").fit_predict(X[sample_idx])

    if key.startswith("louvain_r"):
        import scanpy as sc

        match = re.search(r"louvain_r(.+)$", key)
        if match is None:
            raise ValueError(f"Cannot parse Louvain resolution from {key}")
        resolution = float(match.group(1).replace("_", "."))
        subset = adata[sample_idx].copy()
        pcs = min(n_pcs, subset.obsm["X_pca"].shape[1])
        sc.pp.neighbors(subset, n_neighbors=10, n_pcs=pcs, random_state=random_state)
        return networkx_louvain_labels(
            subset,
            resolution=resolution,
            random_state=random_state,
        ).astype(str)

    raise ValueError(f"Unsupported clustering key for stability: {key}")
