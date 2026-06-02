"""Clustering methods used in the comparison."""

from __future__ import annotations

import pandas as pd

from pbmc_project.config import format_resolution


def get_pca_matrix(adata, n_pcs: int):
    """Return the PCA representation used by clustering and metrics."""
    if "X_pca" not in adata.obsm:
        raise KeyError("AnnData is missing obsm['X_pca']; run preprocessing first.")
    return adata.obsm["X_pca"][:, : min(n_pcs, adata.obsm["X_pca"].shape[1])]


def add_kmeans_labels(adata, k_values: tuple[int, ...], n_pcs: int, random_state: int) -> list[str]:
    """Fit K-means over a range of k values and store labels in adata.obs."""
    from sklearn.cluster import KMeans

    X = get_pca_matrix(adata, n_pcs)
    keys: list[str] = []
    for k in k_values:
        key = f"kmeans_k{k}"
        labels = KMeans(n_clusters=k, n_init=20, random_state=random_state).fit_predict(X)
        adata.obs[key] = pd.Categorical(labels.astype(str))
        keys.append(key)
    return keys


def add_hierarchical_labels(adata, k_values: tuple[int, ...], n_pcs: int) -> list[str]:
    """Fit Ward hierarchical clustering over a range of k values."""
    from sklearn.cluster import AgglomerativeClustering

    X = get_pca_matrix(adata, n_pcs)
    keys: list[str] = []
    for k in k_values:
        key = f"hierarchical_k{k}"
        labels = AgglomerativeClustering(n_clusters=k, linkage="ward").fit_predict(X)
        adata.obs[key] = pd.Categorical(labels.astype(str))
        keys.append(key)
    return keys


def add_louvain_labels(adata, resolutions: tuple[float, ...], random_state: int) -> list[str]:
    """Run Louvain graph clustering for each resolution using the Scanpy KNN graph."""
    keys: list[str] = []
    for resolution in resolutions:
        key = f"louvain_r{format_resolution(resolution)}"
        labels = networkx_louvain_labels(adata, resolution=resolution, random_state=random_state)
        adata.obs[key] = pd.Categorical(labels.astype(str))
        keys.append(key)
    return keys


def networkx_louvain_labels(adata, resolution: float, random_state: int):
    """Compute Louvain labels on adata.obsp['connectivities'] with NetworkX."""
    import networkx as nx
    import numpy as np

    if "connectivities" not in adata.obsp:
        raise KeyError("AnnData is missing obsp['connectivities']; run sc.pp.neighbors first.")

    graph = nx.from_scipy_sparse_array(adata.obsp["connectivities"], edge_attribute="weight")
    communities = nx.community.louvain_communities(
        graph,
        weight="weight",
        resolution=resolution,
        seed=random_state,
    )
    communities = sorted(communities, key=lambda community: min(community))
    labels = np.empty(adata.n_obs, dtype=int)
    for label, community in enumerate(communities):
        labels[list(community)] = label
    return labels


def infer_method_and_param(key: str) -> tuple[str, str]:
    """Parse a cluster label key into method and parameter strings."""
    if key.startswith("kmeans_k"):
        return "K-means", key.removeprefix("kmeans_k")
    if key.startswith("hierarchical_k"):
        return "Hierarchical", key.removeprefix("hierarchical_k")
    if key.startswith("louvain_r"):
        return "Louvain", key.removeprefix("louvain_r").replace("_", ".")
    return "Unknown", key
