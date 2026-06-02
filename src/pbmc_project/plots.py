"""Plotting helpers for the PBMC clustering project."""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def plot_qc(qc_before: pd.DataFrame, qc_after: pd.DataFrame, path: Path) -> None:
    columns = [col for col in ["n_genes_by_counts", "total_counts", "pct_counts_mt"] if col in qc_before]
    fig, axes = plt.subplots(1, len(columns), figsize=(5 * len(columns), 4))
    if len(columns) == 1:
        axes = [axes]
    for ax, col in zip(axes, columns, strict=False):
        ax.hist(qc_before[col], bins=50, alpha=0.45, label="before QC")
        ax.hist(qc_after[col], bins=50, alpha=0.65, label="after QC")
        ax.set_title(col)
        ax.set_ylabel("cells")
        ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)


def plot_pca_elbow(adata, path: Path, max_pcs: int = 50) -> None:
    ratios = np.asarray(adata.uns["pca"]["variance_ratio"])[:max_pcs]
    fig, ax = plt.subplots(figsize=(7, 4))
    ax.plot(np.arange(1, len(ratios) + 1), ratios, marker="o", linewidth=1.5)
    ax.set_xlabel("principal component")
    ax.set_ylabel("explained variance ratio")
    ax.set_title("PCA elbow")
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)


def plot_embedding(adata, basis: str, color_keys: list[str], path: Path) -> None:
    obsm_key = f"X_{basis}"
    if obsm_key not in adata.obsm:
        return

    coords = adata.obsm[obsm_key]
    n_panels = len(color_keys)
    fig, axes = plt.subplots(1, n_panels, figsize=(5 * n_panels, 4.5), squeeze=False)
    for ax, key in zip(axes.ravel(), color_keys, strict=False):
        labels = adata.obs[key].astype("category")
        codes = labels.cat.codes.to_numpy()
        scatter = ax.scatter(coords[:, 0], coords[:, 1], c=codes, s=7, cmap="tab20", alpha=0.82)
        ax.set_title(key)
        ax.set_xlabel(f"{basis.upper()} 1")
        ax.set_ylabel(f"{basis.upper()} 2")
        ax.set_xticks([])
        ax.set_yticks([])
        if labels.cat.categories.size <= 12:
            handles, _ = scatter.legend_elements(num=labels.cat.categories.size)
            ax.legend(
                handles,
                labels.cat.categories.astype(str),
                title="cluster",
                frameon=False,
                fontsize=7,
                loc="best",
            )
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)


def plot_metric_summary(metrics: pd.DataFrame, path: Path) -> None:
    metric_cols = ["silhouette", "calinski_harabasz", "davies_bouldin"]
    fig, axes = plt.subplots(1, len(metric_cols), figsize=(5.4 * len(metric_cols), 4.2))
    for ax, metric in zip(axes, metric_cols, strict=False):
        for method, group in metrics.groupby("method", sort=False):
            group = group.copy()
            x = np.arange(len(group))
            ax.plot(x, group[metric], marker="o", label=method)
            for xpos, label in zip(x, group["param"], strict=False):
                ax.text(xpos, group[metric].iloc[int(xpos)], str(label), fontsize=7, ha="center", va="bottom")
        ax.set_title(metric)
        ax.set_xlabel("parameter scan")
        ax.grid(alpha=0.2)
    axes[0].legend(frameon=False)
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)


def plot_cluster_sizes(adata, keys: list[str], path: Path) -> None:
    fig, axes = plt.subplots(len(keys), 1, figsize=(8, 2.8 * len(keys)), squeeze=False)
    for ax, key in zip(axes.ravel(), keys, strict=False):
        counts = adata.obs[key].astype(str).value_counts().sort_index(key=lambda idx: idx.map(_sort_value))
        ax.bar(counts.index.astype(str), counts.values)
        ax.set_title(key)
        ax.set_ylabel("cells")
        ax.tick_params(axis="x", rotation=45)
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)


def plot_marker_heatmap(scores: pd.DataFrame, path: Path) -> None:
    fig, ax = plt.subplots(figsize=(max(8, scores.shape[1] * 1.1), max(4, scores.shape[0] * 0.45)))
    values = scores.to_numpy(dtype=float)
    image = ax.imshow(values, aspect="auto", cmap="viridis")
    ax.set_xticks(np.arange(scores.shape[1]))
    ax.set_xticklabels(scores.columns, rotation=45, ha="right")
    ax.set_yticks(np.arange(scores.shape[0]))
    ax.set_yticklabels(scores.index.astype(str))
    ax.set_xlabel("marker set")
    ax.set_ylabel("cluster")
    fig.colorbar(image, ax=ax, label="mean log-normalized expression")
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)


def plot_marker_dotplot(means: pd.DataFrame, fractions: pd.DataFrame, path: Path) -> None:
    fig, ax = plt.subplots(figsize=(max(9, means.shape[1] * 0.35), max(4, means.shape[0] * 0.45)))
    x_positions = np.arange(means.shape[1])
    y_positions = np.arange(means.shape[0])
    for y, cluster in enumerate(means.index):
        sizes = 25 + 260 * fractions.loc[cluster].to_numpy(dtype=float)
        colors = means.loc[cluster].to_numpy(dtype=float)
        ax.scatter(x_positions, np.full_like(x_positions, y), s=sizes, c=colors, cmap="magma", alpha=0.82)
    ax.set_xticks(x_positions)
    ax.set_xticklabels(means.columns, rotation=60, ha="right")
    ax.set_yticks(y_positions)
    ax.set_yticklabels(means.index.astype(str))
    ax.set_xlabel("marker gene")
    ax.set_ylabel("cluster")
    ax.set_title("Marker gene dotplot: color = mean expression, size = fraction expressed")
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)


def _sort_value(value: str):
    return int(value) if str(value).isdigit() else str(value)
