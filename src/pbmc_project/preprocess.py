"""Data loading and preprocessing for PBMC3k."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse


@dataclass(frozen=True)
class PreprocessConfig:
    min_genes: int = 200
    min_cells: int = 3
    max_genes_by_counts: int = 2500
    max_pct_mt: float = 5.0
    target_sum: float = 1e4
    n_top_hvgs: int | None = None
    max_scale_value: float = 10.0
    n_pcs: int = 50
    random_state: int = 42
    max_cells: int | None = None


def load_pbmc(source: str, input_h5ad: Path | None = None, max_cells: int | None = None, random_state: int = 42):
    """Load PBMC data from Scanpy or a local h5ad file."""
    import scanpy as sc

    if input_h5ad is not None:
        adata = sc.read_h5ad(input_h5ad)
        source_label = str(input_h5ad)
        is_processed = False
    elif source == "synthetic":
        adata = make_synthetic_pbmc(n_cells=max_cells or 1200, random_state=random_state)
        source_label = "synthetic PBMC-like marker simulation"
        is_processed = False
    elif source == "processed":
        adata = sc.datasets.pbmc3k_processed()
        source_label = "scanpy.datasets.pbmc3k_processed"
        is_processed = True
    elif source == "raw":
        adata = sc.datasets.pbmc3k()
        source_label = "scanpy.datasets.pbmc3k"
        is_processed = False
    else:
        raise ValueError(f"Unsupported source: {source}")

    adata.var_names_make_unique()
    return adata, source_label, is_processed


def make_synthetic_pbmc(n_cells: int = 1200, n_genes: int = 1000, random_state: int = 42):
    """Create a small PBMC-like count matrix for offline pipeline smoke tests."""
    from anndata import AnnData

    rng = np.random.default_rng(random_state)
    marker_by_type = {
        "T cells": ["CD3D", "CD3E", "CD2", "IL7R", "CCR7", "LTB", "CD8A", "CD8B"],
        "B cells": ["MS4A1", "CD79A", "CD79B"],
        "NK cells": ["GNLY", "NKG7", "KLRD1", "GZMK"],
        "Monocytes": ["LYZ", "S100A8", "S100A9", "FCGR3A", "MS4A7", "CST3"],
        "Platelets": ["PPBP", "PF4"],
    }
    mito_genes = ["MT-CO1", "MT-CO2", "MT-ND1", "MT-ND2", "MT-CYB"]
    marker_genes = list(dict.fromkeys(gene for genes in marker_by_type.values() for gene in genes))
    filler_count = max(0, n_genes - len(marker_genes) - len(mito_genes))
    genes = marker_genes + mito_genes + [f"GENE{i:04d}" for i in range(filler_count)]

    cell_types = np.array(list(marker_by_type))
    proportions = np.array([0.32, 0.18, 0.16, 0.27, 0.07])
    labels = rng.choice(cell_types, size=n_cells, p=proportions)

    baseline = rng.gamma(shape=1.6, scale=0.10, size=len(genes))
    X = rng.poisson(lam=np.tile(baseline, (n_cells, 1))).astype(np.float32)
    gene_to_idx = {gene: idx for idx, gene in enumerate(genes)}

    filler_start = len(marker_genes) + len(mito_genes)
    filler_indices = np.arange(filler_start, len(genes))
    module_size = min(70, len(filler_indices) // len(cell_types))
    type_modules = {
        cell_type: filler_indices[i * module_size : (i + 1) * module_size]
        for i, cell_type in enumerate(cell_types)
    }

    for row, cell_type in enumerate(labels):
        for gene in marker_by_type[cell_type]:
            X[row, gene_to_idx[gene]] += rng.poisson(12.0)
        module = type_modules[cell_type]
        X[row, module] += rng.poisson(2.5, size=len(module))
        broad_signal = rng.choice(filler_indices, size=15, replace=False)
        X[row, broad_signal] += rng.poisson(0.8, size=len(broad_signal))

    for gene in mito_genes:
        X[:, gene_to_idx[gene]] += rng.poisson(0.08, size=n_cells)

    obs = pd.DataFrame({"synthetic_cell_type": labels}, index=[f"cell_{i:04d}" for i in range(n_cells)])
    var = pd.DataFrame({"gene_ids": genes}, index=genes)
    return AnnData(X=sparse.csr_matrix(X), obs=obs, var=var)


def calculate_qc(adata) -> pd.DataFrame:
    """Calculate QC metrics and return a detached copy of obs metrics."""
    import scanpy as sc

    adata.var["mt"] = adata.var_names.str.upper().str.startswith("MT-")
    sc.pp.calculate_qc_metrics(
        adata,
        qc_vars=["mt"],
        percent_top=None,
        log1p=False,
        inplace=True,
    )
    cols = ["n_genes_by_counts", "total_counts"]
    if "pct_counts_mt" in adata.obs:
        cols.append("pct_counts_mt")
    return adata.obs[cols].copy()


def preprocess_raw(adata, config: PreprocessConfig):
    """Apply the standard PBMC3k preprocessing workflow."""
    import scanpy as sc

    qc_before = calculate_qc(adata)

    sc.pp.filter_cells(adata, min_genes=config.min_genes)
    sc.pp.filter_genes(adata, min_cells=config.min_cells)
    calculate_qc(adata)
    adata = adata[
        (adata.obs["n_genes_by_counts"] < config.max_genes_by_counts)
        & (adata.obs["pct_counts_mt"] < config.max_pct_mt)
    ].copy()

    if config.max_cells is not None and adata.n_obs > config.max_cells:
        sc.pp.subsample(
            adata,
            n_obs=config.max_cells,
            random_state=config.random_state,
            copy=False,
        )

    qc_after = calculate_qc(adata)

    sc.pp.normalize_total(adata, target_sum=config.target_sum)
    sc.pp.log1p(adata)
    adata.raw = adata

    hvg_kwargs = {"flavor": "seurat"}
    if config.n_top_hvgs is not None:
        hvg_kwargs["n_top_genes"] = config.n_top_hvgs
    else:
        hvg_kwargs.update({"min_mean": 0.0125, "max_mean": 3, "min_disp": 0.5})
    sc.pp.highly_variable_genes(adata, **hvg_kwargs)
    adata = adata[:, adata.var["highly_variable"]].copy()

    sc.pp.scale(adata, max_value=config.max_scale_value)
    n_comps = min(config.n_pcs, adata.n_obs - 1, adata.n_vars - 1)
    sc.tl.pca(adata, svd_solver="arpack", n_comps=n_comps, random_state=config.random_state)

    return adata, qc_before, qc_after


def prepare_processed(adata, config: PreprocessConfig):
    """Prepare a processed Scanpy dataset for the common downstream pipeline."""
    import scanpy as sc

    if config.max_cells is not None and adata.n_obs > config.max_cells:
        sc.pp.subsample(
            adata,
            n_obs=config.max_cells,
            random_state=config.random_state,
            copy=False,
        )

    if adata.raw is None:
        adata.raw = adata

    qc_before = calculate_qc(adata)
    qc_after = qc_before.copy()

    if "X_pca" not in adata.obsm:
        sc.pp.scale(adata, max_value=config.max_scale_value)
        n_comps = min(config.n_pcs, adata.n_obs - 1, adata.n_vars - 1)
        sc.tl.pca(adata, svd_solver="arpack", n_comps=n_comps, random_state=config.random_state)

    return adata, qc_before, qc_after


def compute_embeddings(adata, n_neighbors: int, n_pcs: int, random_state: int, compute_tsne: bool):
    """Compute neighbor graph, UMAP, and optionally t-SNE."""
    import scanpy as sc

    pca_dims = adata.obsm["X_pca"].shape[1]
    n_pcs = min(n_pcs, pca_dims)
    sc.pp.neighbors(adata, n_neighbors=n_neighbors, n_pcs=n_pcs, random_state=random_state)
    sc.tl.umap(adata, random_state=random_state)
    if compute_tsne:
        sc.tl.tsne(adata, n_pcs=n_pcs, random_state=random_state)
    return adata
