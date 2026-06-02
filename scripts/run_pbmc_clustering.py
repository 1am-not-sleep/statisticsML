#!/usr/bin/env python
"""Run the PBMC3k dimensionality reduction and clustering comparison."""

from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
os.environ.setdefault("MPLCONFIGDIR", str(ROOT / ".matplotlib_cache"))
os.environ.setdefault("NUMBA_CACHE_DIR", str(ROOT / ".numba_cache"))
sys.path.insert(0, str(ROOT / "src"))

from pbmc_project.clustering import (  # noqa: E402
    add_hierarchical_labels,
    add_kmeans_labels,
    add_louvain_labels,
)
from pbmc_project.config import HIERARCHICAL_RANGE, KMEANS_RANGE, LOUVAIN_RESOLUTIONS  # noqa: E402
from pbmc_project.markers import marker_gene_tables  # noqa: E402
from pbmc_project.metrics import choose_best_key, compute_stability, evaluate_clusterings  # noqa: E402
from pbmc_project.plots import (  # noqa: E402
    plot_cluster_sizes,
    plot_embedding,
    plot_marker_dotplot,
    plot_marker_heatmap,
    plot_metric_summary,
    plot_pca_elbow,
    plot_qc,
)
from pbmc_project.preprocess import (  # noqa: E402
    PreprocessConfig,
    compute_embeddings,
    load_pbmc,
    prepare_processed,
    preprocess_raw,
)
from pbmc_project.reporting import write_run_summary  # noqa: E402


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-h5ad", type=Path, default=None, help="Optional local AnnData h5ad file.")
    parser.add_argument(
        "--source",
        choices=["raw", "processed", "synthetic"],
        default="raw",
        help="Scanpy PBMC source when --input-h5ad is not supplied.",
    )
    parser.add_argument("--output-dir", type=Path, default=ROOT / "outputs")
    parser.add_argument("--max-cells", type=int, default=None, help="Subsample cells for fast debugging.")
    parser.add_argument("--n-pcs", type=int, default=40)
    parser.add_argument("--n-neighbors", type=int, default=10)
    parser.add_argument("--random-state", type=int, default=42)
    parser.add_argument("--stability-repeats", type=int, default=5)
    parser.add_argument("--stability-fraction", type=float, default=0.8)
    parser.add_argument("--skip-tsne", action="store_true", help="Skip t-SNE for faster runs.")
    parser.add_argument("--skip-save-adata", action="store_true", help="Do not save the labeled h5ad file.")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    output_dir = args.output_dir
    figures_dir = output_dir / "figures"
    figures_dir.mkdir(parents=True, exist_ok=True)

    config = PreprocessConfig(
        n_pcs=max(args.n_pcs, 50),
        random_state=args.random_state,
        max_cells=args.max_cells,
    )

    local_raw = ROOT / "data" / "pbmc3k_raw.h5ad"
    if args.input_h5ad is None and args.source == "raw" and local_raw.exists():
        args.input_h5ad = local_raw

    adata, source_label, is_processed = load_pbmc(
        args.source,
        args.input_h5ad,
        max_cells=args.max_cells,
        random_state=args.random_state,
    )
    if is_processed:
        adata, qc_before, qc_after = prepare_processed(adata, config)
    else:
        adata, qc_before, qc_after = preprocess_raw(adata, config)

    compute_embeddings(
        adata,
        n_neighbors=args.n_neighbors,
        n_pcs=args.n_pcs,
        random_state=args.random_state,
        compute_tsne=not args.skip_tsne,
    )

    cluster_keys: list[str] = []
    cluster_keys.extend(add_kmeans_labels(adata, KMEANS_RANGE, args.n_pcs, args.random_state))
    cluster_keys.extend(add_hierarchical_labels(adata, HIERARCHICAL_RANGE, args.n_pcs))
    cluster_keys.extend(add_louvain_labels(adata, LOUVAIN_RESOLUTIONS, args.random_state))

    metrics = evaluate_clusterings(adata, cluster_keys, args.n_pcs)
    stability = compute_stability(
        adata,
        cluster_keys,
        n_pcs=args.n_pcs,
        repeats=args.stability_repeats,
        sample_fraction=args.stability_fraction,
        random_state=args.random_state,
    )
    if not stability.empty:
        metrics = metrics.merge(stability, on="key", how="left")

    best_kmeans = choose_best_key(metrics, "K-means")
    best_hierarchical = choose_best_key(metrics, "Hierarchical")
    preferred_louvain = "louvain_r0_8" if "louvain_r0_8" in cluster_keys else cluster_keys[-1]
    selected_keys = [key for key in [best_kmeans, best_hierarchical, preferred_louvain] if key is not None]
    annotation_key = preferred_louvain

    means, fractions, marker_scores, annotations = marker_gene_tables(adata, annotation_key)

    metrics.to_csv(output_dir / "metrics.csv", index=False)
    adata.obs[cluster_keys].to_csv(output_dir / "cluster_assignments.csv")
    annotations.to_csv(output_dir / "cluster_annotations.csv", index=False)
    marker_scores.to_csv(output_dir / "marker_scores.csv")
    means.to_csv(output_dir / "marker_gene_means.csv")
    fractions.to_csv(output_dir / "marker_gene_fractions.csv")

    plot_qc(qc_before, qc_after, figures_dir / "qc_histograms.png")
    plot_pca_elbow(adata, figures_dir / "pca_elbow.png", max_pcs=args.n_pcs)
    plot_embedding(adata, "umap", selected_keys, figures_dir / "umap_clusters.png")
    plot_embedding(adata, "tsne", selected_keys, figures_dir / "tsne_clusters.png")
    plot_metric_summary(metrics, figures_dir / "metric_summary.png")
    plot_cluster_sizes(adata, selected_keys, figures_dir / "cluster_sizes.png")
    plot_marker_heatmap(marker_scores, figures_dir / "marker_heatmap.png")
    plot_marker_dotplot(means, fractions, figures_dir / "marker_dotplot.png")

    if not args.skip_save_adata:
        adata.write_h5ad(output_dir / "pbmc3k_clustered.h5ad", compression="gzip")

    summary_path = write_run_summary(
        output_dir=output_dir,
        figures_dir=figures_dir,
        source_label=source_label,
        n_obs=adata.n_obs,
        n_vars=adata.n_vars,
        selected_keys=selected_keys,
        annotation_key=annotation_key,
        metrics=metrics,
        annotations=annotations,
    )

    print(f"Completed PBMC clustering run. Summary: {summary_path}")


if __name__ == "__main__":
    main()
