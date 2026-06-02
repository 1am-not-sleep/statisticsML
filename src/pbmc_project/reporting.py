"""Write lightweight markdown summaries from a run."""

from __future__ import annotations

from pathlib import Path

import pandas as pd


def write_run_summary(
    output_dir: Path,
    figures_dir: Path,
    source_label: str,
    n_obs: int,
    n_vars: int,
    selected_keys: list[str],
    annotation_key: str,
    metrics: pd.DataFrame,
    annotations: pd.DataFrame,
) -> Path:
    """Write a Chinese markdown summary suitable for the report draft."""
    path = output_dir / "run_summary.md"
    best_rows = metrics.sort_values("silhouette", ascending=False).head(8)

    lines = [
        "# PBMC3k 聚类实验运行摘要",
        "",
        f"- 数据来源：`{source_label}`",
        f"- 预处理后规模：`{n_obs}` cells × `{n_vars}` highly-variable genes/features",
        f"- 重点展示聚类结果：{', '.join(f'`{key}`' for key in selected_keys)}",
        f"- marker gene 注释基于：`{annotation_key}`",
        "",
        "## 指标表现 Top 8",
        "",
        best_rows.to_markdown(index=False),
        "",
        "## Marker Gene 注释建议",
        "",
        annotations.to_markdown(index=False),
        "",
        "## 主要图表",
        "",
        f"- QC 分布：`{figures_dir / 'qc_histograms.png'}`",
        f"- PCA elbow：`{figures_dir / 'pca_elbow.png'}`",
        f"- UMAP 聚类图：`{figures_dir / 'umap_clusters.png'}`",
        f"- 指标对比：`{figures_dir / 'metric_summary.png'}`",
        f"- cluster size：`{figures_dir / 'cluster_sizes.png'}`",
        f"- marker heatmap：`{figures_dir / 'marker_heatmap.png'}`",
        f"- marker dotplot：`{figures_dir / 'marker_dotplot.png'}`",
        "",
        "## 报告可写结论",
        "",
        "- PCA 将高维表达矩阵压缩为低维连续表示，便于后续聚类和可视化。",
        "- K-means 与层次聚类主要依赖欧氏距离和簇形状假设；Louvain 在细胞邻接图上寻找社区，更适合单细胞常见的非线性结构。",
        "- 内部指标只能评价几何紧致性，最终解释需要结合 marker genes 和免疫细胞背景。",
    ]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path
