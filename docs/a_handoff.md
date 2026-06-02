# 同学 A 实验交接

## 当前状态

- 依赖已在当前 `python3` 环境中装好：`scanpy/pandas/scikit-learn/scipy/umap/networkx/tabulate` 等。
- 真实 PBMC3k 已下载到 `data/pbmc3k_raw.h5ad`。
- 已用真实 PBMC3k 跑完整 pipeline：QC、HVG、PCA、UMAP、t-SNE、K-means、层次聚类、Louvain、稳定性指标和 marker 注释。

## 已生成输出

- `outputs/metrics.csv`：K-means、Hierarchical、Louvain 的指标和稳定性结果。
- `outputs/cluster_assignments.csv`：每个细胞的聚类标签。
- `outputs/cluster_annotations.csv`：基于 marker genes 的 Louvain cluster 注释建议。
- `outputs/marker_scores.csv`、`outputs/marker_gene_means.csv`、`outputs/marker_gene_fractions.csv`：marker 解释表。
- `outputs/run_summary.md`：可交给同学 B 的自动摘要。
- `outputs/figures/`：QC、PCA elbow、UMAP、t-SNE、指标对比、cluster size、marker heatmap/dotplot。

## 真实 PBMC3k 关键结果

- 预处理后保留 `2638` cells 和 `1838` 个 highly-variable genes/features。
- 内部指标中 `kmeans_k5` 和 `hierarchical_k5` 的 silhouette 最高，分别约为 `0.223` 和 `0.221`。
- Louvain `resolution=0.8` 产生 `8` 个 cluster，适合用于 marker gene 生物解释。
- Louvain marker 注释建议包括 CD4 T cells、B cells、Monocytes、NK cells、T cells、Dendritic cells、Platelets。
- 稳定性指标整体可用：例如 `kmeans_k5` 的 ARI/NMI 均接近 1；Louvain `r=0.8` 的 ARI 约 `0.874`、NMI 约 `0.864`。

## 复现命令

完整复现：

```bash
cd "/Users/olinaaa_/Documents/Statistic ML"
python scripts/run_pbmc_clustering.py --stability-repeats 5
```

如果完整运行太慢，可以先跑：

```bash
python scripts/run_pbmc_clustering.py --max-cells 1200 --skip-tsne --stability-repeats 2
```

把 `outputs/run_summary.md`、`outputs/metrics.csv` 和 `outputs/figures/` 发给同学 B 更新报告。
