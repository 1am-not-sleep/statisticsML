# PBMC3k 聚类实验运行摘要

- 数据来源：`/Users/olinaaa_/Documents/Statistic ML/data/pbmc3k_raw.h5ad`
- 预处理后规模：`2638` cells × `1838` highly-variable genes/features
- 重点展示聚类结果：`kmeans_k5`, `hierarchical_k5`, `louvain_r0_8`
- marker gene 注释基于：`louvain_r0_8`

## 指标表现 Top 8

| key             | method       |   param |   n_clusters |   min_cluster_size |   max_cluster_size |   silhouette |   calinski_harabasz |   davies_bouldin |   ari_mean |     ari_sd |   nmi_mean |     nmi_sd |
|:----------------|:-------------|--------:|-------------:|-------------------:|-------------------:|-------------:|--------------------:|-----------------:|-----------:|-----------:|-----------:|-----------:|
| kmeans_k5       | K-means      |     5   |            5 |                 13 |               1385 |     0.222988 |             410.598 |          1.52024 |   0.995817 | 0.00253841 |   0.991943 | 0.00431578 |
| hierarchical_k5 | Hierarchical |     5   |            5 |                 12 |               1446 |     0.220512 |             395.558 |          1.49547 |   0.984193 | 0.0091336  |   0.971546 | 0.0123831  |
| kmeans_k4       | K-means      |     4   |            4 |                 13 |               1611 |     0.212667 |             412.618 |          1.49868 |   0.882718 | 0.107421   |   0.917436 | 0.0755502  |
| hierarchical_k4 | Hierarchical |     4   |            4 |                 12 |               1606 |     0.209771 |             409.236 |          1.49059 |   0.970561 | 0.0625513  |   0.97505  | 0.047821   |
| kmeans_k6       | K-means      |     6   |            6 |                 13 |               1384 |     0.200651 |             359.594 |          1.8385  |   0.990844 | 0.00679067 |   0.98287  | 0.00810059 |
| louvain_r0_4    | Louvain      |     0.4 |            6 |                 13 |               1185 |     0.198307 |             334.402 |          1.7069  |   0.886134 | 0.115266   |   0.911636 | 0.0640383  |
| hierarchical_k6 | Hierarchical |     6   |            6 |                 12 |               1152 |     0.163671 |             347.596 |          1.97441 |   0.942681 | 0.00527793 |   0.923253 | 0.00346104 |
| kmeans_k8       | K-means      |     8   |            8 |                 13 |               1156 |     0.159567 |             292.298 |          2.01584 |   0.920169 | 0.136226   |   0.945755 | 0.0600439  |

## Marker Gene 注释建议

|   cluster | suggested_cell_type   |   marker_score | supporting_markers                 |
|----------:|:----------------------|---------------:|:-----------------------------------|
|         0 | CD4 T cells           |        1.6714  | IL7R, CCR7, LTB                    |
|         1 | B cells               |        2.53194 | MS4A1, CD79A, CD79B                |
|         2 | Monocytes             |        2.27141 | LYZ, S100A8, S100A9, FCGR3A, MS4A7 |
|         3 | NK cells              |        3.3747  | GNLY, NKG7, KLRD1                  |
|         4 | T cells               |        1.54448 | CD3D, CD3E, CD2, IL7R              |
|         5 | Monocytes             |        2.86696 | LYZ, S100A8, S100A9, FCGR3A, MS4A7 |
|         6 | Dendritic cells       |        3.15556 | FCER1A, CST3                       |
|         7 | Platelets             |        5.46271 | PPBP, PF4                          |

## 主要图表

- QC 分布：`/Users/olinaaa_/Documents/Statistic ML/outputs/figures/qc_histograms.png`
- PCA elbow：`/Users/olinaaa_/Documents/Statistic ML/outputs/figures/pca_elbow.png`
- UMAP 聚类图：`/Users/olinaaa_/Documents/Statistic ML/outputs/figures/umap_clusters.png`
- 指标对比：`/Users/olinaaa_/Documents/Statistic ML/outputs/figures/metric_summary.png`
- cluster size：`/Users/olinaaa_/Documents/Statistic ML/outputs/figures/cluster_sizes.png`
- marker heatmap：`/Users/olinaaa_/Documents/Statistic ML/outputs/figures/marker_heatmap.png`
- marker dotplot：`/Users/olinaaa_/Documents/Statistic ML/outputs/figures/marker_dotplot.png`

## 报告可写结论

- PCA 将高维表达矩阵压缩为低维连续表示，便于后续聚类和可视化。
- K-means 与层次聚类主要依赖欧氏距离和簇形状假设；Louvain 在细胞邻接图上寻找社区，更适合单细胞常见的非线性结构。
- 内部指标只能评价几何紧致性，最终解释需要结合 marker genes 和免疫细胞背景。
