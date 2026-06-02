# PBMC 单细胞 RNA-seq 降维与聚类算法比较

## 摘要

本项目使用公开 PBMC3k 单细胞 RNA-seq 数据，比较 PCA/UMAP/t-SNE 降维与 K-means、层次聚类、Louvain 图聚类在免疫细胞亚群发现中的表现。评价采用内部聚类指标、稳定性分析和 marker gene 生物解释相结合的方式。

## 1. 背景与研究问题

PBMC 指外周血单个核细胞，包含 T cells、B cells、NK cells、monocytes 等免疫细胞群。单细胞 RNA-seq 数据具有高维、稀疏、噪声较强的特点，因此适合检验统计机器学习中的降维和聚类方法。

研究问题：经典统计机器学习方法能否在高维稀疏单细胞表达矩阵中恢复有生物意义的免疫细胞亚群？

## 2. 数据与预处理

- 数据来源：Scanpy PBMC3k / 10x Genomics 公开健康供体 PBMC 数据。
- 质量控制：过滤低基因数细胞、低出现频率基因、高线粒体比例细胞。
- 标准化：total-count normalization、`log1p` 转换、高变基因选择、scale。
- 降维：使用 PCA 得到后续聚类的低维表示。

插入图：`outputs/figures/qc_histograms.png`、`outputs/figures/pca_elbow.png`

## 3. 方法

### 3.1 降维

PCA 用于构建主分析表示；UMAP 和 t-SNE 用于二维可视化。PCA 保留线性方差结构，UMAP/t-SNE 更适合展示局部邻域结构，但二维图本身不作为唯一评价依据。

### 3.2 聚类

- K-means：扫描 `k=4..12`。
- 层次聚类：使用 Ward linkage，扫描 `k=4..12`。
- Louvain：在 KNN 细胞图上做社区发现，扫描 resolution `0.4/0.8/1.2`。

### 3.3 评价指标

使用 silhouette、Calinski-Harabasz、Davies-Bouldin 评价几何结构；使用随机子采样后的 ARI/NMI 评价稳定性；使用 marker gene 表达解释 cluster 的生物含义。

## 4. 实验结果

插入图：`outputs/figures/metric_summary.png`、`outputs/figures/cluster_sizes.png`、`outputs/figures/umap_clusters.png`

填写要点：

- 哪个 K-means `k` 的 silhouette 最高？
- Louvain resolution 改变时 cluster 数如何变化？
- K-means/层次聚类/Louvain 在 UMAP 上呈现的结构有什么差别？

## 5. 细胞群解释

插入图：`outputs/figures/marker_heatmap.png`、`outputs/figures/marker_dotplot.png`

根据 `outputs/cluster_annotations.csv` 填写每个主要 cluster 的注释：

| Cluster | 建议细胞类型 | 支持 marker genes | 解释 |
| --- | --- | --- | --- |
| 0 |  |  |  |

## 6. 讨论与局限

- 内部聚类指标评价的是几何结构，不等价于生物学正确性。
- marker gene 注释是启发式解释，需要结合领域知识或外部标注验证。
- 本项目未涉及 batch correction、trajectory inference 或多数据集泛化。
- Louvain 的 resolution、KNN 参数和 PCA 维数都会影响结果。

## 7. 分工说明

- 同学 A：实验 pipeline、聚类算法、指标与可视化。
- 同学 B：PBMC 背景、marker gene 解释、报告整合与讨论。

## 参考资料

- Scanpy PBMC3k tutorial / dataset documentation
- 课程 Lecture 08 Dimension Reduction
- 课程 Lecture 09 Clustering
- 课程 Lecture 10 Louvain Algorithm
