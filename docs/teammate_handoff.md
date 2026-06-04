# 给同学 B 的交接说明

## 你会收到什么

这个项目已经在真实 PBMC3k 数据上跑完实验，并已经生成正式报告 Markdown 版。你不需要使用同学 A 的电脑，只需要拿到 GitHub 仓库即可查看结果、审核报告，必要时也可以在自己的电脑上复现。

## 最重要的文件

- `outputs/run_summary.md`：实验结果摘要，先读这个。
- `reports/final_report.md`：最终报告 Markdown 版，优先审核这个。
- `outputs/metrics.csv`：K-means、层次聚类、Louvain 的指标对比。
- `outputs/cluster_annotations.csv`：Louvain cluster 的 marker gene 注释建议。
- `outputs/figures/umap_clusters.png`：三类聚类方法在 UMAP 上的结果。
- `outputs/figures/marker_dotplot.png`：marker genes 对 cluster 注释的支持。
- `docs/a_handoff.md`：同学 A 的实验交接。

## 你负责的部分

你主要负责审核和完善报告里的生命科学解释：

- 解释 PBMC 是什么，为什么适合做单细胞聚类案例。
- 根据 marker genes 解释 Louvain `resolution=0.8` 的 8 个 cluster。
- 把 `CD3D/CD3E/IL7R/CCR7/LTB` 解释为 T/CD4 T 相关信号。
- 把 `MS4A1/CD79A/CD79B` 解释为 B cell 相关信号。
- 把 `LYZ/S100A8/S100A9/FCGR3A/MS4A7` 解释为 monocyte 相关信号。
- 把 `GNLY/NKG7/KLRD1` 解释为 NK cell 相关信号。
- 把 `FCER1A/CST3` 解释为 dendritic cell 相关信号。
- 把 `PPBP/PF4` 解释为 platelet 相关信号。

## 报告建议主线

1. 数据经过 QC 后保留 `2638` 个细胞和 `1838` 个 highly-variable genes。
2. K-means 和层次聚类在内部指标上以 `k=5` 表现较好，说明数据中存在若干较清晰的大类细胞群。
3. Louvain `resolution=0.8` 给出 8 个 cluster，更适合做细胞亚群解释。
4. UMAP 图展示了细胞群的几何分离；marker dotplot/heatmap 提供了生物学解释依据。
5. 内部聚类指标不能单独证明生物学正确性，因此最终结论要结合 marker genes。

## 如果你要复现

先安装依赖：

```bash
python -m pip install pandas scipy scikit-learn scanpy anndata umap-learn networkx tabulate seaborn h5py numba certifi
```

然后运行：

```bash
python scripts/run_pbmc_clustering.py --stability-repeats 5
```

如果只是审核报告，不需要复现，直接使用 `reports/final_report.md`、`outputs/` 里的图表和 CSV 即可。

## 协作方式

- 同学 A 后续负责解释代码、补跑实验、调整参数和技术结果核对。
- 同学 B 负责审核报告背景、marker gene 解释、图表文字说明。
- 两人共同确认最终报告里的图编号、结论和局限性。
