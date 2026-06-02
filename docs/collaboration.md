# 两人协作说明

## 角色分工

同学 A：机器学习实验负责人。

- 维护 `scripts/run_pbmc_clustering.py` 和 `src/pbmc_project/`
- 跑通 PBMC3k 下载、QC、PCA、UMAP/t-SNE、K-means、层次聚类、Louvain
- 维护 `outputs/metrics.csv`、图表和运行参数
- 在报告中写“方法”和“实验设置/评价指标”

同学 B：生命科学解释与报告负责人。

- 写 PBMC 和 scRNA-seq 背景
- 根据 `outputs/cluster_annotations.csv`、`outputs/marker_scores.csv` 解释 cluster
- 整理 marker genes 的生物含义和局限性
- 在报告中写“数据背景”“结果解释”“讨论与局限”

## 文件交接规则

- A 每次完整运行后，把输出目录压缩或同步给 B。
- 所有图表文件名保持稳定，不要手动改名。
- B 在报告中引用图时使用脚本生成的文件名，例如 `figures/umap_clusters.png`。
- 任何人工修改的结论都写在报告模板中，不直接改 `outputs/run_summary.md` 的自动生成事实。

## 每日同步模板

```text
日期：
A 今天完成：
B 今天完成：
新的图/表：
目前最重要的发现：
明天需要解决的问题：
```

## 最终检查清单

- 代码能从空输出目录重新运行。
- `metrics.csv` 有 K-means、Hierarchical、Louvain 三类方法。
- 报告中每个主要 cluster 至少有 1-2 个 marker gene 支持。
- 报告明确说明内部聚类指标不能替代生物学验证。
- 图中标题、坐标、图例可读，文件名与正文引用一致。
