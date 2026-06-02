# PBMC 单细胞 RNA-seq 降维与聚类算法比较

这个项目把统计机器学习期末选题落成可复现实验：使用公开 PBMC3k 单细胞 RNA-seq 数据，比较 PCA、UMAP、t-SNE 表示下的 K-means、层次聚类和 Louvain 图聚类，并用免疫细胞 marker genes 做结果解释。

## 课程对应

- 降维：PCA、t-SNE、UMAP，对应 `main_0509.pdf`
- 聚类：K-means、层次聚类，对应 `main_0511.pdf`
- 图聚类：Louvain，对应 `main_0518_1.pdf` 和 `supp_0525.pdf`

## 推荐环境

Scanpy 生态对 Python 版本比较敏感，建议使用 Python 3.11、3.12 或当前已验证可安装的 Python 3.13。

```bash
python3.11 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -e .
```

如果有 conda/mamba：

```bash
conda env create -f environment.yml
conda activate pbmc-sml
```

## 一键运行

```bash
python scripts/run_pbmc_clustering.py
```

常用快速调试：

```bash
python scripts/run_pbmc_clustering.py --max-cells 800 --skip-tsne --stability-repeats 2
```

如果网络暂时无法下载 Scanpy PBMC3k 数据，可以先跑离线模拟数据，验证完整流程和图表生成：

```bash
python scripts/run_pbmc_clustering.py --source synthetic --max-cells 800 --skip-tsne --stability-repeats 2 --skip-save-adata
```

真实数据下载如果在 Codex 沙盒里被网络权限拦住，可以在普通终端运行：

```bash
python scripts/download_pbmc3k.py
python scripts/run_pbmc_clustering.py
```

脚本会优先使用本地 `data/pbmc3k_raw.h5ad`，不再重复下载。

如果已经有本地 `.h5ad` 文件：

```bash
python scripts/run_pbmc_clustering.py --input-h5ad path/to/pbmc3k.h5ad
```

主要输出：

- `outputs/metrics.csv`：聚类指标与稳定性结果
- `outputs/cluster_assignments.csv`：每个细胞的聚类标签
- `outputs/cluster_annotations.csv`：基于 marker genes 的 cluster 注释建议
- `outputs/marker_scores.csv`、`outputs/marker_gene_means.csv`：marker 解释表
- `outputs/figures/`：QC、PCA elbow、UMAP/t-SNE、指标对比、cluster size、marker dotplot/heatmap
- `outputs/run_summary.md`：自动生成的实验摘要，可直接搬进报告初稿

说明：项目里的 Louvain 使用 Scanpy 生成的细胞 KNN 图，并通过 NetworkX 的 Louvain community detection 计算社区，避免安装较难编译的 `louvain` C 扩展。

## 协作建议

同学 A 负责实验 pipeline 和指标，保证脚本从空环境能跑通；同学 B 负责 PBMC 背景、marker gene 注释和报告叙事。每次 A 产出新图时，把图名、参数、它说明什么写进 `outputs/run_summary.md` 或共享文档，B 用同一套文件名整理报告，避免来回找图。

## 报告入口

报告骨架在 `reports/final_report_template.md`，分工细则在 `docs/collaboration.md`。
