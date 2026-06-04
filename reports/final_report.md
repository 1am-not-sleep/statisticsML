# PBMC 单细胞 RNA-seq 降维与聚类算法比较

作者：同学 A、同学 B  
课程：统计机器学习  
项目类型：算法比较与案例应用  
数据集：PBMC3k 单细胞 RNA-seq 公开数据  
版本：最终报告 Markdown 版

## 摘要

单细胞 RNA-seq 数据通常具有高维、稀疏、噪声较强等特点，是检验降维与聚类方法的典型应用场景。本项目使用公开 PBMC3k 数据，围绕“经典统计机器学习方法能否从高维单细胞表达矩阵中发现具有生物意义的免疫细胞亚群”这一问题，比较 PCA、UMAP、t-SNE 降维表示下的 K-means、层次聚类和 Louvain 图聚类。数据经过质量控制、归一化、`log1p` 转换、高变基因筛选、标准化和 PCA 后，保留 `2638` 个细胞和 `1838` 个 highly-variable genes/features。实验采用内部聚类指标、随机子采样稳定性和 marker gene 生物解释相结合的评价方式。结果显示，K-means 与层次聚类在 `k=5` 时取得较好的内部指标和稳定性；Louvain 在 `resolution=0.8` 时得到 8 个 cluster，更适合进行 PBMC 主要免疫细胞亚群解释。基于 marker genes，主要 cluster 可解释为 CD4 T cells、T cells、B cells、Monocytes、NK cells、Dendritic cells 和 Platelets。整体来看，统计机器学习中的降维与聚类方法可以从 PBMC 单细胞表达矩阵中恢复有意义的细胞群结构，但最终解释不能只依赖几何聚类指标，仍需结合生物 marker 和领域知识。

关键词：单细胞 RNA-seq；PBMC；PCA；UMAP；t-SNE；K-means；层次聚类；Louvain；marker genes

## 1. 背景与研究问题

外周血单个核细胞（PBMC, peripheral blood mononuclear cells）是血液中具有单个细胞核的免疫细胞集合，常见细胞类型包括 T cells、B cells、NK cells、monocytes 和 dendritic cells 等。PBMC 是单细胞 RNA-seq 分析中常用的数据对象，因为它包含多个相对明确的免疫细胞群，不同细胞类型通常具有可解释的 marker gene 表达模式。

单细胞 RNA-seq 的原始表达矩阵通常包含数千个细胞和上万个基因。对于每个细胞而言，基因表达向量维度很高，并且存在测序深度差异、dropout、低表达噪声和稀疏性等问题。如果直接在原始基因空间中聚类，距离度量容易受到噪声和维数灾难影响。因此，实际分析中通常先进行质量控制和特征筛选，再使用 PCA 等方法构建低维表示，最后在低维表示或细胞近邻图上进行聚类。

本项目的核心研究问题是：经典统计机器学习方法能否在高维稀疏的 PBMC 单细胞表达矩阵中发现具有生物意义的免疫细胞亚群？为回答该问题，我们比较 K-means、层次聚类和 Louvain 三类聚类方法，并结合 marker genes 对聚类结果进行生物解释。

## 2. 数据与预处理

本项目使用 Scanpy 提供的 PBMC3k 公开数据。该数据来自 10x Genomics 的健康供体 PBMC 单细胞 RNA-seq 实验。项目仓库中不上传 `.h5ad` 原始数据和处理后的 AnnData 文件，因为这类文件体积较大；实验脚本会优先读取本地 `data/pbmc3k_raw.h5ad`，在新环境中也可以通过下载脚本重新获取并复现实验。

预处理流程如下：

1. 读取 PBMC3k 原始表达矩阵。
2. 进行质量控制，过滤低基因数细胞、低出现频率基因和线粒体比例异常细胞。
3. 进行 total-count normalization，使不同细胞之间的测序深度更可比。
4. 对表达值做 `log1p` 转换，降低极端表达值对后续分析的影响。
5. 筛选 highly-variable genes，保留对细胞间差异贡献较大的基因。
6. 对表达矩阵进行标准化，并使用 PCA 得到后续聚类和邻接图构建所需的低维表示。

预处理后，数据规模为 `2638` cells × `1838` highly-variable genes/features。该规模保留了主要细胞群结构，同时减少了原始高维表达矩阵中的噪声。

![QC histograms](../outputs/figures/qc_histograms.png)

**图 1. QC 指标分布。** 该图展示细胞层面的质量控制指标分布，用于检查每个细胞检测到的基因数、总 counts 和线粒体比例等指标是否存在明显异常。质量控制的作用是减少低质量细胞和技术噪声对后续降维与聚类结果的影响。

![PCA elbow](../outputs/figures/pca_elbow.png)

**图 2. PCA elbow 图。** PCA 将高维基因表达矩阵压缩为低维连续表示。该图展示主成分解释方差的下降趋势，用于观察前若干主成分是否保留了主要变异信息。后续聚类和 KNN 图构建均基于 PCA 表示进行。

## 3. 方法

### 3.1 降维方法

PCA 是本项目的主分析表示。PCA 通过线性正交变换寻找数据方差最大的方向，适合在聚类前降低维度、去除部分噪声，并提高算法运行效率。PCA 与课程中的降维内容直接对应，其核心思想是用较少的主成分近似原始高维数据结构。

UMAP 和 t-SNE 主要用于二维可视化。二者都强调局部邻域结构，可以帮助观察细胞群在二维空间中的分离情况。但二维可视化会受到参数、初始化和投影失真的影响，因此本项目不把 UMAP 或 t-SNE 图作为唯一评价依据，而是将其作为解释聚类结果的辅助工具。

### 3.2 聚类方法

K-means 将细胞划分为 `k` 个簇，并最小化样本到簇中心的平方距离。该方法简单、可解释、运行较快，但默认簇结构接近球形，且需要预先指定 `k`。本项目扫描 `k=4..12`。

层次聚类使用 Ward linkage 构造层次结构，再切分为指定数量的 cluster。与 K-means 相比，层次聚类可以体现簇之间的合并关系，但在样本量较大时计算成本更高。本项目同样扫描 `k=4..12`。

Louvain 是图聚类方法。它先基于细胞在 PCA 空间中的近邻关系构造 KNN 图，再通过最大化 modularity 寻找社区结构。单细胞数据中的细胞类型常表现为局部邻域结构或非线性流形，因此 Louvain 常用于细胞群发现。本项目扫描 `resolution=0.4/0.8/1.2`。

### 3.3 评价方式

本项目采用三类评价方式：

1. **内部聚类指标。** 使用 silhouette、Calinski-Harabasz 和 Davies-Bouldin 评价几何结构。silhouette 和 Calinski-Harabasz 越高通常表示簇内更紧密、簇间更分离；Davies-Bouldin 越低通常表示聚类结构更好。
2. **稳定性指标。** 通过随机子采样重复聚类，计算 ARI 和 NMI 的均值与标准差。ARI/NMI 越高、标准差越小，说明结果对采样扰动更稳定。
3. **生物解释。** 使用 marker genes 检查 cluster 是否对应已知免疫细胞类型。该部分用于弥补无监督内部指标无法提供真实生物标签的不足。

## 4. 实验结果

### 4.1 聚类指标比较

主要指标结果如下表所示：

| 方法与参数 | cluster 数 | Silhouette | Calinski-Harabasz | Davies-Bouldin | ARI mean | NMI mean | 解释 |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | --- |
| K-means, `k=5` | 5 | 0.222988 | 410.598 | 1.52024 | 0.995817 | 0.991943 | 内部指标和稳定性表现最好之一，适合概括较大的细胞群结构。 |
| Hierarchical, `k=5` | 5 | 0.220512 | 395.558 | 1.49547 | 0.984193 | 0.971546 | 与 K-means `k=5` 接近，说明 5 个大类结构较稳定。 |
| Louvain, `resolution=0.4` | 6 | 0.198307 | 334.402 | 1.70690 | 0.886134 | 0.911636 | 得到较粗粒度图社区。 |
| Louvain, `resolution=0.8` | 8 | 0.141799 | 285.691 | 2.06665 | 0.873777 | 0.863860 | 本项目用于最终 marker gene 注释，细胞群解释粒度更合适。 |
| Louvain, `resolution=1.2` | 11 | 0.051901 | 208.148 | 3.51947 | 0.604979 | 0.726297 | cluster 更多，但几何紧致性和稳定性下降。 |

从内部指标看，K-means 和层次聚类在 `k=5` 时表现较好。`kmeans_k5` 的 silhouette 为 `0.222988`，ARI 均值为 `0.995817`，NMI 均值为 `0.991943`；`hierarchical_k5` 的 silhouette 为 `0.220512`，ARI 均值为 `0.984193`，NMI 均值为 `0.971546`。这说明当聚成 5 个较大的群时，PBMC3k 数据在 PCA 表示中具有较清晰的几何结构，并且 K-means 与层次聚类的结果相对稳定。

对于 Louvain，`resolution=0.4` 得到 6 个 cluster，`resolution=0.8` 得到 8 个 cluster，`resolution=1.2` 得到 11 个 cluster。随着 resolution 增大，cluster 数量增加，但 silhouette 和稳定性指标下降。这说明更高 resolution 会带来更细的切分，但不一定带来更好的几何聚类质量。综合可解释性和粒度，本项目选择 `louvain_r0_8` 作为最终 marker gene 注释对象。

![Metric summary](../outputs/figures/metric_summary.png)

**图 3. 聚类指标对比。** 该图汇总不同算法和参数下的内部指标与稳定性指标。K-means 和层次聚类在 `k=5` 时表现较好；Louvain 的优势不在于 silhouette 最高，而在于能在细胞近邻图上给出更适合生物解释的社区结构。

![Cluster sizes](../outputs/figures/cluster_sizes.png)

**图 4. Cluster size 分布。** 该图展示不同聚类结果下各 cluster 的细胞数量。cluster size 可用于检查是否存在过多极小 cluster 或单一 cluster 过大的问题。PBMC 数据中某些细胞类型本身数量较少，因此小 cluster 不一定错误，需要结合 marker genes 判断其生物意义。

### 4.2 UMAP 与 t-SNE 可视化

UMAP 和 t-SNE 用于观察不同聚类标签在二维空间中的分布。可视化结果显示，PBMC3k 数据存在若干较明显的细胞群区域。K-means 和层次聚类更倾向于将数据划分为几个较大的几何区域；Louvain 基于近邻图寻找社区，因此能进一步拆分出一些较小但具有 marker gene 支持的细胞群。

![UMAP clusters](../outputs/figures/umap_clusters.png)

**图 5. UMAP 空间中的聚类结果。** 该图展示不同聚类方法的标签在 UMAP 空间中的分布。二维图不能替代定量指标，但可以直观呈现不同算法对细胞群边界的划分差异。

![t-SNE clusters](../outputs/figures/tsne_clusters.png)

**图 6. t-SNE 空间中的聚类结果。** t-SNE 更强调局部邻域，因此适合辅助观察局部群落结构。由于 t-SNE 的全局距离解释有限，本项目仅将其作为可视化补充。

## 5. 细胞群解释

最终生物解释基于 `louvain_r0_8`。该设置得到 8 个 cluster，既比 `resolution=0.4` 更细，也没有像 `resolution=1.2` 那样产生过多小 cluster。marker gene 注释结果如下：

| Cluster | 建议细胞类型 | 支持 marker genes | 生物解释 |
| --- | --- | --- | --- |
| 0 | CD4 T cells | `IL7R`, `CCR7`, `LTB` | `IL7R`、`CCR7` 和 `LTB` 常见于 CD4 T/naive 或 memory T cell 相关群体，说明该 cluster 具有 CD4 T cell 特征。 |
| 1 | B cells | `MS4A1`, `CD79A`, `CD79B` | `MS4A1` 编码 CD20，`CD79A/CD79B` 参与 B cell receptor 复合体，三者共同支持 B cells 注释。 |
| 2 | Monocytes | `LYZ`, `S100A8`, `S100A9`, `FCGR3A`, `MS4A7` | `LYZ`、`S100A8`、`S100A9` 与髓系/炎症相关 monocyte 信号相关，`FCGR3A`、`MS4A7` 进一步支持 monocyte 解释。 |
| 3 | NK cells | `GNLY`, `NKG7`, `KLRD1` | `GNLY`、`NKG7` 和 `KLRD1` 与细胞毒性和 NK cell 功能相关，支持 NK cells 注释。 |
| 4 | T cells | `CD3D`, `CD3E`, `CD2`, `IL7R` | `CD3D/CD3E/CD2` 是 T cell lineage 的典型信号，说明该 cluster 主要为 T cells；与 cluster 0 相比，该 cluster 更适合作为泛 T cell 群解释。 |
| 5 | Monocytes | `LYZ`, `S100A8`, `S100A9`, `FCGR3A`, `MS4A7` | 该 cluster 与 cluster 2 类似，同样具有明显 monocyte marker 表达，可能对应 monocyte 的不同亚群或状态；本项目不进一步强行细分。 |
| 6 | Dendritic cells | `FCER1A`, `CST3` | `FCER1A` 和 `CST3` 支持 dendritic cell 相关解释，说明 Louvain 能识别数量较少但 marker 明确的细胞群。 |
| 7 | Platelets | `PPBP`, `PF4` | `PPBP` 和 `PF4` 是 platelet 相关 marker，支持该 cluster 注释为 platelets。 |

从表中可以看到，T cell、B cell、NK cell、monocyte、dendritic cell 和 platelet 等主要 PBMC 组成均能通过 marker genes 得到解释。这说明 Louvain `resolution=0.8` 的聚类结果不仅具有图结构意义，也具有较明确的生物学含义。需要注意的是，cluster 2 和 cluster 5 都被解释为 monocytes，说明聚类可能捕捉到了 monocyte 内部状态差异，但仅凭当前 marker 表不适合进一步给出过细的生物学命名。

![Marker heatmap](../outputs/figures/marker_heatmap.png)

**图 7. Marker gene heatmap。** 该图展示 marker genes 在不同 Louvain cluster 中的平均表达模式。不同 cluster 对应的 marker 表达具有一定特异性，说明聚类结果不仅反映几何结构，也能对应可解释的免疫细胞类型。

![Marker dotplot](../outputs/figures/marker_dotplot.png)

**图 8. Marker gene dotplot。** 该图同时展示 marker genes 在各 cluster 中的表达强度和表达细胞比例。dotplot 是本报告进行细胞类型注释的主要依据之一。

## 6. 讨论

本项目结果表明，PCA、K-means、层次聚类和 Louvain 等统计机器学习方法可以在 PBMC3k 单细胞数据中发现具有一定生物意义的细胞群结构。K-means 和层次聚类在 `k=5` 时内部指标较好，说明数据中存在若干较明显的大类细胞群。Louvain `resolution=0.8` 虽然 silhouette 不是最高，但能得到 8 个更细的 cluster，并且这些 cluster 可以通过 marker genes 解释为不同免疫细胞类型。因此，在单细胞分析中，内部几何指标和生物解释之间需要平衡。

不同聚类算法的目标函数不同，因此不能只用一个指标决定“最好”的算法。K-means 假设簇围绕中心分布，适合捕捉较规则的大类结构；层次聚类可以提供层次关系；Louvain 依赖细胞邻接图，更适合局部社区发现。对于 PBMC 这类具有已知细胞类型组成的数据，Louvain 的结果在生物解释上更实用。

从课程内容角度看，本项目覆盖了降维、聚类和图社区发现三个核心主题。PCA 解决高维表达矩阵的表示问题；UMAP/t-SNE 辅助可视化局部结构；K-means 和层次聚类代表经典聚类方法；Louvain 则展示了图模型思想在单细胞数据中的应用。该案例能够把课程中的统计机器学习方法连接到真实生命科学数据分析任务。

## 7. 局限性与改进方向

第一，本项目是无监督学习任务，没有把外部人工标注作为唯一真实标签，因此无法用 accuracy 直接评价结果。内部指标只衡量几何结构，不等价于生物学正确性。

第二，marker gene 注释具有启发式性质。虽然主要 marker 支持当前注释，但更严格的生物学结论需要结合人工注释、参考图谱或差异表达检验。

第三，本项目只使用 PBMC3k 单个数据集，没有处理批次效应，也没有验证模型在其他供体或其他测序批次上的泛化能力。

第四，本项目没有扩展到 trajectory inference、batch correction 或多组学整合。这些方向可以作为后续工作，但不适合作为本课程大作业的主要范围，否则容易偏离降维与聚类的核心主题。

后续可以进一步改进的方向包括：引入参考标注计算 ARI/NMI 外部评价；对 marker genes 做差异表达检验；比较 Leiden 与 Louvain；在其他 PBMC 数据集上验证参数选择的稳定性。

## 8. 分工说明

同学 A 主要负责实验 pipeline、环境配置、真实数据下载、预处理、聚类算法实现、指标计算和结果可视化。具体产出包括运行脚本、`outputs/metrics.csv`、`outputs/run_summary.md`、UMAP/t-SNE 图、指标图和 marker 表格。

同学 B 主要负责 PBMC 背景、marker gene 解释、细胞类型注释、图注撰写、讨论部分和报告整合。最终报告中第 1 节背景、第 5 节细胞群解释、第 6 节讨论和图注部分主要对应同学 B 的工作。

## 9. 结论

本项目基于真实 PBMC3k 单细胞 RNA-seq 数据，完成了从数据预处理、降维、聚类比较到 marker gene 生物解释的完整流程。实验结果显示，K-means 和层次聚类在 `k=5` 时具有较好的内部指标和稳定性；Louvain `resolution=0.8` 能获得更适合细胞亚群解释的 8 个 cluster。结合 marker genes，主要 cluster 可以解释为 T cells、CD4 T cells、B cells、Monocytes、NK cells、Dendritic cells 和 Platelets。整体而言，经典统计机器学习方法能够在高维单细胞表达矩阵中恢复有意义的免疫细胞群结构，但最终解释需要同时考虑聚类指标、稳定性和生物 marker。

## 参考资料

- Scanpy PBMC3k dataset documentation and tutorial.
- 10x Genomics PBMC3k public dataset.
- 课程 Lecture 08：Dimension Reduction。
- 课程 Lecture 09：Clustering。
- 课程 Lecture 10：Louvain Algorithm。
- Louvain 补充材料：`supp_0525.pdf`。
