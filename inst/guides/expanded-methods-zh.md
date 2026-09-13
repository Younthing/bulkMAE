# bulkMAE 0.4 方法扩展与选择指南

本文说明新增方法解决什么问题、需要什么数据，以及哪些方法虽然共享后端，
统计含义却不同。所有 assay 级函数仍以 MAE 为公共输入，返回后端原生对象，
不记录运行状态，也不修改 MAE。

## 1. 新增覆盖

| 问题 | 函数 | 后端 | 主要输入 |
|---|---|---|---|
| 一致性分型 | `cluster_consensus()`、`cluster_consensus_diagnostics()`、`cluster_consensus_classes()` | ConsensusClusterPlus | 连续、已变换的基因 × 样本矩阵 |
| 非负分解分型 | `cluster_nmf()`、`cluster_nmf_classes()` | NMF | 非负连续矩阵 |
| 差异共表达（DCA） | `coexpr_differential()` | diffcoexp | 两组连续表达矩阵 |
| 模块保守性 | `coexpr_modules()`、`coexpr_preservation()` | WGCNA | reference/test 两个可比 assay 和参考模块 |
| 连续表达差异 | `de_limma()` | limma | log-expression、芯片或通路分数 |
| 差异剪接/使用 | `dtu_diffsplice()` | limma 或 edgeR | 已拟合的 exon/transcript 模型 |
| 时间/剂量反应 | `de_masigpro_design()`、`de_masigpro()` | maSigPro | 表达矩阵和真实 time/replicate/group 元数据 |
| 去批次可视化 | `adjust_batch()` | limma | 连续表达矩阵和显式保留设计 |
| 非线性样本嵌入 | `reduce_umap()`、`reduce_tsne()` | uwot、Rtsne | 变换后表达矩阵 |
| 偏倚校正 ORA | `enrich_goseq()` | goseq | 全部受检基因的 DE 二元向量和长度/丰度 |
| KEGG/Reactome | `enrich_kegg()`、`enrich_reactome()`，以 `method` 选择 ORA/GSEA | clusterProfiler、ReactomePA | 显式 ID 和物种 |
| 竞争/旋转集合检验 | `enrich_camera()`、`enrich_fry()`、`enrich_roast()` | limma | 表达矩阵、设计和 contrast |
| 稳健样本 signature | `score_singscore()` | singscore | 每个样本的基因秩 |
| 多方法活性推断 | `activity_decouple()` | decoupleR | 表达矩阵和长表网络/集合 |
| 生存曲线/判别 | `surv_km()`、`surv_roc()`、`surv_cox_table()`、`surv_risk_groups()`、`plot_surv_km()`、`plot_surv_forest()`、`plot_surv_risk()` | survival、timeROC | 对齐的随访、结局和 marker；无真实随访时应标明合成诊断结局 |

## 2. 富集和活性方法如何选择

| 科学问题 | 推荐入口 | 单位 | 是否使用表达设计 | 是否使用网络权重/方向 |
|---|---|---|---|---|
| 入选基因是否过多落在某集合 | `enrich_ora()`；`enrich_go()` / `enrich_kegg()` / `enrich_reactome(method = "ora")` | 一个基因列表 | 否 | 否 |
| RNA-seq 入选概率受长度/丰度影响 | `enrich_goseq()` | 一个二元选择向量 | 否 | 偏倚权重，不是调控网络 |
| 全基因排序两端是否富集 | `enrich_fgsea()` / `enrich_gsea()`；数据库函数使用 `method = "gsea"` | 一个排序统计量 | 否 | 否 |
| 在设计矩阵下做竞争集合检验 | `enrich_camera()` | 整个研究 | 是 | 考虑集合内相关性 |
| 在设计矩阵下做旋转集合检验 | `enrich_fry()` / `enrich_roast()` | 整个研究 | 是 | 否 |
| 每个样本得到通路分数 | `score_gsva()` / `score_ssgsea()` | 样本 | 否 | 否 |
| 每个样本得到秩稳健 signature | `score_singscore()` | 样本 | 否 | 可有 up/down 方向 |
| 用调控边推断 TF/通路活性 | `activity_decouple()` | 样本或统计量 | 依方法而定 | 是 |

ORA、GSEA、GSVA 和调控活性不能仅因都输出“pathway score”就视为同一种分析。
它们的零假设、分析单位和所需输入不同。

## 3. decoupleR 的方法边界

`activity_decouple()` 的 `statistics` 直接传给官方 `decouple()`：

| 方法族 | 常见名称 | 适用解释 |
|---|---|---|
| 集合/排序统计 | `aucell`、`fgsea`、`gsva`、`ora` | 目标集合在样本或排序中的富集 |
| 回归/线性活性 | `mlm`、`ulm` | 有方向、可带权的 regulator-target 活性 |
| 调控活性 | `viper` | regulon 模式驱动的活性 |
| 加权汇总 | `wmean`、`wsum` | 用网络权重汇总目标信号 |
| 共识 | `consensus = TRUE` | 汇总所选方法；不是新的生物学数据库 |

用 `activity_methods()` 查看安装版本实际支持的方法。用
`activity_resource()` 获取官方资源时，仍要审查 source、target、weight/
mor 列和物种。调用加权/有符号方法时，用 `mor = "weight"`（或实际列名）
显式映射为 decoupleR 的规范 `mor` 列。`activity_progeny()`、`activity_tf()`
只是常用默认组合，
不是 decoupleR 能力的全部。

## 4. 聚类、DCA 和模块保守性

一致性聚类应先运行 `cluster_consensus()`，再用
`cluster_consensus_diagnostics()` 查看 cluster/item consensus，并结合 CDF、样本稳定性、
外部队列和临床/生物学解释选择 k。`cluster_consensus_classes()` 只提取后端已有的
`consensusClass`，不重新聚类。

本包中的 DCA 明确定义为 differential co-expression analysis：比较两组的
基因相关结构。它不是差异表达，也不是 clinical decision curve analysis。
基因对数量随特征数近似平方增长，因此应先做与结局无关、可复现的特征过滤。

模块保守性使用两个显式 MAE。包装器取 reference assay、test
assay 和已命名的参考模块颜色的三方基因交集，统一顺序后构造 WGCNA
multi-set 输入。`coexpr_modules()` 可从 reference 拟合中提取具名模块颜色；
reference/test 角色和真实独立队列仍由调用者决定。

`coexpr_wgcna()` 遵循 `goodSamplesGenes()` 的质控结果：不合格的样本和特征会
自动剔除，不要求中途交互确认。函数发出警告，并在返回的 WGCNA 原生
list 中附加质控对象及被剔除的样本/特征名，便于事后审核。

## 5. 数据尺度

- 原始整数 counts：DESeq2、edgeR，以及 `counts = TRUE` 的 maSigPro。
- counts 加均值-方差权重：limma-voom、dream。
- 连续/近似同方差表达：PCA、相关性、聚类、DCA、WGCNA、limma、CAMERA、
  GSVA、singscore、decoupleR。
- TPM-like 非 log 表达：部分组织去卷积方法；以相应后端要求为准。
- corrected matrix：优先用于可视化或明确要求校正矩阵的算法；差异分析通常应
  把 batch 放进设计，而不是替换原 counts。

## 6. 尚未统一成单函数的内容

- 样本错标需要表达相似性、基因型/性别标记和实验记录联合判断，不能自动改名。
- 外部验证需要冻结模型、特征、系数和 cutoff，并保持训练/验证队列边界。
- CIBERSORTx 是外部服务/容器；本包只准备输入，不自动上传受保护数据。
- 商业数据库、需要账户授权的资源和数据库快照不由包装器静默下载。

每个实现对应的官方文档、版本审查日期和参数映射见
[official-sources.md](official-sources.md)；发布前操作见
[review-checklist.md](review-checklist.md)。0.3 及更早版本升级时另见
[naming-migration-zh.md](naming-migration-zh.md)。逐函数输入可达性及在线/外部边界见
[input-completeness-zh.md](input-completeness-zh.md)。
