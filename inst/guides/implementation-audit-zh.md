# bulkMAE 实现审查说明

## 1. 核心约束

本包只有一个公共数据契约：分析函数接收 `MultiAssayExperiment`、明确的
`experiment` 和必要时明确的 `assay`。内部使用 `getWithColData()` 解析
`sampleMap`，然后把对齐后的矩阵和样本信息交给后端包。

以下内容没有实现，也不应在后续版本中悄悄加入：

- 不在 MAE 中保存分析成功、失败或历史记录；
- 不设置“当前 assay”或全局活动数据集；
- 不修改输入 MAE；
- 不把 DESeq2、edgeR、limma 等结果统一塞进新结果类；
- 不根据自然语言标签自动猜测设计矩阵或对照方向。

## 2. 文件职责

| 文件 | 职责 | 主要后端 |
|---|---|---|
| `access.R` | MAE 校验、样本映射、assay/metadata 提取、内部断言 | MAE、SE |
| `import.R` | 构建 SE/MAE、导入转录本定量 | tximport |
| `annotation.R` | ID 去版本、注释映射、Ensembl 查询 | AnnotationDbi、biomaRt |
| `qc.R` | 低表达过滤、library 指标、相关性、PCA/MDS/UMAP/t-SNE、离群筛查 | edgeR、limma、uwot、Rtsne |
| `preprocessing.R` | TMM、size factor、VST/rlog/voom、ComBat、SVA、RUV、方差分解 | edgeR、DESeq2、sva、limma、RUVSeq、variancePartition |
| `differential.R` | 常规和复杂设计、连续表达、差异剪接、时间序列/剂量反应 | DESeq2、edgeR、limma、dream、maSigPro |
| `enrichment.R` | ORA/GSEA、GO/KEGG/Reactome、竞争/旋转检验、样本打分、调控活性 | clusterProfiler、fgsea、ReactomePA、limma、GSVA、singscore、decoupleR |
| `network.R` | 一致性聚类、差异共表达、NMF、WGCNA/模块保守性、调控网络、PPI | ConsensusClusterPlus、diffcoexp、NMF、WGCNA、GENIE3、STRINGdb |
| `deconvolution.R` | 免疫去卷积、参考单细胞去卷积、外部输入准备 | immunedeconv、MuSiC、BayesPrism、CIBERSORTx |
| `clinical.R` | signature、KM/Cox、惩罚模型、time ROC、meta-analysis、LINCS | survival、glmnet、timeROC、metafor、signatureSearch |

## 3. 原需求覆盖判断

### 直接实现

- 导入、注释、过滤、library size、PCA/MDS/UMAP/t-SNE、相关性和离群样本筛查；
- TMM、DESeq2 size factor、voom、VST/rlog、ComBat/ComBat-seq、SVA、
  RUVg 和方差分解；
- DESeq2、edgeR、limma/limma-voom、dream、maSigPro、差异剪接；
- ORA/GSEA、GO/KEGG/Reactome、goseq 长度偏倚校正、fgsea、
  CAMERA/FRY/mroast、GSVA/ssGSEA、singscore，以及 decoupleR 的多种集合
  统计与调控活性推断；
- 一致性聚类及 ICL、DCA 差异共表达、NMF、WGCNA/模块保守性、STRING
  PPI、GENIE3；
- immunedeconv 所覆盖的 EPIC、quanTIseq、xCell 等方法，及 MuSiC、
  BayesPrism；
- signature 打分、KM、Cox/惩罚 Cox、time-dependent ROC、meta-analysis、
  LINCS 搜索。

### 由模型公式覆盖

配对、多因素、交互和线性连续变量不需要分别创建函数。
它们由 `design`/`formula` 和显式 contrast 表达；重复测量由 `de_dream()` 的
随机效应公式处理。多组时间序列和剂量反应另有 `de_masigpro()`，其设计表
必须按 maSigPro 官方规范显式提供。

### decoupleR 的归类

`decoupleR` 不是单一“富集包”。`activity_decouple()` 原样暴露其统一引擎：

- `aucell`、`fgsea`、`gsva`、`ora` 更接近集合富集或样本集合评分；
- `mlm`、`ulm`、`viper`、`wmean`、`wsum` 使用带方向/权重的网络推断活性；
- 默认统计量和 consensus 由安装的 decoupleR 版本决定，可用
  `activity_methods()` 审查。

因此包中保留 `activity_progeny()` 和 `activity_tf()` 作为常用快捷入口，
同时提供 `activity_decouple()`，避免把 PROGENy、DoRothEA/CollecTRI 和普通 gene
set 混成同一个概念。

### 只准备输入

`deconv_cibersortx_input()` 只生成 CIBERSORTx 要求的表。它不上传数据，
也不调用外部服务。这是有意的权限和可审查边界。

### 暂不伪实现

- 样本“错标”只能通过表达相似性和实验记录联合判断；当前提供 PCA、相关性和
  异常筛查，不自动改标签。
- DCA 在本包中特指 differential co-expression analysis，由维护中的
  Bioconductor `diffcoexp` 实现。它不自动猜测分组和阈值。
- `coexpr_wgcna()` 自动应用 `goodSamplesGenes()` 的剔除结果，以警告和
  返回 list 中的剔除名单保留审计线索，不进行中途交互。
- `coexpr_preservation()` 负责 reference、test 和模块标签的三方
  基因对齐并调用 WGCNA；reference/test 角色及参考模块标签仍必须
  显式传入。
- 外部验证不是单个拟合函数。当前提供固定权重 signature 和原生模型对象，
  训练/验证队列边界应保留在研究脚本中。

## 4. 依赖策略

`Imports` 只包含 MAE/SE 基础设施和 R 基础包。每个分析后端在调用相应函数时
才检查，因此用户可以按项目安装 DESeq2、edgeR、WGCNA 等，而不必维护一个
包含全部方法的巨大环境。

MuSiC、immunedeconv 和 BayesPrism 已作为可选依赖声明，但不在每个标准
CRAN/Bioconductor 仓库中可用。标准后端与三个官方 GitHub 后端统一按
[pak 依赖安装指南](dependency-installation-zh.md) 安装，并在项目锁文件或容器
中记录 commit/release。正式提交 Bioconductor 前还需确认构建系统能否解析
这些可选依赖，否则应将其拆到 companion package。

`import_tximport()` 会把 `counts`、`abundance`、`length` 和
`countsFromAbundance` 模式保留在 SE 中。DESeq2 与 edgeR 包装器分别重建官方
要求的 tximport 列表并调用 `DESeqDataSetFromTximport()` 或
`DGEListFromTximport()`，不会把 estimated counts 冒充普通整数计数。

## 5. 发布验证要求

本节定义发布门槛，不报告某一次工作目录或临时环境的通过状态。每个候选发布
必须在干净、版本互相兼容的 R/Bioconductor 库中按
[依赖安装指南](dependency-installation-zh.md) 重建环境，并保存可复核的日志。

发布记录至少应包含：

- roxygen 文档和 `NAMESPACE` 的重建差异；
- testthat、`R CMD check --as-cran` 和 `BiocCheck` 的完整日志；
- R、Bioconductor、所有实际运行后端以及 GitHub commit 的版本信息；
- 每个离线后端的小型模拟数据成功路径；
- 使用官方示例数据完成的后端烟雾测试；
- 在线、授权或外部数据库测试的执行日期、资源版本以及运行/跳过状态。

只有发布记录中存在对应日志时，才能声明某个检查或后端已经通过。缺少依赖、
网络、凭据或授权导致的跳过必须明确列出，不能合并计入通过数量。具体项目按
[review-checklist.md](review-checklist.md) 执行。
