# bulkMAE 输入完整性与资源可达性指南

本文回答一个具体问题：从“基因 × 样本”表达矩阵和样本元数据出发，哪些步骤可以由
`bulkMAE` 在本地完成，哪些生物学输入、公共数据库或外部服务仍必须另行取得。
这里的“完整”指公共函数能够构造、校验、按名称对齐并把输入交给相应后端；不表示
软件能够从表达值中创造实验事实、参考知识或临床证据。

本文覆盖最初的 84 个公共 API，并纳入当前 27 个新增的输入、结果和资源桥接函数。函数的
后端包须预先安装；项目约定缺失依赖用 `pak` 安装，具体命令见
[dependency-installation-zh.md](dependency-installation-zh.md)。首次安装软件或下载资源
可能需要联网，因此“可离线”均指软件及所需本地资源已经安装或缓存之后。各函数
采用的官方一手资料及适配边界见 [official-sources.md](official-sources.md)。

## 1. 四类输入与资源边界

| 类别 | 严格含义 | 典型例子 |
|---|---|---|
| **A：表达矩阵 + 样本元数据可离线派生** | 给定有完整 dimnames 的表达矩阵，以及分析所需且真实填写的样本列后，可以在本地计算。A 不会补猜缺失的 `condition`、`batch`、配对、时间、结局等字段。 | MAE 构造、QC、变换、DE、PCA、聚类、SVA、WGCNA 模块、DE 排序和 LINCS 上下调 query。 |
| **B：用户必须提供真实研究或参考输入** | 该输入不能从单一 bulk 矩阵推出；`bulkMAE` 负责格式化、校验、按名称对齐或回填，但不产生输入内容。 | 定量文件及 `tx2gene`、基因长度、本地 gene set/调控网络、冻结 signature、第二队列、真实单细胞参考。 |
| **C：在线或缓存公共资源** | 函数会访问公共服务，或读取一个须先下载/安装到本地缓存的数据库。应记录资源版本、物种、ID 类型、获取日期和缓存位置；有完整缓存时，某些调用可在断网后复用。 | BioMart、OrgDb、MSigDB、KEGG/Reactome、OmniPath/decoupleR 资源、STRING、ExperimentHub LINCS。 |
| **D：外部执行、授权或许可边界** | 包只准备输入，或所选工具/数据库需要在 R 进程之外运行、登录、授权或接受许可。`bulkMAE` 不代替用户同意条款，也不自动上传研究数据。 | FASTQ 定量器、CIBERSORTx 服务/容器、受限 CIBERSORT signature、商业或账户数据库。 |

一个流程可以同时跨越多类。例如 `normalize_tpm()` 的计算本身是 A，但长度来自 B
或 C；`drug_query()` 是 A，而随后 `drug_lincs()` 搜索参考库是 C；
`deconv_cibersortx_input()` 可本地准备表格，但真正执行属于 D。后端包是否已经安装
是软件环境问题，不能把 B、C 或 D 的科学输入变成 A。

## 2. 从何处开始

若已有原始整数 gene-level counts 和样本表，最短入口是
`mae_from_matrix()`。矩阵行、列名必须完整且唯一，样本表行名必须与矩阵列名严格
一致。连续表达、芯片强度或外部 score 也可从该入口构造 MAE，但必须使用与其尺度
匹配的方法，不能把连续值命名为 `counts` 后交给 count 模型。

若起点是转录本定量文件，`import_tximport()` 可直接构建含 `counts`、`abundance`
和 sample-specific `length` 的 experiment leaf；真实文件清单及适用的 `tx2gene`
映射仍是 B。若起点只有 FASTQ，read QC、比对或 Salmon/kallisto 等定量执行仍是 D，
不属于本包的下游分析职责。

`mae_simulate()` 只用于示例、教学、后端 smoke test 和复现缺陷。它生成的条件、
注释、临床结局和表达效应没有真实生物学含义。模拟 bulk 或模拟单细胞数据不能
替代真实研究数据、独立验证队列、细胞参考、gene set、调控网络或临床证据，也不得
被描述为“真实参考资源”。

## 3. 全部 API 的可达性矩阵

表中类别描述完成该分析所需的最高层输入边界；“A（有条件）”表示相应字段必须已
真实存在于样本元数据。所有函数仍要求其后端包已经本地安装。

### 3.1 MAE、导入和注释

| 函数 | 类别 | 输入完整性与边界 |
|---|---|---|
| `mae_create_experiment()`、`mae_create()`、`mae_from_matrix()` | A | 从已命名矩阵、行/列元数据和显式 sample map 构造对象；不按位置猜样本对应关系。 |
| `mae_validate()`、`mae_experiments()`、`mae_assays()`、`mae_pull_experiment()`、`mae_pull_assay()`、`mae_samples()` | A | 本地校验、列举和提取；不会补齐缺失注释。 |
| `mae_add_assay()`、`mae_subset_features()`、`mae_feature_data()`、`mae_add_feature_data()`、`mae_add_sample_data()`、`mae_add_experiment()` | A；内容可来自 B/C | 负责按 feature/sample 名称回填 assay、元数据或新 experiment；不会判断被回填内容的生物学真实性。 |
| `mae_simulate()` | A（仅测试） | 自包含模拟 MAE；禁止作为真实科学输入。 |
| `import_tximport()` | B + D | 本地导入定量器输出并保留 effective length；FASTQ 到定量文件的执行是 D，文件和 `tx2gene` 是 B。 |
| `annotate_ensembl()` | A | 仅离线去除末尾 Ensembl 版本号，不证明 ID 的物种或 assembly。 |
| `annotation_orgdb()`、`annotate_ids()` | C（安装/缓存后本地） | 使用已安装的 OrgDb 做映射；应固定注释包版本并审查一对多映射。 |
| `annotate_rekey()` | A + B/C | 映射表给定后离线重键和显式处理重复；映射表本身必须由用户提供或从注释资源取得。 |
| `annotate_biomart()`、`annotate_gene_lengths()` | C | 在线查询 Ensembl。代表性转录本长度不是样本特异 effective length；有定量结果时优先使用 `import_tximport()` 的 `length` assay。 |

### 3.2 QC、标准化、变换和混杂调整

| 函数 | 类别 | 输入完整性与边界 |
|---|---|---|
| `filter_expr()`、`qc_library()`、`qc_correlation()`、`reduce_pca()`、`reduce_mds()`、`reduce_umap()`、`reduce_tsne()`、`qc_outliers()` | A | 从合适尺度的 assay 离线计算。离群信号不能自动证明样本错标或污染。 |
| `normalize_tmm()`、`normalize_deseq()`、`transform_vst()`、`transform_rlog()`、`transform_voom()` | A | 原始 count assay 和相应设计给定后离线计算；返回后端原生对象，可用 `mae_add_assay()` 回填。 |
| `normalize_tpm()` | A + B/C | 计算离线，但必须有正的基因/effective length。counts 本身不能推出长度。 |
| `adjust_combat()`、`adjust_batch()`、`adjust_combatseq()` | A（有条件） | 必须提供真实 batch 和应保留的生物学设计；不能从矩阵自动决定 batch 含义。统计建模通常优先把 batch 放入设计。 |
| `adjust_sva()`、`adjust_covariates()` | A（有条件） | SVA 可估计并按样本名提取潜变量；`full`/`null` 仍由研究设计决定，潜变量不能自动命名为某种真实混杂。 |
| `adjust_ruv()` | A + B | 计算离线；可靠的 control features 和 `k` 必须由研究者提供/论证，不能由同一个表达矩阵无监督地当作真值。 |
| `de_variance()` | A（有条件） | 在真实公式及样本/受试者元数据存在时离线分解方差。 |

### 3.3 差异表达、复杂设计和结果桥接

| 函数 | 类别 | 输入完整性与边界 |
|---|---|---|
| `de_design()`、`de_contrast()` | A（有条件） | 从对齐样本元数据构建设计和 limma contrast；对照方向、交互项与科学假设必须显式给出。 |
| `de_deseq2()`、`de_deseq2_results()`、`de_edger()`、`de_voom()`、`de_limma()` | A（有条件） | counts 或连续 assay、设计和 contrast 给定后离线拟合；包不猜对照组。 |
| `de_dream()` | A（有条件） | 重复测量/随机效应计算离线；真实 subject、batch 或其他分组 ID 必须在样本元数据中。 |
| `de_masigpro_design()`、`de_masigpro()` | A（有条件） | 可由真实 time/dose、replicate 和 group 元数据构造 maSigPro 设计；这些实验事实不能从表达矩阵推断。 |
| `dtu_diffsplice()` | A + B/C | 检验可离线运行；每个 exon/transcript 对应的 `gene_id`/`feature_id` 分组须来自定量注释或可靠映射。 |
| `de_table()`、`de_ranks()`、`de_selected()` | A | 把支持的原生 DE 对象稳定转换为表、有限命名排序或覆盖全部受检基因的逻辑向量；不会改变统计模型或补算缺失标准误。edgeR 联合检验因无唯一方向而拒绝，DESeq2 LRT 和通用表的 F/LR 非方向 statistic 也不能直接排序（请改用 effect 或明确的方向检验）。 |

### 3.4 Gene set、富集和调控活性

| 函数 | 类别 | 输入完整性与边界 |
|---|---|---|
| `gene_sets_prepare()`、`gene_sets_read_gmt()` | B | 对用户提供的 named list、长表或 GMT 做离线标准化；不会判断集合是否适合当前物种、ID 或问题。 |
| `gene_sets_msigdb()` | C；许可/再分发为 D | `msigdbr` 首次获取时可下载并缓存资源；应固定 MSigDB 版本，并按当前条款使用和再分发。客户端可访问不等于资源没有许可边界。 |
| `enrich_ora()`、`enrich_fgsea()`、`enrich_gsea()` | A + B | DE 列表/排序可由 A 产生，gene sets 是 B 或 C；ID 空间必须一致。 |
| `enrich_go()` | A + C | 计算使用本地安装的 OrgDb/GO 注释；首次安装和版本更新属于 C。 |
| `enrich_goseq()` | A + B/C | 需要覆盖全部受检基因的 `selected`，并需要长度/丰度 bias 及 GO/类别映射；这些不能从选择向量本身推出。 |
| `enrich_kegg()`、`enrich_reactome()` | A + C | DE 输入来自 A，通路和 ID 映射来自公共数据库/本地缓存；记录数据库版本和访问日期。 |
| `enrich_camera()`、`enrich_fry()`、`enrich_roast()` | A + B | 表达、设计和 contrast 来自研究数据；gene sets 另行提供。voom 权重不能静默丢弃。 |
| `score_gsva()`、`score_ssgsea()`、`score_singscore()` | A + B | 样本级计算离线；gene sets 或 up/down signature 必须是真实、版本化输入。 |
| `activity_methods()` | A | 只列出本地安装版本支持的方法，不取得生物学网络。 |
| `activity_resources()`、`activity_resource()` | C；资源特定条款可为 D | `show_resources()`/`get_resource()` 通过 OmniPath 列举或获取资源；需审查物种、source/target、方向、权重、版本及所选上游资源条款。 |
| `activity_decouple()` | A + B | 用户提供长表网络/集合后可离线计算；集合富集与带方向调控活性不是同一统计问题。 |
| `activity_progeny()`、`activity_tf()` | A + C | 快捷函数依赖 PROGENy、CollecTRI/DoRothEA 等公共资源；不能把资源覆盖度当作实验真值。 |
| `activity_matrix()` | A | 将公开 activity 函数的长表按名称变成 source × sample 矩阵，可交给 `mae_add_experiment()`。 |

### 3.5 聚类、共表达和网络

| 函数 | 类别 | 输入完整性与边界 |
|---|---|---|
| `cluster_consensus()`、`cluster_consensus_diagnostics()`、`cluster_consensus_classes()` | A | 在合适连续 assay 上离线聚类、审查稳定性和提取类别；最终 `k` 仍是需论证的研究选择。 |
| `cluster_nmf()`、`cluster_nmf_classes()` | A | 对非负矩阵离线拟合并提取样本类别/基底；不能自动证明分型有临床效度。 |
| `coexpr_differential()` | A（有条件） | 两组标签真实存在时比较相关结构；不是差异表达或 clinical decision curve analysis。 |
| `coexpr_pick_power()`、`coexpr_wgcna()`、`coexpr_modules()` | A | 可离线给出 power 诊断、拟合模块并抽取具名模块标签；power 和网络设定仍需审查。 |
| `coexpr_preservation()` | A + B | 包会对齐 reference/test 和模块标签；真正的第二队列、可比尺度及 reference/test 角色必须由用户提供。 |
| `network_genie3()`、`network_genie3_links()` | A；可含 B | 可从表达矩阵推断并整理候选边；若限制 regulator/target，可靠列表是 B。推断边不是因果证据。 |
| `network_string()` | A + C | 目标基因来自 assay，STRING 映射和互作来自下载/缓存资源；固定 species、版本和 score 阈值。 |

### 3.6 组成去卷积与新增单细胞参考桥接

| 函数 | 类别 | 输入完整性与边界 |
|---|---|---|
| `deconv()` | A + C；部分方法为 D | TPM-like HGNC 矩阵可在长度和映射可用后准备；EPIC、quanTIseq、xCell 等依赖方法自带/缓存参考。xCell 是相对 enrichment score，不是细胞比例。CIBERSORT 类方法另受 signature/许可约束。 |
| `deconv_reference()` | B | 对用户提供的真实 gene × cell 原始整数 counts 和 cell metadata 做严格名称对齐，构造标准 `SingleCellExperiment`，写入 `cell_type`、`sample_id` 和可选 `cell_state`；保留稀疏 Matrix，不提供或模拟生物学参考。 |
| `deconv_music()` | A + B | bulk counts 和真实单细胞参考均必需；默认直接读取 builder 的 `cell_type`/`sample_id`，并识别 `counts` 或唯一自定义 assay。至少需要两个参考生物样本，且每类细胞须跨至少两个样本。 |
| `deconv_bayesprism()` | A + B | 矩阵参考契约保持不变；builder 的 SCE 可直接传入，默认从 `counts` assay、`cell_type` 和可选 `cell_state` 取值。 |
| `deconv_cibersortx_input()` | A/B + D | 仅在本地生成 mixture 表，不写文件、不上传；CIBERSORTx 服务或容器执行、签名矩阵、账户/许可和数据出境审查均在包外。 |

### 3.7 Clinical、meta-analysis 和药物连接性

| 函数 | 类别 | 输入完整性与边界 |
|---|---|---|
| `score_signature()` | A + B | 离线打分；权重、特征、方向、冻结的训练集 center/scale 和阈值必须由真实训练研究提供。对验证队列重新估计标准化量不是锁定验证。 |
| `surv_formula()`、`surv_km()`、`surv_roc()`、`surv_cox()`、`surv_penalized()`、`ml_glmnet()` | A（有条件） | 在真实且对齐的 follow-up、event、outcome 和协变量存在时离线拟合。单队列交叉验证不等于外部验证。 |
| `meta_collect()`、`meta_effect()` | B | 可从多个兼容 DE 结果抽取同一 feature 的 effect/SE 并拟合；独立研究、可比 contrast 和有效标准误不能从单一队列产生。 |
| `drug_query()` | A | 从有限、具名且正负均存在的 DE 统计量离线生成 up/down query；用于 LINCS 时应先映射到参考库要求的 human Entrez ID。 |
| `drug_lincs_databases()` | A（资源目录） | 只返回支持数据库及 ExperimentHub ID 的元数据，不下载数据库，也不表示缓存已存在。 |
| `drug_lincs()` | B/C；受限资源可为 D | 可使用显式本地参考库（B）或让 signatureSearch 从 ExperimentHub 下载/缓存预构建 LINCS 库（C）。连接性反转是假设生成，不是疗效或用药建议。 |
| `drug_lincs_table()` | A | 从已有 `gessResult` 或同列 data frame 抽出排名表，不重算 WTCS/NCS/Tau。 |
| `drug_lincs_example()` | A（诊断玩具表） | 返回与 `gess_lincs()` 同列的合成表，供离线教程和作图；不是 CMap/LINCS 下载，也不能当作湿实验结果。 |

## 4. 只用 bulkMAE 公共 API 的链路示例

以下代码除 base/stats 的对象创建和索引外，只调用 `bulkMAE` 的公共函数。占位对象
如 `counts`、`samples`、`feature_data`、`study_gene_sets`、`locked_*`、
`local_regulatory_network`、`reference_mae`、`test_mae`、`chosen_power` 和单细胞
输入都必须换成真实、已审查的数据；示例不替代样本量、设计和模型诊断。

### 4.1 matrix → MAE → transform → DE → gene sets / GO / drug query

```r
# counts: gene x sample 原始整数矩阵
# samples: 行名与 colnames(counts) 完全一致，含 condition、batch 等真实字段
x <- mae_from_matrix(
  expression = counts,
  samples = samples,
  row_data = feature_data,
  experiment = "rna",
  assay = "counts"
)
mae_validate(x, "rna")

vst <- transform_vst(
  x,
  experiment = "rna",
  assay = "counts",
  blind = FALSE,
  design = ~ batch + condition
)
x <- mae_add_assay(x, "rna", value = vst, name = "vst")

fit <- de_deseq2(
  x,
  experiment = "rna",
  assay = "counts",
  design = ~ batch + condition
)
result <- de_deseq2_results(
  fit,
  contrast = c("condition", "treated", "control")
)

de <- de_table(result)
ranks <- de_ranks(result, column = "statistic")
selected <- de_selected(result, fdr = 0.05, min_abs_effect = 1)

# B：本地、版本化且 ID 空间与结果一致的研究 gene sets
sets <- gene_sets_prepare(study_gene_sets, min_size = 10, max_size = 500)
ora <- enrich_ora(
  genes = names(selected)[selected],
  gene_sets = sets,
  universe = names(selected)
)
gsea <- enrich_fgsea(ranks = ranks, pathways = sets)

# C：已安装/缓存的人类 OrgDb；若结果不是 SYMBOL，先显式映射
orgdb <- annotation_orgdb("human")
go <- enrich_go(
  genes = names(selected)[selected],
  method = "ora",
  org_db = orgdb,
  key_type = "SYMBOL",
  universe = names(selected)
)

id_map <- annotate_ids(
  names(ranks),
  database = orgdb,
  from = "SYMBOL",
  to = "ENTREZID",
  multi_values = "CharacterList"
)
entrez_ranks <- annotate_rekey(ranks, id_map, duplicates = "max_abs")
query <- drug_query(entrez_ranks, n = 150)

# C：只有在资源获取、版本和使用政策确认后才执行搜索
drug_lincs_databases()
# hits <- drug_lincs(query, reference_database = "lincs")
```

`CharacterList` 保留并展开一对多映射；随后必须显式选择 `sum`、`mean`、
`max_abs` 或 `first` 的重复 target 规则。`multi_values = "first"` 是有损选择，
只应在研究方案明确接受该规则时使用。

最后一行被有意保留为可选的 C 阶段：`drug_query()` 已完成离线 query 准备，
但参考库下载/缓存和实际 LINCS 搜索不是 A。

### 4.2 score、activity 和 SVA 结果回填 MAE

```r
# B：locked_* 来自训练研究；外部验证时不要在验证队列重新估计它们
signature_score <- score_signature(
  x,
  experiment = "rna",
  assay = "vst",
  weights = locked_weights,
  center = locked_center,
  scale = locked_scale
)
x <- mae_add_sample_data(
  x,
  experiment = "rna",
  value = signature_score,
  name = "locked_signature_score"
)

# B：sets 是上一节已经版本化的真实 gene sets
pathway_scores <- score_gsva(
  x,
  experiment = "rna",
  assay = "vst",
  gene_sets = sets
)
x <- mae_add_experiment(
  x,
  value = pathway_scores,
  name = "pathway_scores",
  assay_name = "score",
  source_experiment = "rna"
)

# B：本地网络由用户提供；给定后可在 A 阶段计算并通过公开桥接回填
activity_long <- activity_decouple(
  x,
  experiment = "rna",
  assay = "vst",
  network = local_regulatory_network,
  statistics = "ulm",
  source = "regulator",
  target = "target",
  mor = "weight"
)
activity_scores <- activity_matrix(activity_long, statistic = "ulm")
x <- mae_add_experiment(
  x,
  value = activity_scores,
  name = "regulator_activity",
  assay_name = "score",
  source_experiment = "rna"
)

sva_fit <- adjust_sva(
  x,
  experiment = "rna",
  assay = "vst",
  full = ~ batch + condition,
  null = ~ batch
)
sv <- adjust_covariates(sva_fit)
x <- mae_add_sample_data(x, "rna", value = sv)

# 公开 helper 回填后的 SV 列可直接进入后续设计；列数由实际结果决定
sv_formula <- stats::reformulate(c("batch", "condition", names(sv)))
fit_with_sv <- de_limma(
  x,
  experiment = "rna",
  assay = "vst",
  formula = sv_formula
)
```

SVA 的回填保证样本名对齐，不证明任一 `SV` 是某个已知 batch。把 signature、
pathway 或 activity score 回填，也不会自动完成训练/验证隔离或因果解释。

### 4.3 WGCNA modules → preservation

```r
power_diagnostic <- coexpr_pick_power(
  reference_mae,
  experiment = "rna",
  assay = "vst",
  network_type = "signed"
)

# chosen_power 必须根据已保存的诊断、样本规模和预先规定的规则确定
reference_fit <- coexpr_wgcna(
  reference_mae,
  experiment = "rna",
  assay = "vst",
  power = chosen_power,
  network_type = "signed"
)
reference_modules <- coexpr_modules(reference_fit)

preservation <- coexpr_preservation(
  reference = reference_mae,
  reference_experiment = "rna",
  reference_assay = "vst",
  test = test_mae,
  test_experiment = "rna",
  test_assay = "vst",
  module_colors = reference_modules,
  network_type = "signed",
  permutations = 200
)
```

`test_mae` 必须是真实、独立且在技术/组织背景上可比较的队列。用另一次随机模拟
替代 test 只能测试代码路径，不能形成模块保守性证据。

### 4.4 reference builder → MuSiC / BayesPrism

```r
# sc_counts: gene x cell 原始整数矩阵
# sc_cell_data: 行名覆盖全部 cell，含真实 donor、cell type 和可选 state
sc_reference <- deconv_reference(
  counts = sc_counts,
  cell_data = sc_cell_data,
  cell_type = "curated_cell_type",
  sample = "donor_id",
  cell_state = "curated_cell_state"
)

music_result <- deconv_music(
  x,
  experiment = "rna",
  assay = "counts",
  sc_reference = sc_reference
)

bayesprism_result <- deconv_bayesprism(
  x,
  experiment = "rna",
  assay = "counts",
  reference = sc_reference
)
```

builder 只把标准列写入 SCE，并严格对齐 cell 名。它不做 cell-type 注释、双细胞
剔除、供体质控、组织匹配或 reference selection。MuSiC 和 BayesPrism 的结果质量
依赖真实参考的这些属性；小型模拟 reference 只能做参数和 smoke test。

## 5. 单一 bulk 表达矩阵仍无法推断的内容

即使所有 A 类函数均可运行，以下内容仍不能从一个 bulk 矩阵本身恢复：

- 样本的真实条件、对照方向、batch、配对/受试者 ID、时间、剂量、重复及随机化；
- 生存时间、删失/事件、临床结局、治疗、混杂因素，以及独立外部验证队列；
- 物种、基因组版本、可靠 ID 类型、一对多 ID 决策、transcript-to-gene/exon 分组；
- 基因长度，尤其是样本特异的 effective length；
- 问题相关 gene sets、TF/通路调控网络、冻结 signature 权重、训练集中心/尺度、
  cutoff 和适用人群；
- 真实细胞类型、细胞状态、供体结构和细胞比例，以及与 bulk 组织匹配的单细胞参考；
- 第二队列中的模块保守性、一个唯一“正确”的 WGCNA power 或聚类类别数；
- GO、KEGG、Reactome、MSigDB、STRING、OmniPath/decoupleR、LINCS 等数据库的
  内容、版本和许可；
- FASTQ 的定量结果、CIBERSORTx 的外部执行、账户许可及上传受保护数据的授权；
- 样本错标的事实、潜变量的真实来源、调控边的因果方向、药物敏感性、疗效或安全性。

因此“无需其他工具预处理”应理解为：当分析所需的真实矩阵、元数据和参考资源已经
存在时，`bulkMAE` 提供公开的构造/校验桥接，用户不必自己编写隐藏的对象转换、按
位置对齐或结果拆表代码；它不意味着上游测序定量、外部资源获取和科学输入可以省略。

## 6. 在线、缓存和授权测试的声明边界

本指南是能力和输入边界说明，不是一次测试运行日志。这里没有声称 BioMart、KEGG、
Reactome、MSigDB、OmniPath/decoupleR、STRING、ExperimentHub/LINCS 或任何 D 类
服务已经在当前网络环境中跑通。发布验证时应分别记录每个 C/D 测试的执行日期、
资源版本、缓存状态、凭据/许可条件以及实际的运行或跳过结果；模拟数据 smoke test
也不得计作在线资源、真实参考或科学有效性测试。
