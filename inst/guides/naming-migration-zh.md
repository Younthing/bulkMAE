# bulkMAE 0.4 API 命名迁移

0.4.0 将公开 API 统一为 `<family>_<method>()` 或
`<family>_<operation>()`。这是有意的破坏性重命名：旧函数不再导出，也不保留
deprecated aliases，避免自动补全继续出现 `run_*`、`infer_*`、`prepare_*`
和新函数并存的情况。

## 命名边界

- `mae_*`：MAE/SE 构建、校验和取值。
- `qc_*`、`filter_*`、`reduce_*`：质控、过滤和样本降维。
- `normalize_*`、`transform_*`、`adjust_*`：严格区分归一化、尺度变换和批次/
  隐变量校正。
- `de_*`、`dtu_*`：表达模型/差异分析与转录本使用。
- `enrich_*`、`score_*`、`activity_*`：集合检验、样本级评分和调控活性。
- `cluster_*`、`coexpr_*`、`network_*`：分型、共表达和一般网络。
- `deconv_*`、`surv_*`、`ml_*`、`meta_*`、`drug_*`：组织反卷积与临床转化。

本版不额外增加只有 `...` 的总调度器 `de()` 或 `enrich()`。不同后端对
`design`、`formula`、`contrast`、输入尺度和返回对象的要求差异很大；保留
可检查签名的薄函数更符合无状态、不过度设计的边界。`deconv()` 是例外，因为
它直接对应 immunedeconv 已有的统一方法契约。

## 旧名到新名

| 区域 | 0.3 及更早 | 0.4.0 |
|---|---|---|
| MAE | `validate_mae()` | `mae_validate()` |
| MAE | `experiment_names()` | `mae_experiments()` |
| MAE | `assay_names()` | `mae_assays()` |
| MAE | `pull_experiment()` | `mae_pull_experiment()` |
| MAE | `pull_assay()` | `mae_pull_assay()` |
| MAE | `sample_data()` | `mae_samples()` |
| MAE | `make_experiment()` | `mae_create_experiment()` |
| MAE | `make_mae()` | `mae_create()` |
| Annotation | `strip_ensembl_version()` | `annotate_ensembl()` |
| Annotation | `map_feature_ids()` | `annotate_ids()` |
| Annotation | `query_biomart()` | `annotate_biomart()` |
| QC | `filter_low_expression()` | `filter_expr()` |
| QC | `library_metrics()` | `qc_library()` |
| QC | `sample_correlations()` | `qc_correlation()` |
| Reduction | `run_pca()` | `reduce_pca()` |
| Reduction | `run_mds()` | `reduce_mds()` |
| Reduction | `run_umap()` | `reduce_umap()` |
| Reduction | `run_tsne()` | `reduce_tsne()` |
| QC | `flag_pca_outliers()` | `qc_outliers()` |
| Normalize | `estimate_size_factors()` | `normalize_deseq()` |
| Adjust | `remove_batch_effect()` | `adjust_batch()` |
| Adjust | `adjust_combat_seq()` | `adjust_combatseq()` |
| Adjust | `estimate_sva()` | `adjust_sva()` |
| Adjust | `run_ruvg()` | `adjust_ruv()` |
| Model | `partition_variance()` | `de_variance()` |
| Differential | `run_deseq2()` | `de_deseq2()` |
| Differential | `deseq2_results()` | `de_deseq2_results()` |
| Differential | `run_edger()` | `de_edger()` |
| Differential | `run_limma_voom()` | `de_voom()` |
| Differential | `run_limma()` | `de_limma()` |
| DTU | `run_diff_splice()` | `dtu_diffsplice()` |
| Differential | `run_masigpro()` | `de_masigpro()` |
| Differential | `run_dream()` | `de_dream()` |
| Enrichment | `run_ora()` | `enrich_ora()` |
| Enrichment | `run_go_ora()` | `enrich_go(method = "ora")` |
| Enrichment | `run_goseq()` | `enrich_goseq()` |
| Enrichment | `run_kegg_ora()` | `enrich_kegg(method = "ora")` |
| Enrichment | `run_reactome_ora()` | `enrich_reactome(method = "ora")` |
| Enrichment | `run_fgsea()` | `enrich_fgsea()` |
| Enrichment | `run_gene_set_gsea()` | `enrich_gsea()` |
| Enrichment | `run_kegg_gsea()` | `enrich_kegg(method = "gsea")` |
| Enrichment | `run_reactome_gsea()` | `enrich_reactome(method = "gsea")` |
| Enrichment | `run_camera()` | `enrich_camera()` |
| Enrichment | `run_fry()` | `enrich_fry()` |
| Enrichment | `run_mroast()` | `enrich_roast()` |
| Score | `run_gsva()` | `score_gsva()` |
| Score | `run_ssgsea()` | `score_ssgsea()` |
| Score | `run_singscore()` | `score_singscore()` |
| Activity | `decouple_methods()` | `activity_methods()` |
| Activity | `get_decouple_resource()` | `activity_resource()` |
| Activity | `run_decouple()` | `activity_decouple()` |
| Activity | `infer_progeny()` | `activity_progeny()` |
| Activity | `infer_tf_activity()` | `activity_tf()` |
| Cluster | `run_consensus_clustering()` | `cluster_consensus()` |
| Cluster | `run_consensus_icl()` | `cluster_consensus_diagnostics()` |
| Cluster | `consensus_classes()` | `cluster_consensus_classes()` |
| Co-expression | `run_differential_coexpression()` | `coexpr_differential()` |
| Cluster | `run_nmf()` | `cluster_nmf()` |
| Co-expression | `run_wgcna()` | `coexpr_wgcna()` |
| Co-expression | `run_module_preservation()` | `coexpr_preservation()` |
| Network | `run_genie3()` | `network_genie3()` |
| Network | `genie3_links()` | `network_genie3_links()` |
| Network | `run_string_network()` | `network_string()` |
| Deconvolution | `run_immunedeconv()` | `deconv()` |
| Deconvolution | `run_music()` | `deconv_music()` |
| Deconvolution | `prepare_cibersortx_input()` | `deconv_cibersortx_input()` |
| Deconvolution | `run_bayesprism()` | `deconv_bayesprism()` |
| Survival | `run_kaplan_meier()` | `surv_km()` |
| Survival | `run_time_roc()` | `surv_roc()` |
| Survival | `run_cox()` | `surv_cox()` |
| Survival | `run_penalized_cox()` | `surv_penalized()` |
| Machine learning | `run_classifier_glmnet()` | `ml_glmnet()` |
| Meta-analysis | `run_meta_analysis()` | `meta_effect()` |
| Drug | `prepare_lincs_query()` | `drug_query()` |
| Drug | `run_lincs_search()` | `drug_lincs()` |

## 保持不变的名称

`import_tximport()`、`normalize_tmm()`、`transform_voom()`、
`transform_vst()`、`transform_rlog()`、`adjust_combat()` 和
`score_signature()` 已符合新规则，不需重命名。

## 数据库富集的合并

GO、KEGG 和 Reactome 不再为 ORA/GSEA 建立组合式函数名。使用：

```r
enrich_go(genes = selected, method = "ora", org_db = org.Hs.eg.db)
enrich_go(ranks = statistic, method = "gsea", org_db = org.Hs.eg.db)

enrich_kegg(genes = selected, method = "ora")
enrich_kegg(ranks = statistic, method = "gsea")

enrich_reactome(genes = selected, method = "ora")
enrich_reactome(ranks = statistic, method = "gsea")
```

这样将“统计方法”保留为 `method`，将“知识库”保留在函数家族名中，
避免 `ora_*` 和 `gsea_*` 组合随后端增长。
