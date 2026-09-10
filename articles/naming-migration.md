# 0.4 API 命名迁移

0.4 起旧名不再导出。查询替换名：

``` r

library(bulkMAE)
bulkmae_rename("run_gsva")
#> [1] "score_gsva()"
bulkmae_rename("run_go_ora")
#> [1] "enrich_go(method = \"ora\")"
```

完整对照表见随包指南
`system.file("guides", "naming-migration-zh.md", package = "bulkMAE")`。

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

本版不额外增加只有 `...` 的总调度器
[`de()`](https://rdrr.io/r/utils/dataentry.html) 或
`enrich()`。不同后端对
`design`、`formula`、`contrast`、输入尺度和返回对象的要求差异很大；保留
可检查签名的薄函数更符合无状态、不过度设计的边界。[`deconv()`](https://younthing.github.io/bulkMAE/reference/deconv.md)
是例外，因为 它直接对应 immunedeconv 已有的统一方法契约。

## 旧名到新名

| 区域 | 0.3 及更早 | 0.4.0 |
|----|----|----|
| MAE | `validate_mae()` | [`mae_validate()`](https://younthing.github.io/bulkMAE/reference/mae_validate.md) |
| MAE | `experiment_names()` | [`mae_experiments()`](https://younthing.github.io/bulkMAE/reference/mae_experiments.md) |
| MAE | `assay_names()` | [`mae_assays()`](https://younthing.github.io/bulkMAE/reference/mae_assays.md) |
| MAE | `pull_experiment()` | [`mae_pull_experiment()`](https://younthing.github.io/bulkMAE/reference/mae_pull_experiment.md) |
| MAE | `pull_assay()` | [`mae_pull_assay()`](https://younthing.github.io/bulkMAE/reference/mae_pull_assay.md) |
| MAE | `sample_data()` | [`mae_samples()`](https://younthing.github.io/bulkMAE/reference/mae_samples.md) |
| MAE | `make_experiment()` | [`mae_create_experiment()`](https://younthing.github.io/bulkMAE/reference/mae_create_experiment.md) |
| MAE | `make_mae()` | [`mae_create()`](https://younthing.github.io/bulkMAE/reference/mae_create.md) |
| Annotation | `strip_ensembl_version()` | [`annotate_ensembl()`](https://younthing.github.io/bulkMAE/reference/annotate_ensembl.md) |
| Annotation | `map_feature_ids()` | [`annotate_ids()`](https://younthing.github.io/bulkMAE/reference/annotate_ids.md) |
| Annotation | `query_biomart()` | [`annotate_biomart()`](https://younthing.github.io/bulkMAE/reference/annotate_biomart.md) |
| QC | `filter_low_expression()` | [`filter_expr()`](https://younthing.github.io/bulkMAE/reference/filter_expr.md) |
| QC | `library_metrics()` | [`qc_library()`](https://younthing.github.io/bulkMAE/reference/qc_library.md) |
| QC | `sample_correlations()` | [`qc_correlation()`](https://younthing.github.io/bulkMAE/reference/qc_correlation.md) |
| Reduction | `run_pca()` | [`reduce_pca()`](https://younthing.github.io/bulkMAE/reference/reduce_pca.md) |
| Reduction | `run_mds()` | [`reduce_mds()`](https://younthing.github.io/bulkMAE/reference/reduce_mds.md) |
| Reduction | `run_umap()` | [`reduce_umap()`](https://younthing.github.io/bulkMAE/reference/reduce_umap.md) |
| Reduction | `run_tsne()` | [`reduce_tsne()`](https://younthing.github.io/bulkMAE/reference/reduce_tsne.md) |
| QC | `flag_pca_outliers()` | [`qc_outliers()`](https://younthing.github.io/bulkMAE/reference/qc_outliers.md) |
| Normalize | `estimate_size_factors()` | [`normalize_deseq()`](https://younthing.github.io/bulkMAE/reference/normalize_deseq.md) |
| Adjust | `remove_batch_effect()` | [`adjust_batch()`](https://younthing.github.io/bulkMAE/reference/adjust_batch.md) |
| Adjust | `adjust_combat_seq()` | [`adjust_combatseq()`](https://younthing.github.io/bulkMAE/reference/adjust_combatseq.md) |
| Adjust | `estimate_sva()` | [`adjust_sva()`](https://younthing.github.io/bulkMAE/reference/adjust_sva.md) |
| Adjust | `run_ruvg()` | [`adjust_ruv()`](https://younthing.github.io/bulkMAE/reference/adjust_ruv.md) |
| Model | `partition_variance()` | [`de_variance()`](https://younthing.github.io/bulkMAE/reference/de_variance.md) |
| Differential | `run_deseq2()` | [`de_deseq2()`](https://younthing.github.io/bulkMAE/reference/de_deseq2.md) |
| Differential | `deseq2_results()` | [`de_deseq2_results()`](https://younthing.github.io/bulkMAE/reference/de_deseq2_results.md) |
| Differential | `run_edger()` | [`de_edger()`](https://younthing.github.io/bulkMAE/reference/de_edger.md) |
| Differential | `run_limma_voom()` | [`de_voom()`](https://younthing.github.io/bulkMAE/reference/de_voom.md) |
| Differential | `run_limma()` | [`de_limma()`](https://younthing.github.io/bulkMAE/reference/de_limma.md) |
| DTU | `run_diff_splice()` | [`dtu_diffsplice()`](https://younthing.github.io/bulkMAE/reference/dtu_diffsplice.md) |
| Differential | `run_masigpro()` | [`de_masigpro()`](https://younthing.github.io/bulkMAE/reference/de_masigpro.md) |
| Differential | `run_dream()` | [`de_dream()`](https://younthing.github.io/bulkMAE/reference/de_dream.md) |
| Enrichment | `run_ora()` | [`enrich_ora()`](https://younthing.github.io/bulkMAE/reference/enrich_ora.md) |
| Enrichment | `run_go_ora()` | `enrich_go(method = "ora")` |
| Enrichment | `run_goseq()` | [`enrich_goseq()`](https://younthing.github.io/bulkMAE/reference/enrich_goseq.md) |
| Enrichment | `run_kegg_ora()` | `enrich_kegg(method = "ora")` |
| Enrichment | `run_reactome_ora()` | `enrich_reactome(method = "ora")` |
| Enrichment | `run_fgsea()` | [`enrich_fgsea()`](https://younthing.github.io/bulkMAE/reference/enrich_fgsea.md) |
| Enrichment | `run_gene_set_gsea()` | [`enrich_gsea()`](https://younthing.github.io/bulkMAE/reference/enrich_gsea.md) |
| Enrichment | `run_kegg_gsea()` | `enrich_kegg(method = "gsea")` |
| Enrichment | `run_reactome_gsea()` | `enrich_reactome(method = "gsea")` |
| Enrichment | `run_camera()` | [`enrich_camera()`](https://younthing.github.io/bulkMAE/reference/enrich_camera.md) |
| Enrichment | `run_fry()` | [`enrich_fry()`](https://younthing.github.io/bulkMAE/reference/enrich_fry.md) |
| Enrichment | `run_mroast()` | [`enrich_roast()`](https://younthing.github.io/bulkMAE/reference/enrich_roast.md) |
| Score | `run_gsva()` | [`score_gsva()`](https://younthing.github.io/bulkMAE/reference/score_gsva.md) |
| Score | `run_ssgsea()` | [`score_ssgsea()`](https://younthing.github.io/bulkMAE/reference/score_ssgsea.md) |
| Score | `run_singscore()` | [`score_singscore()`](https://younthing.github.io/bulkMAE/reference/score_singscore.md) |
| Activity | `decouple_methods()` | [`activity_methods()`](https://younthing.github.io/bulkMAE/reference/activity_methods.md) |
| Activity | `get_decouple_resource()` | [`activity_resource()`](https://younthing.github.io/bulkMAE/reference/activity_resource.md) |
| Activity | `run_decouple()` | [`activity_decouple()`](https://younthing.github.io/bulkMAE/reference/activity_decouple.md) |
| Activity | `infer_progeny()` | [`activity_progeny()`](https://younthing.github.io/bulkMAE/reference/activity_progeny.md) |
| Activity | `infer_tf_activity()` | [`activity_tf()`](https://younthing.github.io/bulkMAE/reference/activity_tf.md) |
| Cluster | `run_consensus_clustering()` | [`cluster_consensus()`](https://younthing.github.io/bulkMAE/reference/cluster_consensus.md) |
| Cluster | `run_consensus_icl()` | [`cluster_consensus_diagnostics()`](https://younthing.github.io/bulkMAE/reference/cluster_consensus_diagnostics.md) |
| Cluster | `consensus_classes()` | [`cluster_consensus_classes()`](https://younthing.github.io/bulkMAE/reference/cluster_consensus_classes.md) |
| Co-expression | `run_differential_coexpression()` | [`coexpr_differential()`](https://younthing.github.io/bulkMAE/reference/coexpr_differential.md) |
| Cluster | `run_nmf()` | [`cluster_nmf()`](https://younthing.github.io/bulkMAE/reference/cluster_nmf.md) |
| Co-expression | `run_wgcna()` | [`coexpr_wgcna()`](https://younthing.github.io/bulkMAE/reference/coexpr_wgcna.md) |
| Co-expression | `run_module_preservation()` | [`coexpr_preservation()`](https://younthing.github.io/bulkMAE/reference/coexpr_preservation.md) |
| Network | `run_genie3()` | [`network_genie3()`](https://younthing.github.io/bulkMAE/reference/network_genie3.md) |
| Network | `genie3_links()` | [`network_genie3_links()`](https://younthing.github.io/bulkMAE/reference/network_genie3_links.md) |
| Network | `run_string_network()` | [`network_string()`](https://younthing.github.io/bulkMAE/reference/network_string.md) |
| Deconvolution | `run_immunedeconv()` | [`deconv()`](https://younthing.github.io/bulkMAE/reference/deconv.md) |
| Deconvolution | `run_music()` | [`deconv_music()`](https://younthing.github.io/bulkMAE/reference/deconv_music.md) |
| Deconvolution | `prepare_cibersortx_input()` | [`deconv_cibersortx_input()`](https://younthing.github.io/bulkMAE/reference/deconv_cibersortx_input.md) |
| Deconvolution | `run_bayesprism()` | [`deconv_bayesprism()`](https://younthing.github.io/bulkMAE/reference/deconv_bayesprism.md) |
| Survival | `run_kaplan_meier()` | [`surv_km()`](https://younthing.github.io/bulkMAE/reference/surv_km.md) |
| Survival | `run_time_roc()` | [`surv_roc()`](https://younthing.github.io/bulkMAE/reference/surv_roc.md) |
| Survival | `run_cox()` | [`surv_cox()`](https://younthing.github.io/bulkMAE/reference/surv_cox.md) |
| Survival | `run_penalized_cox()` | [`surv_penalized()`](https://younthing.github.io/bulkMAE/reference/surv_penalized.md) |
| Machine learning | `run_classifier_glmnet()` | [`ml_glmnet()`](https://younthing.github.io/bulkMAE/reference/ml_glmnet.md) |
| Meta-analysis | `run_meta_analysis()` | [`meta_effect()`](https://younthing.github.io/bulkMAE/reference/meta_effect.md) |
| Drug | `prepare_lincs_query()` | [`drug_query()`](https://younthing.github.io/bulkMAE/reference/drug_query.md) |
| Drug | `run_lincs_search()` | [`drug_lincs()`](https://younthing.github.io/bulkMAE/reference/drug_lincs.md) |

## 保持不变的名称

[`import_tximport()`](https://younthing.github.io/bulkMAE/reference/import_tximport.md)、[`normalize_tmm()`](https://younthing.github.io/bulkMAE/reference/normalize_tmm.md)、[`transform_voom()`](https://younthing.github.io/bulkMAE/reference/transform_voom.md)、
[`transform_vst()`](https://younthing.github.io/bulkMAE/reference/transform_vst.md)、[`transform_rlog()`](https://younthing.github.io/bulkMAE/reference/transform_rlog.md)、[`adjust_combat()`](https://younthing.github.io/bulkMAE/reference/adjust_combat.md)
和
[`score_signature()`](https://younthing.github.io/bulkMAE/reference/score_signature.md)
已符合新规则，不需重命名。

## 数据库富集的合并

GO、KEGG 和 Reactome 不再为 ORA/GSEA 建立组合式函数名。使用：

``` r

enrich_go(genes = selected, method = "ora", org_db = org.Hs.eg.db)
enrich_go(ranks = statistic, method = "gsea", org_db = org.Hs.eg.db)

enrich_kegg(genes = selected, method = "ora")
enrich_kegg(ranks = statistic, method = "gsea")

enrich_reactome(genes = selected, method = "ora")
enrich_reactome(ranks = statistic, method = "gsea")
```

这样将“统计方法”保留为 `method`，将“知识库”保留在函数家族名中， 避免
`ora_*` 和 `gsea_*` 组合随后端增长。
