test_that("consensus classes are extracted without a new result class", {
  results <- vector("list", 3)
  results[[2]] <- list(consensusClass = c(sample1 = 1L, sample2 = 2L))

  expect_identical(
    cluster_consensus_classes(results, 2),
    c(sample1 = 1L, sample2 = 2L)
  )
  expect_error(cluster_consensus_classes(results, 4), "available consensus solution")
})

test_that("expanded wrappers retain explicit MAE arguments", {
  mae_wrappers <- c(
    "activity_decouple",
    "coexpr_differential",
    "de_masigpro",
    "score_singscore",
    "surv_roc",
    "reduce_umap"
  )

  for (name in mae_wrappers) {
    arguments <- names(formals(getExportedValue("bulkMAE", name)))
    expect_true(all(c("x", "experiment") %in% arguments), info = name)
  }
})

test_that("decouple adapter exposes method and consensus choices", {
  arguments <- names(formals(activity_decouple))
  expect_true(all(c(
    "network",
    "statistics",
    "method_args",
    "consensus",
    "consensus_statistics"
  ) %in% arguments))
})

test_that("public API follows discoverable analysis families", {
  exports <- getNamespaceExports("bulkMAE")
  expected <- c(
    "activity_decouple", "activity_matrix", "activity_methods",
    "activity_progeny", "activity_resource", "activity_resources",
    "activity_tf", "adjust_batch", "adjust_combat", "adjust_combatseq",
    "adjust_covariates", "adjust_ruv", "adjust_sva", "annotate_biomart",
    "annotate_ensembl", "annotate_gene_lengths", "annotate_ids",
    "annotate_rekey", "annotation_orgdb", "bulkmae_rename", "cluster_consensus",
    "cluster_consensus_classes", "cluster_consensus_diagnostics",
    "cluster_nmf", "cluster_nmf_classes", "coexpr_differential",
    "coexpr_hubs", "coexpr_membership", "coexpr_module_trait",
    "coexpr_modules", "coexpr_pick_power", "coexpr_preservation",
    "coexpr_wgcna", "de_contrast", "de_design", "de_deseq2",
    "de_deseq2_results", "de_dream",
    "de_edger", "de_limma", "de_masigpro", "de_masigpro_design",
    "de_ranks", "de_selected", "de_table", "de_variance", "de_voom",
    "deconv", "deconv_bayesprism", "deconv_cibersortx_input",
    "deconv_fractions", "deconv_music", "deconv_reference", "drug_lincs",
    "drug_lincs_databases", "drug_query", "dtu_diffsplice",
    "enrich_camera", "enrich_fgsea", "enrich_fry", "enrich_go",
    "enrich_goseq", "enrich_gsea", "enrich_kegg", "enrich_ora",
    "enrich_reactome", "enrich_roast", "filter_expr", "gene_sets_msigdb",
    "gene_sets_prepare", "gene_sets_read_gmt", "import_tximport",
    "mae_add_assay", "mae_add_experiment", "mae_add_feature_data",
    "mae_add_sample_data", "mae_assays", "mae_create",
    "mae_create_experiment", "mae_experiments",
    "mae_feature_data",
    "mae_from_matrix", "mae_pull_assay", "mae_pull_experiment",
    "mae_samples", "mae_simulate", "mae_subset_features",
    "mae_validate", "mae_variable_features",
    "meta_collect", "meta_effect", "ml_glmnet",
    "network_genie3", "network_genie3_links", "network_string",
    "normalize_deseq", "normalize_tmm", "normalize_tpm", "qc_correlation", "qc_library",
    "qc_outliers", "plot_assay_heatmap", "plot_cluster_cdf",
    "plot_cluster_consensus", "plot_cluster_sizes",
    "plot_coexpr_membership", "plot_coexpr_modules",
    "plot_coexpr_power", "plot_coexpr_trait",
    "plot_de_ma", "plot_de_volcano", "plot_deconv_box",
    "plot_deconv_heatmap", "plot_deconv_stacked",
    "plot_embedding", "plot_gsea_classic",
    "plot_gsea_ridge", "plot_ora_bubble", "plot_ora_network",
    "plot_ora_radial", "plot_qc_correlation", "plot_qc_library",
    "plot_save", "plot_qc_outliers",
    "theme_bulkmae", "reduce_mds", "reduce_pca", "reduce_tsne",
    "reduce_umap", "score_gsva", "score_signature", "score_singscore",
    "score_ssgsea", "surv_cox", "surv_formula", "surv_km",
    "surv_penalized", "surv_roc",
    "transform_rlog", "transform_voom", "transform_vst"
  )

  expect_setequal(exports, expected)
  expect_false(any(grepl("^(run|infer|prepare|estimate|get)_", exports)))
})

test_that("database enrichment uses one method-dispatched API", {
  for (name in c("enrich_go", "enrich_kegg", "enrich_reactome")) {
    arguments <- names(formals(getExportedValue("bulkMAE", name)))
    expect_true(all(c("genes", "ranks", "method") %in% arguments), info = name)
  }
})
