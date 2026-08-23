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
    "activity_decouple", "activity_methods", "activity_progeny",
    "activity_resource", "activity_tf", "adjust_batch", "adjust_combat",
    "adjust_combatseq", "adjust_ruv", "adjust_sva", "annotate_biomart",
    "annotate_ensembl", "annotate_ids", "cluster_consensus",
    "cluster_consensus_classes", "cluster_consensus_diagnostics",
    "cluster_nmf", "coexpr_differential", "coexpr_pick_power",
    "coexpr_preservation", "coexpr_wgcna", "de_deseq2",
    "de_deseq2_results", "de_dream",
    "de_edger", "de_limma", "de_masigpro", "de_variance", "de_voom",
    "deconv", "deconv_bayesprism", "deconv_cibersortx_input",
    "deconv_music", "drug_lincs", "drug_query", "dtu_diffsplice",
    "enrich_camera", "enrich_fgsea", "enrich_fry", "enrich_go",
    "enrich_goseq", "enrich_gsea", "enrich_kegg", "enrich_ora",
    "enrich_reactome", "enrich_roast", "filter_expr", "import_tximport",
    "mae_add_assay", "mae_assays", "mae_create", "mae_create_experiment",
    "mae_experiments", "mae_pull_assay", "mae_pull_experiment",
    "mae_samples", "mae_subset_features", "mae_validate", "meta_effect",
    "ml_glmnet",
    "network_genie3", "network_genie3_links", "network_string",
    "normalize_deseq", "normalize_tmm", "normalize_tpm", "qc_correlation", "qc_library",
    "qc_outliers", "reduce_mds", "reduce_pca", "reduce_tsne",
    "reduce_umap", "score_gsva", "score_signature", "score_singscore",
    "score_ssgsea", "surv_cox", "surv_km", "surv_penalized", "surv_roc",
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
