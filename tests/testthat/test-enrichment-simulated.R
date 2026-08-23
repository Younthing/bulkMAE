.simulated_enrichment_gene_sets <- function() {
  list(
    signal = sprintf("gene%04d", seq_len(80L)),
    background = sprintf("gene%04d", 81:160),
    mixed = c(
      sprintf("gene%04d", seq_len(20L)),
      sprintf("gene%04d", 161:220)
    )
  )
}

.skip_if_enrichment_backend_missing <- function(pkg, minimum_version = NULL) {
  suppressMessages(suppressWarnings(
    testthat::skip_if_not_installed(pkg, minimum_version = minimum_version)
  ))
}

.simulated_enrichment_ranks <- function(mae) {
  expression <- mae_pull_assay(mae, "rna", "log_expression")
  samples <- mae_samples(mae, "rna")
  treated <- samples$condition == "treated"
  ranks <- rowMeans(expression[, treated, drop = FALSE]) -
    rowMeans(expression[, !treated, drop = FALSE])
  stats::setNames(as.numeric(ranks), rownames(expression))
}

test_that("arbitrary gene sets support offline over-representation analysis", {
  .skip_if_enrichment_backend_missing("clusterProfiler")
  mae <- make_simulated_mae(seed = 1701L)
  gene_sets <- .simulated_enrichment_gene_sets()
  universe <- rownames(mae_pull_assay(mae, "rna", "log_expression"))

  result <- suppressMessages(suppressWarnings(enrich_ora(
    genes = sprintf("gene%04d", seq_len(60L)),
    gene_sets = gene_sets,
    universe = universe,
    pvalueCutoff = 1,
    qvalueCutoff = 1,
    minGSSize = 10L,
    maxGSSize = 100L
  )))
  table <- as.data.frame(result)

  expect_s4_class(result, "enrichResult")
  expect_gt(nrow(table), 0L)
  expect_true("signal" %in% table$ID)
  expect_true(all(is.finite(table$pvalue)))
})

test_that("fgsea runs offline on statistics derived from the simulated MAE", {
  .skip_if_enrichment_backend_missing("fgsea")
  mae <- make_simulated_mae(seed = 1701L)
  gene_sets <- .simulated_enrichment_gene_sets()

  result <- suppressWarnings(enrich_fgsea(
    ranks = .simulated_enrichment_ranks(mae),
    pathways = gene_sets,
    min_size = 10L,
    max_size = 100L,
    eps = 0
  ))

  expect_s3_class(result, "data.frame")
  expect_setequal(result$pathway, names(gene_sets))
  expect_true(all(is.finite(result$NES)))
  expect_true(is.list(result$leadingEdge))
})

test_that("clusterProfiler GSEA runs against local arbitrary gene sets", {
  .skip_if_enrichment_backend_missing("clusterProfiler")
  mae <- make_simulated_mae(seed = 1701L)
  gene_sets <- .simulated_enrichment_gene_sets()

  result <- suppressMessages(suppressWarnings(enrich_gsea(
    ranks = .simulated_enrichment_ranks(mae),
    gene_sets = gene_sets,
    pvalueCutoff = 1,
    minGSSize = 10L,
    maxGSSize = 100L,
    eps = 0,
    verbose = FALSE,
    seed = TRUE
  )))
  table <- as.data.frame(result)

  expect_s4_class(result, "gseaResult")
  expect_gt(nrow(table), 0L)
  expect_true("signal" %in% table$ID)
  expect_true(all(is.finite(table$NES)))
})

test_that("CAMERA, FRY, and mroast use local expression and precision weights", {
  .skip_if_enrichment_backend_missing("limma")
  mae <- make_simulated_mae(seed = 1701L)
  gene_sets <- .simulated_enrichment_gene_sets()

  camera <- enrich_camera(
    mae,
    "rna",
    gene_sets = gene_sets,
    formula = ~ condition,
    contrast = 2L,
    assay = "log_expression",
    weights = "weights"
  )
  fry <- enrich_fry(
    mae,
    "rna",
    gene_sets = gene_sets,
    formula = ~ condition,
    contrast = 2L,
    assay = "log_expression",
    weights = "weights"
  )

  set.seed(983)
  caller_seed <- .Random.seed
  roast <- enrich_roast(
    mae,
    "rna",
    gene_sets = gene_sets,
    formula = ~ condition,
    contrast = 2L,
    assay = "log_expression",
    weights = "weights",
    rotations = 99L,
    seed = 33L
  )

  for (result in list(camera, fry, roast)) {
    expect_s3_class(result, "data.frame")
    expect_setequal(rownames(result), names(gene_sets))
    expect_true(all(is.finite(result$PValue)))
  }
  expect_identical(.Random.seed, caller_seed)
})

test_that("goseq runs offline with simulated gene lengths and local categories", {
  .skip_if_enrichment_backend_missing("goseq")
  mae <- make_simulated_mae(seed = 1701L)
  experiment <- mae_pull_experiment(mae, "rna")
  genes <- rownames(experiment)
  selected <- stats::setNames(seq_along(genes) <= 80L, genes)
  bias <- stats::setNames(
    as.numeric(SummarizedExperiment::rowData(experiment)$gene_length),
    genes
  )
  gene_to_category <- data.frame(
    gene = rep(genes[seq_len(360L)], 2L),
    category = rep(paste0("set_", letters[1:6]), each = 120L),
    stringsAsFactors = FALSE
  )

  result <- suppressMessages(suppressWarnings(enrich_goseq(
    selected = selected,
    genome = "custom",
    id = "custom",
    bias = bias,
    gene_to_category = gene_to_category,
    categories = NULL,
    plot_fit = FALSE
  )))

  expect_s3_class(result, "data.frame")
  expect_gt(nrow(result), 0L)
  expect_true("set_a" %in% result$category)
  expect_true(all(is.finite(result$over_represented_pvalue)))
})

test_that("GSVA and ssGSEA score every simulated sample offline", {
  .skip_if_enrichment_backend_missing("GSVA", minimum_version = "2.0.0")
  mae <- make_simulated_mae(seed = 1701L)
  gene_sets <- .simulated_enrichment_gene_sets()
  samples <- mae_samples(mae, "rna")

  gsva <- score_gsva(
    mae,
    "rna",
    gene_sets = gene_sets,
    assay = "log_expression",
    kcdf = "Gaussian",
    min_size = 10L,
    max_size = 100L,
    verbose = FALSE
  )
  ssgsea <- score_ssgsea(
    mae,
    "rna",
    gene_sets = gene_sets,
    assay = "log_expression",
    min_size = 10L,
    max_size = 100L,
    verbose = FALSE
  )

  for (result in list(gsva, ssgsea)) {
    expect_true(is.matrix(result))
    expect_identical(dim(result), c(3L, 24L))
    expect_identical(colnames(result), rownames(samples))
    expect_true(all(is.finite(result)))
    expect_gt(
      mean(result["signal", samples$condition == "treated"]),
      mean(result["signal", samples$condition == "control"])
    )
  }
})

test_that("singscore returns finite rank-based scores for every sample", {
  .skip_if_enrichment_backend_missing("singscore")
  mae <- make_simulated_mae(seed = 1701L)
  samples <- mae_samples(mae, "rna")

  result <- suppressMessages(score_singscore(
    mae,
    "rna",
    up_set = sprintf("gene%04d", seq_len(30L)),
    down_set = sprintf("gene%04d", 81:110),
    assay = "log_expression"
  ))

  expect_s3_class(result, "data.frame")
  expect_identical(rownames(result), rownames(samples))
  expect_true(all(is.finite(result$TotalScore)))
  expect_gt(
    mean(result$TotalScore[samples$condition == "treated"]),
    mean(result$TotalScore[samples$condition == "control"])
  )
})

test_that("GO enrichment uses a locally installed human OrgDb", {
  .skip_if_enrichment_backend_missing("clusterProfiler")
  .skip_if_enrichment_backend_missing("AnnotationDbi")
  .skip_if_enrichment_backend_missing("org.Hs.eg.db")

  org_db <- suppressWarnings(suppressPackageStartupMessages(
    getExportedValue("org.Hs.eg.db", "org.Hs.eg.db")
  ))
  identifiers <- AnnotationDbi::keys(org_db, keytype = "ENTREZID")
  mapping <- suppressMessages(AnnotationDbi::select(
    org_db,
    keys = utils::head(identifiers, 2000L),
    columns = c("GO", "ONTOLOGY"),
    keytype = "ENTREZID"
  ))
  mapping <- mapping[
    !is.na(mapping$GO) & mapping$ONTOLOGY == "BP",
    ,
    drop = FALSE
  ]
  sizes <- sort(table(mapping$GO), decreasing = TRUE)
  eligible <- names(sizes)[sizes >= 20L & sizes <= 500L]
  skip_if(!length(eligible), "The installed OrgDb has no suitable local BP term.")
  selected_term <- eligible[[1L]]
  selected_genes <- unique(mapping$ENTREZID[mapping$GO == selected_term])
  universe <- unique(mapping$ENTREZID)

  result <- suppressMessages(suppressWarnings(enrich_go(
    genes = selected_genes,
    method = "ora",
    org_db = org_db,
    key_type = "ENTREZID",
    ontology = "BP",
    universe = universe,
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.2,
    minGSSize = 5L,
    maxGSSize = 1000L,
    readable = FALSE
  )))
  table <- as.data.frame(result)

  expect_s4_class(result, "enrichResult")
  expect_gt(nrow(table), 0L)
  expect_true(selected_term %in% table$ID)
})

test_that("Reactome enrichment can use its installed local pathway data", {
  .skip_if_enrichment_backend_missing("ReactomePA")
  .skip_if_enrichment_backend_missing("AnnotationDbi")
  .skip_if_enrichment_backend_missing("org.Hs.eg.db")

  org_db <- suppressWarnings(suppressPackageStartupMessages(
    getExportedValue("org.Hs.eg.db", "org.Hs.eg.db")
  ))
  selected_genes <- c(
    "7157", "4609", "1029", "1956", "672", "675", "5925", "894",
    "891", "890", "983", "991", "9133", "1062", "4171", "4172",
    "4173", "4174", "4175", "4176"
  )
  universe <- unique(c(
    selected_genes,
    utils::head(AnnotationDbi::keys(org_db, keytype = "ENTREZID"), 5000L)
  ))

  result <- suppressMessages(suppressWarnings(enrich_reactome(
    genes = selected_genes,
    method = "ora",
    organism = "human",
    universe = universe,
    pvalueCutoff = 1,
    qvalueCutoff = 1,
    minGSSize = 5L,
    maxGSSize = 500L,
    readable = FALSE
  )))
  table <- as.data.frame(result)

  expect_s4_class(result, "enrichResult")
  expect_gt(nrow(table), 0L)
  expect_true(all(is.finite(table$pvalue)))
})

test_that("KEGG enrichment is an explicit opt-in online test", {
  run_online <- identical(
    tolower(Sys.getenv("BULKMAE_RUN_ONLINE_TESTS", "false")),
    "true"
  )
  skip_if(!run_online, "KEGG uses a remote service; set BULKMAE_RUN_ONLINE_TESTS=true.")
  .skip_if_enrichment_backend_missing("clusterProfiler")

  result <- enrich_kegg(
    genes = c("7157", "4609", "1029", "1956", "672", "675"),
    method = "ora",
    organism = "hsa",
    pvalueCutoff = 1,
    qvalueCutoff = 1
  )
  expect_s4_class(result, "enrichResult")
})
