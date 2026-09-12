network_api <- function(name) {
  get(name, envir = asNamespace("bulkMAE"), inherits = FALSE)
}

network_mae_with_assay <- function(x, value, assay_name = "expression") {
  samples <- mae_samples(x, "rna")
  experiment <- mae_create_experiment(
    assays = stats::setNames(list(value), assay_name),
    col_data = samples
  )
  mae_create(list(rna = experiment), samples)
}

test_that("consensus clustering and its extractors run offline", {
  skip_if_not_installed("ConsensusClusterPlus")
  mae <- make_model_mae(n_features = 30L, n_samples = 12L, seed = 101L)
  output_directory <- tempfile("bulkmae-consensus-")
  dir.create(output_directory)

  fit <- cluster_consensus(
    mae,
    "rna",
    assay = "log_expression",
    max_k = 3L,
    repetitions = 5L,
    item_fraction = 1,
    feature_fraction = 1,
    seed = 7L,
    title = output_directory,
    plot = NULL
  )
  classes <- cluster_consensus_classes(fit, 2L)
  diagnostics <- cluster_consensus_diagnostics(
    fit,
    title = output_directory,
    plot = NULL
  )
  unlink(output_directory, recursive = TRUE, force = TRUE)

  expect_type(fit, "list")
  expect_length(fit, 3L)
  expect_length(classes, 12L)
  expect_setequal(names(classes), paste0("sample", seq_len(12L)))
  expect_type(diagnostics, "list")
  expect_true(all(c("clusterConsensus", "itemConsensus") %in% names(diagnostics)))
  pac <- cluster_consensus_pac(fit)
  delta <- cluster_consensus_delta_area(fit)
  expect_named(pac, c("2", "3"))
  expect_true(all(is.finite(pac) & pac >= 0 & pac <= 1))
  expect_identical(delta$k, c(2L, 3L))
  expect_true(all(is.finite(delta$area)))
})

test_that("consensus clustering records removed constant features", {
  skip_if_not_installed("ConsensusClusterPlus")
  mae <- make_model_mae(n_features = 20L, n_samples = 10L, seed = 102L)
  expression <- mae_pull_assay(mae, "rna", "log_expression")
  expression[1L, ] <- 3
  mae <- network_mae_with_assay(mae, expression)
  output_directory <- tempfile("bulkmae-consensus-filter-")
  dir.create(output_directory)

  expect_warning(
    fit <- cluster_consensus(
      mae,
      "rna",
      assay = "expression",
      max_k = 3L,
      repetitions = 5L,
      item_fraction = 1,
      feature_fraction = 1,
      title = output_directory,
      plot = NULL
    ),
    "removed 1 non-variable"
  )
  unlink(output_directory, recursive = TRUE, force = TRUE)

  expect_identical(attr(fit, "bulkMAERemovedFeatures"), "gene1")
})

test_that("consensus clustering rejects malformed requests before fitting", {
  skip_if_not_installed("ConsensusClusterPlus")
  mae <- make_model_mae(n_features = 8L, n_samples = 6L)

  expect_error(
    cluster_consensus(mae, "rna", "log_expression", max_k = 7L),
    "cannot exceed"
  )
  expect_error(
    cluster_consensus(
      mae,
      "rna",
      "log_expression",
      max_k = 3L,
      item_fraction = 0.2
    ),
    "too few items"
  )
  expect_error(
    cluster_consensus_diagnostics(list(list()), plot = NULL),
    "result list"
  )
  expect_error(cluster_consensus_classes(vector("list", 2L), 3L), "available")
  expect_error(
    cluster_consensus_classes(list(NULL, list()), 2L),
    "consensusClass"
  )
})

test_that("differential co-expression runs and audits filtered features", {
  skip_if_not_installed("diffcoexp")
  mae <- make_model_mae(n_features = 12L, n_samples = 20L, seed = 103L)
  expression <- mae_pull_assay(mae, "rna", "log_expression")
  expression[1L, ] <- 2
  mae <- network_mae_with_assay(mae, expression)

  expect_warning(
    invisible(utils::capture.output(
      fit <- coexpr_differential(
        mae,
        "rna",
        group = "condition",
        contrast = c("control", "treated"),
        assay = "expression",
        max_pairs = 100L
      )
    )),
    "removed 1 feature"
  )

  expect_type(fit, "list")
  expect_true(all(c("DCGs", "DCLs") %in% names(fit)))
  expect_identical(attr(fit, "bulkMAERemovedFeatures"), "gene1")
})

test_that("differential co-expression validates groups and analysis size", {
  skip_if_not_installed("diffcoexp")
  mae <- make_model_mae(n_features = 8L, n_samples = 8L)

  expect_error(
    coexpr_differential(
      mae,
      "rna",
      "condition",
      c("control", "missing"),
      "log_expression"
    ),
    "observed group labels"
  )
  expect_error(
    coexpr_differential(
      mae,
      "rna",
      "condition",
      c("control", "treated"),
      "log_expression",
      max_pairs = 2L
    ),
    "exceeding"
  )
  small <- make_model_mae(n_features = 4L, n_samples = 6L)
  expect_error(
    coexpr_differential(
      small,
      "rna",
      "condition",
      c("control", "treated"),
      "log_expression"
    ),
    "at least four samples per group"
  )
})

test_that("NMF wrapper passes a non-negative assay to the backend", {
  skip_if_not_installed("NMF")
  mae <- make_model_mae(n_features = 12L, n_samples = 8L, seed = 104L)
  observed <- new.env(parent = emptyenv())

  local_mocked_bindings(
    nmf = function(x, rank, method, nrun, seed, ...) {
      observed$matrix <- x
      observed$rank <- rank
      observed$method <- method
      observed$runs <- nrun
      observed$seed <- seed
      structure(list(rank = rank), class = "bulkMAETestNMF")
    },
    .package = "NMF"
  )
  fit <- cluster_nmf(
    mae,
    "rna",
    assay = "counts",
    rank = 2L,
    method = "brunet",
    runs = 2L,
    seed = 11L
  )

  expect_s3_class(fit, "bulkMAETestNMF")
  expect_identical(dim(observed$matrix), c(12L, 8L))
  expect_identical(observed$rank, 2L)
  expect_identical(observed$method, "brunet")
  expect_identical(observed$runs, 2L)
  expect_identical(observed$seed, 11L)
})

test_that("NMF has a real-backend smoke test", {
  skip_if_not_installed("NMF")
  mae <- make_model_mae(n_features = 12L, n_samples = 8L, seed = 105L)
  fit <- suppressWarnings(cluster_nmf(
    mae,
    "rna",
    assay = "counts",
    rank = 2L,
    runs = 1L,
    seed = 16L,
    maxIter = 50L,
    .pbackend = NA
  ))
  expect_s4_class(fit, "NMFfit")
})

test_that("NMF rejects negative assays and invalid factorization settings", {
  skip_if_not_installed("NMF")
  mae <- make_model_mae(n_features = 6L, n_samples = 4L)
  counts <- mae_pull_assay(mae, "rna", "counts")
  counts[1L, 1L] <- -1
  negative <- network_mae_with_assay(mae, counts, "counts")

  expect_error(cluster_nmf(negative, "rna", "counts", rank = 2L), "non-negative")
  expect_error(cluster_nmf(mae, "rna", "counts", rank = 5L), "no larger")
  expect_error(cluster_nmf(mae, "rna", "counts", rank = 2L, runs = 0L), "runs")
})

test_that("WGCNA module wrapper keeps orientation and quality audit", {
  skip_if_not_installed("WGCNA")
  mae <- make_model_mae(n_features = 24L, n_samples = 12L, seed = 106L)
  observed <- new.env(parent = emptyenv())

  local_mocked_bindings(
    blockwiseModules = function(
        datExpr,
        power,
        networkType,
        TOMType,
        minModuleSize,
        mergeCutHeight,
        numericLabels,
        randomSeed,
        verbose,
        ...
    ) {
      observed$data <- datExpr
      observed$power <- power
      observed$network_type <- networkType
      list(colors = rep(1L, ncol(datExpr)))
    },
    .package = "WGCNA"
  )
  fit <- coexpr_wgcna(
    mae,
    "rna",
    assay = "log_expression",
    power = 2,
    top_n = 20L,
    min_module_size = 4L,
    seed = 12L,
    verbose = 0
  )

  expect_identical(dim(observed$data), c(12L, 20L))
  expect_identical(observed$power, 2)
  expect_identical(observed$network_type, "signed")
  expect_true(all(c(
    "bulkMAEQuality",
    "bulkMAERemovedSamples",
    "bulkMAERemovedFeatures",
    "bulkMAENetworkType"
  ) %in% names(fit)))
})

test_that("WGCNA module detection has a real-backend smoke test", {
  skip_if_not_installed("WGCNA")
  mae <- make_model_mae(n_features = 8L, n_samples = 8L, seed = 107L)
  fit <- suppressWarnings(coexpr_wgcna(
    mae,
    "rna",
    assay = "log_expression",
    power = 2,
    min_module_size = 2L,
    seed = 13L,
    verbose = 0,
    maxBlockSize = 10L,
    pamRespectsDendro = FALSE
  ))
  expect_true("colors" %in% names(fit))
})

test_that("WGCNA module detection validates core settings", {
  skip_if_not_installed("WGCNA")
  mae <- make_model_mae(n_features = 8L, n_samples = 8L)

  expect_error(coexpr_wgcna(mae, "rna", "log_expression", power = 0), "positive")
  expect_error(
    coexpr_wgcna(
      mae,
      "rna",
      "log_expression",
      power = 2,
      min_module_size = 9L
    ),
    "cannot exceed"
  )
  expect_error(
    coexpr_wgcna(
      mae,
      "rna",
      "log_expression",
      power = 2,
      numeric_labels = NA
    ),
    "TRUE or FALSE"
  )
})

test_that("WGCNA power selection runs on an MAE assay", {
  skip_if_not_installed("WGCNA")
  pick_power <- network_api("coexpr_pick_power")
  mae <- make_model_mae(n_features = 30L, n_samples = 12L, seed = 108L)

  fit <- suppressWarnings(pick_power(
    mae,
    "rna",
    assay = "log_expression",
    powers = c(1, 2),
    top_n = 20L,
    network_type = "signed",
    verbose = 0
  ))

  expect_type(fit, "list")
  expect_identical(as.numeric(fit$fitIndices$Power), c(1, 2))
  expect_identical(fit$bulkMAENetworkType, "signed")
  expect_true("bulkMAEQuality" %in% names(fit))
})

test_that("WGCNA power selection validates candidates", {
  skip_if_not_installed("WGCNA")
  pick_power <- network_api("coexpr_pick_power")
  mae <- make_model_mae(n_features = 8L, n_samples = 8L)

  expect_true(is.function(pick_power))
  expect_error(
    pick_power(mae, "rna", "log_expression", powers = c(1, 1)),
    "unique"
  )
  expect_error(
    pick_power(mae, "rna", "log_expression", powers = 1:2, r_squared = 1.1),
    "r_squared"
  )
})

test_that("WGCNA module preservation runs on aligned simulated cohorts", {
  skip_if_not_installed("WGCNA")
  reference <- make_model_mae(n_features = 30L, n_samples = 12L, seed = 109L)
  test <- make_model_mae(n_features = 30L, n_samples = 12L, seed = 110L)
  module_colors <- stats::setNames(
    rep(c("blue", "brown"), each = 10L),
    paste0("gene", seq_len(20L))
  )

  fit <- coexpr_preservation(
    reference,
    "rna",
    test,
    "rna",
    module_colors = module_colors,
    reference_assay = "log_expression",
    test_assay = "log_expression",
    permutations = 1L,
    seed = 14L,
    verbose = 0
  )

  expect_type(fit, "list")
  expect_true("preservation" %in% names(fit))
  expect_identical(
    fit$bulkMAEInputQuality$analyzedFeatures,
    paste0("gene", seq_len(20L))
  )
})

test_that("WGCNA preservation rejects unaligned module declarations", {
  skip_if_not_installed("WGCNA")
  reference <- make_model_mae(n_features = 8L, n_samples = 8L, seed = 111L)
  test <- make_model_mae(n_features = 8L, n_samples = 8L, seed = 112L)

  expect_error(
    coexpr_preservation(
      reference,
      "rna",
      test,
      "rna",
      module_colors = rep("blue", 8L),
      reference_assay = "log_expression",
      test_assay = "log_expression"
    ),
    "named"
  )
  expect_error(
    coexpr_preservation(
      reference,
      "rna",
      test,
      "rna",
      module_colors = c(other1 = "blue", other2 = "brown"),
      reference_assay = "log_expression",
      test_assay = "log_expression"
    ),
    "at least three common"
  )
  expect_error(
    coexpr_preservation(
      reference,
      "rna",
      test,
      "rna",
      module_colors = stats::setNames(rep("blue", 8L), paste0("gene", 1:8)),
      reference_assay = "log_expression",
      test_assay = "log_expression",
      permutations = 0L
    ),
    "permutations"
  )
})

test_that("GENIE3 inference and link ranking run offline", {
  skip_if_not_installed("GENIE3")
  mae <- make_model_mae(n_features = 16L, n_samples = 10L, seed = 113L)

  weights <- network_genie3(
    mae,
    "rna",
    assay = "log_expression",
    regulators = paste0("gene", 1:3),
    targets = paste0("gene", 4:10),
    trees = 10L,
    cores = 1L,
    seed = 15L
  )
  links <- network_genie3_links(weights, threshold = 0, top = 5L)

  expect_identical(dim(weights), c(3L, 7L))
  expect_true(all(is.finite(weights)))
  expect_s3_class(links, "data.frame")
  expect_lte(nrow(links), 5L)
  expect_true(all(c("regulatoryGene", "targetGene", "weight") %in% names(links)))
})

test_that("GENIE3 validates computation and link filters", {
  skip_if_not_installed("GENIE3")
  mae <- make_model_mae(n_features = 8L, n_samples = 8L)

  expect_error(
    network_genie3(mae, "rna", "log_expression", trees = 0L),
    "trees"
  )
  expect_error(
    network_genie3(mae, "rna", "log_expression", cores = 1.5),
    "cores"
  )
  expect_error(network_genie3_links(matrix(1, 2, 2), threshold = -1), "threshold")
  expect_error(network_genie3_links(matrix(1, 2, 2), top = 0L), "top")
})

test_that("STRING wrapper performs all offline preflight checks", {
  skip_if_not_installed("STRINGdb")
  mae <- make_model_mae(n_features = 8L, n_samples = 6L)

  expect_error(
    network_string(mae, "rna", "log_expression", genes = character()),
    "genes"
  )
  expect_error(
    network_string(mae, "rna", "log_expression", species = 0L),
    "species"
  )
  expect_error(
    network_string(mae, "rna", "log_expression", version = ""),
    "version"
  )
  expect_error(
    network_string(mae, "rna", "log_expression", score_threshold = 1001),
    "score_threshold"
  )
  expect_error(
    network_string(mae, "rna", "log_expression", input_directory = NA_character_),
    "input_directory"
  )
  expect_error(
    network_string(mae, "rna", "log_expression", remove_unmapped = NA),
    "remove_unmapped"
  )
})

test_that("STRING mapping integration is explicitly online", {
  skip_if_not_installed("STRINGdb")
  skip_if(
    !identical(
      tolower(Sys.getenv("BULKMAE_RUN_ONLINE_TESTS", "false")),
      "true"
    ),
    paste(
      "STRING identifier mapping may download species resources; set",
      "BULKMAE_RUN_ONLINE_TESTS=true to opt in."
    )
  )
  mae <- make_model_mae(n_features = 6L, n_samples = 6L)
  expression <- mae_pull_assay(mae, "rna", "log_expression")
  rownames(expression) <- c("TP53", "EGFR", "MYC", "AKT1", "MAPK1", "STAT3")
  mae <- network_mae_with_assay(mae, expression)

  result <- network_string(
    mae,
    "rna",
    assay = "expression",
    genes = rownames(expression)
  )

  expect_s3_class(result, "data.frame")
  expect_true(all(c("from", "to") %in% names(result)))
})
