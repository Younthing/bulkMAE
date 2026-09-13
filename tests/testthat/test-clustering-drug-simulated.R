.make_small_clustering_mae <- function(seed = 501L) {
  source <- make_simulated_mae(seed = seed)
  source_samples <- mae_samples(source, "rna")
  selected_samples <- c(seq_len(4L), 13:16)
  sample_names <- rownames(source_samples)[selected_samples]
  feature_names <- sprintf("gene%04d", seq_len(30L))
  tpm <- mae_pull_assay(source, "rna", "tpm")[
    feature_names,
    sample_names,
    drop = FALSE
  ]
  sample_data <- source_samples[sample_names, , drop = FALSE]
  experiment <- mae_create_experiment(
    assays = list(tpm = tpm),
    col_data = data.frame(row.names = sample_names)
  )
  mae_create(list(rna = experiment), sample_data)
}

test_that("cluster_nmf fits a reproducible non-negative factorization", {
  skip_if_not_installed("NMF")
  mae <- .make_small_clustering_mae()

  set.seed(840)
  caller_seed <- .Random.seed
  search_path <- search()
  fit <- cluster_nmf(
    mae,
    "rna",
    assay = "tpm",
    rank = 2L,
    runs = 1L,
    seed = 41L,
    maxIter = 30L
  )

  expect_identical(.Random.seed, caller_seed)
  expect_identical(search(), search_path)
  expect_true(methods::is(fit, "NMFfit"))
  expect_identical(dim(NMF::basis(fit)), c(30L, 2L))
  expect_identical(dim(NMF::coef(fit)), c(2L, 8L))
  expect_true(all(is.finite(NMF::basis(fit))))
  expect_true(all(is.finite(NMF::coef(fit))))
  expect_true(all(NMF::basis(fit) >= 0))
  expect_true(all(NMF::coef(fit) >= 0))
  expect_error(
    cluster_nmf(mae, "rna", assay = "tpm", rank = 0L, runs = 1L),
    "rank.*positive integers"
  )
})

test_that("cluster_nmf supports repeated runs through a namespace-only backend", {
  skip_if_not_installed("NMF")
  mae <- .make_small_clustering_mae()

  fit <- cluster_nmf(
    mae,
    "rna",
    assay = "tpm",
    rank = 2L,
    runs = 2L,
    seed = 43L,
    maxIter = 10L,
    .pbackend = NA
  )

  expect_true(methods::is(fit, "NMFfit"))
})

test_that("cluster_nmf rejects negative expression before backend fitting", {
  skip_if_not_installed("NMF")
  mae <- .make_small_clustering_mae()
  sample_data <- mae_samples(mae, "rna")
  expression <- mae_pull_assay(mae, "rna", "tpm")
  expression[1L, 1L] <- -1
  experiment <- mae_create_experiment(
    assays = list(expression = expression),
    col_data = data.frame(row.names = rownames(sample_data))
  )
  invalid <- mae_create(list(rna = experiment), sample_data)

  expect_error(
    cluster_nmf(invalid, "rna", assay = "expression", rank = 2L, runs = 1L),
    "finite, non-negative"
  )
})

test_that("drug_query creates ranked disjoint LINCS inputs", {
  statistics <- c(
    TP53 = 4.2,
    EGFR = 2.1,
    AKT1 = 0.5,
    MYC = -3.8,
    MAPK1 = -1.7,
    STAT3 = -0.2
  )

  query <- drug_query(statistics, n = 2L)

  expect_identical(query$upset, c("TP53", "EGFR"))
  expect_identical(query$downset, c("MYC", "MAPK1"))
  expect_length(intersect(query$upset, query$downset), 0L)
})

test_that("drug_lincs validates local query options without opening a database", {
  suppressWarnings(skip_if_not_installed("signatureSearch"))
  valid_query <- list(upset = c("TP53", "EGFR"), downset = c("MYC", "MAPK1"))

  expect_error(
    drug_lincs(list(upset = "TP53"), reference_database = "unused"),
    "upset.*downset"
  )
  expect_error(
    drug_lincs(
      list(upset = c("TP53", "EGFR"), downset = c("TP53", "MYC")),
      reference_database = "unused"
    ),
    "must be disjoint"
  )
  expect_error(
    drug_lincs(valid_query, reference_database = "unused", workers = 0L),
    "workers.*positive integer"
  )
  expect_error(
    drug_lincs(valid_query, reference_database = "unused", tau = NA),
    "tau.*annotations"
  )
})

test_that("drug_lincs_table ranks native-style columns without rescoring", {
  raw <- drug_lincs_example()
  expect_identical(attr(raw, "bulkmae_lincs_source", exact = TRUE), "diagnostic_toy")
  expect_identical(
    names(raw),
    c(
      "pert", "cell", "type", "trend", "WTCS", "WTCS_Pval", "WTCS_FDR",
      "NCS", "NCSct", "Tau", "N_upset", "N_downset"
    )
  )

  table <- drug_lincs_table(raw)
  expect_identical(attr(table, "score_column", exact = TRUE), "NCS")
  expect_identical(
    attr(table, "score_label", exact = TRUE),
    "Normalized connectivity score (NCS)"
  )
  expect_identical(attr(table, "bulkmae_lincs_source", exact = TRUE), "diagnostic_toy")
  expect_identical(table$NCS, sort(raw$NCS, decreasing = TRUE))
  expect_identical(
    as.character(table$direction[table$NCS < 0]),
    rep("Reverse", sum(table$NCS < 0))
  )
  expect_identical(
    as.character(table$direction[table$NCS > 0]),
    rep("Mimic", sum(table$NCS > 0))
  )

  tau <- drug_lincs_table(raw, score = "Tau")
  expect_identical(attr(tau, "score_column", exact = TRUE), "Tau")
  expect_identical(tau$Tau, sort(raw$Tau, decreasing = TRUE))
  expect_error(drug_lincs_table(raw, score = "missing"), "not in `result`")
  expect_error(drug_lincs_table(raw[, setdiff(names(raw), c("NCS", "Tau", "WTCS", "NCSct")), drop = FALSE]),
               "finite NCS, Tau, WTCS, or NCSct")
})

test_that("drug_lincs_table reads a gessResult slot when the class exists", {
  raw <- drug_lincs_example()
  expect_error(drug_lincs_table(list()), "gessResult or a LINCS result")
  skip_if_not_installed("signatureSearch")
  if (!methods::isClass("gessResult")) {
    skip("Installed signatureSearch does not define gessResult.")
  }
  fake <- methods::new("gessResult")
  methods::slot(fake, "result") <- raw
  table <- drug_lincs_table(fake)
  expect_identical(table$NCS, sort(raw$NCS, decreasing = TRUE))
})

test_that("drug_lincs execution requires an explicit local reference database", {
  suppressWarnings(skip_if_not_installed("signatureSearch"))
  reference_database <- Sys.getenv("BULKMAE_LINCS_DB")
  skip_if(
    !nzchar(reference_database),
    "Set BULKMAE_LINCS_DB to run the local LINCS database integration test."
  )

  result <- drug_lincs(
    list(upset = c("TP53", "EGFR"), downset = c("MYC", "MAPK1")),
    reference_database = reference_database,
    workers = 1L,
    tau = FALSE,
    annotations = FALSE
  )

  expect_true(methods::is(result, "gessResult") || inherits(result, "gessResult"))
})
