test_that("simulated MAE construction and access preserve the public contract", {
  mae <- make_simulated_mae()

  expect_s4_class(mae, "MultiAssayExperiment")
  expect_identical(mae_validate(mae, "rna"), mae)
  expect_identical(mae_experiments(mae), "rna")
  expect_identical(
    mae_assays(mae, "rna"),
    c("counts", "log_expression", "tpm", "weights")
  )

  experiment <- mae_pull_experiment(mae, "rna")
  counts <- mae_pull_assay(mae, "rna", "counts")
  expect_s4_class(experiment, "SummarizedExperiment")
  expect_identical(counts, mae_pull_assay(mae, "rna", 1L))
  expect_identical(dim(counts), c(1200L, 24L))
  expect_identical(colnames(counts), rownames(mae_samples(mae, "rna")))
  expect_identical(
    rownames(SummarizedExperiment::rowData(experiment)),
    rownames(counts)
  )
})

test_that("MAE accessors reject objects, experiments, and assays outside the contract", {
  mae <- make_simulated_mae()

  expect_error(mae_validate(list(), "rna"), "MultiAssayExperiment")
  expect_error(mae_validate(mae, "missing"), "Unknown experiment")
  expect_error(mae_pull_experiment(mae, "missing"), "Unknown experiment")
  expect_error(mae_pull_assay(mae, "rna", "missing"), "Unknown assay")
  expect_error(mae_pull_assay(mae, "rna", 0L), "name or index")
  expect_error(mae_pull_assay(mae, "rna", 99L), "name or index")
  expect_error(mae_pull_assay(mae, "rna", c(1L, 2L)), "name or index")
})

test_that("mae_create_experiment constructs aligned multi-assay leaves", {
  sample_names <- paste0("sample", 1:4)
  feature_names <- paste0("gene", 1:6)
  counts <- matrix(
    seq_len(24L),
    nrow = 6L,
    dimnames = list(feature_names, sample_names)
  )
  samples <- data.frame(
    group = rep(c("a", "b"), each = 2L),
    row.names = sample_names
  )
  features <- data.frame(
    length = seq(500, 1000, length.out = 6L),
    row.names = feature_names
  )

  experiment <- mae_create_experiment(
    assays = list(counts = counts, log_expression = log2(counts + 1)),
    col_data = samples,
    row_data = features
  )
  expect_s4_class(experiment, "SummarizedExperiment")
  expect_identical(
    SummarizedExperiment::assayNames(experiment),
    c("counts", "log_expression")
  )
  expect_identical(colnames(experiment), sample_names)
  expect_identical(rownames(experiment), feature_names)
})

test_that("MAE constructors reject misaligned or malformed inputs", {
  samples <- data.frame(group = c("a", "b"), row.names = c("s1", "s2"))
  counts <- matrix(
    1:6,
    nrow = 3L,
    dimnames = list(paste0("g", 1:3), c("s1", "s2"))
  )

  expect_error(mae_create_experiment(list(), samples), "at least one assay")
  expect_error(
    mae_create_experiment(unname(list(counts)), samples),
    "unique, non-empty names"
  )

  mismatched <- counts[, c("s2", "s1"), drop = FALSE]
  expect_error(
    mae_create_experiment(
      list(counts = counts, second = mismatched),
      samples
    ),
    "identical dimensions and dimnames"
  )

  duplicated <- counts
  rownames(duplicated) <- c("g1", "g1", "g3")
  expect_error(
    mae_create_experiment(duplicated, samples),
    "feature names must be non-empty and unique"
  )

  wrong_samples <- samples[c("s2", "s1"), , drop = FALSE]
  expect_error(
    mae_create_experiment(counts, wrong_samples),
    "exactly match"
  )

  experiment <- mae_create_experiment(counts, samples)
  expect_error(mae_create(list(), samples), "at least one experiment")
  expect_error(mae_create(unname(list(experiment)), samples), "unique, non-empty")
  expect_error(mae_create(list(rna = counts), samples), "SummarizedExperiment")
})

test_that("annotate_ensembl removes only terminal numeric version suffixes", {
  identifiers <- c(
    "ENSG00000141510.18",
    "ENST00000269305.9",
    "TP53",
    "gene.2.extra"
  )

  expect_identical(
    annotate_ensembl(identifiers),
    c(
      "ENSG00000141510",
      "ENST00000269305",
      "TP53",
      "gene.2.extra"
    )
  )
  expect_identical(annotate_ensembl(c(101L, 202L)), c("101", "202"))
  expect_error(annotate_ensembl(character()), "non-missing identifiers")
  expect_error(annotate_ensembl(c("ENSG1.1", NA_character_)), "non-missing")
  expect_error(annotate_ensembl(list("ENSG1.1")), "atomic vector")
})

test_that("import_tximport imports minimal Salmon quantification files", {
  skip_if_not_installed("tximport")

  directory <- tempfile("bulkMAE-salmon-")
  dir.create(directory)
  on.exit(unlink(directory, recursive = TRUE), add = TRUE)

  sample_names <- c("sample_a", "sample_b")
  transcript_names <- paste0("tx", 1:4)
  files <- character(length(sample_names))
  for (index in seq_along(sample_names)) {
    sample_directory <- file.path(directory, sample_names[[index]])
    dir.create(sample_directory)
    files[[index]] <- file.path(sample_directory, "quant.sf")
    quantification <- data.frame(
      Name = transcript_names,
      Length = c(1000, 1500, 800, 2000),
      EffectiveLength = c(800, 1300, 600, 1800),
      TPM = c(100000, 200000, 300000, 400000),
      NumReads = c(20, 40, 60, 80) * index,
      check.names = FALSE
    )
    utils::write.table(
      quantification,
      file = files[[index]],
      sep = "\t",
      quote = FALSE,
      row.names = FALSE
    )
  }
  names(files) <- sample_names
  samples <- data.frame(
    condition = c("control", "treated"),
    row.names = sample_names
  )

  experiment <- import_tximport(
    files,
    samples,
    type = "salmon",
    counts_from_abundance = "no",
    txOut = TRUE
  )
  expect_s4_class(experiment, "SummarizedExperiment")
  expect_setequal(
    SummarizedExperiment::assayNames(experiment),
    c("counts", "abundance", "length")
  )
  expect_identical(dim(experiment), c(4L, 2L))
  expect_identical(rownames(experiment), transcript_names)
  expect_identical(colnames(experiment), sample_names)
  expect_identical(
    S4Vectors::metadata(experiment)$bulkMAE$counts_from_abundance,
    "no"
  )
})

test_that("import_tximport rejects unnamed and misaligned sample inputs", {
  skip_if_not_installed("tximport")

  directory <- tempfile("bulkMAE-tximport-contract-")
  dir.create(directory)
  on.exit(unlink(directory, recursive = TRUE), add = TRUE)
  samples <- data.frame(row.names = c("sample_a", "sample_b"))
  files <- stats::setNames(
    file.path(directory, c("sample_a.sf", "sample_b.sf")),
    rownames(samples)
  )
  file.create(files)

  expect_error(
    import_tximport(unname(files), samples, type = "salmon"),
    "Sample names in `files`"
  )
  expect_error(
    import_tximport(files, samples[c("sample_b", "sample_a"), , drop = FALSE]),
    "exactly match"
  )
  invalid <- files
  invalid[[1L]] <- ""
  expect_error(import_tximport(invalid, samples), "non-empty path")
})

test_that("mae_add_assay accepts matrix, SummarizedExperiment, and EList values", {
  if (!exists("mae_add_assay", envir = asNamespace("bulkMAE"), inherits = FALSE)) {
    skip("mae_add_assay() is not implemented yet")
  }
  mae_add_assay <- getFromNamespace("mae_add_assay", "bulkMAE")
  mae <- make_simulated_mae()
  counts <- mae_pull_assay(mae, "rna", "counts")
  transformed <- log2(counts + 1)
  shuffled <- transformed[
    rev(rownames(transformed)),
    rev(colnames(transformed)),
    drop = FALSE
  ]

  expect_silent(
    with_matrix <- mae_add_assay(
      mae,
      "rna",
      value = shuffled,
      name = "log2_counts"
    )
  )
  expect_identical(
    mae_pull_assay(with_matrix, "rna", "log2_counts"),
    transformed
  )
  expect_false("log2_counts" %in% mae_assays(mae, "rna"))

  value_experiment <- SummarizedExperiment::SummarizedExperiment(
    assays = list(value = transformed)
  )
  expect_silent(
    with_experiment <- mae_add_assay(
      mae,
      "rna",
      value = value_experiment,
      name = "from_experiment"
    )
  )
  expect_identical(
    mae_pull_assay(with_experiment, "rna", "from_experiment"),
    transformed
  )

  value_elist <- structure(list(E = transformed), class = "EList")
  expect_silent(
    with_elist <- mae_add_assay(
      mae,
      "rna",
      value = value_elist,
      name = "from_elist"
    )
  )
  expect_identical(
    mae_pull_assay(with_elist, "rna", "from_elist"),
    transformed
  )

  wrong <- transformed[-1L, , drop = FALSE]
  expect_error(
    mae_add_assay(mae, "rna", value = wrong, name = "wrong"),
    "dimension|feature"
  )
  expect_error(
    mae_add_assay(mae, "rna", value = transformed, name = "counts"),
    "already exists"
  )
  multiple <- SummarizedExperiment::SummarizedExperiment(
    assays = list(first = transformed, second = transformed)
  )
  expect_error(
    mae_add_assay(mae, "rna", value = multiple, name = "multiple"),
    "exactly one assay"
  )
  expect_error(
    mae_add_assay(
      mae,
      "rna",
      value = structure(list(), class = "EList"),
      name = "empty_elist"
    ),
    "contain an `E` matrix"
  )
})

test_that("mae_subset_features subsets every assay and preserves requested order", {
  if (
    !exists("mae_subset_features", envir = asNamespace("bulkMAE"), inherits = FALSE)
  ) {
    skip("mae_subset_features() is not implemented yet")
  }
  mae_subset_features <- getFromNamespace("mae_subset_features", "bulkMAE")
  mae <- make_simulated_mae()
  features <- rev(sprintf("gene%04d", 1:20))

  expect_silent(subset <- mae_subset_features(mae, "rna", features))
  experiment <- mae_pull_experiment(subset, "rna")
  expect_identical(rownames(experiment), features)
  for (assay in mae_assays(subset, "rna")) {
    expect_identical(rownames(mae_pull_assay(subset, "rna", assay)), features)
  }
  expect_identical(
    rownames(SummarizedExperiment::rowData(experiment)),
    features
  )
  expect_identical(nrow(mae_pull_experiment(mae, "rna")), 1200L)
  expect_error(
    mae_subset_features(mae, "rna", c(features, "missing_gene")),
    "Unknown features|missing_gene"
  )
  expect_error(
    mae_subset_features(mae, "rna", c(features[[1L]], features[[1L]])),
    "unique"
  )
})

test_that("normalize_tpm aligns named feature lengths and closes each library", {
  if (!exists("normalize_tpm", envir = asNamespace("bulkMAE"), inherits = FALSE)) {
    skip("normalize_tpm() is not implemented yet")
  }
  normalize_tpm <- getFromNamespace("normalize_tpm", "bulkMAE")
  mae <- make_simulated_mae()
  experiment <- mae_pull_experiment(mae, "rna")
  lengths <- SummarizedExperiment::rowData(experiment)$gene_length
  names(lengths) <- rownames(experiment)
  lengths <- lengths[rev(names(lengths))]

  observed <- normalize_tpm(mae, "rna", lengths = lengths, assay = "counts")
  expected <- mae_pull_assay(mae, "rna", "tpm")
  expect_identical(dim(observed), dim(expected))
  expect_identical(dimnames(observed), dimnames(expected))
  expect_equal(
    unname(colSums(observed)),
    rep(1e6, ncol(observed)),
    tolerance = 1e-7
  )
  expect_equal(observed, expected, tolerance = 1e-7)
  expect_equal(
    normalize_tpm(mae, "rna", lengths = "gene_length", assay = "counts"),
    expected,
    tolerance = 1e-7
  )

  invalid <- lengths
  invalid[[1L]] <- 0
  expect_error(
    normalize_tpm(mae, "rna", lengths = invalid, assay = "counts"),
    "length"
  )
  expect_error(
    normalize_tpm(mae, "rna", lengths = lengths[-1L], assay = "counts"),
    "length|feature"
  )
})
