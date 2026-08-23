make_deconvolution_mae <- function(
    n_genes = 120L,
    n_samples = 4L,
    seed = 201L
) {
  set.seed(seed)
  genes <- paste0("gene", seq_len(n_genes))
  samples <- paste0("sample", seq_len(n_samples))
  counts <- matrix(
    stats::rpois(n_genes * n_samples, lambda = 80),
    nrow = n_genes,
    dimnames = list(genes, samples)
  )
  storage.mode(counts) <- "integer"
  tpm <- sweep(counts, 2L, colSums(counts), FUN = "/") * 1e6
  sample_data <- data.frame(
    condition = rep(c("control", "treated"), length.out = n_samples),
    row.names = samples
  )
  experiment <- mae_create_experiment(
    assays = list(counts = counts, tpm = tpm),
    col_data = sample_data
  )
  mae_create(list(rna = experiment), sample_data)
}

deconvolution_mae_with_assay <- function(x, value, assay_name) {
  samples <- mae_samples(x, "rna")
  experiment <- mae_create_experiment(
    assays = stats::setNames(list(value), assay_name),
    col_data = samples
  )
  mae_create(list(rna = experiment), samples)
}

make_music_reference <- function(
    genes,
    n_donors = 4L,
    cells_per_type = 2L,
    seed = 202L
) {
  set.seed(seed)
  cell_types <- c("B_cell", "T_cell")
  cells_per_donor <- length(cell_types) * cells_per_type
  donors <- paste0("donor", seq_len(n_donors))
  cell_names <- paste0("cell", seq_len(n_donors * cells_per_donor))
  cell_type <- rep(rep(cell_types, each = cells_per_type), n_donors)
  donor <- rep(donors, each = cells_per_donor)
  counts <- matrix(
    stats::rpois(length(genes) * length(cell_names), lambda = 20),
    nrow = length(genes),
    dimnames = list(genes, cell_names)
  )
  SingleCellExperiment::SingleCellExperiment(
    assays = list(counts = counts),
    colData = S4Vectors::DataFrame(
      cell_type = cell_type,
      donor = donor,
      row.names = cell_names
    )
  )
}

make_bayesprism_reference <- function(genes, n_cells = 8L, seed = 203L) {
  set.seed(seed)
  matrix(
    stats::rpois(length(genes) * n_cells, lambda = 20),
    nrow = length(genes),
    dimnames = list(genes, paste0("cell", seq_len(n_cells)))
  )
}

test_that("immunedeconv wrapper passes a named non-negative expression matrix", {
  suppressWarnings(skip_if_not_installed("immunedeconv"))
  mae <- make_deconvolution_mae(seed = 204L)
  observed <- new.env(parent = emptyenv())

  local_mocked_bindings(
    .call_backend = function(package, function_name, arguments) {
      observed$package <- package
      observed$function_name <- function_name
      observed$arguments <- arguments
      data.frame(cell_type = "B_cell", sample1 = 0.5, check.names = FALSE)
    },
    .package = "bulkMAE"
  )
  result <- deconv(
    mae,
    "rna",
    method = "quantiseq",
    assay = "tpm",
    tumor = FALSE
  )

  expect_s3_class(result, "data.frame")
  expect_identical(observed$package, "immunedeconv")
  expect_identical(observed$function_name, "deconvolute")
  expect_identical(observed$arguments$method, "quantiseq")
  expect_identical(dim(observed$arguments$gene_expression), c(120L, 4L))
  expect_identical(rownames(observed$arguments$gene_expression), paste0("gene", 1:120))
  expect_false(observed$arguments$tumor)
})

test_that("xCell runs from a namespace-only session and restores the search path", {
  suppressWarnings(skip_if_not_installed("immunedeconv"))
  resource <- suppressMessages(suppressWarnings(
    getExportedValue("immunedeconv", "xCell.data")
  ))
  genes <- resource$genes
  samples <- paste0("sample", seq_len(4L))
  tpm <- matrix(
    (seq_len(length(genes) * length(samples)) %% 101L) + 1,
    nrow = length(genes),
    dimnames = list(genes, samples)
  )
  tpm <- sweep(tpm, 2L, colSums(tpm), FUN = "/") * 1e6
  sample_data <- data.frame(row.names = samples)
  mae <- mae_create(
    list(rna = mae_create_experiment(list(tpm = tpm), sample_data)),
    sample_data
  )
  search_path <- search()

  result <- suppressMessages(suppressWarnings(
    deconv(mae, "rna", method = "xcell", assay = "tpm")
  ))

  expect_identical(search(), search_path)
  expect_s3_class(result, "data.frame")
  expect_identical(names(result), c("cell_type", samples))
  expect_gt(nrow(result), 0L)
  expect_true(all(is.finite(as.matrix(result[-1L]))))
})

test_that("immunedeconv wrapper documents xCell score semantics", {
  suppressWarnings(skip_if_not_installed("immunedeconv"))
  mae <- make_deconvolution_mae(n_genes = 20L, seed = 205L)
  local_mocked_bindings(
    .call_backend = function(package, function_name, arguments) list(),
    .package = "bulkMAE"
  )

  expect_message(
    deconv(mae, "rna", method = "xcell", assay = "tpm"),
    "enrichment score"
  )
})

test_that("immunedeconv wrapper rejects invalid methods and expression values", {
  suppressWarnings(skip_if_not_installed("immunedeconv"))
  mae <- make_deconvolution_mae(n_genes = 20L, seed = 206L)
  tpm <- mae_pull_assay(mae, "rna", "tpm")
  tpm[1L, 1L] <- -1
  negative <- deconvolution_mae_with_assay(mae, tpm, "tpm")

  expect_error(deconv(mae, "rna", method = "", assay = "tpm"), "method")
  expect_error(
    deconv(negative, "rna", method = "quantiseq", assay = "tpm"),
    "non-negative"
  )
})

test_that("MuSiC adapter passes aligned bulk and single-cell inputs", {
  skip_if_not_installed("MuSiC")
  skip_if_not_installed("SingleCellExperiment")
  mae <- make_deconvolution_mae(seed = 207L)
  genes <- rownames(mae_pull_assay(mae, "rna", "counts"))
  reference <- make_music_reference(genes, seed = 208L)
  observed <- new.env(parent = emptyenv())

  local_mocked_bindings(
    .call_backend = function(package, function_name, arguments) {
      observed$package <- package
      observed$function_name <- function_name
      observed$arguments <- arguments
      list(
        Est.prop.weighted = matrix(
          0.5,
          nrow = ncol(arguments$bulk.mtx),
          ncol = 2L
        )
      )
    },
    .package = "bulkMAE"
  )
  result <- deconv_music(
    mae,
    "rna",
    sc_reference = reference,
    clusters = "cell_type",
    samples = "donor",
    assay = "counts",
    iter.max = 50L
  )

  expect_identical(observed$package, "MuSiC")
  expect_identical(observed$function_name, "music_prop")
  expect_identical(dim(observed$arguments$bulk.mtx), c(120L, 4L))
  expect_identical(observed$arguments$sc.sce, reference)
  expect_identical(observed$arguments$clusters, "cell_type")
  expect_identical(observed$arguments$samples, "donor")
  expect_identical(observed$arguments$iter.max, 50L)
  expect_identical(dim(result$Est.prop.weighted), c(4L, 2L))
})

test_that("MuSiC has a real-backend smoke test", {
  skip_if_not_installed("MuSiC")
  skip_if_not_installed("SingleCellExperiment")
  mae <- make_deconvolution_mae(seed = 209L)
  genes <- rownames(mae_pull_assay(mae, "rna", "counts"))
  reference <- make_music_reference(genes, seed = 210L)

  invisible(utils::capture.output(
    result <- suppressMessages(suppressWarnings(deconv_music(
      mae,
      "rna",
      sc_reference = reference,
      clusters = "cell_type",
      samples = "donor",
      assay = "counts"
    )))
  ))
  expect_identical(dim(result$Est.prop.weighted), c(4L, 2L))
})

test_that("MuSiC validates reference metadata and identifier overlap", {
  skip_if_not_installed("MuSiC")
  skip_if_not_installed("SingleCellExperiment")
  mae <- make_deconvolution_mae(seed = 211L)
  genes <- rownames(mae_pull_assay(mae, "rna", "counts"))
  reference <- make_music_reference(genes, seed = 212L)

  expect_error(
    deconv_music(mae, "rna", matrix(1, 2, 2), "cell_type", "donor"),
    "SingleCellExperiment"
  )
  expect_error(
    deconv_music(mae, "rna", reference, "missing", "donor"),
    "metadata is missing"
  )

  one_type <- reference
  SummarizedExperiment::colData(one_type)$cell_type <- "only_type"
  expect_error(
    deconv_music(mae, "rna", one_type, "cell_type", "donor"),
    "at least two reference cell types"
  )

  no_overlap <- reference
  rownames(no_overlap) <- paste0("other", seq_len(nrow(no_overlap)))
  expect_error(
    deconv_music(mae, "rna", no_overlap, "cell_type", "donor"),
    "no common feature"
  )
})

test_that("CIBERSORTx mixture preparation is complete and offline", {
  mae <- make_deconvolution_mae(n_genes = 25L, n_samples = 4L, seed = 213L)

  mixture <- deconv_cibersortx_input(
    mae,
    "rna",
    assay = "tpm",
    gene_column = "GeneSymbol"
  )

  expect_s3_class(mixture, "data.frame")
  expect_identical(dim(mixture), c(25L, 5L))
  expect_identical(names(mixture), c("GeneSymbol", paste0("sample", 1:4)))
  expect_identical(mixture$GeneSymbol, paste0("gene", 1:25))
  expect_true(all(vapply(mixture[-1L], is.numeric, logical(1))))
})

test_that("CIBERSORTx preparation validates values and identifier-column name", {
  mae <- make_deconvolution_mae(n_genes = 20L, seed = 214L)
  tpm <- mae_pull_assay(mae, "rna", "tpm")
  tpm[1L, 1L] <- Inf
  non_finite <- deconvolution_mae_with_assay(mae, tpm, "tpm")

  expect_error(
    deconv_cibersortx_input(mae, "rna", "tpm", gene_column = ""),
    "gene_column"
  )
  expect_error(
    deconv_cibersortx_input(non_finite, "rna", "tpm"),
    "finite, non-negative"
  )
})

test_that("BayesPrism adapter aligns both inputs as observations by genes", {
  suppressWarnings(skip_if_not_installed("BayesPrism"))
  mae <- make_deconvolution_mae(seed = 215L)
  genes <- rownames(mae_pull_assay(mae, "rna", "counts"))
  reference <- make_bayesprism_reference(genes, seed = 216L)
  cells <- colnames(reference)
  labels <- stats::setNames(
    rep(c("B_cell", "T_cell"), length.out = length(cells)),
    rev(cells)
  )
  observed <- new.env(parent = emptyenv())
  observed$calls <- list()

  local_mocked_bindings(
    .call_backend = function(package, function_name, arguments) {
      observed$calls[[length(observed$calls) + 1L]] <- list(
        package = package,
        function_name = function_name,
        arguments = arguments
      )
      if (identical(function_name, "new.prism")) {
        return(structure(list(id = "prism"), class = "bulkMAETestPrism"))
      }
      structure(list(status = "complete"), class = "bulkMAETestBayesResult")
    },
    .package = "bulkMAE"
  )
  result <- deconv_bayesprism(
    mae,
    "rna",
    reference = reference,
    cell_type_labels = labels,
    assay = "counts",
    reference_orientation = "genes_by_cells",
    cores = 2L,
    outlier.cut = 0.02
  )

  constructor <- observed$calls[[1L]]
  runner <- observed$calls[[2L]]
  expect_s3_class(result, "bulkMAETestBayesResult")
  expect_identical(constructor$package, "BayesPrism")
  expect_identical(constructor$function_name, "new.prism")
  expect_identical(dim(constructor$arguments$mixture), c(4L, 120L))
  expect_identical(dim(constructor$arguments$reference), c(8L, 120L))
  expect_identical(rownames(constructor$arguments$mixture), paste0("sample", 1:4))
  expect_identical(rownames(constructor$arguments$reference), cells)
  expect_identical(
    constructor$arguments$cell.type.labels,
    unname(labels[cells])
  )
  expect_identical(constructor$arguments$outlier.cut, 0.02)
  expect_identical(runner$function_name, "run.prism")
  expect_s3_class(runner$arguments$prism, "bulkMAETestPrism")
  expect_identical(runner$arguments$n.cores, 2L)
  expect_true(runner$arguments$update.gibbs)
  expect_identical(runner$arguments$gibbs.control, list())
  expect_identical(runner$arguments$opt.control, list())
})

test_that("BayesPrism accepts an explicitly cells-by-genes reference", {
  suppressWarnings(skip_if_not_installed("BayesPrism"))
  mae <- make_deconvolution_mae(seed = 217L)
  genes <- rownames(mae_pull_assay(mae, "rna", "counts"))
  reference <- t(make_bayesprism_reference(genes, seed = 218L))
  observed <- new.env(parent = emptyenv())

  local_mocked_bindings(
    .call_backend = function(package, function_name, arguments) {
      if (identical(function_name, "new.prism")) observed$reference <- arguments$reference
      list()
    },
    .package = "bulkMAE"
  )
  deconv_bayesprism(
    mae,
    "rna",
    reference = reference,
    cell_type_labels = rep(c("B_cell", "T_cell"), 4L),
    reference_orientation = "cells_by_genes"
  )

  expect_equal(observed$reference, reference)
  expect_identical(dimnames(observed$reference), dimnames(reference))
})

test_that("BayesPrism validates raw counts, labels, overlap, and workers", {
  suppressWarnings(skip_if_not_installed("BayesPrism"))
  mae <- make_deconvolution_mae(seed = 219L)
  genes <- rownames(mae_pull_assay(mae, "rna", "counts"))
  reference <- make_bayesprism_reference(genes, seed = 220L)
  labels <- rep(c("B_cell", "T_cell"), 4L)

  expect_error(
    deconv_bayesprism(mae, "rna", reference, labels, assay = "tpm"),
    "raw integer counts"
  )

  non_integer_reference <- reference
  non_integer_reference[1L, 1L] <- non_integer_reference[1L, 1L] + 0.5
  expect_error(
    deconv_bayesprism(mae, "rna", non_integer_reference, labels),
    "raw integer counts"
  )
  expect_error(
    deconv_bayesprism(mae, "rna", reference, labels[-1L]),
    "one value per reference cell"
  )

  named_labels <- stats::setNames(labels, colnames(reference))
  names(named_labels)[1L] <- "unknown_cell"
  expect_error(
    deconv_bayesprism(mae, "rna", reference, named_labels),
    "every reference cell name"
  )

  unrelated <- reference
  rownames(unrelated) <- paste0("other", seq_len(nrow(unrelated)))
  expect_error(
    deconv_bayesprism(mae, "rna", unrelated, labels),
    "no common gene"
  )
  expect_error(
    deconv_bayesprism(mae, "rna", reference, labels, cores = 0L),
    "cores"
  )
  expect_error(
    deconv_bayesprism(mae, "rna", reference, labels, input_type = ""),
    "input_type"
  )
  expect_error(
    deconv_bayesprism(mae, "rna", reference, labels, update_gibbs = NA),
    "update_gibbs"
  )
  expect_error(
    deconv_bayesprism(mae, "rna", reference, labels, gibbs_control = 1),
    "gibbs_control"
  )
})

test_that("BayesPrism runs an offline short-chain integration", {
  suppressWarnings(skip_if_not_installed("BayesPrism"))
  mae <- make_deconvolution_mae(seed = 221L)
  genes <- rownames(mae_pull_assay(mae, "rna", "counts"))
  reference <- make_bayesprism_reference(genes, seed = 222L)
  labels <- rep(c("B_cell", "T_cell"), 4L)

  invisible(utils::capture.output(
    result <- suppressMessages(suppressWarnings(deconv_bayesprism(
      mae,
      "rna",
      reference = reference,
      cell_type_labels = labels,
      cores = 1L,
      update_gibbs = FALSE,
      gibbs_control = list(
        chain.length = 20L,
        burn.in = 10L,
        thinning = 2L,
        seed = 123L
      ),
      opt_control = list(maxit = 20L)
    )))
  ))

  expect_s4_class(result, "BayesPrism")
  expect_false(result@control_param$update.gibbs)
  expect_identical(result@control_param$gibbs.control$chain.length, 20L)
})
