.reference_bridge_api <- function(name) {
  get(name, envir = asNamespace("bulkMAE"), inherits = FALSE)
}

.make_reference_bridge_fixture <- function(
    n_genes = 120L,
    n_bulk_samples = 4L,
    n_donors = 4L,
    cells_per_type = 2L,
    seed = 301L
) {
  set.seed(seed)
  genes <- sprintf("gene%03d", seq_len(n_genes))
  bulk_samples <- sprintf("bulk%02d", seq_len(n_bulk_samples))
  cell_types <- c("B_cell", "T_cell")
  n_cells <- n_donors * length(cell_types) * cells_per_type
  cells <- sprintf("cell%03d", seq_len(n_cells))
  donor <- rep(
    sprintf("donor%02d", seq_len(n_donors)),
    each = length(cell_types) * cells_per_type
  )
  annotation <- rep(rep(cell_types, each = cells_per_type), n_donors)
  state <- ifelse(annotation == "B_cell", "resting_B", "activated_T")

  reference_mean <- matrix(20, nrow = n_genes, ncol = n_cells)
  reference_mean[seq_len(20L), annotation == "B_cell"] <- 55
  reference_mean[21:40, annotation == "T_cell"] <- 55
  reference_counts <- matrix(
    stats::rpois(length(reference_mean), lambda = as.vector(reference_mean)),
    nrow = n_genes,
    dimnames = list(genes, cells)
  )
  storage.mode(reference_counts) <- "integer"

  cell_data <- data.frame(
    annotation = annotation,
    donor = donor,
    state = state,
    library_batch = rep(c("batch_a", "batch_b"), length.out = n_cells),
    quality = seq_len(n_cells) / n_cells,
    row.names = cells,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  bulk_mean <- matrix(80, nrow = n_genes, ncol = n_bulk_samples)
  bulk_mean[seq_len(20L), seq_len(n_bulk_samples) %% 2L == 1L] <- 140
  bulk_mean[21:40, seq_len(n_bulk_samples) %% 2L == 0L] <- 140
  bulk_counts <- matrix(
    stats::rpois(length(bulk_mean), lambda = as.vector(bulk_mean)),
    nrow = n_genes,
    dimnames = list(genes, bulk_samples)
  )
  storage.mode(bulk_counts) <- "integer"
  bulk_data <- data.frame(
    condition = rep(c("control", "treated"), length.out = n_bulk_samples),
    row.names = bulk_samples,
    stringsAsFactors = FALSE
  )
  experiment <- mae_create_experiment(
    assays = list(counts = bulk_counts),
    col_data = bulk_data
  )

  list(
    mae = mae_create(list(rna = experiment), bulk_data),
    reference_counts = reference_counts,
    cell_data = cell_data
  )
}

test_that("deconv_reference aligns canonical and additional cell metadata", {
  suppressWarnings(skip_if_not_installed("SingleCellExperiment"))
  builder <- .reference_bridge_api("deconv_reference")
  fixture <- .make_reference_bridge_fixture(seed = 302L)
  cells <- colnames(fixture$reference_counts)
  reversed <- rev(cells)
  sample_labels <- stats::setNames(
    fixture$cell_data$donor,
    rownames(fixture$cell_data)
  )[reversed]

  reference <- suppressWarnings(builder(
    counts = fixture$reference_counts,
    cell_data = fixture$cell_data[reversed, , drop = FALSE],
    cell_type = "annotation",
    sample = sample_labels,
    cell_state = "state",
    assay_name = "raw_counts"
  ))
  metadata <- as.data.frame(
    SummarizedExperiment::colData(reference),
    optional = TRUE
  )

  expect_s4_class(reference, "SingleCellExperiment")
  expect_identical(SummarizedExperiment::assayNames(reference), "raw_counts")
  expect_equal(
    SummarizedExperiment::assay(reference, "raw_counts"),
    fixture$reference_counts
  )
  expect_identical(colnames(reference), cells)
  expect_identical(rownames(metadata), cells)
  expect_identical(metadata$cell_type, fixture$cell_data[cells, "annotation"])
  expect_identical(metadata$sample_id, fixture$cell_data[cells, "donor"])
  expect_identical(metadata$cell_state, fixture$cell_data[cells, "state"])
  expect_identical(metadata$library_batch, fixture$cell_data[cells, "library_batch"])
  expect_identical(metadata$quality, fixture$cell_data[cells, "quality"])
})

test_that("deconv_reference validates raw counts, cell alignment, and labels", {
  suppressWarnings(skip_if_not_installed("SingleCellExperiment"))
  builder <- .reference_bridge_api("deconv_reference")
  fixture <- .make_reference_bridge_fixture(seed = 303L)
  counts <- fixture$reference_counts
  cell_data <- fixture$cell_data
  cells <- colnames(counts)

  non_integer <- counts
  non_integer[1L, 1L] <- non_integer[1L, 1L] + 0.5
  expect_error(
    builder(non_integer, cell_data, "annotation", "donor"),
    "raw integer counts"
  )

  unnamed <- counts
  colnames(unnamed) <- NULL
  expect_error(
    builder(unnamed, cell_data, "annotation", "donor"),
    "row and column names"
  )

  duplicated <- counts
  colnames(duplicated)[2L] <- colnames(duplicated)[1L]
  expect_error(
    builder(duplicated, cell_data, "annotation", "donor"),
    "cell identifiers.*unique"
  )

  missing_cell <- cell_data[-1L, , drop = FALSE]
  expect_error(
    builder(counts, missing_cell, "annotation", "donor"),
    "every count-matrix cell name"
  )
  expect_error(
    builder(counts, unname(as.list(cell_data)), "annotation", "donor"),
    "cell_data.*data frame"
  )
  expect_error(
    builder(counts, cell_data, "unknown_column", "donor"),
    "cell_type.*named vector"
  )
  expect_error(
    builder(counts, cell_data, cell_data$annotation, "donor"),
    "cell_type.*named vector"
  )
  empty_names <- stats::setNames(cell_data$annotation, rep("", length(cells)))
  expect_error(
    builder(counts, cell_data, empty_names, "donor"),
    "cell_type.*named vector"
  )

  bad_samples <- stats::setNames(cell_data$donor[-1L], cells[-1L])
  expect_error(
    builder(counts, cell_data, "annotation", bad_samples),
    "one value per reference cell"
  )

  invalid_labels <- cell_data
  invalid_labels$donor[1L] <- NA_character_
  expect_error(
    builder(counts, invalid_labels, "annotation", "donor"),
    "sample.*non-missing label"
  )
  expect_error(
    builder(counts, cell_data, "annotation", "donor", cell_state = "missing"),
    "cell_state.*named vector"
  )
  expect_error(
    builder(counts, cell_data, "annotation", "donor", assay_name = ""),
    "assay_name"
  )
})

test_that("deconv_music uses canonical builder columns by default", {
  suppressWarnings(skip_if_not_installed("SingleCellExperiment"))
  suppressWarnings(skip_if_not_installed("MuSiC"))
  builder <- .reference_bridge_api("deconv_reference")
  fixture <- .make_reference_bridge_fixture(seed = 304L)
  reference <- builder(
    fixture$reference_counts,
    fixture$cell_data,
    cell_type = "annotation",
    sample = "donor"
  )
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
    fixture$mae,
    "rna",
    sc_reference = reference,
    assay = "counts"
  )

  expect_identical(formals(deconv_music)$clusters, "cell_type")
  expect_identical(formals(deconv_music)$samples, "sample_id")
  expect_identical(observed$package, "MuSiC")
  expect_identical(observed$function_name, "music_prop")
  expect_identical(observed$arguments$clusters, "cell_type")
  expect_identical(observed$arguments$samples, "sample_id")
  expect_identical(observed$arguments$sc.sce, reference)
  expect_identical(dim(result$Est.prop.weighted), c(4L, 2L))
})

test_that("deconv_music validates biological replication in the reference", {
  suppressWarnings(skip_if_not_installed("SingleCellExperiment"))
  suppressWarnings(skip_if_not_installed("MuSiC"))
  builder <- .reference_bridge_api("deconv_reference")
  fixture <- .make_reference_bridge_fixture(seed = 3041L)
  one_donor <- fixture$cell_data
  one_donor$donor <- "donor01"
  reference <- builder(
    fixture$reference_counts,
    one_donor,
    cell_type = "annotation",
    sample = "donor"
  )
  expect_error(
    deconv_music(fixture$mae, "rna", reference, assay = "counts"),
    "at least two reference biological samples"
  )

  one_type_one_donor <- fixture$cell_data
  one_type_one_donor$donor[
    one_type_one_donor$annotation == "B_cell"
  ] <- "donor01"
  reference <- builder(
    fixture$reference_counts,
    one_type_one_donor,
    cell_type = "annotation",
    sample = "donor"
  )
  expect_error(
    deconv_music(fixture$mae, "rna", reference, assay = "counts"),
    "Every MuSiC reference cell type.*B_cell"
  )
})

test_that("deconv_reference reaches the real MuSiC backend", {
  suppressWarnings(skip_if_not_installed("SingleCellExperiment"))
  suppressWarnings(skip_if_not_installed("MuSiC"))
  builder <- .reference_bridge_api("deconv_reference")
  fixture <- .make_reference_bridge_fixture(seed = 305L)
  reference <- builder(
    fixture$reference_counts,
    fixture$cell_data,
    cell_type = "annotation",
    sample = "donor"
  )

  invisible(utils::capture.output(
    result <- suppressMessages(suppressWarnings(deconv_music(
      fixture$mae,
      "rna",
      sc_reference = reference,
      assay = "counts"
    )))
  ))

  expect_identical(dim(result$Est.prop.weighted), c(4L, 2L))
  expect_true(all(is.finite(result$Est.prop.weighted)))
})

test_that("custom and sparse builder assays reach MuSiC without densification", {
  suppressWarnings(skip_if_not_installed("SingleCellExperiment"))
  suppressWarnings(skip_if_not_installed("MuSiC"))
  skip_if_not_installed("Matrix")
  builder <- .reference_bridge_api("deconv_reference")
  fixture <- .make_reference_bridge_fixture(seed = 3051L)
  sparse_counts <- methods::as(
    Matrix::Matrix(fixture$reference_counts, sparse = TRUE),
    "dgCMatrix"
  )
  reference <- builder(
    sparse_counts,
    fixture$cell_data,
    cell_type = "annotation",
    sample = "donor",
    assay_name = "raw_counts"
  )

  expect_true(methods::is(
    SummarizedExperiment::assay(reference, "raw_counts"),
    "sparseMatrix"
  ))
  invisible(utils::capture.output(
    result <- suppressMessages(suppressWarnings(deconv_music(
      fixture$mae,
      "rna",
      sc_reference = reference,
      assay = "counts"
    )))
  ))
  expect_identical(dim(result$Est.prop.weighted), c(4L, 2L))
  expect_true(methods::is(
    SummarizedExperiment::assay(reference, "raw_counts"),
    "sparseMatrix"
  ))
  expect_identical(SummarizedExperiment::assayNames(reference), "raw_counts")
})

test_that("BayesPrism resolves SCE assays and canonical or named metadata columns", {
  suppressWarnings(skip_if_not_installed("SingleCellExperiment"))
  suppressWarnings(skip_if_not_installed("BayesPrism"))
  builder <- .reference_bridge_api("deconv_reference")
  fixture <- .make_reference_bridge_fixture(seed = 306L)
  reference <- builder(
    fixture$reference_counts,
    fixture$cell_data,
    cell_type = "annotation",
    sample = "donor",
    cell_state = "state",
    assay_name = "raw_counts"
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
    fixture$mae,
    "rna",
    reference = reference,
    update_gibbs = FALSE
  )

  constructor <- observed$calls[[1L]]
  expect_s3_class(result, "bulkMAETestBayesResult")
  expect_identical(constructor$function_name, "new.prism")
  expect_identical(dim(constructor$arguments$mixture), c(4L, 120L))
  expect_identical(
    dim(constructor$arguments$reference),
    c(ncol(fixture$reference_counts), 120L)
  )
  expect_identical(
    constructor$arguments$cell.type.labels,
    fixture$cell_data[colnames(fixture$reference_counts), "annotation"]
  )
  expect_identical(
    constructor$arguments$cell.state.labels,
    fixture$cell_data[colnames(fixture$reference_counts), "state"]
  )

  observed$calls <- list()
  deconv_bayesprism(
    fixture$mae,
    "rna",
    reference = reference,
    cell_type_labels = "annotation",
    cell_state_labels = "state",
    reference_assay = "raw_counts",
    update_gibbs = FALSE
  )
  constructor <- observed$calls[[1L]]
  expect_identical(
    constructor$arguments$cell.type.labels,
    fixture$cell_data[colnames(fixture$reference_counts), "annotation"]
  )
  expect_identical(
    constructor$arguments$cell.state.labels,
    fixture$cell_data[colnames(fixture$reference_counts), "state"]
  )
})

test_that("BayesPrism validates SCE assay and metadata contracts", {
  suppressWarnings(skip_if_not_installed("SingleCellExperiment"))
  suppressWarnings(skip_if_not_installed("BayesPrism"))
  builder <- .reference_bridge_api("deconv_reference")
  fixture <- .make_reference_bridge_fixture(seed = 307L)
  reference <- builder(
    fixture$reference_counts,
    fixture$cell_data,
    cell_type = "annotation",
    sample = "donor"
  )

  expect_error(
    deconv_bayesprism(
      fixture$mae,
      "rna",
      reference,
      reference_assay = "missing"
    ),
    "Unknown reference assay"
  )

  missing_type <- reference
  SummarizedExperiment::colData(missing_type)$cell_type <- NULL
  expect_error(
    deconv_bayesprism(fixture$mae, "rna", missing_type),
    "missing the default `cell_type`"
  )
  expect_error(
    deconv_bayesprism(
      fixture$mae,
      "rna",
      reference,
      cell_type_labels = "missing"
    ),
    "unknown reference metadata column"
  )
  expect_error(
    deconv_bayesprism(
      fixture$mae,
      "rna",
      reference,
      cell_state_labels = "missing"
    ),
    "unknown reference metadata column"
  )
  expect_error(
    deconv_bayesprism(
      fixture$mae,
      "rna",
      reference,
      reference_orientation = "cells_by_genes"
    ),
    "always uses genes-by-cells"
  )

  non_integer <- reference
  SummarizedExperiment::assay(non_integer, "counts")[1L, 1L] <-
    SummarizedExperiment::assay(non_integer, "counts")[1L, 1L] + 0.5
  expect_error(
    deconv_bayesprism(fixture$mae, "rna", non_integer),
    "raw integer counts"
  )
})

test_that("builder SCE reaches the real BayesPrism short chain", {
  suppressWarnings(skip_if_not_installed("SingleCellExperiment"))
  suppressWarnings(skip_if_not_installed("BayesPrism"))
  builder <- .reference_bridge_api("deconv_reference")
  fixture <- .make_reference_bridge_fixture(
    cells_per_type = 1L,
    seed = 308L
  )
  reference <- builder(
    fixture$reference_counts,
    fixture$cell_data,
    cell_type = "annotation",
    sample = "donor"
  )

  invisible(utils::capture.output(
    result <- suppressMessages(suppressWarnings(deconv_bayesprism(
      fixture$mae,
      "rna",
      reference = reference,
      cores = 1L,
      update_gibbs = FALSE,
      gibbs_control = list(
        chain.length = 20L,
        burn.in = 10L,
        thinning = 2L,
        seed = 321L
      ),
      opt_control = list(maxit = 20L)
    )))
  ))

  expect_s4_class(result, "BayesPrism")
  expect_false(result@control_param$update.gibbs)
  expect_identical(result@control_param$gibbs.control$chain.length, 20L)
})
