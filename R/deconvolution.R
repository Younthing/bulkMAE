#' Deconvolve bulk expression with immunedeconv
#'
#' The selected assay must contain non-log TPM-like expression with HGNC gene
#' symbols as row names, as required by immunedeconv.
#' xCell returns relative enrichment scores rather than cell fractions. EPIC
#' and quanTIseq produce proportion-like estimates, but their output scales and
#' supported cell types remain method-specific and should not be pooled as if
#' they were interchangeable measurements.
#'
#' @inheritParams mae_pull_assay
#' @param method Method supported by `immunedeconv::deconvolute()`, such as
#'   `epic`, `quantiseq`, or `xcell`.
#' @param ... Additional arguments passed to `immunedeconv::deconvolute()`.
#'
#' @return The native immunedeconv result table.
#' @export
deconv <- function(x, experiment, method, assay = "tpm", ...) {
  .require_backend("immunedeconv", "to run immune deconvolution")
  if (
    !is.character(method) || length(method) != 1L || is.na(method) ||
      !nzchar(method)
  ) {
    stop("`method` must name one immunedeconv method.", call. = FALSE)
  }
  matrix <- .pull_matrix(x, experiment, assay)
  .assert_nonnegative_matrix(matrix, "The selected deconvolution assay")
  .assert_deconvolution_ids(matrix, "Deconvolution input")
  if (tolower(method) == "xcell") {
    message("xCell output is an enrichment score, not an estimated cell fraction.")
  }
  arguments <- c(
    list(gene_expression = matrix, method = method),
    list(...)
  )
  .with_attached_namespaces(
    "immunedeconv",
    .call_backend("immunedeconv", "deconvolute", arguments)
  )
}

#' Construct a single-cell deconvolution reference
#'
#' Builds the canonical `SingleCellExperiment` reference accepted directly by
#' [deconv_music()] and [deconv_bayesprism()]. This helper only standardizes a
#' user-supplied single-cell count matrix and its metadata; it does not provide
#' or simulate a biological reference.
#'
#' @param counts Raw integer count matrix with genes in rows and cells in
#'   columns. Dense matrices and sparse `Matrix` objects are accepted without
#'   densifying the latter. Complete, unique gene and cell names are required.
#' @param cell_data Cell-level metadata with row names containing every count
#'   matrix cell name exactly once. Rows are aligned by cell name.
#' @param cell_type,sample A column name in `cell_data`, or a named vector with
#'   one non-missing label per cell. Values are stored in the canonical
#'   `cell_type` and `sample_id` columns.
#' @param cell_state Optional cell-state column name or named vector. When
#'   supplied, values are stored in the canonical `cell_state` column.
#' @param assay_name Name used for the raw-count assay.
#'
#' @return A `SingleCellExperiment` with aligned cell metadata.
#' @export
deconv_reference <- function(
    counts,
    cell_data,
    cell_type,
    sample,
    cell_state = NULL,
    assay_name = "counts"
) {
  .require_backend(
    "SingleCellExperiment",
    "to construct a single-cell deconvolution reference"
  )
  .assert_scalar_character(assay_name, "assay_name")
  counts <- .as_reference_count_matrix(
    counts,
    "Single-cell `counts`",
    preserve_sparse = TRUE
  )
  .assert_identifier_vector(
    rownames(counts),
    "Single-cell reference gene identifiers"
  )
  cells <- colnames(counts)
  .assert_identifier_vector(cells, "Single-cell reference cell identifiers")
  cell_data <- .align_reference_cell_data(cell_data, cells)
  cell_data$cell_type <- .resolve_builder_annotation(
    cell_type,
    cell_data,
    cells,
    "cell_type"
  )
  cell_data$sample_id <- .resolve_builder_annotation(
    sample,
    cell_data,
    cells,
    "sample"
  )
  if (!is.null(cell_state)) {
    cell_data$cell_state <- .resolve_builder_annotation(
      cell_state,
      cell_data,
      cells,
      "cell_state"
    )
  }

  SingleCellExperiment::SingleCellExperiment(
    assays = stats::setNames(list(counts), assay_name),
    colData = S4Vectors::DataFrame(cell_data, check.names = FALSE)
  )
}

#' Deconvolve bulk counts with MuSiC
#'
#' @inheritParams mae_pull_assay
#' @param sc_reference Single-cell reference accepted by `MuSiC::music_prop()`.
#' @param clusters Cell-type annotation column in the reference. The default
#'   matches the canonical column created by [deconv_reference()].
#' @param samples Sample identifier column in the reference. The default
#'   matches the canonical column created by [deconv_reference()].
#' @param reference_assay Raw-count assay in `sc_reference`. `NULL` selects
#'   `counts`, or the sole assay when no `counts` assay is present.
#' @param ... Additional arguments passed to `MuSiC::music_prop()`.
#'
#' @return The native MuSiC result list. Its weighted estimates are cell-type
#'   proportions; they are not on the same scale as xCell enrichment scores.
#' @export
deconv_music <- function(
    x,
    experiment,
    sc_reference,
    clusters = "cell_type",
    samples = "sample_id",
    assay = "counts",
    ...,
    reference_assay = NULL
) {
  .require_backend("MuSiC", "to run MuSiC deconvolution")
  if (
    !is.character(clusters) || length(clusters) != 1L || is.na(clusters) ||
      !nzchar(clusters) || !is.character(samples) || length(samples) != 1L ||
      is.na(samples) || !nzchar(samples)
  ) {
    stop("`clusters` and `samples` must each name one reference column.", call. = FALSE)
  }
  if (!methods::is(sc_reference, "SingleCellExperiment")) {
    stop("`sc_reference` must be a SingleCellExperiment.", call. = FALSE)
  }
  reference_assay <- .resolve_reference_assay(
    sc_reference,
    reference_assay,
    "sc_reference"
  )
  reference_counts <- .as_reference_count_matrix(
    SummarizedExperiment::assay(sc_reference, reference_assay),
    "MuSiC reference assay",
    preserve_sparse = TRUE
  )
  .assert_deconvolution_ids(reference_counts, "MuSiC reference")
  music_reference <- sc_reference
  if (!identical(reference_assay, "counts")) {
    SummarizedExperiment::assay(
      music_reference,
      "counts",
      withDimnames = FALSE
    ) <- reference_counts
  }
  reference_data <- as.data.frame(
    SummarizedExperiment::colData(sc_reference),
    optional = TRUE
  )
  missing_columns <- setdiff(c(clusters, samples), names(reference_data))
  if (length(missing_columns)) {
    stop(
      "Reference metadata is missing: ", paste(missing_columns, collapse = ", "),
      call. = FALSE
    )
  }
  if (
    anyNA(reference_data[[clusters]]) || anyNA(reference_data[[samples]]) ||
      any(!nzchar(as.character(reference_data[[clusters]]))) ||
      any(!nzchar(as.character(reference_data[[samples]])))
  ) {
    stop("Reference cluster and sample labels cannot be missing or empty.", call. = FALSE)
  }
  if (length(unique(reference_data[[clusters]])) < 2L) {
    stop("MuSiC requires at least two reference cell types.", call. = FALSE)
  }
  sample_values <- as.character(reference_data[[samples]])
  cluster_values <- as.character(reference_data[[clusters]])
  if (length(unique(sample_values)) < 2L) {
    stop("MuSiC requires at least two reference biological samples.", call. = FALSE)
  }
  samples_per_cluster <- vapply(
    split(sample_values, cluster_values),
    function(value) length(unique(value)),
    integer(1)
  )
  if (any(samples_per_cluster < 2L)) {
    stop(
      "Every MuSiC reference cell type must occur in at least two biological ",
      "samples; insufficient: ",
      paste(names(samples_per_cluster)[samples_per_cluster < 2L], collapse = ", "),
      ".",
      call. = FALSE
    )
  }
  matrix <- .as_count_matrix(x, experiment, assay)
  .assert_deconvolution_ids(matrix, "MuSiC bulk input")
  reference_ids <- rownames(reference_counts)
  .assert_identifier_vector(reference_ids, "MuSiC reference feature names")
  overlap <- intersect(rownames(matrix), reference_ids)
  if (!length(overlap)) {
    stop(
      "MuSiC bulk and reference inputs have no common feature identifiers.",
      call. = FALSE
    )
  }
  if (length(overlap) < 100L) {
    warning(
      "MuSiC inputs share only ", length(overlap),
      " features; check identifier types and reference coverage.",
      call. = FALSE
    )
  }
  arguments <- c(
    list(
      bulk.mtx = matrix,
      sc.sce = music_reference,
      clusters = clusters,
      samples = samples
    ),
    list(...)
  )
  .with_attached_namespaces(
    c("SummarizedExperiment", "SingleCellExperiment"),
    .call_backend("MuSiC", "music_prop", arguments)
  )
}

#' Prepare a CIBERSORTx mixture table
#'
#' CIBERSORTx is an external service or container rather than an R backend.
#' This helper prepares its required gene-by-sample table without uploading or
#' writing data.
#'
#' @inheritParams mae_pull_assay
#' @param gene_column Name of the first identifier column.
#'
#' @return A data frame ready to write as a tab-delimited mixture file.
#' @export
deconv_cibersortx_input <- function(
    x,
    experiment,
    assay = "tpm",
    gene_column = "GeneSymbol"
) {
  matrix <- .pull_matrix(x, experiment, assay)
  .assert_nonnegative_matrix(matrix, "CIBERSORTx input")
  .assert_deconvolution_ids(matrix, "CIBERSORTx input")
  if (
    !is.character(gene_column) || length(gene_column) != 1L ||
      is.na(gene_column) || !nzchar(gene_column)
  ) {
    stop("`gene_column` must be one non-empty column name.", call. = FALSE)
  }
  result <- data.frame(rownames(matrix), matrix, check.names = FALSE)
  names(result)[1L] <- gene_column
  rownames(result) <- NULL
  result
}

#' Deconvolve bulk counts with BayesPrism
#'
#' Both bulk and single-cell reference inputs must be raw integer counts. The
#' MAE assay is unambiguously gene-by-sample and is converted to BayesPrism's
#' sample-by-gene layout. `reference_orientation` explicitly declares the
#' supplied reference layout; no orientation is inferred from dimensions.
#'
#' @inheritParams mae_pull_assay
#' @param reference Raw reference count matrix, or a `SingleCellExperiment`
#'   created by [deconv_reference()].
#' @param cell_type_labels Cell-type label for each reference cell. For a
#'   `SingleCellExperiment`, `NULL` uses its `cell_type` column and a single
#'   string names another `colData` column.
#' @param cell_state_labels Optional cell-state label for each reference cell.
#'   For a `SingleCellExperiment`, `NULL` uses `cell_state` when that column is
#'   present and a single string names another `colData` column.
#' @param reference_orientation Layout of `reference`. The default preserves
#'   the package's documented gene-by-cell input, while `cells_by_genes`
#'   accepts BayesPrism's native layout directly.
#' @param key Optional malignant/tumour key accepted by BayesPrism.
#' @param input_type BayesPrism reference input type.
#' @param cores Number of worker cores passed to `run.prism()`.
#' @param update_gibbs Run the final Gibbs update.
#' @param gibbs_control,opt_control Named control lists passed to
#'   [BayesPrism::run.prism()].
#' @param reference_assay Raw-count assay used when `reference` is a
#'   `SingleCellExperiment`. `NULL` selects `counts`, or the sole assay when no
#'   `counts` assay is present.
#' @param ... Additional arguments passed to `BayesPrism::new.prism()`.
#'
#' @return The native object returned by `BayesPrism::run.prism()`.
#' @export
deconv_bayesprism <- function(
    x,
    experiment,
    reference,
    cell_type_labels = NULL,
    assay = "counts",
    reference_orientation = c("genes_by_cells", "cells_by_genes"),
    cell_state_labels = NULL,
    key = NULL,
    input_type = "count.matrix",
    cores = 1,
    update_gibbs = TRUE,
    gibbs_control = list(),
    opt_control = list(),
    ...,
    reference_assay = NULL
) {
  .require_backend("BayesPrism", "to run BayesPrism deconvolution")
  reference_orientation <- match.arg(reference_orientation)
  mixture <- .as_count_matrix(x, experiment, assay)
  .assert_integerish_counts(mixture, "BayesPrism mixture")
  .assert_deconvolution_ids(mixture, "BayesPrism mixture")
  mixture_sample_names <- colnames(mixture)

  reference_is_sce <- methods::is(reference, "SingleCellExperiment")
  if (reference_is_sce) {
    if (!identical(reference_orientation, "genes_by_cells")) {
      stop(
        "A `SingleCellExperiment` reference always uses genes-by-cells ",
        "orientation; do not set `reference_orientation = \"cells_by_genes\"`.",
        call. = FALSE
      )
    }
    reference_assay <- .resolve_reference_assay(
      reference,
      reference_assay,
      "reference"
    )
    reference_data <- as.data.frame(
      SummarizedExperiment::colData(reference),
      optional = TRUE
    )
    reference <- .as_reference_count_matrix(
      SummarizedExperiment::assay(reference, reference_assay)
    )
    .assert_deconvolution_ids(reference, "BayesPrism reference")
    cell_names <- colnames(reference)
    cell_type_labels <- .resolve_sce_reference_labels(
      cell_type_labels,
      reference_data,
      cell_names,
      "cell_type_labels",
      default_column = "cell_type",
      optional = FALSE
    )
    cell_state_labels <- .resolve_sce_reference_labels(
      cell_state_labels,
      reference_data,
      cell_names,
      "cell_state_labels",
      default_column = "cell_state",
      optional = TRUE
    )
    reference <- t(reference)
  } else {
    reference <- .as_reference_count_matrix(reference)
    if (reference_orientation == "genes_by_cells") {
      .assert_deconvolution_ids(reference, "BayesPrism reference")
      cell_names <- colnames(reference)
      reference <- t(reference)
    } else {
      .assert_identifier_vector(colnames(reference), "BayesPrism reference genes")
      if (anyDuplicated(rownames(reference))) {
        stop("BayesPrism reference cell names must be unique.", call. = FALSE)
      }
      cell_names <- rownames(reference)
    }
    .assert_identifier_vector(cell_names, "BayesPrism reference cell names")
    cell_type_labels <- .align_reference_labels(
      cell_type_labels,
      cell_names,
      "cell_type_labels"
    )
    if (!is.null(cell_state_labels)) {
      cell_state_labels <- .align_reference_labels(
        cell_state_labels,
        cell_names,
        "cell_state_labels"
      )
    }
  }
  .assert_nonnegative_matrix(reference, "BayesPrism reference")
  .assert_identifier_vector(cell_names, "BayesPrism reference cell names")

  mixture <- t(mixture)
  common_genes <- intersect(colnames(mixture), colnames(reference))
  if (!length(common_genes)) {
    stop(
      "BayesPrism mixture and reference have no common gene identifiers.",
      call. = FALSE
    )
  }
  if (length(common_genes) < 100L) {
    warning(
      "BayesPrism inputs share only ", length(common_genes),
      " genes; check identifier types and reference coverage.",
      call. = FALSE
    )
  }
  mixture <- mixture[, common_genes, drop = FALSE]
  reference <- reference[, common_genes, drop = FALSE]
  if (
    !identical(rownames(mixture), mixture_sample_names) ||
      !identical(rownames(reference), cell_names) ||
      !identical(colnames(mixture), colnames(reference))
  ) {
    stop(
      "Internal BayesPrism alignment failed: expected sample/cell-by-gene inputs.",
      call. = FALSE
    )
  }
  if (
    !is.numeric(cores) || length(cores) != 1L || !is.finite(cores) ||
      cores < 1L || cores != as.integer(cores)
  ) {
    stop("`cores` must be a positive integer.", call. = FALSE)
  }
  if (
    !is.character(input_type) || length(input_type) != 1L ||
      is.na(input_type) || !nzchar(input_type)
  ) {
    stop("`input_type` must be one non-empty string.", call. = FALSE)
  }
  if (
    !is.logical(update_gibbs) || length(update_gibbs) != 1L ||
      is.na(update_gibbs)
  ) {
    stop("`update_gibbs` must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.list(gibbs_control) || is.data.frame(gibbs_control)) {
    stop("`gibbs_control` must be a list.", call. = FALSE)
  }
  if (!is.list(opt_control) || is.data.frame(opt_control)) {
    stop("`opt_control` must be a list.", call. = FALSE)
  }
  arguments <- c(
    list(
      reference = reference,
      mixture = mixture,
      input.type = input_type,
      cell.type.labels = cell_type_labels,
      cell.state.labels = cell_state_labels,
      key = key
    ),
    list(...)
  )
  prism <- .call_backend("BayesPrism", "new.prism", arguments)
  .call_backend(
    "BayesPrism",
    "run.prism",
    list(
      prism = prism,
      n.cores = cores,
      update.gibbs = update_gibbs,
      gibbs.control = gibbs_control,
      opt.control = opt_control
    )
  )
}

.call_backend <- function(package, function_name, arguments) {
  function_object <- getExportedValue(package, function_name)
  do.call(function_object, arguments)
}

.assert_nonnegative_matrix <- function(x, label) {
  if (!is.matrix(x) || !is.numeric(x) || any(!is.finite(x)) || any(x < 0)) {
    stop(label, " must contain finite, non-negative values.", call. = FALSE)
  }
  invisible(x)
}

.assert_identifier_vector <- function(value, label) {
  if (
    is.null(value) || !length(value) || anyNA(value) ||
      any(!nzchar(as.character(value))) || anyDuplicated(value)
  ) {
    stop(label, " must be present, non-empty, and unique.", call. = FALSE)
  }
  invisible(value)
}

.assert_deconvolution_ids <- function(matrix, label) {
  .assert_identifier_vector(rownames(matrix), paste(label, "gene identifiers"))
  .assert_identifier_vector(colnames(matrix), paste(label, "sample identifiers"))
  invisible(matrix)
}

.as_reference_count_matrix <- function(
    reference,
    label = "BayesPrism `reference`",
    preserve_sparse = FALSE
) {
  is_sparse_matrix <- methods::is(reference, "sparseMatrix")
  numeric_sparse_matrix <- is_sparse_matrix && any(vapply(
    c("dMatrix", "iMatrix", "lMatrix", "nMatrix"),
    function(class_name) methods::is(reference, class_name),
    logical(1)
  ))
  is_numeric_matrix <- is.matrix(reference) && is.numeric(reference)
  if (is_sparse_matrix && !numeric_sparse_matrix) {
    stop(label, " must contain numeric raw counts.", call. = FALSE)
  }
  if (!is_sparse_matrix && !is_numeric_matrix) {
    reference <- as.matrix(reference)
    is_numeric_matrix <- is.numeric(reference)
  }
  if (!is_sparse_matrix && !is_numeric_matrix) {
    stop(label, " must contain numeric raw counts.", call. = FALSE)
  }
  if (is_sparse_matrix && !preserve_sparse) {
    reference <- as.matrix(reference)
    is_sparse_matrix <- FALSE
  }
  if (!is_sparse_matrix) storage.mode(reference) <- "double"
  if (is.null(rownames(reference)) || is.null(colnames(reference))) {
    stop(label, " must have row and column names.", call. = FALSE)
  }
  values <- if (is_sparse_matrix) {
    if ("x" %in% methods::slotNames(reference)) {
      methods::slot(reference, "x")
    } else {
      # Pattern sparse matrices store only positions; every stored value is 1.
      1
    }
  } else {
    reference
  }
  if (
    any(!is.finite(values)) || any(values < 0) ||
      any(abs(values - round(values)) > sqrt(.Machine$double.eps)) ||
      any(values > .Machine$integer.max)
  ) {
    stop(label, " must contain raw integer counts.", call. = FALSE)
  }
  reference
}

.resolve_reference_assay <- function(reference, requested, argument) {
  available <- SummarizedExperiment::assayNames(reference)
  if (!length(available)) {
    stop("`", argument, "` contains no assays.", call. = FALSE)
  }
  if (is.null(requested)) {
    if ("counts" %in% available) return("counts")
    if (length(available) == 1L) return(available[[1L]])
    stop(
      "`reference_assay` is required when `", argument,
      "` has multiple assays and none is named `counts`.",
      call. = FALSE
    )
  }
  .assert_scalar_character(requested, "reference_assay")
  if (!requested %in% available) {
    stop(
      "Unknown reference assay `", requested, "`. Available assays: ",
      paste(available, collapse = ", "),
      call. = FALSE
    )
  }
  requested
}

.assert_integerish_counts <- function(value, label) {
  if (any(abs(value - round(value)) > sqrt(.Machine$double.eps))) {
    stop(label, " must contain raw integer counts.", call. = FALSE)
  }
  invisible(value)
}

.align_reference_cell_data <- function(cell_data, cells) {
  if (!is.data.frame(cell_data) && !methods::is(cell_data, "DataFrame")) {
    stop("`cell_data` must be a data frame with cell names as row names.", call. = FALSE)
  }
  cell_data <- as.data.frame(cell_data, optional = TRUE)
  cell_names <- rownames(cell_data)
  if (
    is.null(cell_names) || length(cell_names) != nrow(cell_data) ||
      anyNA(cell_names) || any(!nzchar(cell_names)) || anyDuplicated(cell_names) ||
      !setequal(cell_names, cells)
  ) {
    stop(
      "`cell_data` row names must contain every count-matrix cell name exactly once.",
      call. = FALSE
    )
  }
  cell_data[cells, , drop = FALSE]
}

.resolve_builder_annotation <- function(value, cell_data, cells, argument) {
  if (
    is.character(value) && length(value) == 1L && !is.na(value) &&
      nzchar(value) && value %in% names(cell_data)
  ) {
    labels <- cell_data[[value]]
  } else {
    value_names <- names(value)
    if (
      (!is.atomic(value) && !is.factor(value)) || is.null(value_names) ||
        length(value_names) != length(value) || anyNA(value_names) ||
        any(!nzchar(value_names))
    ) {
      stop(
        "`", argument,
        "` must name one `cell_data` column or be a named vector.",
        call. = FALSE
      )
    }
    labels <- .align_reference_labels(value, cells, argument)
  }
  .validate_reference_annotation(labels, cells, argument)
}

.resolve_sce_reference_labels <- function(
    labels,
    reference_data,
    cells,
    argument,
    default_column,
    optional
) {
  if (is.null(labels)) {
    if (default_column %in% names(reference_data)) {
      labels <- reference_data[[default_column]]
    } else if (optional) {
      return(NULL)
    } else {
      stop(
        "SingleCellExperiment reference metadata is missing the default `",
        default_column, "` column; supply `", argument, "` explicitly.",
        call. = FALSE
      )
    }
  } else if (
    is.character(labels) && length(labels) == 1L && !is.na(labels) &&
      nzchar(labels)
  ) {
    if (!labels %in% names(reference_data)) {
      stop(
        "`", argument, "` names an unknown reference metadata column: ",
        labels, ".",
        call. = FALSE
      )
    }
    labels <- reference_data[[labels]]
  }
  .align_reference_labels(labels, cells, argument)
}

.validate_reference_annotation <- function(labels, cells, argument) {
  if (
    (!is.atomic(labels) && !is.factor(labels)) || length(labels) != length(cells) ||
      anyNA(labels) || any(!nzchar(as.character(labels)))
  ) {
    stop(
      "`", argument, "` must provide one non-missing label per reference cell.",
      call. = FALSE
    )
  }
  unname(labels)
}

.align_reference_labels <- function(labels, cells, argument) {
  if (length(labels) != length(cells)) {
    stop("`", argument, "` must have one value per reference cell.", call. = FALSE)
  }
  label_names <- names(labels)
  if (!is.null(label_names) && any(nzchar(label_names))) {
    if (
      any(!nzchar(label_names)) || anyDuplicated(label_names) ||
        !setequal(label_names, cells)
    ) {
      stop(
        "Named `", argument,
        "` must contain every reference cell name exactly once.",
        call. = FALSE
      )
    }
    labels <- labels[cells]
  }
  if (anyNA(labels) || any(!nzchar(as.character(labels)))) {
    stop("`", argument, "` must not contain missing or empty labels.", call. = FALSE)
  }
  unname(labels)
}
