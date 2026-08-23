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

#' Deconvolve bulk counts with MuSiC
#'
#' @inheritParams mae_pull_assay
#' @param sc_reference Single-cell reference accepted by `MuSiC::music_prop()`.
#' @param clusters Cell-type annotation column in the reference.
#' @param samples Sample identifier column in the reference.
#' @param ... Additional arguments passed to `MuSiC::music_prop()`.
#'
#' @return The native MuSiC result list. Its weighted estimates are cell-type
#'   proportions; they are not on the same scale as xCell enrichment scores.
#' @export
deconv_music <- function(
    x,
    experiment,
    sc_reference,
    clusters,
    samples,
    assay = "counts",
    ...
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
  matrix <- .as_count_matrix(x, experiment, assay)
  .assert_deconvolution_ids(matrix, "MuSiC bulk input")
  reference_ids <- rownames(sc_reference)
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
      sc.sce = sc_reference,
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
#' @param reference Raw reference count matrix.
#' @param cell_type_labels Cell-type label for each reference cell.
#' @param cell_state_labels Optional cell-state label for each reference cell.
#' @param reference_orientation Layout of `reference`. The default preserves
#'   the package's documented gene-by-cell input, while `cells_by_genes`
#'   accepts BayesPrism's native layout directly.
#' @param key Optional malignant/tumour key accepted by BayesPrism.
#' @param input_type BayesPrism reference input type.
#' @param cores Number of worker cores passed to `run.prism()`.
#' @param update_gibbs Run the final Gibbs update.
#' @param gibbs_control,opt_control Named control lists passed to
#'   [BayesPrism::run.prism()].
#' @param ... Additional arguments passed to `BayesPrism::new.prism()`.
#'
#' @return The native object returned by `BayesPrism::run.prism()`.
#' @export
deconv_bayesprism <- function(
    x,
    experiment,
    reference,
    cell_type_labels,
    assay = "counts",
    reference_orientation = c("genes_by_cells", "cells_by_genes"),
    cell_state_labels = NULL,
    key = NULL,
    input_type = "count.matrix",
    cores = 1,
    update_gibbs = TRUE,
    gibbs_control = list(),
    opt_control = list(),
    ...
) {
  .require_backend("BayesPrism", "to run BayesPrism deconvolution")
  reference_orientation <- match.arg(reference_orientation)
  mixture <- .as_count_matrix(x, experiment, assay)
  .assert_integerish_counts(mixture, "BayesPrism mixture")
  .assert_deconvolution_ids(mixture, "BayesPrism mixture")
  mixture_sample_names <- colnames(mixture)
  reference <- .as_reference_count_matrix(reference)
  .assert_nonnegative_matrix(reference, "BayesPrism reference")
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

.as_reference_count_matrix <- function(reference) {
  reference <- as.matrix(reference)
  if (!is.numeric(reference)) {
    stop("BayesPrism `reference` must contain numeric raw counts.", call. = FALSE)
  }
  storage.mode(reference) <- "double"
  if (is.null(rownames(reference)) || is.null(colnames(reference))) {
    stop("BayesPrism `reference` must have row and column names.", call. = FALSE)
  }
  if (
    any(!is.finite(reference)) || any(reference < 0) ||
      any(abs(reference - round(reference)) > sqrt(.Machine$double.eps)) ||
      any(reference > .Machine$integer.max)
  ) {
    stop("BayesPrism `reference` must contain raw integer counts.", call. = FALSE)
  }
  reference
}

.assert_integerish_counts <- function(value, label) {
  if (any(abs(value - round(value)) > sqrt(.Machine$double.eps))) {
    stop(label, " must contain raw integer counts.", call. = FALSE)
  }
  invisible(value)
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
