#' Construct a SummarizedExperiment leaf
#'
#' @param assays A matrix, matrix-like object, or named list of assays.
#' @param col_data Sample metadata coercible to `S4Vectors::DataFrame`.
#' @param row_data Optional feature metadata coercible to `S4Vectors::DataFrame`.
#' @param assay_name Name used when `assays` is a single matrix.
#'
#' @return A [SummarizedExperiment::SummarizedExperiment].
#'
#' @examples
#' counts <- matrix(
#'   1:6,
#'   nrow = 3,
#'   dimnames = list(paste0("gene", 1:3), c("sample1", "sample2"))
#' )
#' samples <- data.frame(
#'   condition = c("control", "treated"),
#'   row.names = colnames(counts)
#' )
#' mae_create_experiment(counts, samples)
#' @export
mae_create_experiment <- function(
    assays,
    col_data,
    row_data = NULL,
    assay_name = "counts") {
  if (methods::is(assays, "List") && !methods::is(assays, "DataFrame")) {
    assays <- as.list(assays)
  }
  if (
    !is.list(assays) || is.data.frame(assays) ||
      methods::is(assays, "DataFrame")
  ) {
    .assert_scalar_character(assay_name, "assay_name")
    assays <- stats::setNames(list(assays), assay_name)
  }
  if (!length(assays)) {
    stop("`assays` must contain at least one assay.", call. = FALSE)
  }
  if (
    is.null(names(assays)) || anyNA(names(assays)) ||
      any(!nzchar(names(assays))) || anyDuplicated(names(assays))
  ) {
    stop("`assays` must have unique, non-empty names.", call. = FALSE)
  }

  dimensions <- lapply(assays, dim)
  not_matrix <- vapply(
    dimensions,
    function(value) length(value) != 2L,
    logical(1)
  )
  if (any(not_matrix)) {
    stop("Every assay must be a two-dimensional matrix-like object.", call. = FALSE)
  }
  if (any(vapply(dimensions, function(value) any(value == 0L), logical(1)))) {
    stop("Assays must contain features and samples.", call. = FALSE)
  }

  first <- assays[[1L]]
  .assert_ids(rownames(first), "Assay feature names")
  .assert_ids(colnames(first), "Assay sample names")
  if (length(assays) > 1L) {
    for (index in seq.int(2L, length(assays))) {
      current <- assays[[index]]
      .assert_ids(rownames(current), "Assay feature names")
      .assert_ids(colnames(current), "Assay sample names")
      if (
        !identical(dim(current), dim(first)) ||
          !identical(rownames(current), rownames(first)) ||
          !identical(colnames(current), colnames(first))
      ) {
        stop(
          paste(
            "All assays in one experiment must have identical dimensions",
            "and dimnames."
          ),
          call. = FALSE
        )
      }
    }
  }

  col_data <- S4Vectors::DataFrame(col_data)
  .assert_ids(rownames(col_data), "`col_data` row names")
  .assert_ids(names(col_data), "`col_data` column names")
  if (!identical(colnames(first), rownames(col_data))) {
    stop("Assay columns must exactly match `col_data` row names.", call. = FALSE)
  }

  if (is.null(row_data)) {
    row_data <- S4Vectors::DataFrame(row.names = rownames(first))
  } else {
    row_data <- S4Vectors::DataFrame(row_data)
    .assert_ids(rownames(row_data), "`row_data` row names")
    .assert_ids(names(row_data), "`row_data` column names")
    if (!identical(rownames(first), rownames(row_data))) {
      stop("Assay rows must exactly match `row_data` row names.", call. = FALSE)
    }
  }

  SummarizedExperiment::SummarizedExperiment(
    assays = assays,
    rowData = row_data,
    colData = col_data
  )
}

#' Construct a MultiAssayExperiment
#'
#' @param experiments A named list or `ExperimentList` of experiment leaves.
#' @param col_data Primary sample metadata with unique row names.
#' @param sample_map Optional explicit MAE sample map. When omitted, experiment
#'   column names must map directly to primary sample names.
#'
#' @return A [MultiAssayExperiment::MultiAssayExperiment].
#'
#' @examples
#' counts <- matrix(
#'   1:6,
#'   nrow = 3,
#'   dimnames = list(paste0("gene", 1:3), c("sample1", "sample2"))
#' )
#' samples <- data.frame(
#'   condition = c("control", "treated"),
#'   row.names = colnames(counts)
#' )
#' rna <- mae_create_experiment(counts, samples)
#' mae_create(list(rna = rna), samples)
#' @export
mae_create <- function(experiments, col_data, sample_map = NULL) {
  if (!length(experiments)) {
    stop("`experiments` must contain at least one experiment.", call. = FALSE)
  }
  if (
    is.null(names(experiments)) || anyNA(names(experiments)) ||
      any(!nzchar(names(experiments))) || anyDuplicated(names(experiments))
  ) {
    stop("`experiments` must have unique, non-empty names.", call. = FALSE)
  }
  is_se <- vapply(
    experiments,
    methods::is,
    logical(1),
    class2 = "SummarizedExperiment"
  )
  if (!all(is_se)) {
    stop(
      "Every experiment must inherit from SummarizedExperiment.",
      call. = FALSE
    )
  }
  col_data <- S4Vectors::DataFrame(col_data)
  .assert_ids(rownames(col_data), "`col_data` primary sample row names")
  .assert_ids(names(col_data), "`col_data` column names")

  arguments <- list(
    experiments = MultiAssayExperiment::ExperimentList(experiments),
    colData = col_data
  )
  if (!is.null(sample_map)) {
    sample_map <- as.data.frame(sample_map, stringsAsFactors = FALSE)
    required_columns <- c("assay", "primary", "colname")
    if (!all(required_columns %in% names(sample_map))) {
      stop(
        "`sample_map` must contain `assay`, `primary`, and `colname` columns.",
        call. = FALSE
      )
    }
    sample_map <- sample_map[, required_columns, drop = FALSE]
    invalid_value <- vapply(
      sample_map,
      function(value) {
        !is.atomic(value) || anyNA(value) || any(!nzchar(as.character(value)))
      },
      logical(1)
    )
    if (any(invalid_value)) {
      stop("`sample_map` values must be non-missing identifiers.", call. = FALSE)
    }
    sample_map[] <- lapply(sample_map, as.character)
    if (anyDuplicated(sample_map[c("assay", "colname")])) {
      stop(
        "Each `sample_map` assay/colname pair must be unique.",
        call. = FALSE
      )
    }
    if (any(!sample_map$assay %in% names(experiments))) {
      stop("`sample_map$assay` contains an unknown experiment.", call. = FALSE)
    }
    if (any(!sample_map$primary %in% rownames(col_data))) {
      stop("`sample_map$primary` contains an unknown primary sample.", call. = FALSE)
    }
    valid_colname <- vapply(
      seq_len(nrow(sample_map)),
      function(index) {
        sample_map$colname[[index]] %in%
          colnames(experiments[[sample_map$assay[[index]]]])
      },
      logical(1)
    )
    if (any(!valid_colname)) {
      stop(
        "`sample_map$colname` contains an unknown experiment column.",
        call. = FALSE
      )
    }
    sample_map$assay <- factor(sample_map$assay, levels = names(experiments))
    arguments$sampleMap <- S4Vectors::DataFrame(sample_map)
  }

  result <- do.call(MultiAssayExperiment::MultiAssayExperiment, arguments)
  mae_validate(result)
  result
}

#' Import transcript abundance estimates with tximport
#'
#' The returned object is an experiment leaf, not a stateful analysis object.
#' Add it to a `MultiAssayExperiment` with [mae_create()].
#'
#' @param files Named character vector of quantification files.
#' @param col_data Sample metadata; row names must equal `names(files)`.
#' @param type Quantifier type accepted by [tximport::tximport()].
#' @param tx2gene Optional transcript-to-gene mapping.
#' @param counts_from_abundance Value passed to `countsFromAbundance`.
#' @param ... Additional arguments passed to [tximport::tximport()].
#'
#' @return A `SummarizedExperiment` with available `counts`, `abundance`, and
#'   `length` assays.
#' @export
import_tximport <- function(
    files,
    col_data,
    type = "salmon",
    tx2gene = NULL,
    counts_from_abundance = "no",
    ...) {
  .require_backend("tximport", "to import transcript abundance estimates")
  if (
    !is.character(files) || !length(files) || anyNA(files) ||
      any(!nzchar(files))
  ) {
    stop("`files` must contain one non-empty path per sample.", call. = FALSE)
  }
  .assert_ids(names(files), "Sample names in `files`")
  missing_files <- files[!file.exists(files)]
  if (length(missing_files)) {
    stop(
      "Quantification files do not exist: ",
      paste(utils::head(missing_files, 10L), collapse = ", "),
      call. = FALSE
    )
  }

  col_data <- as.data.frame(col_data, optional = TRUE)
  .assert_ids(rownames(col_data), "`col_data` row names")
  if (!identical(names(files), rownames(col_data))) {
    stop("`names(files)` must exactly match `col_data` row names.", call. = FALSE)
  }

  txi <- tximport::tximport(
    files = files,
    type = type,
    tx2gene = tx2gene,
    countsFromAbundance = counts_from_abundance,
    ...
  )
  assay_keys <- intersect(c("counts", "abundance", "length"), names(txi))
  assays <- txi[assay_keys]

  se <- mae_create_experiment(assays = assays, col_data = col_data)
  object_metadata <- S4Vectors::metadata(se)
  object_metadata$bulkMAE <- list(
    importer = "tximport",
    counts_from_abundance = txi$countsFromAbundance %||% counts_from_abundance
  )
  S4Vectors::metadata(se) <- object_metadata
  se
}

.tximport_list <- function(se, assay = "counts") {
  selected_assay <- if (is.character(assay)) {
    assay
  } else {
    SummarizedExperiment::assayNames(se)[assay]
  }
  info <- S4Vectors::metadata(se)$bulkMAE
  required <- c("counts", "abundance", "length")

  if (
    !identical(selected_assay, "counts") ||
      is.null(info) ||
      !identical(info$importer, "tximport") ||
      !all(required %in% SummarizedExperiment::assayNames(se))
  ) {
    return(NULL)
  }

  matrices <- lapply(required, function(name) {
    as.matrix(SummarizedExperiment::assay(se, name))
  })
  names(matrices) <- required
  same_shape <- vapply(
    matrices,
    function(value) {
      identical(dim(value), dim(matrices$counts)) &&
        identical(dimnames(value), dimnames(matrices$counts))
    },
    logical(1)
  )
  if (!all(same_shape)) {
    stop(
      "The tximport counts, abundance, and length assays must align exactly.",
      call. = FALSE
    )
  }
  valid_numeric <- vapply(
    matrices,
    function(value) {
      is.numeric(value) && all(is.finite(value)) && all(value >= 0)
    },
    logical(1)
  )
  if (!all(valid_numeric)) {
    stop(
      "The tximport assays must contain finite, non-negative numeric values.",
      call. = FALSE
    )
  }
  counts <- matrices$counts
  library_sizes <- colSums(counts)
  if (any(!is.finite(library_sizes)) || any(library_sizes <= 0)) {
    stop(
      "The tximport counts assay must be non-negative with positive libraries.",
      call. = FALSE
    )
  }
  if (any(matrices$length <= 0)) {
    stop("The tximport length assay must be strictly positive.", call. = FALSE)
  }

  list(
    counts = counts,
    abundance = matrices$abundance,
    length = matrices$length,
    countsFromAbundance = info$counts_from_abundance
  )
}

`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}
