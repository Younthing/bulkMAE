#' Validate the package's MultiAssayExperiment contract
#'
#' Checks only the structural assumptions shared by the stateless wrappers:
#' `x` is a [MultiAssayExperiment::MultiAssayExperiment], the selected
#' experiment exists, and its columns can be aligned to primary samples.
#'
#' @param x A `MultiAssayExperiment` object.
#' @param experiment Optional experiment name.
#'
#' @return `x`, invisibly.
#' @export
mae_validate <- function(x, experiment = NULL) {
  .assert_mae(x)

  available <- mae_experiments(x)
  if (!length(available)) {
    stop("`x` must contain at least one experiment.", call. = FALSE)
  }

  if (is.null(experiment)) {
    selected <- available
  } else {
    .assert_scalar_character(experiment, "experiment")
    if (!experiment %in% available) {
      stop("Unknown experiment: ", experiment, call. = FALSE)
    }
    selected <- experiment
  }
  for (name in selected) .pull_se(x, name)

  invisible(x)
}

#' List experiments in a MultiAssayExperiment
#'
#' @inheritParams mae_validate
#'
#' @return A character vector.
#' @export
mae_experiments <- function(x) {
  .assert_mae(x)
  names(MultiAssayExperiment::experiments(x))
}

#' List assays in one experiment
#'
#' @inheritParams mae_validate
#'
#' @return A character vector.
#' @export
mae_assays <- function(x, experiment) {
  se <- .pull_se(x, experiment)
  SummarizedExperiment::assayNames(se)
}

#' Extract an experiment aligned to primary sample data
#'
#' Uses [MultiAssayExperiment::getWithColData()] so replicated or differently
#' named assay columns are mapped through `sampleMap` before analysis.
#'
#' @inheritParams mae_validate
#'
#' @return A `SummarizedExperiment` with aligned `colData`.
#' @export
mae_pull_experiment <- function(x, experiment) {
  .pull_se(x, experiment)
}

#' Extract an assay matrix
#'
#' @inheritParams mae_validate
#' @param assay Assay name or one-based assay index.
#'
#' @return A matrix-like assay object.
#' @export
mae_pull_assay <- function(x, experiment, assay = 1L) {
  se <- .pull_se(x, experiment)
  .extract_assay(se, assay)
}

#' Extract aligned sample metadata
#'
#' @inheritParams mae_validate
#'
#' @return A base `data.frame`, with row names matching assay columns.
#' @export
mae_samples <- function(x, experiment) {
  se <- .pull_se(x, experiment)
  data <- as.data.frame(SummarizedExperiment::colData(se), optional = TRUE)
  .assert_ids(rownames(data), "Sample metadata row names")
  .assert_ids(names(data), "Sample metadata column names")
  if (!identical(rownames(data), colnames(se))) {
    stop(
      "Aligned sample metadata row names must match experiment columns.",
      call. = FALSE
    )
  }
  data
}

#' Add an aligned assay to one experiment
#'
#' This pure helper returns a modified copy of `x`. It accepts a matrix, a
#' single-assay `SummarizedExperiment` such as a DESeq2 transformation, or a
#' limma `EList` containing an `E` matrix. Feature and sample identifiers are
#' aligned by name before insertion.
#'
#' @inheritParams mae_pull_assay
#' @param value Matrix-like values, a single-assay `SummarizedExperiment`, or
#'   an `EList` with an `E` matrix.
#' @param name Name for the new assay.
#' @param overwrite Replace an existing assay of the same name.
#'
#' @return A `MultiAssayExperiment` copy containing the added assay.
#' @export
mae_add_assay <- function(
    x,
    experiment,
    value,
    name,
    overwrite = FALSE
) {
  mae_validate(x, experiment)
  se <- MultiAssayExperiment::experiments(x)[[experiment]]
  .assert_scalar_character(name, "name")
  if (!is.logical(overwrite) || length(overwrite) != 1L || is.na(overwrite)) {
    stop("`overwrite` must be TRUE or FALSE.", call. = FALSE)
  }
  if (name %in% SummarizedExperiment::assayNames(se) && !overwrite) {
    stop(
      "Assay `", name, "` already exists; set `overwrite = TRUE` to replace it.",
      call. = FALSE
    )
  }

  if (methods::is(value, "SummarizedExperiment")) {
    available <- SummarizedExperiment::assayNames(value)
    if (length(available) != 1L) {
      stop(
        "A `SummarizedExperiment` value must contain exactly one assay.",
        call. = FALSE
      )
    }
    value <- SummarizedExperiment::assay(value, 1L)
  } else if (inherits(value, "EList")) {
    if (is.null(value$E)) {
      stop("An `EList` value must contain an `E` matrix.", call. = FALSE)
    }
    value <- value$E
  }
  value <- as.matrix(value)
  if (!is.numeric(value) && !is.logical(value)) {
    stop("`value` must contain numeric or logical assay values.", call. = FALSE)
  }
  if (is.null(rownames(value)) || is.null(colnames(value))) {
    stop("`value` must have feature and sample names.", call. = FALSE)
  }
  .assert_ids(rownames(value), "`value` feature names")
  .assert_ids(colnames(value), "`value` sample names")
  if (
    !setequal(rownames(value), rownames(se)) ||
      !setequal(colnames(value), colnames(se))
  ) {
    stop(
      "`value` must describe exactly the experiment features and samples.",
      call. = FALSE
    )
  }
  value <- value[rownames(se), colnames(se), drop = FALSE]
  SummarizedExperiment::assay(se, name, withDimnames = FALSE) <- value

  experiment_list <- MultiAssayExperiment::experiments(x)
  experiment_list[[experiment]] <- se
  MultiAssayExperiment::experiments(x) <- experiment_list
  mae_validate(x, experiment)
  x
}

#' Subset features in one experiment
#'
#' This pure helper returns a modified copy of `x`, leaving all other
#' experiments and the `sampleMap` unchanged.
#'
#' @inheritParams mae_pull_assay
#' @param features Unique feature names to retain, in the requested order.
#'
#' @return A `MultiAssayExperiment` copy with one experiment subset by row.
#' @export
mae_subset_features <- function(x, experiment, features) {
  mae_validate(x, experiment)
  se <- MultiAssayExperiment::experiments(x)[[experiment]]
  if (is.null(features)) {
    stop("`features` must contain feature names to retain.", call. = FALSE)
  }
  .select_features(
    as.matrix(SummarizedExperiment::assay(se, 1L)),
    features
  )
  se <- se[features, , drop = FALSE]

  experiment_list <- MultiAssayExperiment::experiments(x)
  experiment_list[[experiment]] <- se
  MultiAssayExperiment::experiments(x) <- experiment_list
  mae_validate(x, experiment)
  x
}

.assert_mae <- function(x) {
  if (!methods::is(x, "MultiAssayExperiment")) {
    stop("`x` must be a MultiAssayExperiment.", call. = FALSE)
  }
  invisible(x)
}

.assert_scalar_character <- function(x, argument) {
  if (!is.character(x) || length(x) != 1L || is.na(x) || !nzchar(x)) {
    stop("`", argument, "` must be one non-empty string.", call. = FALSE)
  }
  invisible(x)
}

.require_backend <- function(package, reason = NULL) {
  if (!requireNamespace(package, quietly = TRUE)) {
    message <- paste0("Package `", package, "` is required")
    if (!is.null(reason)) {
      message <- paste0(message, " ", reason)
    }
    stop(message, ".", call. = FALSE)
  }
  invisible(TRUE)
}

.with_attached_namespaces <- function(packages, code) {
  packages <- unique(packages)
  before <- search()
  on.exit({
    added <- setdiff(search(), before)
    for (name in added[startsWith(added, "package:")]) {
      if (name %in% search()) {
        try(detach(name, character.only = TRUE), silent = TRUE)
      }
    }
  }, add = TRUE)

  for (package in packages) {
    label <- paste0("package:", package)
    if (!label %in% search()) {
      suppressPackageStartupMessages(base::attachNamespace(package))
    }
  }
  force(code)
}

.pull_se <- function(x, experiment) {
  .assert_mae(x)
  .assert_scalar_character(experiment, "experiment")

  if (!experiment %in% mae_experiments(x)) {
    stop("Unknown experiment: ", experiment, call. = FALSE)
  }

  se <- withCallingHandlers(
    MultiAssayExperiment::getWithColData(x, experiment),
    warning = function(condition) {
      message <- conditionMessage(condition)
      if (
        startsWith(message, "Ignoring redundant column names in 'colData(x)'") ||
          identical(message, "'experiments' dropped; see 'drops()'")
      ) {
        invokeRestart("muffleWarning")
      }
    }
  )
  if (!methods::is(se, "SummarizedExperiment")) {
    stop(
      "Experiment `", experiment,
      "` must resolve to a SummarizedExperiment.",
      call. = FALSE
    )
  }

  if (!nrow(se) || !ncol(se)) {
    stop(
      "Experiment `", experiment, "` must contain features and samples.",
      call. = FALSE
    )
  }
  .assert_ids(rownames(se), "Feature names")
  .assert_ids(colnames(se), "Sample names")

  data_names <- rownames(SummarizedExperiment::colData(se))
  .assert_ids(data_names, "Sample metadata row names")
  if (!identical(colnames(se), data_names)) {
    stop(
      "Experiment columns and aligned sample metadata are out of order.",
      call. = FALSE
    )
  }

  se
}

.extract_assay <- function(se, assay) {
  available <- SummarizedExperiment::assayNames(se)

  if (is.character(assay)) {
    .assert_scalar_character(assay, "assay")
    if (!assay %in% available) {
      stop(
        "Unknown assay `", assay, "`. Available assays: ",
        paste(available, collapse = ", "),
        call. = FALSE
      )
    }
  } else if (
    !is.numeric(assay) || length(assay) != 1L || is.na(assay) ||
      !is.finite(assay) || assay != trunc(assay) || assay < 1L ||
      assay > length(available)
  ) {
    stop("`assay` must name or index one assay.", call. = FALSE)
  }

  SummarizedExperiment::assay(se, assay)
}

.pull_matrix <- function(x, experiment, assay = 1L, storage = "double") {
  value <- mae_pull_assay(x, experiment, assay)
  value <- as.matrix(value)

  if (!is.numeric(value) && !is.logical(value)) {
    stop("The selected assay must contain numeric values.", call. = FALSE)
  }
  storage.mode(value) <- storage

  if (!nrow(value) || !ncol(value)) {
    stop("The selected assay must contain features and samples.", call. = FALSE)
  }
  .assert_ids(rownames(value), "Feature names")
  .assert_ids(colnames(value), "Sample names")

  value
}

.as_count_matrix <- function(x, experiment, assay = "counts", integer = FALSE) {
  counts <- .pull_matrix(x, experiment, assay)

  if (any(!is.finite(counts)) || any(counts < 0)) {
    stop("Counts must be finite and non-negative.", call. = FALSE)
  }
  library_sizes <- colSums(counts)
  if (any(!is.finite(library_sizes)) || any(library_sizes <= 0)) {
    stop("Each sample must have a positive finite library size.", call. = FALSE)
  }
  if (integer && any(abs(counts - round(counts)) > sqrt(.Machine$double.eps))) {
    stop("This backend requires integer counts.", call. = FALSE)
  }
  if (integer && any(counts > .Machine$integer.max)) {
    stop("Counts exceed R's integer range for this backend.", call. = FALSE)
  }
  if (integer) {
    storage.mode(counts) <- "integer"
  }

  counts
}

.model_matrix <- function(formula, data) {
  if (!inherits(formula, "formula")) {
    stop("`formula` must be a formula.", call. = FALSE)
  }
  data <- as.data.frame(data, optional = TRUE)
  .assert_ids(rownames(data), "Sample metadata row names")
  .assert_ids(names(data), "Sample metadata column names")

  frame <- tryCatch(
    stats::model.frame(
      formula,
      data = data,
      na.action = stats::na.fail,
      drop.unused.levels = TRUE
    ),
    error = function(error) {
      stop(
        "Could not construct the model frame: ", conditionMessage(error),
        call. = FALSE
      )
    }
  )
  design <- stats::model.matrix(formula, data = frame)
  design <- .validate_design_matrix(design, data, "design")
  design
}

.column_or_vector <- function(value, data, argument) {
  data <- as.data.frame(data, optional = TRUE)
  sample_names <- rownames(data)
  .assert_ids(sample_names, "Sample metadata row names")
  .assert_ids(names(data), "Sample metadata column names")

  if (is.character(value) && length(value) == 1L && value %in% names(data)) {
    value <- data[[value]]
    names(value) <- sample_names
  } else if (!is.null(names(value))) {
    .assert_ids(names(value), paste0("Names of `", argument, "`"))
    missing <- setdiff(sample_names, names(value))
    extra <- setdiff(names(value), sample_names)
    if (length(missing) || length(extra)) {
      stop(
        "Named `", argument, "` must contain each sample exactly once.",
        call. = FALSE
      )
    }
    value <- value[sample_names]
  }

  if (!is.null(dim(value)) || (!is.atomic(value) && !is.factor(value))) {
    stop("`", argument, "` must be an atomic vector.", call. = FALSE)
  }
  if (length(value) != nrow(data)) {
    stop(
      "`", argument, "` must be a metadata column or one value per sample.",
      call. = FALSE
    )
  }
  if (anyNA(value)) {
    stop("`", argument, "` must not contain missing values.", call. = FALSE)
  }
  names(value) <- sample_names
  value
}

.validate_design_matrix <- function(
    design,
    data,
    argument = "design",
    full_rank = TRUE
) {
  data <- as.data.frame(data, optional = TRUE)
  sample_names <- rownames(data)
  .assert_ids(sample_names, "Sample metadata row names")

  if (is.null(dim(design)) || length(dim(design)) != 2L) {
    stop("`", argument, "` must be a matrix-like object.", call. = FALSE)
  }
  design <- as.matrix(design)
  if (!is.numeric(design)) {
    stop("`", argument, "` must be numeric.", call. = FALSE)
  }

  if (!is.null(rownames(design))) {
    .assert_ids(rownames(design), paste0("Row names of `", argument, "`"))
    missing <- setdiff(sample_names, rownames(design))
    extra <- setdiff(rownames(design), sample_names)
    if (length(missing) || length(extra)) {
      stop(
        "Row names of `", argument, "` must contain each sample exactly once.",
        call. = FALSE
      )
    }
    design <- design[sample_names, , drop = FALSE]
  } else if (nrow(design) != length(sample_names)) {
    stop("`", argument, "` must have one row per sample.", call. = FALSE)
  } else {
    rownames(design) <- sample_names
  }

  if (!ncol(design)) {
    stop("`", argument, "` must contain at least one column.", call. = FALSE)
  }
  if (any(!is.finite(design))) {
    stop("`", argument, "` must contain only finite values.", call. = FALSE)
  }
  if (is.null(colnames(design))) {
    colnames(design) <- paste0("V", seq_len(ncol(design)))
  }
  .assert_ids(colnames(design), paste0("Column names of `", argument, "`"))

  if (full_rank) {
    decomposition <- qr(design)
    if (decomposition$rank < ncol(design)) {
      dependent <- colnames(design)[
        decomposition$pivot[seq.int(decomposition$rank + 1L, ncol(design))]
      ]
      stop(
        "`", argument, "` is not full rank; aliased columns include: ",
        paste(dependent, collapse = ", "),
        call. = FALSE
      )
    }
  }

  design
}

.sample_matrix <- function(value, data, argument, numeric = FALSE) {
  data <- as.data.frame(data, optional = TRUE)
  sample_names <- rownames(data)
  .assert_ids(sample_names, "Sample metadata row names")

  if (is.atomic(value) && is.null(dim(value))) {
    value <- .column_or_vector(value, data, argument)
    value <- matrix(value, ncol = 1L, dimnames = list(sample_names, argument))
  } else {
    value <- as.matrix(value)
    if (!is.null(rownames(value))) {
      .assert_ids(rownames(value), paste0("Row names of `", argument, "`"))
      missing <- setdiff(sample_names, rownames(value))
      extra <- setdiff(rownames(value), sample_names)
      if (length(missing) || length(extra)) {
        stop(
          "Row names of `", argument, "` must contain each sample exactly once.",
          call. = FALSE
        )
      }
      value <- value[sample_names, , drop = FALSE]
    } else if (nrow(value) != length(sample_names)) {
      stop("`", argument, "` must have one row per sample.", call. = FALSE)
    } else {
      rownames(value) <- sample_names
    }
  }

  if (!ncol(value)) {
    stop("`", argument, "` must contain at least one column.", call. = FALSE)
  }
  if (is.null(colnames(value))) {
    colnames(value) <- paste0(argument, seq_len(ncol(value)))
  }
  .assert_ids(colnames(value), paste0("Column names of `", argument, "`"))
  if (numeric && !is.numeric(value)) {
    stop("`", argument, "` must contain only numeric columns.", call. = FALSE)
  }
  if (numeric && any(!is.finite(value))) {
    stop("`", argument, "` must contain only finite values.", call. = FALSE)
  }
  value
}

.with_seed <- function(seed, code) {
  if (is.null(seed)) {
    return(force(code))
  }
  if (
    !is.numeric(seed) || length(seed) != 1L || is.na(seed) ||
      !is.finite(seed) || seed < -.Machine$integer.max ||
      seed > .Machine$integer.max || seed != trunc(seed)
  ) {
    stop("`seed` must be NULL or one finite integer.", call. = FALSE)
  }

  has_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  old_kind <- RNGkind()
  if (has_seed) {
    old_seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  }
  on.exit({
    do.call(RNGkind, as.list(old_kind))
    if (has_seed) {
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)

  set.seed(as.integer(seed))
  force(code)
}

.assert_ids <- function(ids, description) {
  if (is.null(ids)) {
    stop(description, " must be present.", call. = FALSE)
  }
  if (
    !is.character(ids) || anyNA(ids) || any(!nzchar(ids)) ||
      anyDuplicated(ids)
  ) {
    stop(description, " must be non-empty and unique.", call. = FALSE)
  }
  invisible(ids)
}

.select_features <- function(matrix, features = NULL) {
  if (is.null(features)) {
    return(matrix)
  }
  if (
    !is.character(features) || !length(features) || anyNA(features) ||
      any(!nzchar(features)) || anyDuplicated(features)
  ) {
    stop("`features` must contain unique, non-empty feature names.", call. = FALSE)
  }
  missing <- setdiff(features, rownames(matrix))
  if (length(missing)) {
    stop(
      "Unknown features: ", paste(utils::head(missing, 10L), collapse = ", "),
      call. = FALSE
    )
  }
  matrix[features, , drop = FALSE]
}

#' bulkMAE: stateless bulk-transcriptomics adapters
#'
#' `bulkMAE` uses `MultiAssayExperiment` as its public data contract, resolves
#' one `SummarizedExperiment` leaf and its aligned sample metadata, then returns
#' the native result produced by an established analysis backend. It does not
#' mutate the input object or keep an analysis-status registry.
#'
#' Start with [mae_create()] or [import_tximport()], inspect alignment with
#' [mae_validate()] and [mae_samples()], and then call an analysis family such
#' as `qc_*`, `de_*`, `enrich_*`, `coexpr_*`, or `surv_*`. Statistical engines
#' are optional dependencies and are checked only when their adapter is called.
#'
#' @seealso `vignette("getting-started", package = "bulkMAE")`
#' @keywords internal
#' @aliases bulkMAE-package
"_PACKAGE"
