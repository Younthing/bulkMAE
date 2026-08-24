#' Construct a one-experiment MAE from a matrix
#'
#' @param expression A feature-by-sample matrix-like expression object with
#'   unique feature and sample names.
#' @param samples Sample metadata whose row names exactly match the expression
#'   columns.
#' @param row_data Optional feature metadata whose row names exactly match the
#'   expression rows.
#' @param experiment Name of the experiment to create.
#' @param assay Name of the expression assay.
#'
#' @return A [MultiAssayExperiment::MultiAssayExperiment].
#' @export
mae_from_matrix <- function(
    expression,
    samples,
    row_data = NULL,
    experiment = "rna",
    assay = "counts"
) {
  .assert_scalar_character(experiment, "experiment")
  .assert_scalar_character(assay, "assay")
  samples <- as.data.frame(samples, optional = TRUE)
  .assert_ids(rownames(samples), "`samples` row names")
  leaf_samples <- data.frame(row.names = rownames(samples))
  leaf <- mae_create_experiment(
    assays = expression,
    col_data = leaf_samples,
    row_data = row_data,
    assay_name = assay
  )
  mae_create(stats::setNames(list(leaf), experiment), samples)
}

#' Extract feature metadata from one experiment
#'
#' @inheritParams mae_validate
#'
#' @return A base data frame aligned to the raw experiment rows.
#' @export
mae_feature_data <- function(x, experiment) {
  mae_validate(x, experiment)
  leaf <- MultiAssayExperiment::experiments(x)[[experiment]]
  result <- as.data.frame(
    SummarizedExperiment::rowData(leaf),
    optional = TRUE
  )
  if (!identical(rownames(result), rownames(leaf))) {
    stop("Experiment feature metadata are out of order.", call. = FALSE)
  }
  result
}

#' Add aligned feature metadata to one experiment
#'
#' This pure helper updates the raw experiment leaf and returns a modified copy
#' of `x`. It never changes feature identifiers or assay values.
#'
#' @inheritParams mae_validate
#' @param value A feature-named atomic vector, or a feature-by-covariate data
#'   frame or matrix with feature row names.
#' @param name Optional name for a single metadata column. It is required for
#'   vector `value` and may rename a one-column data frame or matrix.
#' @param overwrite Replace metadata columns already present in the experiment.
#'
#' @return A `MultiAssayExperiment` copy containing the feature metadata.
#' @export
mae_add_feature_data <- function(
    x,
    experiment,
    value,
    name = NULL,
    overwrite = FALSE
) {
  mae_validate(x, experiment)
  .workflow_assert_flag(overwrite, "overwrite")
  leaf <- MultiAssayExperiment::experiments(x)[[experiment]]
  value <- .workflow_prepare_metadata(
    value = value,
    ids = rownames(leaf),
    axis = "feature",
    name = name
  )
  current <- SummarizedExperiment::rowData(leaf)
  .workflow_assert_columns_available(
    incoming = names(value),
    existing = names(current),
    overwrite = overwrite,
    label = "feature metadata"
  )
  for (column in names(value)) {
    current[[column]] <- value[[column]]
  }
  SummarizedExperiment::rowData(leaf) <- current
  .workflow_replace_experiment(x, experiment, leaf)
}

#' Add aligned sample metadata to one experiment
#'
#' This pure helper updates `colData` on the raw experiment leaf and returns a
#' modified copy of `x`. New columns may not duplicate primary MAE metadata;
#' `overwrite` applies only to columns already stored on the raw leaf.
#'
#' @inheritParams mae_validate
#' @param value A sample-named atomic vector, or a sample-by-covariate data
#'   frame or matrix with assay-sample row names.
#' @param name Optional name for a single metadata column. It is required for
#'   vector `value` and may rename a one-column data frame or matrix.
#' @param overwrite Replace metadata columns already present on the raw
#'   experiment leaf.
#'
#' @return A `MultiAssayExperiment` copy containing the sample metadata.
#' @export
mae_add_sample_data <- function(
    x,
    experiment,
    value,
    name = NULL,
    overwrite = FALSE
) {
  mae_validate(x, experiment)
  .workflow_assert_flag(overwrite, "overwrite")
  leaf <- MultiAssayExperiment::experiments(x)[[experiment]]
  value <- .workflow_prepare_metadata(
    value = value,
    ids = colnames(leaf),
    axis = "sample",
    name = name
  )

  primary_names <- names(MultiAssayExperiment::colData(x))
  primary_collision <- intersect(names(value), primary_names)
  if (length(primary_collision)) {
    stop(
      "Sample metadata columns duplicate primary MAE metadata: ",
      paste(primary_collision, collapse = ", "),
      ". Update primary metadata explicitly instead of shadowing it on an ",
      "experiment leaf.",
      call. = FALSE
    )
  }

  current <- SummarizedExperiment::colData(leaf)
  .workflow_assert_columns_available(
    incoming = names(value),
    existing = names(current),
    overwrite = overwrite,
    label = "experiment sample metadata"
  )
  for (column in names(value)) {
    current[[column]] <- value[[column]]
  }
  SummarizedExperiment::colData(leaf) <- current
  .workflow_replace_experiment(x, experiment, leaf)
}

#' Add an experiment while preserving the MAE sample map
#'
#' A matrix is wrapped in a `SummarizedExperiment`. With neither
#' `source_experiment` nor `sample_map`, its column names must be primary MAE
#' sample identifiers. `source_experiment` copies an existing experiment's map
#' when the new value has the same assay-column identifiers.
#'
#' @param x A `MultiAssayExperiment` object.
#' @param value A feature-by-sample matrix-like object or a
#'   `SummarizedExperiment`.
#' @param name Name for the new experiment.
#' @param assay_name Assay name used when `value` is a matrix-like object.
#' @param source_experiment Optional existing experiment whose sample map is
#'   copied for the new experiment.
#' @param sample_map Optional map for the new experiment, with `primary` and
#'   `colname` columns and an optional `assay` column equal to `name`.
#' @param overwrite Replace an existing experiment of the same name.
#'
#' @return A `MultiAssayExperiment` copy containing the new experiment.
#' @export
mae_add_experiment <- function(
    x,
    value,
    name,
    assay_name = "value",
    source_experiment = NULL,
    sample_map = NULL,
    overwrite = FALSE
) {
  mae_validate(x)
  .assert_scalar_character(name, "name")
  .assert_scalar_character(assay_name, "assay_name")
  .workflow_assert_flag(overwrite, "overwrite")
  if (!is.null(source_experiment) && !is.null(sample_map)) {
    stop(
      "Supply only one of `source_experiment` and `sample_map`.",
      call. = FALSE
    )
  }
  if (!is.null(source_experiment)) {
    .assert_scalar_character(source_experiment, "source_experiment")
    if (!source_experiment %in% mae_experiments(x)) {
      stop("Unknown source experiment: ", source_experiment, call. = FALSE)
    }
  }

  exists <- name %in% mae_experiments(x)
  if (exists && !overwrite) {
    stop(
      "Experiment `", name, "` already exists; set `overwrite = TRUE` to replace it.",
      call. = FALSE
    )
  }

  if (methods::is(value, "SummarizedExperiment")) {
    leaf <- value
    if (!nrow(leaf) || !ncol(leaf)) {
      stop("`value` must contain features and samples.", call. = FALSE)
    }
    .assert_ids(rownames(leaf), "New experiment feature names")
    .assert_ids(colnames(leaf), "New experiment sample names")
  } else {
    if (is.null(dim(value)) || length(dim(value)) != 2L) {
      stop(
        "`value` must be a matrix-like object or SummarizedExperiment.",
        call. = FALSE
      )
    }
    .assert_ids(colnames(value), "New experiment sample names")
    leaf <- mae_create_experiment(
      assays = value,
      col_data = data.frame(row.names = colnames(value)),
      assay_name = assay_name
    )
  }

  old_map <- as.data.frame(
    MultiAssayExperiment::sampleMap(x),
    optional = TRUE,
    stringsAsFactors = FALSE
  )
  old_map[] <- lapply(old_map, as.character)
  new_map <- .workflow_new_experiment_map(
    x = x,
    leaf = leaf,
    name = name,
    source_experiment = source_experiment,
    sample_map = sample_map,
    old_map = old_map
  )

  experiment_list <- MultiAssayExperiment::experiments(x)
  experiments <- stats::setNames(
    lapply(seq_along(experiment_list), function(index) experiment_list[[index]]),
    names(experiment_list)
  )
  experiments[[name]] <- leaf
  old_map <- old_map[old_map$assay != name, , drop = FALSE]
  combined_map <- rbind(old_map, new_map)
  rownames(combined_map) <- NULL
  primary_data <- as.data.frame(
    MultiAssayExperiment::colData(x),
    optional = TRUE
  )
  result <- mae_create(
    experiments = experiments,
    col_data = primary_data,
    sample_map = combined_map
  )
  S4Vectors::metadata(result) <- S4Vectors::metadata(x)
  MultiAssayExperiment::drops(result) <- MultiAssayExperiment::drops(x)
  result
}

#' Extract aligned adjustment covariates
#'
#' Supports an SVA result containing `sv`, an RUV matrix result containing `W`,
#' an RUVSeq `SeqExpressionSet` with `W_` phenotype columns, or an already
#' sample-by-covariate matrix/data frame. Sample names are never inferred from
#' row positions.
#'
#' @param result An adjustment result or sample-by-covariate object.
#' @param columns Optional unique column names or indices to retain.
#'
#' @return A finite numeric data frame with unique sample row names.
#' @export
adjust_covariates <- function(result, columns = NULL) {
  sample_names <- NULL
  default_columns <- NULL

  if (methods::is(result, "SeqExpressionSet")) {
    .require_backend("Biobase", "to extract RUVSeq phenotype covariates")
    phenotype <- getExportedValue("Biobase", "pData")(result)
    default_columns <- names(phenotype)[startsWith(names(phenotype), "W_")]
    if (!length(default_columns)) {
      stop(
        "The SeqExpressionSet contains no phenotype columns beginning with `W_`.",
        call. = FALSE
      )
    }
    value <- phenotype[, default_columns, drop = FALSE]
  } else if (is.list(result) && !is.null(result$sv)) {
    value <- result$sv
    sample_names <- attr(result, "bulkMAE_sample_names", exact = TRUE)
    default_columns <- "SV"
  } else if (is.list(result) && !is.null(result$W)) {
    value <- result$W
    sample_names <- attr(result, "bulkMAE_sample_names", exact = TRUE)
    if (
      is.null(sample_names) && !is.null(result$normalizedCounts) &&
        !is.null(dim(result$normalizedCounts))
    ) {
      sample_names <- colnames(result$normalizedCounts)
    }
    default_columns <- "W_"
  } else if (is.matrix(result) || is.data.frame(result)) {
    value <- result
  } else {
    stop(
      "`result` must contain SVA `sv`, RUV `W`, be a SeqExpressionSet, or be ",
      "a sample-by-covariate matrix/data frame.",
      call. = FALSE
    )
  }

  value <- .workflow_covariate_frame(
    value,
    sample_names = sample_names,
    default_columns = default_columns
  )
  .workflow_select_covariates(value, columns)
}

.workflow_prepare_metadata <- function(value, ids, axis, name = NULL) {
  .assert_ids(ids, paste0("Experiment ", axis, " names"))
  if (!is.null(name)) {
    .assert_scalar_character(name, "name")
  }

  if (is.atomic(value) && is.null(dim(value))) {
    if (is.null(name)) {
      stop("`name` is required when `value` is a vector.", call. = FALSE)
    }
    value_ids <- names(value)
    .assert_ids(value_ids, paste0("Names of ", axis, " metadata vector"))
    result <- data.frame(value, check.names = FALSE)
    names(result) <- name
    rownames(result) <- value_ids
  } else if (is.data.frame(value) || is.matrix(value)) {
    result <- as.data.frame(value, optional = TRUE)
    if (!ncol(result)) {
      stop("`value` must contain at least one metadata column.", call. = FALSE)
    }
    .assert_ids(
      rownames(result),
      paste0("Row names of ", axis, " metadata")
    )
    .assert_ids(
      names(result),
      paste0("Column names of ", axis, " metadata")
    )
    if (!is.null(name)) {
      if (ncol(result) != 1L) {
        stop("`name` can rename only one metadata column.", call. = FALSE)
      }
      names(result) <- name
    }
  } else {
    stop(
      "`value` must be a named vector or a data frame/matrix with row names.",
      call. = FALSE
    )
  }

  missing <- setdiff(ids, rownames(result))
  extra <- setdiff(rownames(result), ids)
  if (length(missing) || length(extra)) {
    stop(
      "`value` row names must contain every ", axis, " exactly once. Missing: ",
      .workflow_format_ids(missing), "; extra: ",
      .workflow_format_ids(extra), ".",
      call. = FALSE
    )
  }
  result[ids, , drop = FALSE]
}

.workflow_assert_columns_available <- function(
    incoming,
    existing,
    overwrite,
    label
) {
  collision <- intersect(incoming, existing)
  if (length(collision) && !overwrite) {
    stop(
      "Columns already exist in ", label, ": ",
      paste(collision, collapse = ", "),
      ". Set `overwrite = TRUE` to replace them.",
      call. = FALSE
    )
  }
  invisible(incoming)
}

.workflow_replace_experiment <- function(x, experiment, leaf) {
  experiments <- MultiAssayExperiment::experiments(x)
  experiments[[experiment]] <- leaf
  MultiAssayExperiment::experiments(x) <- experiments
  mae_validate(x, experiment)
  x
}

.workflow_new_experiment_map <- function(
    x,
    leaf,
    name,
    source_experiment,
    sample_map,
    old_map
) {
  sample_ids <- colnames(leaf)
  primary_ids <- rownames(MultiAssayExperiment::colData(x))

  if (!is.null(source_experiment)) {
    source_map <- old_map[old_map$assay == source_experiment, , drop = FALSE]
    if (!setequal(source_map$colname, sample_ids)) {
      stop(
        "New experiment sample names must equal the source experiment's ",
        "mapped column names.",
        call. = FALSE
      )
    }
    source_map$assay <- name
    result <- source_map[, c("assay", "primary", "colname"), drop = FALSE]
  } else if (!is.null(sample_map)) {
    result <- as.data.frame(sample_map, stringsAsFactors = FALSE)
    required <- c("primary", "colname")
    if (!all(required %in% names(result))) {
      stop("`sample_map` must contain `primary` and `colname` columns.", call. = FALSE)
    }
    if ("assay" %in% names(result)) {
      if (anyNA(result$assay) || any(as.character(result$assay) != name)) {
        stop("`sample_map$assay` must equal the new experiment name.", call. = FALSE)
      }
    } else {
      result$assay <- name
    }
    result <- result[, c("assay", "primary", "colname"), drop = FALSE]
  } else {
    unknown <- setdiff(sample_ids, primary_ids)
    if (length(unknown)) {
      stop(
        "Without `source_experiment` or `sample_map`, new experiment columns ",
        "must be primary sample IDs. Unknown columns: ",
        paste(unknown, collapse = ", "),
        ".",
        call. = FALSE
      )
    }
    result <- data.frame(
      assay = name,
      primary = sample_ids,
      colname = sample_ids,
      stringsAsFactors = FALSE
    )
  }

  result[] <- lapply(result, as.character)
  invalid <- vapply(
    result,
    function(column) anyNA(column) || any(!nzchar(column)),
    logical(1)
  )
  if (any(invalid)) {
    stop("`sample_map` values must be non-missing identifiers.", call. = FALSE)
  }
  if (anyDuplicated(result$colname)) {
    stop("`sample_map$colname` must map each assay column once.", call. = FALSE)
  }
  if (!setequal(result$colname, sample_ids)) {
    stop(
      "`sample_map$colname` must contain every new experiment column exactly once.",
      call. = FALSE
    )
  }
  unknown_primary <- setdiff(result$primary, primary_ids)
  if (length(unknown_primary)) {
    stop(
      "`sample_map$primary` contains unknown primary samples: ",
      paste(unknown_primary, collapse = ", "),
      ".",
      call. = FALSE
    )
  }
  result
}

.workflow_covariate_frame <- function(
    value,
    sample_names = NULL,
    default_columns = NULL
) {
  if (is.atomic(value) && is.null(dim(value))) {
    if (is.null(sample_names)) sample_names <- names(value)
    value <- matrix(value, ncol = 1L)
  } else if (!is.matrix(value) && !is.data.frame(value)) {
    stop("Adjustment covariates must be a vector, matrix, or data frame.", call. = FALSE)
  }
  result <- as.data.frame(value, optional = TRUE)
  if (!nrow(result) || !ncol(result)) {
    stop("Adjustment covariates must contain samples and columns.", call. = FALSE)
  }

  if (is.null(rownames(result)) || .workflow_default_rownames(rownames(result))) {
    if (!is.null(sample_names)) rownames(result) <- sample_names
  }
  .assert_ids(rownames(result), "Adjustment covariate sample names")
  if (.workflow_default_rownames(rownames(result))) {
    stop(
      "Adjustment covariates must carry real sample row names; row positions ",
      "are not treated as sample identifiers.",
      call. = FALSE
    )
  }

  if (is.null(names(result)) || any(!nzchar(names(result)))) {
    prefix <- default_columns %||% "covariate"
    if (length(prefix) != 1L) prefix <- "covariate"
    names(result) <- paste0(prefix, seq_len(ncol(result)))
  } else if (
    identical(default_columns, "SV") &&
      identical(names(result), paste0("V", seq_len(ncol(result))))
  ) {
    names(result) <- paste0("SV", seq_len(ncol(result)))
  }
  .assert_ids(names(result), "Adjustment covariate column names")
  numeric_columns <- vapply(result, is.numeric, logical(1))
  if (!all(numeric_columns)) {
    stop("Adjustment covariates must contain only numeric columns.", call. = FALSE)
  }
  if (any(!is.finite(as.matrix(result)))) {
    stop("Adjustment covariates must contain only finite values.", call. = FALSE)
  }
  result
}

.workflow_select_covariates <- function(value, columns) {
  if (is.null(columns)) return(value)
  if (is.character(columns)) {
    if (
      !length(columns) || anyNA(columns) || any(!nzchar(columns)) ||
        anyDuplicated(columns)
    ) {
      stop("Character `columns` must contain unique column names.", call. = FALSE)
    }
    missing <- setdiff(columns, names(value))
    if (length(missing)) {
      stop(
        "Unknown adjustment covariate columns: ",
        paste(missing, collapse = ", "),
        ".",
        call. = FALSE
      )
    }
  } else if (
    !is.numeric(columns) || !length(columns) || anyNA(columns) ||
      any(!is.finite(columns)) || any(columns != trunc(columns)) ||
      any(columns < 1L) || any(columns > ncol(value)) ||
      anyDuplicated(columns)
  ) {
    stop("`columns` must contain unique valid names or indices.", call. = FALSE)
  }
  value[, columns, drop = FALSE]
}

.workflow_default_rownames <- function(value) {
  identical(value, as.character(seq_along(value)))
}

.workflow_assert_flag <- function(value, argument) {
  if (!is.logical(value) || length(value) != 1L || is.na(value)) {
    stop("`", argument, "` must be TRUE or FALSE.", call. = FALSE)
  }
  invisible(value)
}

.workflow_format_ids <- function(value, n = 10L) {
  if (!length(value)) return("none")
  shown <- paste(utils::head(value, n), collapse = ", ")
  if (length(value) > n) paste0(shown, ", ...") else shown
}
