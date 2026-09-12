#' Construct a model design matrix from MAE sample metadata
#'
#' @inheritParams mae_samples
#' @param formula Model formula evaluated against aligned sample metadata.
#'
#' @return A numeric design matrix with rows aligned to assay samples.
#' @export
de_design <- function(x, experiment, formula) {
  .model_matrix(formula, mae_samples(x, experiment))
}

#' Construct limma-compatible contrasts from MAE sample metadata
#'
#' This helper keeps design-column construction and contrast parsing inside
#' `bulkMAE`. The returned matrix can be passed directly to [de_voom()] or
#' [de_limma()].
#'
#' @inheritParams de_design
#' @param contrasts Character contrast expressions accepted by
#'   [limma::makeContrasts()].
#'
#' @return A numeric contrast matrix whose rows match [de_design()] columns.
#' @export
de_contrast <- function(x, experiment, formula, contrasts) {
  .require_backend("limma", "to construct model contrasts")
  if (
    !is.character(contrasts) || !length(contrasts) || anyNA(contrasts) ||
      any(!nzchar(contrasts))
  ) {
    stop("`contrasts` must contain non-empty contrast expressions.", call. = FALSE)
  }
  design <- de_design(x, experiment, formula)
  result <- limma::makeContrasts(contrasts = contrasts, levels = design)
  # limma makes `(Intercept)` syntactic internally and returns `Intercept`.
  # The coefficient order is unchanged, so restore the design names expected
  # by contrasts.fit() and promised by this wrapper.
  rownames(result) <- colnames(design)
  result
}

#' Build a maSigPro experimental-design table from sample metadata
#'
#' @inheritParams mae_samples
#' @param time Metadata column containing numeric time or dose.
#' @param replicate Metadata column identifying biological replicates.
#' @param group One categorical metadata column, or multiple existing binary
#'   indicator-column names.
#'
#' @return A sample-aligned data frame containing `Time`, `Replicates`, and
#'   one binary column per group.
#' @export
de_masigpro_design <- function(x, experiment, time, replicate, group) {
  data <- mae_samples(x, experiment)
  .assert_metadata_column(data, time, "time")
  .assert_metadata_column(data, replicate, "replicate")
  if (identical(time, replicate)) {
    stop("`time` and `replicate` must name different columns.", call. = FALSE)
  }

  time_values <- data[[time]]
  replicate_values <- data[[replicate]]
  if (!is.numeric(time_values) || any(!is.finite(time_values))) {
    stop("The selected `time` column must be finite and numeric.", call. = FALSE)
  }
  if (anyNA(replicate_values) || any(!nzchar(as.character(replicate_values)))) {
    stop("The selected `replicate` column cannot contain missing or empty values.", call. = FALSE)
  }

  if (
    !is.character(group) || !length(group) || anyNA(group) ||
      any(!nzchar(group)) || anyDuplicated(group) ||
      any(!group %in% names(data))
  ) {
    stop("`group` must name one categorical column or binary metadata columns.", call. = FALSE)
  }
  if (any(group %in% c(time, replicate))) {
    stop("Group columns must differ from `time` and `replicate`.", call. = FALSE)
  }

  if (length(group) == 1L) {
    group_values <- data[[group]]
    if (anyNA(group_values) || any(!nzchar(as.character(group_values)))) {
      stop("The selected group column cannot contain missing or empty values.", call. = FALSE)
    }
    levels <- unique(as.character(group_values))
    indicator_names <- make.unique(c(
      "Time",
      "Replicates",
      make.names(levels, unique = TRUE)
    ))[-c(1L, 2L)]
    indicators <- vapply(
      levels,
      function(level) as.integer(as.character(group_values) == level),
      integer(nrow(data))
    )
    if (is.null(dim(indicators))) {
      indicators <- matrix(indicators, ncol = 1L)
    }
    colnames(indicators) <- indicator_names
    rownames(indicators) <- rownames(data)
    group_map <- stats::setNames(levels, indicator_names)
  } else {
    invalid <- vapply(
      data[, group, drop = FALSE],
      function(value) {
        (!is.numeric(value) && !is.logical(value)) ||
          anyNA(value) || !all(value %in% c(0, 1))
      },
      logical(1)
    )
    if (any(invalid)) {
      stop(
        "Multiple `group` columns must be binary 0/1 indicators: ",
        paste(group[invalid], collapse = ", "),
        ".",
        call. = FALSE
      )
    }
    indicators <- as.matrix(data[, group, drop = FALSE])
    storage.mode(indicators) <- "integer"
    indicator_names <- make.unique(c(
      "Time",
      "Replicates",
      make.names(group, unique = TRUE)
    ))[-c(1L, 2L)]
    colnames(indicators) <- indicator_names
    group_map <- stats::setNames(group, indicator_names)
  }
  if (any(rowSums(indicators) < 1L)) {
    stop("Every sample must belong to at least one group.", call. = FALSE)
  }

  result <- data.frame(
    Time = as.numeric(time_values),
    Replicates = replicate_values,
    indicators,
    row.names = rownames(data),
    check.names = FALSE
  )
  attr(result, "group_levels") <- group_map
  result
}

#' Convert native differential-expression results to a stable table
#'
#' The helper does not replace native result objects. It extracts the small set
#' of fields needed by downstream enrichment, drug-query, and meta-analysis
#' functions while retaining feature identifiers explicitly.
#' For edgeR likelihood-ratio and quasi-likelihood F tests, `statistic` is a
#' directional signed square-root score for ranking only. Such tests do not
#' expose a coefficient standard error, so `standard_error` remains `NA` and
#' [meta_collect()] rejects them instead of manufacturing an invalid value.
#' A DESeq2 likelihood-ratio statistic is retained in the table but marked as
#' non-directional; [de_ranks()] requires `column = "effect"` for that result,
#' or a Wald result when a directional test statistic is desired.
#' Generic data-frame `F` and `LR` columns are likewise retained as
#' non-directional because their degrees of freedom cannot be inferred safely.
#'
#' @param result A `DESeqResults`, edgeR `DGELRT`/`DGEExact`, limma/dream
#'   `MArrayLM`, or data-frame-like result.
#' @param coef Coefficient name or index for an `MArrayLM` with more than one
#'   coefficient.
#' @param feature_column,effect_column,standard_error_column,statistic_column
#'   Optional feature, effect, standard-error, and statistic column names for a
#'   generic data frame. Common backend names are detected automatically.
#' @param p_value_column,adjusted_p_value_column Optional raw and adjusted
#'   p-value column names for a generic data frame. Common backend names are
#'   detected automatically.
#'
#' @return A data frame with `feature_id`, `effect`, `standard_error`,
#'   `statistic`, `p_value`, and `adjusted_p_value` columns.
#' @export
de_table <- function(
    result,
    coef = NULL,
    feature_column = NULL,
    effect_column = NULL,
    standard_error_column = NULL,
    statistic_column = NULL,
    p_value_column = NULL,
    adjusted_p_value_column = NULL
) {
  statistic_type <- NULL
  if (inherits(result, "MArrayLM")) {
    table <- .de_table_limma(result, coef)
  } else if (inherits(result, c("DGELRT", "DGEExact"))) {
    table <- .de_table_edger(result)
  } else if (inherits(result, "TopTags") && !is.null(result$table)) {
    table <- .de_table_edger(result)
  } else if (methods::is(result, "DESeqResults")) {
    statistic_type <- .de_deseq_statistic_type(result)
    table <- as.data.frame(result, optional = TRUE)
  } else if (is.data.frame(result) || is.matrix(result)) {
    table <- as.data.frame(result, optional = TRUE)
  } else {
    stop(
      "`result` must be a DESeqResults, edgeR test, MArrayLM, or data frame.",
      call. = FALSE
    )
  }
  standardized <- .de_standardize_table(
    table,
    feature_column = feature_column,
    effect_column = effect_column,
    standard_error_column = standard_error_column,
    statistic_column = statistic_column,
    p_value_column = p_value_column,
    adjusted_p_value_column = adjusted_p_value_column
  )
  if (!is.null(statistic_type)) {
    attr(standardized, "statistic_type") <- statistic_type
  }
  standardized
}

#' Extract finite named ranks from differential-expression results
#'
#' @param result Input accepted by [de_table()].
#' @param column Normalized [de_table()] column used for ranking.
#' @param decreasing Sort ranks from high to low.
#' @param ... Additional arguments passed to [de_table()].
#'
#' @return A finite named numeric vector suitable for `enrich_*()` and
#'   [drug_query()].
#' @export
de_ranks <- function(
    result,
    column = c("statistic", "effect"),
    decreasing = TRUE,
    ...
) {
  column <- match.arg(column)
  if (!is.logical(decreasing) || length(decreasing) != 1L || is.na(decreasing)) {
    stop("`decreasing` must be TRUE or FALSE.", call. = FALSE)
  }
  table <- de_table(result, ...)
  statistic_type <- attr(table, "statistic_type", exact = TRUE)
  if (identical(column, "statistic") && !is.null(statistic_type)) {
    if (identical(statistic_type, "nondirectional_lrt")) {
      stop(
        "DESeq2 LRT statistics are non-directional. Use `column = \"effect\"` ",
        "or fit a Wald test for directional statistic ranks.",
        call. = FALSE
      )
    }
    if (identical(statistic_type, "nondirectional_f_or_lr")) {
      stop(
        "Generic F/LR statistics may be multi-degree-of-freedom and are ",
        "non-directional. Use `column = \"effect\"` or supply a separately ",
        "defined directional statistic.",
        call. = FALSE
      )
    }
  }
  values <- table[[column]]
  keep <- is.finite(values)
  if (!any(keep)) {
    stop("No finite `", column, "` values are available for ranking.", call. = FALSE)
  }
  values <- stats::setNames(as.numeric(values[keep]), table$feature_id[keep])
  sort(values, decreasing = decreasing)
}

#' Select differential features as a named logical vector
#'
#' The returned vector covers every tested feature and can be passed directly
#' to [enrich_goseq()]. Its names also define the tested universe for ORA.
#'
#' @param result Input accepted by [de_table()].
#' @param fdr Maximum adjusted p-value.
#' @param min_abs_effect Minimum absolute effect size.
#' @param direction Select both directions, positive effects, or negative
#'   effects.
#' @param ... Additional arguments passed to [de_table()].
#'
#' @return A named logical vector over all result features.
#' @export
de_selected <- function(
    result,
    fdr = 0.05,
    min_abs_effect = 0,
    direction = c("both", "up", "down"),
    ...
) {
  direction <- match.arg(direction)
  .differential_assert_probability(fdr, "fdr")
  if (
    !is.numeric(min_abs_effect) || length(min_abs_effect) != 1L ||
      is.na(min_abs_effect) || !is.finite(min_abs_effect) || min_abs_effect < 0
  ) {
    stop("`min_abs_effect` must be one finite, non-negative number.", call. = FALSE)
  }
  table <- de_table(result, ...)
  if (all(is.na(table$adjusted_p_value))) {
    stop("Adjusted p-values are unavailable in `result`.", call. = FALSE)
  }
  selected <- !is.na(table$adjusted_p_value) &
    table$adjusted_p_value <= fdr &
    is.finite(table$effect) & table$effect != 0 &
    abs(table$effect) >= min_abs_effect
  if (direction == "up") selected <- selected & table$effect > 0
  if (direction == "down") selected <- selected & table$effect < 0
  stats::setNames(selected, table$feature_id)
}

#' Convert long-format activity results to a source-by-sample matrix
#'
#' @param result A data-frame-like result returned by [activity_decouple()],
#'   [activity_progeny()], or [activity_tf()].
#' @param value,source,sample Column names containing activity values, sources,
#'   and sample identifiers.
#' @param statistic Optional statistic/method value to retain when a result
#'   contains multiple methods.
#' @param statistic_column Column containing statistic/method names.
#'
#' @return A numeric source-by-sample matrix suitable for
#'   [mae_add_experiment()].
#' @export
activity_matrix <- function(
    result,
    value = "score",
    source = "source",
    sample = "condition",
    statistic = NULL,
    statistic_column = "statistic"
) {
  result <- as.data.frame(result, optional = TRUE)
  for (argument in c("value", "source", "sample")) {
    column <- get(argument)
    .de_assert_column(result, column, argument)
  }
  if (!is.null(statistic)) {
    .assert_scalar_character(statistic, "statistic")
    .de_assert_column(result, statistic_column, "statistic_column")
    result <- result[as.character(result[[statistic_column]]) == statistic, , drop = FALSE]
    if (!nrow(result)) {
      stop("No rows match the requested `statistic`.", call. = FALSE)
    }
  } else if (statistic_column %in% names(result)) {
    available <- unique(as.character(result[[statistic_column]]))
    available <- available[!is.na(available)]
    if (length(available) > 1L) {
      stop(
        "`statistic` is required when results contain multiple statistics: ",
        paste(available, collapse = ", "),
        ".",
        call. = FALSE
      )
    }
  }
  source_ids <- as.character(result[[source]])
  sample_ids <- as.character(result[[sample]])
  values <- result[[value]]
  if (
    anyNA(source_ids) || any(!nzchar(source_ids)) || anyNA(sample_ids) ||
      any(!nzchar(sample_ids))
  ) {
    stop("Activity source and sample identifiers cannot be missing or empty.", call. = FALSE)
  }
  if (!is.numeric(values) || any(!is.finite(values))) {
    stop("The activity `value` column must contain finite numeric values.", call. = FALSE)
  }
  keys <- paste(source_ids, sample_ids, sep = "\r")
  if (anyDuplicated(keys)) {
    stop("Activity results contain duplicate source/sample rows.", call. = FALSE)
  }
  sources <- unique(source_ids)
  samples <- unique(sample_ids)
  matrix <- matrix(
    NA_real_,
    nrow = length(sources),
    ncol = length(samples),
    dimnames = list(sources, samples)
  )
  matrix[cbind(match(source_ids, sources), match(sample_ids, samples))] <- values
  if (anyNA(matrix)) {
    stop("Activity results do not contain a complete source-by-sample grid.", call. = FALSE)
  }
  matrix
}

#' Extract named WGCNA module labels
#'
#' @param fit A native result returned by [coexpr_wgcna()].
#'
#' @return A feature-named vector accepted by [coexpr_preservation()].
#' @export
coexpr_modules <- function(fit) {
  if (!is.list(fit) || is.null(fit$colors)) {
    stop("`fit` must be a WGCNA result containing `colors`.", call. = FALSE)
  }
  colors <- fit$colors
  if (
    (!is.atomic(colors) && !is.factor(colors)) || is.null(names(colors)) ||
      anyNA(names(colors)) || any(!nzchar(names(colors))) ||
      anyDuplicated(names(colors)) || anyNA(colors) ||
      any(!nzchar(as.character(colors)))
  ) {
    stop("WGCNA `colors` must be named by unique, non-missing features.", call. = FALSE)
  }
  colors
}

#' Correlate WGCNA module eigengenes with sample traits
#'
#' Uses the `MEs` matrix stored on a [coexpr_wgcna()] result. Trait columns
#' are aligned by sample name. Factor and character columns are converted to
#' integer codes; logical columns become 0/1. This is an extractor, not a new
#' result class.
#'
#' @param fit A native result returned by [coexpr_wgcna()] containing `MEs`.
#' @param traits A sample-named atomic vector, or a sample-by-trait data frame
#'   or matrix whose row names cover every eigengene sample.
#'
#' @return A list with `correlation` and `p_value` matrices (modules by
#'   traits).
#' @export
coexpr_module_trait <- function(fit, traits) {
  mes <- .coexpr_module_eigengenes(fit)
  traits <- .coexpr_trait_frame(traits, rownames(mes))
  modules <- colnames(mes)
  trait_names <- names(traits)
  correlation <- matrix(
    NA_real_,
    nrow = length(modules),
    ncol = length(trait_names),
    dimnames = list(modules, trait_names)
  )
  p_value <- correlation
  for (module in modules) {
    for (trait in trait_names) {
      left <- mes[, module]
      right <- traits[[trait]]
      if (stats::sd(left) == 0 || stats::sd(right) == 0) {
        next
      }
      test <- stats::cor.test(left, right, method = "pearson")
      correlation[module, trait] <- unname(as.numeric(test$estimate))
      p_value[module, trait] <- test$p.value
    }
  }
  list(correlation = correlation, p_value = p_value)
}

#' Compute feature-to-module membership from a WGCNA fit
#'
#' Pearson-correlates each retained feature with each module eigengene. The
#' expression assay is aligned to `fit$MEs` by sample name and to
#' [coexpr_modules()] by feature name.
#'
#' @param fit A native result returned by [coexpr_wgcna()] containing `MEs`
#'   and `colors`.
#' @inheritParams mae_pull_assay
#'
#' @return A feature-by-module numeric matrix of membership values.
#' @export
coexpr_membership <- function(fit, x, experiment, assay) {
  mes <- .coexpr_module_eigengenes(fit)
  modules <- coexpr_modules(fit)
  matrix <- as.matrix(.pull_matrix(x, experiment, assay))
  features <- intersect(names(modules), rownames(matrix))
  samples <- intersect(rownames(mes), colnames(matrix))
  if (length(features) < 2L) {
    stop(
      "The WGCNA fit and assay must share at least two named features.",
      call. = FALSE
    )
  }
  if (length(samples) < 3L) {
    stop(
      "The WGCNA eigengenes and assay must share at least three samples.",
      call. = FALSE
    )
  }
  expression <- t(matrix[features, samples, drop = FALSE])
  mes <- mes[samples, , drop = FALSE]
  if (any(!is.finite(expression)) || any(!is.finite(mes))) {
    stop("Module membership requires finite expression and eigengenes.", call. = FALSE)
  }
  stats::cor(expression, mes)
}

#' Rank intramodular hub features from a membership matrix
#'
#' @param membership A feature-by-module matrix from [coexpr_membership()].
#' @param modules A feature-named module vector from [coexpr_modules()].
#' @param n Positive number of hubs to retain per module.
#'
#' @return A data frame with `feature`, `module`, and `membership` columns,
#'   ordered by module and decreasing absolute membership.
#' @export
coexpr_hubs <- function(membership, modules, n = 5) {
  membership <- .coexpr_membership_matrix(membership)
  if (
    (!is.atomic(modules) && !is.factor(modules)) || is.null(names(modules)) ||
      anyNA(names(modules)) || any(!nzchar(names(modules))) ||
      anyDuplicated(names(modules)) || anyNA(modules) ||
      any(!nzchar(as.character(modules)))
  ) {
    stop("`modules` must be named by unique, non-missing features.", call. = FALSE)
  }
  if (
    !is.numeric(n) || length(n) != 1L || is.na(n) || !is.finite(n) ||
      n != trunc(n) || n < 1L
  ) {
    stop("`n` must be a positive integer.", call. = FALSE)
  }
  features <- intersect(rownames(membership), names(modules))
  if (!length(features)) {
    stop("`membership` and `modules` must share feature names.", call. = FALSE)
  }
  rows <- lapply(unique(as.character(modules[features])), function(module) {
    members <- features[as.character(modules[features]) == module]
    column <- .coexpr_match_module_column(module, colnames(membership))
    values <- membership[members, column]
    keep <- order(abs(values), decreasing = TRUE)
    keep <- keep[seq_len(min(as.integer(n), length(keep)))]
    data.frame(
      feature = members[keep],
      module = module,
      membership = unname(values[keep]),
      row.names = NULL,
      check.names = FALSE
    )
  })
  do.call(rbind, rows)
}

.coexpr_module_eigengenes <- function(fit) {
  if (!is.list(fit) || is.null(fit$MEs)) {
    stop("`fit` must be a WGCNA result containing `MEs`.", call. = FALSE)
  }
  mes <- as.data.frame(fit$MEs, optional = TRUE)
  if (!nrow(mes) || !ncol(mes)) {
    stop("WGCNA `MEs` must contain samples and modules.", call. = FALSE)
  }
  .assert_ids(rownames(mes), "Module eigengene sample names")
  .assert_ids(names(mes), "Module eigengene names")
  numeric_columns <- vapply(mes, is.numeric, logical(1))
  if (!all(numeric_columns) || any(!is.finite(as.matrix(mes)))) {
    stop("WGCNA `MEs` must contain finite numeric values.", call. = FALSE)
  }
  as.matrix(mes)
}

.coexpr_trait_frame <- function(traits, samples) {
  .assert_ids(samples, "Module eigengene sample names")
  if (is.atomic(traits) && is.null(dim(traits))) {
    .assert_ids(names(traits), "`traits` names")
    traits <- data.frame(
      trait = traits,
      row.names = names(traits),
      check.names = FALSE
    )
  } else if (is.data.frame(traits) || is.matrix(traits)) {
    traits <- as.data.frame(traits, optional = TRUE)
    if (!ncol(traits)) {
      stop("`traits` must contain at least one column.", call. = FALSE)
    }
    .assert_ids(rownames(traits), "`traits` row names")
    .assert_ids(names(traits), "`traits` column names")
  } else {
    stop(
      "`traits` must be a named vector or a data frame/matrix with row names.",
      call. = FALSE
    )
  }
  missing <- setdiff(samples, rownames(traits))
  if (length(missing)) {
    stop(
      "`traits` is missing eigengene samples: ",
      paste(utils::head(missing, 10L), collapse = ", "),
      ".",
      call. = FALSE
    )
  }
  traits <- traits[samples, , drop = FALSE]
  converted <- lapply(names(traits), function(name) {
    column <- traits[[name]]
    if (is.logical(column)) {
      column <- as.integer(column)
    } else if (is.factor(column) || is.character(column)) {
      column <- as.integer(factor(column))
    }
    if (!is.numeric(column) || anyNA(column) || any(!is.finite(column))) {
      stop(
        "Trait `", name, "` must be finite after conversion to numeric.",
        call. = FALSE
      )
    }
    column
  })
  names(converted) <- names(traits)
  as.data.frame(converted, row.names = samples, optional = TRUE)
}

.coexpr_membership_matrix <- function(membership) {
  membership <- as.matrix(membership)
  if (!is.numeric(membership) || !nrow(membership) || !ncol(membership)) {
    stop("`membership` must be a non-empty numeric matrix.", call. = FALSE)
  }
  .assert_ids(rownames(membership), "Membership feature names")
  .assert_ids(colnames(membership), "Membership module names")
  if (any(!is.finite(membership))) {
    stop("`membership` must contain finite values.", call. = FALSE)
  }
  membership
}

.coexpr_match_module_column <- function(module, columns) {
  module <- as.character(module)
  if (!length(module) || is.na(module) || !nzchar(module)) {
    stop("`module` must be one non-empty module label.", call. = FALSE)
  }
  stripped <- sub("^ME", "", columns)
  candidates <- unique(c(module, paste0("ME", module), sub("^ME", "", module)))
  hits <- columns[columns %in% candidates | stripped %in% candidates]
  hits <- unique(hits)
  if (length(hits) == 1L) {
    return(hits)
  }
  if (length(hits) > 1L) {
    stop(
      "Module `", module, "` matches multiple membership columns: ",
      paste(hits, collapse = ", "),
      ".",
      call. = FALSE
    )
  }
  stop(
    "No membership column matches module `", module, "`.",
    call. = FALSE
  )
}

#' Extract sample or feature classes from an NMF fit
#'
#' @param fit A native result returned by [cluster_nmf()].
#' @param what Return sample or feature classes.
#'
#' @return A named class vector from the NMF backend.
#' @export
cluster_nmf_classes <- function(
    fit,
    what = c("samples", "features")
) {
  .require_backend("NMF", "to extract NMF classes")
  what <- match.arg(what)
  classes <- .with_attached_namespaces(
    "NMF",
    NMF::predict(fit, what = what)
  )
  if (
    is.null(names(classes)) || anyNA(names(classes)) ||
      any(!nzchar(names(classes))) || anyDuplicated(names(classes))
  ) {
    stop("NMF classes must be named by unique samples or features.", call. = FALSE)
  }
  classes
}

.de_table_limma <- function(result, coef) {
  coefficients <- as.matrix(result$coefficients)
  if (is.null(coefficients) || !nrow(coefficients) || !ncol(coefficients)) {
    stop("The MArrayLM object has no coefficient matrix.", call. = FALSE)
  }
  if (is.null(coef)) {
    if (ncol(coefficients) != 1L) {
      stop(
        "`coef` is required for an MArrayLM with multiple coefficients: ",
        paste(colnames(coefficients), collapse = ", "),
        ".",
        call. = FALSE
      )
    }
    coef <- 1L
  }
  if (is.character(coef)) {
    if (length(coef) != 1L || is.na(coef) || !coef %in% colnames(coefficients)) {
      stop("Unknown limma coefficient in `coef`.", call. = FALSE)
    }
    coefficient_index <- match(coef, colnames(coefficients))
  } else if (
    !is.numeric(coef) || length(coef) != 1L || is.na(coef) ||
      !is.finite(coef) || coef != trunc(coef) || coef < 1L ||
      coef > ncol(coefficients)
  ) {
    stop("`coef` must name or index one limma coefficient.", call. = FALSE)
  } else {
    coefficient_index <- as.integer(coef)
  }
  feature_ids <- rownames(coefficients)
  .assert_ids(feature_ids, "MArrayLM feature names")
  effect <- coefficients[, coefficient_index]
  statistic <- if (!is.null(result$t)) {
    as.matrix(result$t)[, coefficient_index]
  } else {
    rep(NA_real_, length(effect))
  }
  p_value <- if (!is.null(result$p.value)) {
    as.matrix(result$p.value)[, coefficient_index]
  } else {
    rep(NA_real_, length(effect))
  }
  standard_error <- rep(NA_real_, length(effect))
  if (!is.null(result$stdev.unscaled) && !is.null(result$s2.post)) {
    standard_error <- as.matrix(result$stdev.unscaled)[, coefficient_index] *
      sqrt(result$s2.post)
  } else {
    valid <- is.finite(effect) & is.finite(statistic) & statistic != 0
    standard_error[valid] <- abs(effect[valid] / statistic[valid])
  }
  data.frame(
    feature_id = feature_ids,
    effect = as.numeric(effect),
    standard_error = as.numeric(standard_error),
    statistic = as.numeric(statistic),
    p_value = as.numeric(p_value),
    adjusted_p_value = stats::p.adjust(p_value, method = "BH"),
    row.names = feature_ids,
    check.names = FALSE
  )
}

.de_table_edger <- function(result) {
  table <- as.data.frame(result$table, optional = TRUE)
  feature_ids <- rownames(table)
  .assert_ids(feature_ids, "edgeR feature names")
  .de_reject_multicoefficient_edger(table)
  effect_name <- .de_find_column(table, NULL, c("logFC"), required = TRUE, "effect")
  test_name <- .de_find_column(
    table,
    NULL,
    c("t", "F", "LR"),
    required = FALSE,
    "statistic"
  )
  effect <- as.numeric(table[[effect_name]])
  statistic <- rep(NA_real_, nrow(table))
  if (!is.null(test_name)) {
    raw <- as.numeric(table[[test_name]])
    if (test_name %in% c("F", "LR") && any(raw < 0, na.rm = TRUE)) {
      stop("edgeR F/LR statistics cannot be negative.", call. = FALSE)
    }
    statistic <- if (test_name %in% c("F", "LR")) {
      sign(effect) * sqrt(pmax(raw, 0))
    } else {
      raw
    }
  }
  p_name <- .de_find_column(
    table,
    NULL,
    c("PValue", "P.Value", "pvalue"),
    required = FALSE,
    "p-value"
  )
  p_value <- if (is.null(p_name)) rep(NA_real_, nrow(table)) else table[[p_name]]
  standard_error <- rep(NA_real_, nrow(table))
  data.frame(
    feature_id = feature_ids,
    effect = effect,
    standard_error = standard_error,
    statistic = statistic,
    p_value = as.numeric(p_value),
    adjusted_p_value = stats::p.adjust(p_value, method = "BH"),
    row.names = feature_ids,
    check.names = FALSE
  )
}

.de_reject_multicoefficient_edger <- function(table) {
  if (!"logFC" %in% names(table)) {
    coefficient_effects <- grep("^logFC\\.", names(table), value = TRUE)
    if (length(coefficient_effects) > 1L) {
      stop(
        "Multi-coefficient edgeR tests do not have one effect direction. ",
        "Run a one-degree-of-freedom test or construct a separately defined ",
        "effect table.",
        call. = FALSE
      )
    }
  }
  invisible(table)
}

.de_deseq_statistic_type <- function(result) {
  columns <- S4Vectors::mcols(result)
  if (!"stat" %in% rownames(columns) || !"description" %in% names(columns)) {
    return(NULL)
  }
  description <- as.character(columns["stat", "description"])
  if (length(description) == 1L && !is.na(description) &&
      grepl("^LRT statistic", description)) {
    return("nondirectional_lrt")
  }
  NULL
}

.de_standardize_table <- function(
    table,
    feature_column,
    effect_column,
    standard_error_column,
    statistic_column,
    p_value_column,
    adjusted_p_value_column
) {
  if (!nrow(table)) stop("`result` contains no features.", call. = FALSE)
  if (identical(names(table), c(
    "feature_id", "effect", "standard_error", "statistic", "p_value",
    "adjusted_p_value"
  ))) {
    table$feature_id <- as.character(table$feature_id)
    .assert_ids(table$feature_id, "Differential-result feature identifiers")
    for (column in setdiff(names(table), "feature_id")) {
      if (!is.numeric(table[[column]])) {
        stop("`", column, "` must be numeric.", call. = FALSE)
      }
      invalid <- !is.na(table[[column]]) & !is.finite(table[[column]])
      if (any(invalid)) {
        stop("`", column, "` contains infinite values.", call. = FALSE)
      }
    }
    for (column in c("p_value", "adjusted_p_value")) {
      invalid <- !is.na(table[[column]]) &
        (table[[column]] < 0 | table[[column]] > 1)
      if (any(invalid)) {
        stop("`", column, "` must lie between zero and one.", call. = FALSE)
      }
    }
    invalid_se <- !is.na(table$standard_error) & table$standard_error <= 0
    if (any(invalid_se)) {
      stop("`standard_error` must be strictly positive when present.", call. = FALSE)
    }
    rownames(table) <- table$feature_id
    return(table)
  }
  feature_id <- if (is.null(feature_column)) {
    candidate <- intersect(c("feature_id", "feature", "gene_id", "gene"), names(table))
    if (length(candidate)) as.character(table[[candidate[[1L]]]]) else rownames(table)
  } else {
    .de_assert_column(table, feature_column, "feature_column")
    as.character(table[[feature_column]])
  }
  .assert_ids(feature_id, "Differential-result feature identifiers")

  effect_column <- .de_find_column(
    table,
    effect_column,
    c("effect", "log2FoldChange", "logFC", "coef"),
    required = TRUE,
    "effect"
  )
  standard_error_column <- .de_find_column(
    table,
    standard_error_column,
    c("standard_error", "lfcSE", "SE", "Std. Error"),
    required = FALSE,
    "standard error"
  )
  statistic_column <- .de_find_column(
    table,
    statistic_column,
    c("statistic", "stat", "t", "WaldStatistic", "F", "LR"),
    required = FALSE,
    "statistic"
  )
  p_value_column <- .de_find_column(
    table,
    p_value_column,
    c("p_value", "pvalue", "PValue", "P.Value", "p.value"),
    required = FALSE,
    "p-value"
  )
  adjusted_p_value_column <- .de_find_column(
    table,
    adjusted_p_value_column,
    c("adjusted_p_value", "padj", "FDR", "adj.P.Val", "adj.p.value"),
    required = FALSE,
    "adjusted p-value"
  )

  effect <- .de_numeric_column(table, effect_column, "effect")
  standard_error <- .de_optional_numeric_column(table, standard_error_column)
  statistic <- .de_optional_numeric_column(table, statistic_column)
  statistic_type <- NULL
  if (
    !is.null(statistic_column) && statistic_column %in% c("F", "LR") &&
      !all(is.na(statistic))
  ) {
    if (any(statistic < 0, na.rm = TRUE)) {
      stop("F/LR statistics cannot be negative.", call. = FALSE)
    }
    statistic_type <- "nondirectional_f_or_lr"
  }
  p_value <- .de_optional_numeric_column(table, p_value_column)
  adjusted_p_value <- .de_optional_numeric_column(table, adjusted_p_value_column)
  if (all(is.na(adjusted_p_value)) && !all(is.na(p_value))) {
    adjusted_p_value <- stats::p.adjust(p_value, method = "BH")
  }
  result <- data.frame(
    feature_id = feature_id,
    effect = effect,
    standard_error = standard_error,
    statistic = statistic,
    p_value = p_value,
    adjusted_p_value = adjusted_p_value,
    row.names = feature_id,
    check.names = FALSE
  )
  invalid_se <- !is.na(result$standard_error) & result$standard_error <= 0
  if (any(invalid_se)) {
    stop("`standard_error` must be strictly positive when present.", call. = FALSE)
  }
  for (column in setdiff(names(result), "feature_id")) {
    value <- result[[column]]
    invalid <- !is.na(value) & !is.finite(value)
    if (any(invalid)) {
      stop("`", column, "` contains infinite values.", call. = FALSE)
    }
  }
  for (column in c("p_value", "adjusted_p_value")) {
    invalid <- !is.na(result[[column]]) &
      (result[[column]] < 0 | result[[column]] > 1)
    if (any(invalid)) {
      stop("`", column, "` must lie between zero and one.", call. = FALSE)
    }
  }
  if (!is.null(statistic_type)) {
    attr(result, "statistic_type") <- statistic_type
  }
  result
}

.de_find_column <- function(table, supplied, candidates, required, label) {
  if (!is.null(supplied)) {
    .de_assert_column(table, supplied, paste0(label, " column"))
    return(supplied)
  }
  found <- intersect(candidates, names(table))
  if (length(found)) return(found[[1L]])
  if (required) {
    stop(
      "Could not identify the ", label, " column; supply it explicitly.",
      call. = FALSE
    )
  }
  NULL
}

.de_assert_column <- function(table, column, argument) {
  if (
    !is.character(column) || length(column) != 1L || is.na(column) ||
      !nzchar(column) || !column %in% names(table)
  ) {
    stop("`", argument, "` must name one result column.", call. = FALSE)
  }
  invisible(column)
}

.de_numeric_column <- function(table, column, label) {
  value <- table[[column]]
  if (!is.numeric(value)) {
    stop("The selected ", label, " column must be numeric.", call. = FALSE)
  }
  as.numeric(value)
}

.de_optional_numeric_column <- function(table, column) {
  if (is.null(column)) return(rep(NA_real_, nrow(table)))
  .de_numeric_column(table, column, column)
}
