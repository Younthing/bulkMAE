#' Fit a DESeq2 differential-expression model
#'
#' @inheritParams mae_pull_assay
#' @param design DESeq2 design formula.
#' @param test Wald or likelihood-ratio test.
#' @param reduced Reduced formula required for a likelihood-ratio test.
#' @param ... Additional arguments passed to [DESeq2::DESeq()].
#'
#' @return A fitted native `DESeqDataSet`.
#' @export
de_deseq2 <- function(
    x,
    experiment,
    design,
    assay = "counts",
    test = c("Wald", "LRT"),
    reduced = NULL,
    ...
) {
  .require_backend("DESeq2", "to fit a differential-expression model")
  test <- match.arg(test)
  if (test == "LRT" && is.null(reduced)) {
    stop("`reduced` is required when `test = \"LRT\"`.", call. = FALSE)
  }

  dds <- .make_deseq_dataset(x, experiment, assay, design)
  arguments <- list(object = dds, test = test, ...)
  if (!is.null(reduced)) {
    arguments$reduced <- reduced
  }
  do.call(DESeq2::DESeq, arguments)
}

#' Extract DESeq2 results
#'
#' @param fit A fitted `DESeqDataSet`, usually from [de_deseq2()].
#' @param contrast Optional DESeq2 contrast.
#' @param name Optional coefficient name. Do not supply together with
#'   `contrast`.
#' @param alpha Independent-filtering false discovery rate threshold.
#' @param independent_filtering Use DESeq2 independent filtering.
#' @param ... Additional arguments passed to [DESeq2::results()].
#'
#' @return A native `DESeqResults` object.
#' @export
de_deseq2_results <- function(
    fit,
    contrast = NULL,
    name = NULL,
    alpha = 0.05,
    independent_filtering = TRUE,
    ...
) {
  .require_backend("DESeq2", "to extract DESeq2 results")
  if (!is.null(contrast) && !is.null(name)) {
    stop("Supply only one of `contrast` and `name`.", call. = FALSE)
  }

  arguments <- list(
    object = fit,
    alpha = alpha,
    independentFiltering = independent_filtering,
    ...
  )
  if (!is.null(contrast)) arguments$contrast <- contrast
  if (!is.null(name)) arguments$name <- name
  do.call(DESeq2::results, arguments)
}

#' Fit an edgeR quasi-likelihood model
#'
#' @inheritParams mae_pull_assay
#' @param formula Fixed-effect model formula.
#' @param filter Remove low-expression features with
#'   [edgeR::filterByExpr()].
#' @param normalize_method Library normalization method.
#' @param robust Use robust empirical Bayes dispersion estimation.
#' @param coef Coefficient index or name for [edgeR::glmQLFTest()].
#' @param contrast Optional numeric contrast for [edgeR::glmQLFTest()].
#' @param ... Additional arguments passed to [edgeR::glmQLFit()].
#'
#' @return A native `DGEGLM` fit when neither `coef` nor `contrast` is given;
#'   otherwise a native `DGELRT` test object.
#' @export
de_edger <- function(
    x,
    experiment,
    formula,
    assay = "counts",
    filter = TRUE,
    normalize_method = "TMM",
    robust = TRUE,
    coef = NULL,
    contrast = NULL,
    ...
) {
  .require_backend("edgeR", "to fit a quasi-likelihood model")
  if (!is.null(coef) && !is.null(contrast)) {
    stop("Supply only one of `coef` and `contrast`.", call. = FALSE)
  }

  data <- mae_samples(x, experiment)
  design <- .model_matrix(formula, data)
  y <- .make_dge_list(x, experiment, assay)
  if (filter) {
    y <- y[edgeR::filterByExpr(y, design = design), , keep.lib.sizes = FALSE]
  }
  if (!nrow(y)) {
    stop("No features remain after edgeR expression filtering.", call. = FALSE)
  }
  y <- edgeR::normLibSizes(y, method = normalize_method)
  y <- edgeR::estimateDisp(y, design = design, robust = robust)
  fit <- edgeR::glmQLFit(y, design = design, robust = robust, ...)

  if (is.null(coef) && is.null(contrast)) {
    return(fit)
  }
  edgeR::glmQLFTest(fit, coef = coef, contrast = contrast)
}

#' Fit a limma-voom differential-expression model
#'
#' @inheritParams de_edger
#' @param contrasts Optional contrast matrix passed to
#'   [limma::contrasts.fit()].
#' @param trend Passed to [limma::eBayes()].
#' @param voom_plot Draw the voom mean-variance diagnostic plot.
#' @param ... Additional arguments passed to [limma::eBayes()].
#'
#' @return A native `MArrayLM` object.
#' @export
de_voom <- function(
    x,
    experiment,
    formula,
    assay = "counts",
    filter = TRUE,
    normalize_method = "TMM",
    contrasts = NULL,
    trend = FALSE,
    voom_plot = FALSE,
    ...
) {
  .require_backend("edgeR", "to prepare counts for voom")
  .require_backend("limma", "to fit a voom model")
  data <- mae_samples(x, experiment)
  design <- .model_matrix(formula, data)
  y <- .make_dge_list(x, experiment, assay)
  if (filter) {
    y <- y[edgeR::filterByExpr(y, design = design), , keep.lib.sizes = FALSE]
  }
  if (!nrow(y)) {
    stop("No features remain after voom expression filtering.", call. = FALSE)
  }
  y <- edgeR::normLibSizes(y, method = normalize_method)
  voom <- limma::voom(y, design = design, plot = voom_plot)
  fit <- limma::lmFit(voom, design = design)
  if (!is.null(contrasts)) {
    fit <- limma::contrasts.fit(fit, contrasts = contrasts)
  }
  limma::eBayes(fit, trend = trend, ...)
}

#' Fit a limma model to continuous expression or score data
#'
#' Use this adapter for microarray intensities, log-expression assays, pathway
#' scores, or other approximately continuous outcomes. For counts, use
#' [de_voom()] instead.
#'
#' @inheritParams mae_pull_assay
#' @param formula Fixed-effect model formula.
#' @param contrasts Optional contrast matrix passed to
#'   [limma::contrasts.fit()].
#' @param trend Passed to [limma::eBayes()].
#' @param ... Additional arguments passed to [limma::eBayes()].
#'
#' @return A native `MArrayLM` object.
#' @export
de_limma <- function(
    x,
    experiment,
    formula,
    assay,
    contrasts = NULL,
    trend = FALSE,
    ...
) {
  .require_backend("limma", "to fit a continuous-expression model")
  matrix <- .pull_matrix(x, experiment, assay)
  if (any(!is.finite(matrix))) {
    stop("limma requires a finite continuous-expression matrix.", call. = FALSE)
  }
  design <- .model_matrix(formula, mae_samples(x, experiment))
  fit <- limma::lmFit(matrix, design = design)
  if (!is.null(contrasts)) {
    fit <- limma::contrasts.fit(fit, contrasts = contrasts)
  }
  limma::eBayes(fit, trend = trend, ...)
}

#' Test differential exon or transcript usage from a fitted limma/edgeR model
#'
#' @param fit A native `MArrayLM` or `DGEGLM` model fit.
#' @param gene_id Gene identifier for every fitted feature.
#' @param feature_id Optional exon or transcript identifier.
#' @param robust Use robust empirical Bayes estimation.
#' @param coef,contrast Coefficient or contrast for the edgeR method. Supply at
#'   most one. For limma, apply contrasts to `fit` before calling this function;
#'   limma's `diffSplice()` calculates statistics for the fitted coefficients.
#' @param ... Additional arguments passed to the selected backend's
#'   `diffSplice()` method.
#'
#' @details `gene_id` and `feature_id` are aligned to the fitted rows when
#' they are named. This is particularly useful after expression filtering.
#' Supply an untested `DGEGLM` to the edgeR method; choose its test with `coef`
#' or `contrast`. For limma,
#' apply the desired contrast with [limma::contrasts.fit()] before calling this
#' function when a contrast, rather than a coefficient, is required.
#'
#' @return A native limma or edgeR differential-splicing result.
#' @export
dtu_diffsplice <- function(
    fit,
    gene_id,
    feature_id = NULL,
    robust = FALSE,
    coef = NULL,
    contrast = NULL,
    ...
) {
  if (!inherits(fit, "MArrayLM") && !inherits(fit, "DGEGLM")) {
    stop("`fit` must inherit from MArrayLM or DGEGLM.", call. = FALSE)
  }
  if (!is.logical(robust) || length(robust) != 1L || is.na(robust)) {
    stop("`robust` must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.null(coef) && !is.null(contrast)) {
    stop("Supply only one of `coef` and `contrast`.", call. = FALSE)
  }
  coefficients <- fit$coefficients
  if (is.null(coefficients) || is.null(dim(coefficients))) {
    stop("`fit` must contain a coefficient matrix.", call. = FALSE)
  }
  fitted_features <- rownames(coefficients)
  if (
    !is.null(fitted_features) &&
      (anyNA(fitted_features) || any(!nzchar(fitted_features)) ||
        anyDuplicated(fitted_features))
  ) {
    stop("Fitted feature names must be unique and non-missing.", call. = FALSE)
  }
  gene_id <- .align_splice_id(
    gene_id,
    fitted_features = fitted_features,
    n_features = nrow(coefficients),
    argument = "gene_id"
  )
  if (!is.null(feature_id)) {
    feature_id <- .align_splice_id(
      feature_id,
      fitted_features = fitted_features,
      n_features = nrow(coefficients),
      argument = "feature_id"
    )
  }

  if (inherits(fit, "MArrayLM")) {
    .require_backend("limma", "to test differential exon usage")
    if (!is.null(coef) || !is.null(contrast)) {
      stop(
        "For limma, apply the desired contrast to `fit` before calling ",
        "`dtu_diffsplice()`.",
        call. = FALSE
      )
    }
    return(limma::diffSplice(
      fit,
      geneid = gene_id,
      exonid = feature_id,
      robust = robust,
      ...
    ))
  }
  if (inherits(fit, "DGEGLM")) {
    .require_backend("edgeR", "to test differential exon usage")
    arguments <- list(
      glmfit = fit,
      geneid = gene_id,
      exonid = feature_id,
      robust = robust,
      ...
    )
    if (!is.null(coef)) arguments$coef <- coef
    if (!is.null(contrast)) arguments$contrast <- contrast
    return(do.call(edgeR::diffSplice, arguments))
  }
  stop("No differential-splicing method is available for `fit`.", call. = FALSE)
}

.align_splice_id <- function(value, fitted_features, n_features, argument) {
  if (
    is.data.frame(value) || is.matrix(value) ||
      (!is.atomic(value) && !is.factor(value))
  ) {
    stop("`", argument, "` must be an atomic vector.", call. = FALSE)
  }

  value_names <- names(value)
  if (!is.null(value_names)) {
    if (anyNA(value_names) || any(!nzchar(value_names)) || anyDuplicated(value_names)) {
      stop(
        "Names on `", argument, "` must be unique, non-missing feature IDs.",
        call. = FALSE
      )
    }
    if (!is.null(fitted_features)) {
      missing <- setdiff(fitted_features, value_names)
      if (length(missing)) {
        stop(
          "`", argument, "` is missing fitted features: ",
          paste(utils::head(missing, 10L), collapse = ", "),
          call. = FALSE
        )
      }
      value <- value[fitted_features]
    }
  }

  if (length(value) != n_features) {
    stop(
      "`", argument, "` must contain one identifier for every fitted row ",
      "(", n_features, "). Name it by feature to align after filtering.",
      call. = FALSE
    )
  }
  if (anyNA(value) || any(!nzchar(as.character(value)))) {
    stop("`", argument, "` cannot contain missing or empty IDs.", call. = FALSE)
  }
  unname(value)
}

#' Fit a maSigPro time-course or dose-response model
#'
#' Follows the official maSigPro sequence: `make.design.matrix()`, `p.vector()`,
#' `T.fit()`, and `get.siggenes()`.
#'
#' @inheritParams mae_pull_assay
#' @param edesign Experimental design data frame required by maSigPro. Rows
#'   must be row-named by assay sample ID and include `Time`, `Replicates`, and
#'   at least one binary group-indicator column. Rows are reordered to assay
#'   sample order before fitting.
#' @param degree Polynomial degree for the regression model.
#' @param q,p_adjust Significance threshold and adjustment method for
#'   [maSigPro::p.vector()].
#' @param min_observations Optional minimum observations per gene.
#' @param counts Whether the assay contains counts.
#' @param step_method Variable-selection method used by [maSigPro::T.fit()].
#' @param alpha Significance threshold for variable selection.
#' @param r_squared Minimum model R-squared.
#' @param variables maSigPro variable grouping passed as `vars`.
#' @param ... Additional arguments passed to [maSigPro::p.vector()].
#'
#' @return The native object returned by [maSigPro::get.siggenes()].
#' @export
de_masigpro <- function(
    x,
    experiment,
    edesign,
    assay,
    degree = 2,
    q = 0.05,
    p_adjust = "BH",
    min_observations = NULL,
    counts = FALSE,
    step_method = "backward",
    alpha = 0.05,
    r_squared = 0.6,
    variables = "groups",
    ...
) {
  .require_backend("maSigPro", "to fit a time-course model")
  if (!is.logical(counts) || length(counts) != 1L || is.na(counts)) {
    stop("`counts` must be TRUE or FALSE.", call. = FALSE)
  }
  matrix <- if (counts) {
    .as_count_matrix(x, experiment, assay, integer = TRUE)
  } else {
    .pull_matrix(x, experiment, assay)
  }
  edesign <- as.data.frame(edesign, check.names = FALSE)
  if (nrow(edesign) != ncol(matrix)) {
    stop("`edesign` must have one row per assay sample.", call. = FALSE)
  }
  sample_ids <- colnames(matrix)
  design_ids <- rownames(edesign)
  if (
    is.null(design_ids) || anyNA(design_ids) || any(!nzchar(design_ids)) ||
      anyDuplicated(design_ids)
  ) {
    stop("`edesign` must have unique, non-missing sample row names.", call. = FALSE)
  }
  missing_samples <- setdiff(sample_ids, design_ids)
  extra_samples <- setdiff(design_ids, sample_ids)
  if (length(missing_samples) || length(extra_samples)) {
    stop(
      "`edesign` row names must equal assay sample names. Missing: ",
      .differential_format_ids(missing_samples), "; extra: ",
      .differential_format_ids(extra_samples),
      ".",
      call. = FALSE
    )
  }
  edesign <- edesign[sample_ids, , drop = FALSE]

  required <- c("Time", "Replicates")
  if (!all(required %in% names(edesign))) {
    stop("`edesign` must contain `Time` and `Replicates` columns.", call. = FALSE)
  }
  group_columns <- setdiff(names(edesign), required)
  if (!length(group_columns)) {
    stop("`edesign` must contain at least one group-indicator column.", call. = FALSE)
  }
  if (!is.numeric(edesign$Time) || any(!is.finite(edesign$Time))) {
    stop("`edesign$Time` must be finite and numeric.", call. = FALSE)
  }
  if (anyNA(edesign$Replicates)) {
    stop("`edesign$Replicates` cannot contain missing values.", call. = FALSE)
  }
  invalid_group <- vapply(
    edesign[group_columns],
    function(value) {
      (!is.numeric(value) && !is.logical(value)) ||
        anyNA(value) || !all(value %in% c(0, 1))
    },
    logical(1)
  )
  if (any(invalid_group)) {
    stop(
      "Group-indicator columns must contain only 0/1 values: ",
      paste(names(invalid_group)[invalid_group], collapse = ", "),
      call. = FALSE
    )
  }
  if (any(rowSums(as.matrix(edesign[group_columns])) < 1L)) {
    stop("Every sample must belong to at least one indicated group.", call. = FALSE)
  }
  edesign <- edesign[, c(required, group_columns), drop = FALSE]
  if (
    length(degree) != 1L || is.na(degree) || !is.finite(degree) ||
      degree < 1L || degree != trunc(degree)
  ) {
    stop("`degree` must be one positive integer.", call. = FALSE)
  }
  if (length(unique(edesign$Time)) <= degree) {
    stop("The number of distinct time/dose values must exceed `degree`.", call. = FALSE)
  }
  insufficient_time_groups <- vapply(
    group_columns,
    function(group_column) {
      active <- edesign[[group_column]] == 1
      length(unique(edesign$Time[active])) <= degree
    },
    logical(1)
  )
  if (any(insufficient_time_groups)) {
    stop(
      "Each group needs more distinct time/dose values than `degree`: ",
      paste(group_columns[insufficient_time_groups], collapse = ", "),
      call. = FALSE
    )
  }
  .differential_assert_probability(q, "q")
  .differential_assert_probability(alpha, "alpha")
  .differential_assert_probability(r_squared, "r_squared")
  p_adjust <- match.arg(p_adjust, stats::p.adjust.methods)
  if (
    !is.null(min_observations) &&
      (!is.numeric(min_observations) || length(min_observations) != 1L ||
        is.na(min_observations) || !is.finite(min_observations) ||
        min_observations != trunc(min_observations) || min_observations < 1L ||
        min_observations > ncol(matrix))
  ) {
    stop(
      "`min_observations` must be NULL or a positive integer no larger than ",
      "the sample count.",
      call. = FALSE
    )
  }
  if (any(!is.finite(matrix))) {
    stop("The maSigPro assay must contain only finite values.", call. = FALSE)
  }
  design <- maSigPro::make.design.matrix(edesign, degree = degree)
  design_matrix <- if (is.list(design) && !is.null(design$dis)) {
    as.matrix(design$dis)
  } else {
    as.matrix(design)
  }
  if (nrow(design_matrix) != ncol(matrix)) {
    stop("maSigPro produced a design with the wrong number of samples.", call. = FALSE)
  }
  design_rank <- qr(design_matrix)$rank
  if (design_rank < ncol(design_matrix)) {
    stop("The maSigPro design matrix is rank deficient.", call. = FALSE)
  }
  if (nrow(design_matrix) <= design_rank) {
    stop("The maSigPro design has no residual degrees of freedom.", call. = FALSE)
  }
  arguments <- list(
    data = matrix,
    design = design,
    Q = q,
    MT.adjust = p_adjust,
    counts = counts,
    ...
  )
  if (!is.null(min_observations)) arguments$min.obs <- min_observations
  first_fit <- do.call(maSigPro::p.vector, arguments)
  selected_fit <- maSigPro::T.fit(
    first_fit,
    step.method = step_method,
    alfa = alpha
  )
  maSigPro::get.siggenes(selected_fit, rsq = r_squared, vars = variables)
}

.differential_format_ids <- function(value, n = 10L) {
  if (!length(value)) return("none")
  shown <- paste(utils::head(value, n), collapse = ", ")
  if (length(value) > n) paste0(shown, ", ...") else shown
}

.differential_assert_probability <- function(value, argument) {
  if (
    !is.numeric(value) || length(value) != 1L || is.na(value) ||
      !is.finite(value) || value < 0 || value > 1
  ) {
    stop("`", argument, "` must be one number in [0, 1].", call. = FALSE)
  }
  invisible(value)
}

#' Fit a repeated-measures model with dream
#'
#' @inheritParams de_edger
#' @param formula Formula containing fixed and optional random effects.
#' @param contrast Optional contrast matrix passed as `L` to
#'   [variancePartition::dream()].
#' @param ... Additional arguments passed to [variancePartition::dream()].
#'
#' @return A native, empirical-Bayes moderated `MArrayLM` object.
#' @export
de_dream <- function(
    x,
    experiment,
    formula,
    assay = "counts",
    filter = TRUE,
    normalize_method = "TMM",
    contrast = NULL,
    ...
) {
  .require_backend("edgeR", "to prepare counts for dream")
  .require_backend("lme4", "to identify fixed effects in a dream formula")
  .require_backend("limma", "to moderate dream results")
  .require_backend("variancePartition", "to fit a dream model")
  data <- mae_samples(x, experiment)
  fixed_formula <- lme4::nobars(formula)
  design <- .model_matrix(fixed_formula, data)
  y <- .make_dge_list(x, experiment, assay)
  if (filter) {
    y <- y[edgeR::filterByExpr(y, design = design), , keep.lib.sizes = FALSE]
  }
  y <- edgeR::normLibSizes(y, method = normalize_method)
  voom <- variancePartition::voomWithDreamWeights(y, formula, data)

  arguments <- list(exprObj = voom, formula = formula, data = data, ...)
  if (!is.null(contrast)) arguments$L <- contrast
  fit <- do.call(variancePartition::dream, arguments)
  limma::eBayes(fit)
}
