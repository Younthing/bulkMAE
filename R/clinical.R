#' Score a gene-expression signature
#'
#' @inheritParams mae_pull_assay
#' @param weights Named numeric feature weights, or an unweighted character
#'   vector of feature names.
#' @param center,scale Either a logical value, or a named numeric vector of
#'   training-set feature centers/scales. `TRUE` estimates values from the
#'   current cohort and is therefore unsuitable for locked external
#'   validation; supply frozen training values instead.
#' @param na_rm Remove missing values when summing scores.
#'
#' @return A named numeric vector with one score per sample.
#' @export
score_signature <- function(
    x,
    experiment,
    weights,
    assay,
    center = FALSE,
    scale = FALSE,
    na_rm = FALSE
) {
  matrix <- .pull_matrix(x, experiment, assay)
  if (!is.logical(na_rm) || length(na_rm) != 1L || is.na(na_rm)) {
    stop("`na_rm` must be `TRUE` or `FALSE`.", call. = FALSE)
  }
  if (is.character(weights)) {
    features <- unique(weights)
    if (!length(features) || anyNA(features) || any(!nzchar(features))) {
      stop("`weights` must contain non-missing feature identifiers.", call. = FALSE)
    }
    weights <- stats::setNames(rep(1 / length(features), length(features)), features)
  } else {
    .assert_named_numeric(weights, "weights")
    features <- names(weights)
  }
  matrix <- .select_features(matrix, features)
  if (any(is.infinite(matrix)) || (!na_rm && anyNA(matrix))) {
    stop(
      "Signature features must be finite; set `na_rm = TRUE` only to ignore NA values.",
      call. = FALSE
    )
  }
  matrix <- .standardize_signature(matrix, center, scale, na_rm)
  score <- colSums(matrix * weights[rownames(matrix)], na.rm = na_rm)
  if (any(!is.finite(score))) {
    stop("Signature scoring produced non-finite values.", call. = FALSE)
  }
  stats::setNames(score, colnames(matrix))
}

#' Fit Kaplan-Meier survival curves
#'
#' @inheritParams mae_samples
#' @param formula Survival formula evaluated against aligned sample metadata.
#' @param ... Additional arguments passed to [survival::survfit()].
#'
#' @return A native `survfit` object.
#' @export
surv_km <- function(x, experiment, formula, ...) {
  .require_backend("survival", "to fit Kaplan-Meier curves")
  data <- mae_samples(x, experiment)
  .validate_survival_formula(formula, data)
  survival::survfit(formula = formula, data = data, ...)
}

#' Construct a namespace-safe survival formula
#'
#' This convenience helper lets users name MAE metadata columns without first
#' attaching or calling the `survival` package. More complex time-dependent or
#' stratified formulas can still be supplied directly to [surv_km()] and
#' [surv_cox()].
#'
#' @param time,event Metadata column names containing follow-up time and event.
#' @param predictors Optional right-hand-side variable names. These must be
#'   sample-metadata columns for [surv_km()]. For [surv_cox()], expression
#'   features are also available when the same names are supplied through its
#'   `features` argument.
#'
#' @return A formula with a `survival::Surv()` response.
#' @export
surv_formula <- function(time, event, predictors = NULL) {
  .assert_scalar_character(time, "time")
  .assert_scalar_character(event, "event")
  if (identical(time, event)) {
    stop("`time` and `event` must name different columns.", call. = FALSE)
  }
  if (!is.null(predictors)) {
    if (
      !is.character(predictors) || !length(predictors) || anyNA(predictors) ||
        any(!nzchar(predictors)) || anyDuplicated(predictors)
    ) {
      stop("`predictors` must contain unique, non-empty names.", call. = FALSE)
    }
    if (any(predictors %in% c(time, event))) {
      stop("Predictors must differ from `time` and `event`.", call. = FALSE)
    }
    if (any(predictors == ".")) {
      stop(
        "`.` is a formula expansion operator, not a supported predictor name; ",
        "supply explicit metadata columns or a formula directly.",
        call. = FALSE
      )
    }
  }
  survival_call <- call(
    "::",
    as.name("survival"),
    as.name("Surv")
  )
  response <- as.call(c(list(survival_call), as.name(time), as.name(event)))
  right_hand_side <- if (is.null(predictors)) {
    1
  } else {
    Reduce(
      function(left, right) call("+", left, as.name(right)),
      predictors[-1L],
      init = as.name(predictors[[1L]])
    )
  }
  stats::as.formula(call("~", response, right_hand_side), env = parent.frame())
}

#' Calculate time-dependent ROC curves
#'
#' @inheritParams mae_samples
#' @param time,event Metadata columns containing follow-up time and event.
#' @param marker Metadata column or numeric marker vector.
#' @param times Evaluation times.
#' @param cause Event code treated as the outcome of interest.
#' @param ... Additional arguments passed to [timeROC::timeROC()].
#'
#' @return A native `timeROC` object.
#' @export
surv_roc <- function(
    x,
    experiment,
    time,
    event,
    marker,
    times,
    cause = 1,
    ...
) {
  .require_backend("timeROC", "to calculate time-dependent ROC curves")
  data <- mae_samples(x, experiment)
  .assert_metadata_column(data, time, "time")
  .assert_metadata_column(data, event, "event")
  if (identical(time, event)) {
    stop("`time` and `event` must name different metadata columns.", call. = FALSE)
  }
  marker <- .column_or_vector(marker, data, "marker")
  .assert_finite_numeric(data[[time]], "time")
  .assert_finite_numeric(data[[event]], "event")
  .assert_finite_numeric(marker, "marker")
  .assert_finite_numeric(times, "times")
  if (any(data[[time]] < 0) || any(times <= 0)) {
    stop("Follow-up `time` must be non-negative and `times` must be positive.", call. = FALSE)
  }
  if (length(cause) != 1L || is.na(cause) || !cause %in% data[[event]]) {
    stop("`cause` must be one observed event code.", call. = FALSE)
  }
  if (!0 %in% data[[event]]) {
    stop("`event` must contain zero-coded censored observations.", call. = FALSE)
  }
  if (!any(data[[event]] == cause) || !any(data[[event]] == 0)) {
    stop("Time-dependent ROC needs both events of interest and censoring.", call. = FALSE)
  }
  if (any(times > max(data[[time]]))) {
    warning("Some evaluation `times` exceed the observed follow-up range.", call. = FALSE)
  }
  backend <- getExportedValue("timeROC", "timeROC")
  backend_environment <- new.env(parent = environment(backend))
  backend_environment$Surv <- survival::Surv
  environment(backend) <- backend_environment
  do.call(
    backend,
    c(
      list(
        T = data[[time]],
        delta = data[[event]],
        marker = marker,
        cause = cause,
        times = times
      ),
      list(...)
    )
  )
}

#' Fit a Cox proportional-hazards model
#'
#' @inheritParams mae_pull_assay
#' @param formula Cox model formula. Clinical variables are read from MAE
#'   sample metadata.
#' @param features Optional expression features appended to the model data.
#' @param ... Additional arguments passed to [survival::coxph()].
#'
#' @return A native `coxph` fit.
#' @export
surv_cox <- function(x, experiment, formula, assay = NULL, features = NULL, ...) {
  .require_backend("survival", "to fit a Cox model")
  data <- mae_samples(x, experiment)
  if (!is.null(features)) {
    if (is.null(assay)) {
      stop("`assay` is required when `features` are supplied.", call. = FALSE)
    }
    matrix <- .select_features(.pull_matrix(x, experiment, assay), features)
    expression <- as.data.frame(t(matrix), check.names = FALSE)
    duplicates <- intersect(names(data), names(expression))
    if (length(duplicates)) {
      stop(
        "Feature names duplicate clinical columns: ",
        paste(duplicates, collapse = ", "),
        call. = FALSE
      )
    }
    data <- cbind(data, expression)
  }
  .validate_survival_formula(formula, data)
  survival::coxph(formula = formula, data = data, ...)
}

#' Extract a hazard-ratio table from a Cox model
#'
#' Converts an already-fitted [survival::coxph()] object into a data frame of
#' coefficients, hazard ratios, and confidence intervals. The function does
#' not refit the model or change the native `coxph` return type.
#'
#' @param fit A `coxph` object, usually from [surv_cox()].
#' @param conf_level Confidence level for the hazard-ratio interval.
#'
#' @return A data frame with `term`, `coefficient`, `hazard_ratio`,
#'   `conf_low`, `conf_high`, `statistic`, `p_value`, `n`, and `events`.
#' @export
surv_cox_table <- function(fit, conf_level = 0.95) {
  .require_backend("survival", "to summarise a Cox model")
  if (!inherits(fit, "coxph")) {
    stop("`fit` must be a coxph object.", call. = FALSE)
  }
  .assert_confidence_level(conf_level)
  summarised <- summary(fit, conf.int = conf_level)
  coefficients <- summarised$coefficients
  intervals <- summarised$conf.int
  if (is.null(coefficients) || !nrow(coefficients) || is.null(intervals)) {
    stop("The Cox model has no coefficient table to extract.", call. = FALSE)
  }
  lower_column <- grep("^lower", colnames(intervals), value = TRUE)
  upper_column <- grep("^upper", colnames(intervals), value = TRUE)
  if (length(lower_column) != 1L || length(upper_column) != 1L) {
    stop("The Cox summary does not contain a unique confidence interval.", call. = FALSE)
  }
  table <- data.frame(
    term = rownames(coefficients),
    coefficient = as.numeric(coefficients[, "coef"]),
    hazard_ratio = as.numeric(intervals[, "exp(coef)"]),
    conf_low = as.numeric(intervals[, lower_column]),
    conf_high = as.numeric(intervals[, upper_column]),
    statistic = as.numeric(coefficients[, "z"]),
    p_value = as.numeric(coefficients[, "Pr(>|z|)"]),
    n = as.integer(fit$n),
    events = as.integer(fit$nevent),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  rownames(table) <- NULL
  if (any(!is.finite(c(table$hazard_ratio, table$conf_low, table$conf_high)))) {
    stop("Cox hazard ratios and confidence limits must be finite.", call. = FALSE)
  }
  table
}

#' Fit univariable Cox models for selected predictors
#'
#' Fits one [surv_cox()] model per predictor and stacks the coefficient tables.
#' Clinical predictors are read from aligned sample metadata. Expression
#' features are appended only when they are named in `features` and `assay` is
#' supplied. This is a convenience loop around the existing Cox wrapper, not a
#' new estimator.
#'
#' @inheritParams surv_cox
#' @param time,event Metadata columns containing follow-up time and event.
#' @param predictors Unique metadata column names and/or expression feature
#'   identifiers. Feature names must also appear in `features`.
#' @inheritParams surv_cox_table
#'
#' @return A data frame in the [surv_cox_table()] schema, with one or more
#'   rows per predictor.
#' @export
surv_cox_univariable <- function(
    x,
    experiment,
    time,
    event,
    predictors,
    assay = NULL,
    features = NULL,
    conf_level = 0.95,
    ...
) {
  .require_backend("survival", "to fit univariable Cox models")
  if (
    !is.character(predictors) || !length(predictors) || anyNA(predictors) ||
      any(!nzchar(predictors)) || anyDuplicated(predictors)
  ) {
    stop("`predictors` must contain unique, non-empty names.", call. = FALSE)
  }
  if (any(predictors %in% c(time, event))) {
    stop("Predictors must differ from `time` and `event`.", call. = FALSE)
  }
  if (!is.null(features)) {
    if (
      !is.character(features) || !length(features) || anyNA(features) ||
        any(!nzchar(features)) || anyDuplicated(features)
    ) {
      stop("`features` must contain unique, non-empty names.", call. = FALSE)
    }
  }
  sample_data <- mae_samples(x, experiment)
  sample_predictors <- intersect(predictors, names(sample_data))
  feature_predictors <- if (is.null(features)) {
    character()
  } else {
    intersect(predictors, features)
  }
  unknown <- setdiff(predictors, c(sample_predictors, feature_predictors))
  if (length(unknown)) {
    stop(
      "Unknown univariable predictors: ",
      paste(unknown, collapse = ", "),
      ". Clinical names must exist in sample metadata; expression names must ",
      "also be supplied in `features`.",
      call. = FALSE
    )
  }
  collision <- intersect(sample_predictors, feature_predictors)
  if (length(collision)) {
    stop(
      "Predictors collide with both sample metadata and `features`: ",
      paste(collision, collapse = ", "),
      ".",
      call. = FALSE
    )
  }
  if (length(feature_predictors) && is.null(assay)) {
    stop("`assay` is required when expression `features` are supplied.", call. = FALSE)
  }

  rows <- lapply(predictors, function(predictor) {
    is_feature <- predictor %in% feature_predictors
    fit <- surv_cox(
      x,
      experiment,
      formula = surv_formula(time, event, predictors = predictor),
      assay = if (is_feature) assay else NULL,
      features = if (is_feature) predictor else NULL,
      ...
    )
    table <- surv_cox_table(fit, conf_level = conf_level)
    table$predictor <- predictor
    table
  })
  result <- do.call(rbind, rows)
  rownames(result) <- NULL
  result
}

#' Stratify a named risk score into groups
#'
#' Cuts an already-computed, sample-named score into discrete risk groups.
#' Quantile probabilities in `(0, 1)` are the default; supply score-scale
#' cutpoints to use explicit thresholds. The helper does not fit a survival
#' model.
#'
#' @param score A uniquely named finite numeric vector, one value per sample.
#' @param cuts Quantile probabilities in `(0, 1)`, or explicit finite cutpoints
#'   on the score scale. The default `0.5` is a median split.
#' @param labels Optional group labels. When omitted, two groups are labelled
#'   `Low` and `High`; more groups are labelled `G1`, `G2`, ...
#'
#' @return A named factor aligned to `names(score)`.
#' @export
surv_risk_groups <- function(score, cuts = 0.5, labels = NULL) {
  .assert_named_numeric(score, "score")
  if (!is.numeric(cuts) || !length(cuts) || any(!is.finite(cuts))) {
    stop("`cuts` must contain finite numeric values.", call. = FALSE)
  }
  if (anyDuplicated(cuts)) {
    stop("`cuts` must be unique.", call. = FALSE)
  }
  quantile_cuts <- all(cuts > 0 & cuts < 1)
  if (quantile_cuts) {
    probs <- sort(unique(c(0, cuts, 1)))
    breaks <- unname(stats::quantile(score, probs = probs, names = FALSE, type = 7L))
  } else {
    if (any(cuts > 0 & cuts < 1)) {
      stop(
        "`cuts` must be either quantile probabilities in (0, 1) or explicit ",
        "score-scale cutpoints, not a mixture.",
        call. = FALSE
      )
    }
    breaks <- c(-Inf, sort(cuts), Inf)
  }
  breaks <- .surv_unique_breaks(breaks, quantile_cuts)
  n_groups <- length(breaks) - 1L
  if (is.null(labels)) {
    labels <- if (identical(n_groups, 2L)) {
      c("Low", "High")
    } else {
      paste0("G", seq_len(n_groups))
    }
  }
  if (
    !is.character(labels) || length(labels) != n_groups || anyNA(labels) ||
      any(!nzchar(labels)) || anyDuplicated(labels)
  ) {
    stop(
      "`labels` must contain one unique, non-empty name per risk group.",
      call. = FALSE
    )
  }
  groups <- cut(
    score,
    breaks = breaks,
    labels = labels,
    include.lowest = TRUE,
    right = TRUE
  )
  if (anyNA(groups)) {
    stop("Risk-group cuts left samples unclassified.", call. = FALSE)
  }
  stats::setNames(groups, names(score))
}

#' Fit a cross-validated penalized Cox model
#'
#' @inheritParams mae_pull_assay
#' @param time,event Metadata columns containing follow-up time and event.
#' @param features Optional expression feature subset.
#' @param alpha Elastic-net mixing parameter.
#' @param folds Number of cross-validation folds.
#' @param seed Optional local random seed used while constructing folds.
#' @param ... Additional arguments passed to [glmnet::cv.glmnet()].
#'
#' @return A native `cv.glmnet` object.
#' @export
surv_penalized <- function(
    x,
    experiment,
    time,
    event,
    assay,
    features = NULL,
    alpha = 1,
    folds = 10,
    seed = NULL,
    ...
) {
  .require_backend("survival", "to construct a survival response")
  .require_backend("glmnet", "to fit a penalized Cox model")
  data <- mae_samples(x, experiment)
  .assert_metadata_column(data, time, "time")
  .assert_metadata_column(data, event, "event")
  if (identical(time, event)) {
    stop("`time` and `event` must name different metadata columns.", call. = FALSE)
  }
  predictors <- t(.select_features(.pull_matrix(x, experiment, assay), features))
  .assert_finite_matrix(predictors, "The predictor matrix")
  .assert_finite_numeric(data[[time]], "time")
  .assert_finite_numeric(data[[event]], "event")
  if (any(data[[time]] < 0)) {
    stop("`time` must be non-negative.", call. = FALSE)
  }
  event_count <- .count_survival_events(data[[event]])
  if (event_count < 2L) {
    stop("At least two events are required for penalized Cox modelling.", call. = FALSE)
  }
  .assert_alpha(alpha)
  .assert_folds(folds, nrow(predictors), event_count)
  response <- survival::Surv(data[[time]], data[[event]])
  .with_seed(
    seed,
    glmnet::cv.glmnet(
      x = predictors,
      y = response,
      family = "cox",
      alpha = alpha,
      nfolds = folds,
      ...
    )
  )
}

#' Fit a cross-validated glmnet prediction model
#'
#' @inheritParams surv_penalized
#' @param outcome Metadata column containing the outcome.
#' @param family glmnet model family.
#'
#' @return A native `cv.glmnet` object.
#' @export
ml_glmnet <- function(
    x,
    experiment,
    outcome,
    assay,
    features = NULL,
    family = "binomial",
    alpha = 1,
    folds = 10,
    seed = NULL,
    ...
) {
  .require_backend("glmnet", "to fit a prediction model")
  data <- mae_samples(x, experiment)
  .assert_metadata_column(data, outcome, "outcome")
  predictors <- t(.select_features(.pull_matrix(x, experiment, assay), features))
  response <- data[[outcome]]
  .assert_finite_matrix(predictors, "The predictor matrix")
  if (anyNA(response)) {
    stop("`outcome` must not contain missing values.", call. = FALSE)
  }
  if (!is.character(family) || length(family) != 1L || is.na(family)) {
    stop("`family` must name one glmnet model family.", call. = FALSE)
  }
  supported_families <- c("gaussian", "binomial", "poisson", "multinomial")
  if (!family %in% supported_families) {
    stop(
      "`family` must be one of: ", paste(supported_families, collapse = ", "),
      ". Use `surv_penalized()` for Cox models.",
      call. = FALSE
    )
  }
  if (is.numeric(response)) {
    .assert_finite_numeric(response, "outcome")
  }
  .assert_alpha(alpha)
  fold_limit <- nrow(predictors)
  if (family %in% c("binomial", "multinomial")) {
    class_sizes <- table(response)
    if (length(class_sizes) < 2L) {
      stop("Classification requires at least two outcome classes.", call. = FALSE)
    }
    if (family == "binomial" && length(class_sizes) != 2L) {
      stop("Binomial classification requires exactly two outcome classes.", call. = FALSE)
    }
    fold_limit <- min(class_sizes)
  }
  if (family == "poisson" && (!is.numeric(response) || any(response < 0))) {
    stop("Poisson outcomes must be non-negative numeric values.", call. = FALSE)
  }
  .assert_folds(folds, nrow(predictors), fold_limit)
  .with_seed(
    seed,
    glmnet::cv.glmnet(
      x = predictors,
      y = response,
      family = family,
      alpha = alpha,
      nfolds = folds,
      ...
    )
  )
}

#' Fit a univariate or meta-regression model
#'
#' @param effects Numeric effect estimates, or a data frame containing
#'   `effect` and `standard_error` columns, such as the output of
#'   [meta_collect()].
#' @param standard_errors Numeric standard errors. Leave `NULL` when `effects`
#'   is a data frame.
#' @param moderators Optional moderator matrix or formula accepted by
#'   [metafor::rma.uni()].
#' @param method Meta-analysis estimator.
#' @param ... Additional arguments passed to [metafor::rma.uni()].
#'
#' @return A native `rma.uni` fit.
#' @export
meta_effect <- function(
    effects,
    standard_errors = NULL,
    moderators = NULL,
    method = "REML",
    ...
) {
  .require_backend("metafor", "to fit a meta-analysis model")
  model_data <- NULL
  if (is.data.frame(effects)) {
    if (!is.null(standard_errors)) {
      stop(
        "Do not supply `standard_errors` when `effects` is a data frame.",
        call. = FALSE
      )
    }
    required <- c("effect", "standard_error")
    if (!all(required %in% names(effects))) {
      stop(
        "A data-frame `effects` input must contain `effect` and ",
        "`standard_error` columns.",
        call. = FALSE
      )
    }
    model_data <- effects
    standard_errors <- model_data$standard_error
    effects <- model_data$effect
  }
  .assert_finite_numeric(effects, "effects")
  .assert_finite_numeric(standard_errors, "standard_errors")
  if (length(effects) != length(standard_errors)) {
    stop("`effects` and `standard_errors` must have equal lengths.", call. = FALSE)
  }
  if (any(standard_errors <= 0)) {
    stop("`standard_errors` must be strictly positive.", call. = FALSE)
  }
  if (is.matrix(moderators)) {
    if (nrow(moderators) != length(effects)) {
      stop("`moderators` must have one row per effect estimate.", call. = FALSE)
    }
    .assert_finite_matrix(moderators, "Numeric moderators")
  } else if (is.data.frame(moderators)) {
    if (nrow(moderators) != length(effects)) {
      stop("`moderators` must have one row per effect estimate.", call. = FALSE)
    }
    numeric_columns <- vapply(moderators, is.numeric, logical(1))
    if (any(numeric_columns)) {
      values <- as.matrix(moderators[, numeric_columns, drop = FALSE])
      .assert_finite_matrix(values, "Numeric moderators")
    }
  }
  .assert_scalar_character(method, "method")
  arguments <- list(
    yi = effects,
    sei = standard_errors,
    method = method,
    ...
  )
  if (!is.null(moderators)) arguments$mods <- moderators
  if (!is.null(model_data) && inherits(moderators, "formula")) {
    missing_moderators <- setdiff(all.vars(moderators), names(model_data))
    if (length(missing_moderators)) {
      stop(
        "Moderator formula refers to missing data-frame columns: ",
        paste(missing_moderators, collapse = ", "),
        ".",
        call. = FALSE
      )
    }
    arguments$data <- model_data
  }
  do.call(metafor::rma.uni, arguments)
}

#' Collect one feature's effect estimates across differential analyses
#'
#' @param results A named list of inputs accepted by [de_table()].
#' @param feature One feature identifier present in every result.
#' @param coef Optional coefficient shared by all results, or a list with one
#'   coefficient per result.
#'
#' @return A data frame with `study`, `feature_id`, `effect`, and
#'   `standard_error`, accepted directly by [meta_effect()].
#' @export
meta_collect <- function(results, feature, coef = NULL) {
  if (
    !is.list(results) || !length(results) || is.null(names(results)) ||
      anyNA(names(results)) || any(!nzchar(names(results))) ||
      anyDuplicated(names(results))
  ) {
    stop("`results` must be a uniquely named, non-empty list.", call. = FALSE)
  }
  .assert_scalar_character(feature, "feature")
  coefficients <- if (is.list(coef)) {
    if (length(coef) != length(results)) {
      stop("A list `coef` must contain one value per result.", call. = FALSE)
    }
    coefficient_names <- names(coef)
    if (!is.null(coefficient_names) && any(nzchar(coefficient_names))) {
      if (
        anyNA(coefficient_names) || any(!nzchar(coefficient_names)) ||
          anyDuplicated(coefficient_names) ||
          !setequal(coefficient_names, names(results))
      ) {
        stop(
          "A named `coef` list must have exactly the same unique names as ",
          "`results`.",
          call. = FALSE
        )
      }
      coef[names(results)]
    } else {
      coef
    }
  } else {
    rep(list(coef), length(results))
  }

  rows <- lapply(seq_along(results), function(index) {
    table <- de_table(results[[index]], coef = coefficients[[index]])
    match_index <- match(feature, table$feature_id)
    if (is.na(match_index)) {
      stop(
        "Feature `", feature, "` is absent from result `", names(results)[[index]],
        "`.",
        call. = FALSE
      )
    }
    effect <- table$effect[[match_index]]
    standard_error <- table$standard_error[[match_index]]
    if (
      !is.finite(effect) || !is.finite(standard_error) || standard_error <= 0
    ) {
      stop(
        "Result `", names(results)[[index]], "` lacks a finite positive ",
        "standard error for `", feature, "`.",
        call. = FALSE
      )
    }
    data.frame(
      study = names(results)[[index]],
      feature_id = feature,
      effect = effect,
      standard_error = standard_error,
      stringsAsFactors = FALSE
    )
  })
  result <- do.call(rbind, rows)
  rownames(result) <- NULL
  result
}

#' Prepare an up/down query for LINCS or CMap
#'
#' @param statistics Named numeric differential-expression statistics.
#' @param n Maximum number of genes selected separately from the positive and
#'   negative statistics.
#'
#' @return A list with `upset` and `downset` gene identifiers.
#'
#' @examples
#' statistic <- c(gene1 = 3, gene2 = 2, gene3 = -1, gene4 = -4)
#' drug_query(statistic, n = 2)
#' @export
drug_query <- function(statistics, n = 150) {
  .assert_named_numeric(statistics, "statistics")
  if (
    length(n) != 1L || !is.numeric(n) || !is.finite(n) || n < 1L ||
      n != as.integer(n)
  ) {
    stop("`n` must be a positive integer.", call. = FALSE)
  }
  up <- names(sort(statistics[statistics > 0], decreasing = TRUE))
  down <- names(sort(statistics[statistics < 0], decreasing = FALSE))
  if (!length(up) || !length(down)) {
    stop(
      "LINCS queries require at least one positive and one negative statistic.",
      call. = FALSE
    )
  }
  n <- as.integer(n)
  list(
    upset = utils::head(up, n),
    downset = utils::head(down, n)
  )
}

#' Search a LINCS reference database for connected signatures
#'
#' @param query A list with `upset` and `downset`, usually from
#'   [drug_query()].
#' @param reference_database Path or identifier accepted by
#'   [signatureSearch::qSig()].
#' @param workers Number of workers.
#' @param sort_by LINCS score used to rank results.
#' @param tau Calculate the standardized Tau score.
#' @param annotations Add compound annotations when available.
#' @param ... Additional arguments passed to [signatureSearch::gess_lincs()].
#'
#' @return A native `gessResult` object.
#' @export
drug_lincs <- function(
    query,
    reference_database,
    workers = 1,
    sort_by = "NCS",
    tau = FALSE,
    annotations = TRUE,
    ...
) {
  .require_backend("signatureSearch", "to search a LINCS reference database")
  if (!is.list(query) || !all(c("upset", "downset") %in% names(query))) {
    stop("`query` must contain `upset` and `downset`.", call. = FALSE)
  }
  upset <- .clean_query_set(query$upset, "upset")
  downset <- .clean_query_set(query$downset, "downset")
  if (length(intersect(upset, downset))) {
    stop("LINCS `upset` and `downset` must be disjoint.", call. = FALSE)
  }
  query <- list(upset = upset, downset = downset)
  if (
    !is.numeric(workers) || length(workers) != 1L || !is.finite(workers) ||
      workers < 1L || workers != as.integer(workers)
  ) {
    stop("`workers` must be a positive integer.", call. = FALSE)
  }
  if (
    !is.logical(tau) || length(tau) != 1L || is.na(tau) ||
      !is.logical(annotations) || length(annotations) != 1L ||
      is.na(annotations)
  ) {
    stop("`tau` and `annotations` must be `TRUE` or `FALSE`.", call. = FALSE)
  }
  qsig <- signatureSearch::qSig(
    query = query,
    gess_method = "LINCS",
    refdb = reference_database
  )
  signatureSearch::gess_lincs(
    qSig = qsig,
    sortby = sort_by,
    tau = tau,
    workers = workers,
    addAnnotations = annotations,
    ...
  )
}

.assert_metadata_column <- function(data, column, argument) {
  if (
    !is.character(column) || length(column) != 1L || is.na(column) ||
      !nzchar(column) || !column %in% names(data)
  ) {
    stop(
      "`", argument, "` must name exactly one metadata column.",
      call. = FALSE
    )
  }
  invisible(column)
}

.standardize_signature <- function(matrix, center, scale, na_rm) {
  centers <- .signature_parameter(
    center,
    matrix,
    "center",
    function(value) rowMeans(value, na.rm = na_rm)
  )
  scales <- .signature_parameter(
    scale,
    matrix,
    "scale",
    function(value) apply(value, 1L, stats::sd, na.rm = na_rm)
  )
  if (!is.null(centers)) {
    matrix <- sweep(matrix, 1L, centers, FUN = "-")
  }
  if (!is.null(scales)) {
    if (any(scales <= 0)) {
      stop("Signature scales must be strictly positive.", call. = FALSE)
    }
    matrix <- sweep(matrix, 1L, scales, FUN = "/")
  }
  matrix
}

.signature_parameter <- function(value, matrix, argument, estimate) {
  if (is.logical(value) && length(value) == 1L && !is.na(value)) {
    if (!value) {
      return(NULL)
    }
    warning(
      "`", argument, " = TRUE` estimates feature ", argument,
      " values in the current cohort. Use named values frozen from the ",
      "training cohort for external validation.",
      call. = FALSE
    )
    result <- estimate(matrix)
    names(result) <- rownames(matrix)
  } else {
    .assert_named_numeric(value, argument)
    missing <- setdiff(rownames(matrix), names(value))
    if (length(missing)) {
      stop(
        "`", argument, "` is missing signature features: ",
        paste(utils::head(missing, 10L), collapse = ", "),
        call. = FALSE
      )
    }
    result <- value[rownames(matrix)]
  }
  if (any(!is.finite(result))) {
    stop("Signature `", argument, "` values must be finite.", call. = FALSE)
  }
  result
}

.validate_survival_formula <- function(formula, data) {
  if (!inherits(formula, "formula")) {
    stop("`formula` must be a formula.", call. = FALSE)
  }
  frame <- stats::model.frame(
    formula,
    data = data,
    na.action = stats::na.fail
  )
  if (!inherits(stats::model.response(frame), "Surv")) {
    stop("The left side of `formula` must construct a survival::Surv response.", call. = FALSE)
  }
  for (column in frame) {
    if (is.numeric(column) && any(!is.finite(column))) {
      stop("Numeric survival-model variables must be finite.", call. = FALSE)
    }
  }
  invisible(frame)
}

.assert_finite_numeric <- function(value, argument) {
  if (!is.numeric(value) || !length(value) || any(!is.finite(value))) {
    stop("`", argument, "` must contain finite numeric values.", call. = FALSE)
  }
  invisible(value)
}

.assert_finite_matrix <- function(value, label) {
  if (
    !is.matrix(value) || !is.numeric(value) || !length(value) ||
      any(!is.finite(value))
  ) {
    stop(label, " must contain finite numeric values.", call. = FALSE)
  }
  invisible(value)
}

.assert_confidence_level <- function(conf_level) {
  if (
    !is.numeric(conf_level) || length(conf_level) != 1L || !is.finite(conf_level) ||
      conf_level <= 0 || conf_level >= 1
  ) {
    stop("`conf_level` must be one finite number between zero and one.", call. = FALSE)
  }
  invisible(conf_level)
}

.surv_unique_breaks <- function(breaks, from_quantiles) {
  collapsed <- unique(breaks)
  if (length(collapsed) < 3L) {
    stop(
      if (from_quantiles) {
        "Quantile `cuts` must produce at least two distinct risk groups."
      } else {
        "Score-scale `cuts` must produce at least two distinct risk groups."
      },
      call. = FALSE
    )
  }
  if (from_quantiles && length(collapsed) != length(breaks)) {
    stop(
      "Quantile `cuts` produced tied break-points; choose different probabilities ",
      "or supply explicit score-scale cutpoints.",
      call. = FALSE
    )
  }
  if (is.unsorted(collapsed, strictly = TRUE)) {
    stop("Risk-group break-points must be strictly increasing.", call. = FALSE)
  }
  collapsed
}

.assert_alpha <- function(alpha) {
  if (
    !is.numeric(alpha) || length(alpha) != 1L || !is.finite(alpha) ||
      alpha < 0 || alpha > 1
  ) {
    stop("`alpha` must be one finite number between zero and one.", call. = FALSE)
  }
  invisible(alpha)
}

.assert_folds <- function(folds, n_samples, effective_limit = n_samples) {
  limit <- min(n_samples, effective_limit)
  if (
    !is.numeric(folds) || length(folds) != 1L || !is.finite(folds) ||
      folds != as.integer(folds) || folds < 3L || folds > limit
  ) {
    stop(
      "`folds` must be an integer from 3 through ", limit,
      " for the available samples/events/classes.",
      call. = FALSE
    )
  }
  invisible(as.integer(folds))
}

.count_survival_events <- function(status) {
  observed <- sort(unique(status))
  if (all(observed %in% c(0, 1))) {
    return(sum(status == 1))
  }
  if (all(observed %in% c(1, 2))) {
    return(sum(status == 2))
  }
  stop("Cox `event` must use either 0/1 or 1/2 status coding.", call. = FALSE)
}

.clean_query_set <- function(value, argument) {
  if (!is.atomic(value) || is.list(value)) {
    stop("`", argument, "` must be a vector of gene identifiers.", call. = FALSE)
  }
  value <- unique(as.character(value))
  if (!length(value) || anyNA(value) || any(!nzchar(value))) {
    stop("`", argument, "` must contain non-missing identifiers.", call. = FALSE)
  }
  value
}
