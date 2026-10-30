#' Calibrate a Cox model against Kaplan-Meier survival
#'
#' Predicts each sample's survival at the requested times from an
#' already-fitted Cox model, splits samples into quantile groups of that
#' prediction, and compares the mean prediction in each group with the
#' group's Kaplan-Meier estimate. Use a cohort the model was not fitted on
#' for external calibration.
#'
#' @inheritParams mae_samples
#' @param fit A `coxph` fit, usually from [surv_cox()] called with
#'   `model = TRUE`. Its covariates are read from the sample metadata of `x`.
#' @param time,event Metadata columns containing follow-up time and a zero/one
#'   event indicator.
#' @param times Positive evaluation times on the scale of `time`.
#' @param groups Number of quantile groups of predicted survival.
#'
#' @return A data frame with one row per time and group: `time`, `group`,
#'   `n`, `predicted` (mean predicted survival), `observed` (Kaplan-Meier
#'   survival), and its `lower` and `upper` confidence limits.
#' @export
surv_calibration <- function(
    x,
    experiment,
    fit,
    time,
    event,
    times,
    groups = 4L
) {
  .require_backend("survival", "to predict survival from a Cox model")
  data <- .survival_evaluation_data(x, experiment, list(fit), time, event)
  .assert_finite_numeric(times, "times")
  if (any(times <= 0) || anyDuplicated(times)) {
    stop("`times` must be unique and positive.", call. = FALSE)
  }
  if (!is.numeric(groups) || length(groups) != 1L || groups < 2 ||
      groups != trunc(groups) || groups > nrow(data)) {
    stop("`groups` must be one whole number between 2 and the sample count.", call. = FALSE)
  }
  predicted <- .survival_predict(fit, data, times)
  response <- survival::Surv(data[[time]], data[[event]])
  rows <- lapply(seq_along(times), function(index) {
    value <- predicted[index, ]
    breaks <- unique(stats::quantile(value, seq(0, 1, length.out = groups + 1L)))
    if (length(breaks) < 3L) {
      stop("Predicted survival at time ", times[[index]], " is too uniform to group.",
           call. = FALSE)
    }
    group <- cut(value, breaks, include.lowest = TRUE, labels = FALSE)
    do.call(rbind, lapply(sort(unique(group)), function(level) {
      members <- group == level
      km <- survival::survfit(response[members] ~ 1)
      at <- summary(km, times = times[[index]], extend = TRUE)
      data.frame(
        time = times[[index]], group = level, n = sum(members),
        predicted = mean(value[members]), observed = at$surv,
        lower = at$lower, upper = at$upper
      )
    }))
  })
  result <- do.call(rbind, rows)
  rownames(result) <- NULL
  result
}

#' Decision-curve analysis for Cox models
#'
#' Converts each already-fitted Cox model into a predicted probability of the
#' event by `horizon` and passes those probabilities to [dcurves::dca()].
#'
#' @inheritParams surv_calibration
#' @param models A named list of `coxph` fits, usually from [surv_cox()]
#'   called with `model = TRUE`. Names become the curve labels.
#' @param horizon One positive time at which the event risk is evaluated.
#' @param thresholds Threshold probabilities passed to [dcurves::dca()].
#'
#' @return A native `dca` object.
#' @export
surv_dca <- function(
    x,
    experiment,
    models,
    time,
    event,
    horizon,
    thresholds = seq(0.01, 0.6, by = 0.01)
) {
  .require_backend("survival", "to predict survival from a Cox model")
  .require_backend("dcurves", "to run decision-curve analysis")
  if (!is.list(models) || !length(models) || inherits(models, "coxph")) {
    stop("`models` must be a named list of coxph fits.", call. = FALSE)
  }
  .assert_ids(names(models), "Names of `models`")
  data <- .survival_evaluation_data(x, experiment, models, time, event)
  .assert_finite_numeric(horizon, "horizon")
  if (length(horizon) != 1L || horizon <= 0) {
    stop("`horizon` must be one positive time.", call. = FALSE)
  }
  .assert_finite_numeric(thresholds, "thresholds")
  columns <- paste0(".bulkmae_model_", seq_along(models))
  frame <- data.frame(.bulkmae_time = data[[time]], .bulkmae_event = data[[event]])
  for (index in seq_along(models)) {
    frame[[columns[[index]]]] <- 1 - .survival_predict(models[[index]], data, horizon)[1L, ]
  }
  formula <- stats::as.formula(paste(
    "survival::Surv(.bulkmae_time, .bulkmae_event) ~",
    paste(columns, collapse = " + ")
  ))
  dcurves::dca(
    formula, data = frame, time = horizon, thresholds = thresholds,
    label = stats::setNames(as.list(names(models)), columns)
  )
}

.survival_evaluation_data <- function(x, experiment, fits, time, event) {
  data <- mae_samples(x, experiment)
  .assert_metadata_column(data, time, "time")
  .assert_metadata_column(data, event, "event")
  .assert_finite_numeric(data[[time]], "time")
  if (any(data[[time]] < 0) || !all(data[[event]] %in% c(0, 1))) {
    stop("`time` must be non-negative and `event` coded 0 or 1.", call. = FALSE)
  }
  for (fit in fits) {
    if (!inherits(fit, "coxph")) {
      stop("Survival models must be coxph fits.", call. = FALSE)
    }
    variables <- all.vars(stats::delete.response(stats::terms(fit)))
    missing <- setdiff(variables, names(data))
    if (length(missing)) {
      stop("Sample metadata lack model variables: ", paste(missing, collapse = ", "), ".",
           call. = FALSE)
    }
    if (anyNA(data[variables])) {
      stop("Model variables must not contain missing values; subset complete samples first.",
           call. = FALSE)
    }
  }
  data
}

.survival_predict <- function(fit, data, times) {
  curves <- tryCatch(
    survival::survfit(fit, newdata = data, se.fit = FALSE),
    error = function(error) {
      stop("Could not predict survival from `fit`; refit it with `model = TRUE`. ",
           conditionMessage(error), call. = FALSE)
    }
  )
  at <- summary(curves, times = times, extend = TRUE)$surv
  matrix(at, nrow = length(times), dimnames = list(NULL, rownames(data)))
}
