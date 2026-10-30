#' Plot Kaplan-Meier survival curves
#'
#' Draws the step curves and censoring marks stored in an already-fitted
#' `survfit` object, with the legend inside the panel and an optional
#' annotation of the hazard ratio and log-rank p-value. The constructor does
#' not run a test or fit a model: pass the `survdiff` and `coxph` results the
#' annotation should report. A number-at-risk panel can be added below the
#' curves as a second facet of the same ggplot.
#'
#' @param fit A single-state `survfit` object, usually from [surv_km()].
#' @param test Optional `survival::survdiff()` result for the same strata.
#' @param cox Optional `coxph` fit whose first coefficient is reported as the
#'   hazard ratio, for example a model of the same two-level group.
#' @param hr_label Text printed before the hazard ratio.
#' @param conf_int Draw the pointwise confidence band when the fit has one.
#' @param risk_table Draw a number-at-risk panel below the curves.
#' @param risk_times Optional non-negative times for the number-at-risk panel.
#'   Defaults to pretty breaks across the observed follow-up.
#'
#' @return An unprinted standard ggplot object carrying recommended physical
#'   dimensions for [plot_save()].
#' @family plotting
#' @export
plot_surv_km <- function(
    fit,
    test = NULL,
    cox = NULL,
    hr_label = "HR",
    conf_int = FALSE,
    risk_table = FALSE,
    risk_times = NULL
) {
  .require_backend("survival", "to read Kaplan-Meier fits")
  .plot_assert_flag(conf_int, "conf_int")
  .plot_assert_flag(risk_table, "risk_table")
  .assert_scalar_character(hr_label, "hr_label")
  extracted <- .plot_km_data(fit)
  curves <- extracted$curves
  strata <- extracted$strata
  panels <- c("Survival", "Number at risk")
  curves$panel <- factor(panels[[1L]], levels = panels)
  colours <- .plot_survival_values(strata)
  conf_int <- conf_int && !is.null(fit$upper)

  plot <- ggplot2::ggplot(curves, ggplot2::aes(
    x = .data[["time"]], y = .data[["surv"]], colour = .data[["stratum"]]
  ))
  if (conf_int) {
    plot <- plot + ggplot2::geom_ribbon(
      data = .plot_km_stairs(curves),
      mapping = ggplot2::aes(
        x = .data[["time"]], ymin = .data[["lower"]], ymax = .data[["upper"]],
        fill = .data[["stratum"]]
      ),
      inherit.aes = FALSE, alpha = 0.12, colour = NA, na.rm = TRUE
    )
  }
  plot <- plot +
    ggplot2::geom_step(direction = "hv", linewidth = 0.45) +
    ggplot2::geom_point(
      data = curves[curves$n_censor > 0, , drop = FALSE],
      shape = 3, size = 0.5, stroke = 0.3, show.legend = FALSE
    )
  breaks <- seq(0, 1, by = 0.25)
  labels <- paste0(100 * breaks, "%")
  annotation <- .plot_km_annotation(test, cox, hr_label, strata)
  if (!is.null(annotation)) {
    plot <- plot + ggplot2::geom_text(
      data = data.frame(
        time = 0, surv = 0.03, label = annotation,
        panel = factor(panels[[1L]], levels = panels)
      ),
      mapping = ggplot2::aes(x = .data[["time"]], y = .data[["surv"]], label = .data[["label"]]),
      inherit.aes = FALSE, hjust = 0, vjust = 0, lineheight = 0.9,
      size = .bulkmae_survival_text_pt, size.unit = "pt", colour = .bulkmae_survival[["ink"]]
    )
  }
  if (risk_table) {
    table <- .plot_km_risk_table(fit, extracted, risk_times)
    table$panel <- factor(panels[[2L]], levels = panels)
    plot <- plot +
      ggplot2::geom_text(
        data = table,
        mapping = ggplot2::aes(
          x = .data[["time"]], y = .data[["y"]], label = .data[["n_risk"]],
          colour = .data[["stratum"]]
        ),
        size = .bulkmae_survival_text_pt, size.unit = "pt", show.legend = FALSE
      ) +
      ggplot2::facet_grid(
        rows = ggplot2::vars(.data[["panel"]]), scales = "free_y", space = "free_y"
      )
    rows <- unique(table[c("stratum", "y")])
    breaks <- c(breaks, rows$y)
    labels <- c(labels, as.character(rows$stratum))
  }
  plot <- plot +
    ggplot2::scale_colour_manual(values = colours, breaks = strata, drop = FALSE, name = NULL) +
    ggplot2::scale_fill_manual(values = colours, breaks = strata, drop = FALSE, guide = "none") +
    ggplot2::scale_y_continuous(breaks = breaks, labels = labels,
                                expand = ggplot2::expansion(mult = 0.03, add = c(0.04, 0))) +
    ggplot2::scale_x_continuous(
      expand = ggplot2::expansion(mult = c(if (risk_table) 0.06 else 0.03, 0.03))
    ) +
    ggplot2::coord_cartesian(clip = "off") +
    ggplot2::labs(
      x = "Time", y = "Survival probability",
      alt = paste(
        "Kaplan-Meier step curves by stratum with censoring marks",
        if (risk_table) "and a number-at-risk panel below." else "."
      )
    ) +
    .plot_survival_theme() +
    ggplot2::theme(
      strip.text = ggplot2::element_blank(),
      legend.position = "inside",
      legend.position.inside = c(0.98, 0.98),
      legend.justification = c(1, 1),
      panel.spacing.y = grid::unit(0.5, "lines")
    )
  if (risk_table) {
    plot <- plot + ggplot2::theme(axis.ticks.y = ggplot2::element_line(
      colour = .bulkmae_survival[["ink"]], linewidth = 0.3
    ))
  }
  .plot_with_dimensions(
    plot,
    width = 6,
    height = if (risk_table) 6 + 0.35 * length(strata) else 5.5
  )
}

#' Plot hazard ratios from Cox models as a forest plot
#'
#' Draws already-estimated hazard ratios with capped confidence intervals on
#' a log axis. The right-hand axis prints each estimate, interval, and Wald
#' p-value. Estimates whose interval excludes one are coloured by direction.
#' Pass one multivariable `coxph` fit, or a named list of univariate fits.
#'
#' @param fit A `coxph` fit, usually from [surv_cox()], or a named list of
#'   `coxph` fits.
#' @param terms Optional coefficient names (or, for a list, element names) to
#'   draw, in display order from top to bottom.
#' @param term_labels Optional term-named character display labels.
#' @param conf_level Confidence level of the drawn intervals.
#'
#' @return An unprinted standard ggplot object carrying recommended physical
#'   dimensions for [plot_save()].
#' @family plotting
#' @export
plot_surv_forest <- function(
    fit,
    terms = NULL,
    term_labels = NULL,
    conf_level = 0.95
) {
  .require_backend("survival", "to read Cox model fits")
  .differential_assert_probability(conf_level, "conf_level")
  data <- .plot_cox_table(fit, conf_level)
  if (!is.null(terms)) {
    .plot_assert_ids(as.character(terms), "`terms`")
    unknown <- setdiff(terms, data$term_id)
    if (length(unknown)) {
      stop("Unknown `terms`: ", paste(unknown, collapse = ", "), ".", call. = FALSE)
    }
    data <- data[match(terms, data$term_id), , drop = FALSE]
  }
  for (column in c("hazard_ratio", "lower", "upper")) {
    if (any(!is.finite(data[[column]])) || any(data[[column]] <= 0)) {
      stop("Hazard ratios and confidence limits must be finite and positive.", call. = FALSE)
    }
  }
  labels <- data$term_id
  if (!is.null(term_labels)) {
    data$term_label <- data$term_id
    labels <- .plot_enrichment_labels(data, term_labels)
  }
  data$term <- factor(make.unique(labels), levels = rev(make.unique(labels)))
  data$effect <- factor(
    ifelse(data$lower > 1, "Hazard", ifelse(data$upper < 1, "Protective", "Not significant")),
    levels = c("Protective", "Not significant", "Hazard")
  )
  table_labels <- sprintf(
    paste0("%.2f (%.2f", intToUtf8(8211L), "%.2f)  %s"),
    data$hazard_ratio, data$lower, data$upper, .plot_format_p(data$p_value)
  )
  names(table_labels) <- as.character(data$term)

  plot <- ggplot2::ggplot(data, ggplot2::aes(
    x = .data[["hazard_ratio"]], y = .data[["term"]], colour = .data[["effect"]]
  )) +
    ggplot2::geom_vline(xintercept = 1, colour = .bulkmae_survival[["neutral"]],
                        linetype = "dashed", linewidth = 0.3) +
    ggplot2::geom_errorbar(
      mapping = ggplot2::aes(xmin = .data[["lower"]], xmax = .data[["upper"]]),
      width = 0.25, linewidth = 0.4, orientation = "y"
    ) +
    ggplot2::geom_point(size = 1.4) +
    ggplot2::scale_colour_manual(
      values = c(Protective = .bulkmae_survival[["low"]],
                 `Not significant` = .bulkmae_survival[["neutral"]],
                 Hazard = .bulkmae_survival[["high"]]),
      drop = FALSE, guide = "none"
    ) +
    ggplot2::scale_x_log10(breaks = .plot_hr_breaks(c(data$lower, data$upper))) +
    ggplot2::scale_y_discrete(
      sec.axis = ggplot2::dup_axis(
        labels = unname(table_labels[levels(data$term)]),
        name = NULL
      )
    ) +
    ggplot2::labs(
      x = sprintf("Hazard ratio (%g%% CI)", 100 * conf_level),
      y = NULL,
      alt = paste(
        "Forest plot of Cox hazard ratios with capped confidence intervals on",
        "a log axis; estimates and p-values are printed on the right."
      )
    ) +
    .plot_survival_theme() +
    ggplot2::theme(
      axis.line.y.right = ggplot2::element_blank(),
      axis.ticks.y.right = ggplot2::element_blank(),
      axis.text.y.right = ggplot2::element_text(hjust = 0)
    )
  .plot_with_dimensions(
    plot,
    width = 10,
    height = .plot_clamped_dimension(nrow(data), base = 1.4, per_item = 0.4,
                                     minimum = 3.5, maximum = 16)
  )
}

#' Plot the cross-validation curve of a glmnet model
#'
#' Shows the cross-validated error (mean and one standard error) along the
#' penalty path stored in a `cv.glmnet` object, with dashed lines at
#' `lambda.min` and `lambda.1se`. The top axis gives the number of non-zero
#' coefficients. Works for [surv_penalized()] and [ml_glmnet()] results.
#'
#' @param fit A `cv.glmnet` object.
#'
#' @return An unprinted standard ggplot object carrying recommended physical
#'   dimensions for [plot_save()].
#' @family plotting
#' @export
plot_glmnet_cv <- function(fit) {
  .plot_assert_cv_glmnet(fit)
  data <- data.frame(
    log_lambda = log(fit$lambda), cvm = fit$cvm, cvlo = fit$cvlo,
    cvup = fit$cvup, nzero = as.integer(fit$nzero)
  )
  .plot_assert_finite_numeric(data$log_lambda, "`fit$lambda` (log)")
  .plot_assert_finite_numeric(data$cvm, "`fit$cvm`")
  index <- unique(round(seq(1L, nrow(data), length.out = min(7L, nrow(data)))))
  chosen <- data.frame(
    log_lambda = log(c(fit$lambda.min, fit$lambda.1se)),
    rule = c("lambda.min", "lambda.1se")
  )
  plot <- ggplot2::ggplot(data, ggplot2::aes(x = .data[["log_lambda"]], y = .data[["cvm"]])) +
    ggplot2::geom_vline(
      data = chosen, mapping = ggplot2::aes(xintercept = .data[["log_lambda"]]),
      colour = .bulkmae_survival[["neutral"]], linetype = "dashed", linewidth = 0.3
    ) +
    ggplot2::geom_errorbar(
      mapping = ggplot2::aes(ymin = .data[["cvlo"]], ymax = .data[["cvup"]]),
      colour = .bulkmae_survival[["miss"]], linewidth = 0.3, width = 0
    ) +
    ggplot2::geom_point(colour = .bulkmae_survival[["high"]], size = 0.9) +
    ggplot2::scale_x_continuous(
      sec.axis = ggplot2::sec_axis(
        transform = ~ ., breaks = data$log_lambda[index], labels = data$nzero[index],
        name = "Non-zero coefficients"
      )
    ) +
    ggplot2::labs(
      x = expression(log(lambda)), y = unname(as.character(fit$name)),
      alt = paste(
        "Cross-validated error with one-standard-error bars along the glmnet",
        "penalty path, with dashed lines at lambda.min and lambda.1se."
      )
    ) +
    .plot_survival_theme()
  .plot_with_dimensions(plot, width = 6.5, height = 5.5)
}

#' Plot glmnet coefficient paths
#'
#' Draws each feature's coefficient across the penalty path of a `cv.glmnet`
#' fit. Features that are non-zero at the chosen `lambda` are coloured by the
#' sign of that coefficient; all other paths are light grey. Text labels are
#' added only for features named in `label_features`.
#'
#' @param fit A single-response `cv.glmnet` object.
#' @param lambda Penalty rule marked on the plot and used to colour paths.
#' @param label_features Optional feature names to label at the left edge.
#'
#' @return An unprinted standard ggplot object carrying recommended physical
#'   dimensions for [plot_save()].
#' @family plotting
#' @export
plot_glmnet_path <- function(
    fit,
    lambda = c("lambda.min", "lambda.1se"),
    label_features = NULL
) {
  .plot_assert_cv_glmnet(fit)
  lambda <- match.arg(lambda)
  beta <- fit$glmnet.fit$beta
  if (is.list(beta)) {
    stop("Coefficient paths are only supported for single-response glmnet fits.", call. = FALSE)
  }
  .require_backend("Matrix", "to read sparse glmnet coefficients")
  beta <- as.matrix(beta)
  features <- rownames(beta)
  .plot_assert_ids(features, "glmnet feature names")
  path_lambda <- fit$glmnet.fit$lambda
  chosen <- which.min(abs(path_lambda - fit[[lambda]]))
  at_chosen <- beta[, chosen]
  status <- ifelse(at_chosen > 0, "Positive", ifelse(at_chosen < 0, "Negative", "Zero"))
  data <- data.frame(
    feature = rep(features, times = ncol(beta)),
    log_lambda = rep(log(path_lambda), each = nrow(beta)),
    coefficient = as.vector(beta),
    status = factor(rep(status, times = ncol(beta)),
                    levels = c("Negative", "Zero", "Positive"))
  )
  data <- data[order(data$status == "Zero", decreasing = TRUE), , drop = FALSE]
  palette <- c(Negative = .bulkmae_survival[["low"]], Zero = .bulkmae_survival[["miss"]],
               Positive = .bulkmae_survival[["high"]])
  plot <- ggplot2::ggplot(data, ggplot2::aes(
    x = .data[["log_lambda"]], y = .data[["coefficient"]],
    group = .data[["feature"]], colour = .data[["status"]]
  )) +
    ggplot2::geom_hline(yintercept = 0, colour = .bulkmae_survival[["neutral"]], linewidth = 0.25) +
    ggplot2::geom_vline(xintercept = log(path_lambda[[chosen]]),
                        colour = .bulkmae_survival[["neutral"]], linetype = "dashed", linewidth = 0.3) +
    ggplot2::geom_line(linewidth = 0.35) +
    ggplot2::scale_colour_manual(values = palette, breaks = names(palette), drop = FALSE,
                                 name = paste("At", lambda)) +
    ggplot2::labs(
      x = expression(log(lambda)), y = "Coefficient",
      alt = paste(
        "Coefficient paths along the glmnet penalty, coloured by the sign of",
        "each coefficient at the marked lambda."
      )
    ) +
    .plot_survival_theme() +
    ggplot2::theme(
      legend.position = "inside",
      legend.position.inside = c(0.98, 0.98),
      legend.justification = c(1, 1)
    )
  if (!is.null(label_features)) {
    .plot_assert_ids(as.character(label_features), "`label_features`")
    unknown <- setdiff(label_features, features)
    if (length(unknown)) {
      stop("Unknown `label_features`: ", paste(unknown, collapse = ", "), ".", call. = FALSE)
    }
    first <- which.min(path_lambda)
    labelled <- data.frame(
      feature = label_features,
      log_lambda = log(path_lambda[[first]]),
      coefficient = beta[label_features, first],
      status = factor(status[match(label_features, features)], levels = levels(data$status))
    )
    plot <- plot +
      ggplot2::geom_text(
        data = labelled, mapping = ggplot2::aes(label = .data[["feature"]]),
        hjust = 1, nudge_x = -0.05, check_overlap = TRUE, show.legend = FALSE,
        size = .bulkmae_survival_text_pt, size.unit = "pt"
      ) +
      ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0.2, 0.03)))
  }
  .plot_with_dimensions(plot, width = 6.5, height = 5.5)
}

#' Plot a risk-score distribution with follow-up status
#'
#' Ranks samples by an already-computed risk score and draws two aligned
#' panels: the sorted score split at `cutoff`, and each sample's follow-up
#' time coloured by event status. Samples are joined to `sample_data` by name.
#'
#' @param score A sample-named numeric risk score, for example from
#'   [score_signature()] or `predict()` on a Cox fit.
#' @param sample_data Sample metadata whose row names match `score` exactly.
#' @param time,event Columns in `sample_data` with follow-up time and a
#'   zero/one event indicator.
#' @param cutoff One finite score threshold. Samples above it are `High`.
#'   Fix it before looking at outcomes, for example the training-set median.
#' @param time_label Strip label for the follow-up panel.
#'
#' @return An unprinted standard ggplot object carrying recommended physical
#'   dimensions for [plot_save()].
#' @family plotting
#' @export
plot_surv_risk <- function(
    score,
    sample_data,
    time,
    event,
    cutoff,
    time_label = "Follow-up time"
) {
  if (!is.numeric(score) || is.null(names(score))) {
    stop("`score` must be a sample-named numeric vector.", call. = FALSE)
  }
  .plot_assert_ids(names(score), "Names of `score`")
  .plot_assert_finite_numeric(score, "`score`")
  .plot_assert_finite_numeric(cutoff, "`cutoff`")
  if (length(cutoff) != 1L) stop("`cutoff` must be one finite number.", call. = FALSE)
  .assert_scalar_character(time_label, "time_label")
  sample_data <- as.data.frame(sample_data, optional = TRUE)
  .plot_align_frame(sample_data, names(score), "sample_data", "`score` names")
  .plot_require_columns(sample_data, c(time, event), "sample_data")
  times <- sample_data[names(score), time, drop = TRUE]
  events <- sample_data[names(score), event, drop = TRUE]
  .plot_assert_finite_numeric(times, paste0("`sample_data$", time, "`"))
  if (any(times < 0)) stop("Follow-up time must be non-negative.", call. = FALSE)
  if (!is.numeric(events) || anyNA(events) || !all(events %in% c(0, 1))) {
    stop("`sample_data$", event, "` must be coded 0 (censored) or 1 (event).", call. = FALSE)
  }
  ordered <- order(score)
  panels <- c("Risk score", time_label)
  group <- factor(ifelse(score > cutoff, "High", "Low"), levels = c("Low", "High"))
  status <- factor(ifelse(events == 1, "Event", "Censored"), levels = c("Censored", "Event"))
  rank <- integer(length(score))
  rank[ordered] <- seq_along(score)
  top <- data.frame(rank = rank, value = unname(score), group = group,
                    panel = factor(panels[[1L]], levels = panels))
  bottom <- data.frame(rank = rank, value = times, status = status,
                       panel = factor(panels[[2L]], levels = panels))
  bottom <- bottom[order(bottom$status), , drop = FALSE]
  boundary <- sum(group == "Low") + 0.5
  plot <- ggplot2::ggplot(mapping = ggplot2::aes(x = .data[["rank"]], y = .data[["value"]])) +
    ggplot2::geom_vline(xintercept = boundary, colour = .bulkmae_survival[["neutral"]],
                        linetype = "dashed", linewidth = 0.3) +
    ggplot2::geom_hline(
      data = data.frame(value = cutoff, panel = factor(panels[[1L]], levels = panels)),
      mapping = ggplot2::aes(yintercept = .data[["value"]]),
      colour = .bulkmae_survival[["neutral"]], linetype = "dashed", linewidth = 0.3
    ) +
    ggplot2::geom_point(data = top, mapping = ggplot2::aes(colour = .data[["group"]]),
                        size = 0.6) +
    ggplot2::geom_point(data = bottom, mapping = ggplot2::aes(fill = .data[["status"]]),
                        shape = 21, colour = "transparent", stroke = 0, size = 0.9) +
    ggplot2::facet_grid(rows = ggplot2::vars(.data[["panel"]]), scales = "free_y",
                        switch = "y") +
    ggplot2::scale_colour_manual(
      values = c(Low = .bulkmae_survival[["low"]], High = .bulkmae_survival[["high"]]),
      drop = FALSE, name = "Risk"
    ) +
    ggplot2::scale_fill_manual(
      values = c(Censored = .bulkmae_survival[["miss"]], Event = .bulkmae_survival[["ink"]]),
      drop = FALSE, name = "Status"
    ) +
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = 0.01)) +
    ggplot2::guides(
      colour = ggplot2::guide_legend(override.aes = list(size = 1.4)),
      fill = ggplot2::guide_legend(override.aes = list(size = 1.4))
    ) +
    ggplot2::labs(
      x = "Patients ranked by risk score", y = NULL,
      alt = paste(
        "Two aligned panels: samples sorted by risk score and split at the",
        "cutoff, and each sample's follow-up time coloured by event status."
      )
    ) +
    .plot_survival_theme() +
    ggplot2::theme(
      strip.placement = "outside",
      strip.text.y.left = ggplot2::element_text(angle = 90, face = "plain"),
      legend.position = "right",
      legend.box.just = "left",
      panel.spacing.y = grid::unit(0.6, "lines")
    )
  .plot_with_dimensions(plot, width = 8, height = 6.5)
}

#' Plot time-dependent ROC curves
#'
#' Draws the true- and false-positive rates stored in a `timeROC` result,
#' one curve per evaluation time, with each time's AUC in the legend.
#'
#' @param roc A `timeROC` result, usually from [surv_roc()].
#' @param time_labels Optional character labels, one per evaluation time, such
#'   as `c("3 years", "5 years")`.
#'
#' @return An unprinted standard ggplot object carrying recommended physical
#'   dimensions for [plot_save()].
#' @family plotting
#' @export
plot_surv_roc <- function(roc, time_labels = NULL) {
  parts <- .plot_timeroc_parts(roc)
  time_labels <- .plot_time_labels(parts$times, time_labels)
  keep <- is.finite(parts$auc)
  if (!any(keep)) stop("No evaluation time has a finite AUC.", call. = FALSE)
  legend <- sprintf("%s  AUC = %.3f", time_labels, parts$auc)
  data <- do.call(rbind, lapply(which(keep), function(index) {
    curve <- data.frame(fp = parts$fp[, index], tp = parts$tp[, index])
    curve <- curve[is.finite(curve$fp) & is.finite(curve$tp), , drop = FALSE]
    curve <- curve[order(curve$fp, curve$tp), , drop = FALSE]
    curve$time <- legend[[index]]
    curve
  }))
  data$time <- factor(data$time, levels = legend[keep])
  plot <- ggplot2::ggplot(data, ggplot2::aes(
    x = .data[["fp"]], y = .data[["tp"]], colour = .data[["time"]]
  )) +
    ggplot2::geom_abline(slope = 1, intercept = 0, colour = .bulkmae_survival[["neutral"]],
                         linetype = "dashed", linewidth = 0.3) +
    ggplot2::geom_path(linewidth = 0.5) +
    ggplot2::scale_colour_manual(values = .plot_survival_values(levels(data$time)), name = NULL) +
    ggplot2::scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.25),
                                labels = .plot_unit_labels,
                                expand = ggplot2::expansion(mult = 0.01)) +
    ggplot2::scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.25),
                                labels = .plot_unit_labels,
                                expand = ggplot2::expansion(mult = 0.01)) +
    ggplot2::coord_equal() +
    ggplot2::labs(
      x = "1 - Specificity", y = "Sensitivity",
      alt = "Time-dependent ROC curves, one per evaluation time, with AUC in the legend."
    ) +
    .plot_survival_theme() +
    ggplot2::theme(
      legend.position = "inside",
      legend.position.inside = c(0.98, 0.02),
      legend.justification = c(1, 0)
    )
  .plot_with_dimensions(plot, width = 5.5, height = 5.5)
}

#' Compare time-dependent AUC across predictors
#'
#' Draws one line per predictor across the evaluation times of several
#' `timeROC` results, for example a risk score, a clinical stage, and a
#' combined model evaluated in the same cohort.
#'
#' @param rocs A named list of `timeROC` results that share evaluation times.
#'   Names are used as predictor labels, in legend order.
#' @param time_labels Optional character labels, one per evaluation time.
#'
#' @return An unprinted standard ggplot object carrying recommended physical
#'   dimensions for [plot_save()].
#' @family plotting
#' @export
plot_surv_auc <- function(rocs, time_labels = NULL) {
  if (!is.list(rocs) || !length(rocs) || inherits(rocs, "ipcwsurvivalROC")) {
    stop("`rocs` must be a named list of timeROC results.", call. = FALSE)
  }
  .plot_assert_ids(names(rocs), "Names of `rocs`")
  parts <- lapply(rocs, .plot_timeroc_parts)
  times <- parts[[1L]]$times
  if (!all(vapply(parts, function(part) isTRUE(all.equal(part$times, times)), logical(1)))) {
    stop("All `rocs` must use the same evaluation times.", call. = FALSE)
  }
  time_labels <- .plot_time_labels(times, time_labels)
  data <- do.call(rbind, Map(function(part, name) {
    data.frame(predictor = name, time = time_labels, auc = unname(part$auc))
  }, parts, names(rocs)))
  data <- data[is.finite(data$auc), , drop = FALSE]
  data$predictor <- factor(data$predictor, levels = names(rocs))
  data$time <- factor(data$time, levels = time_labels)
  plot <- ggplot2::ggplot(data, ggplot2::aes(
    x = .data[["time"]], y = .data[["auc"]],
    colour = .data[["predictor"]], group = .data[["predictor"]]
  )) +
    ggplot2::geom_hline(yintercept = 0.5, colour = .bulkmae_survival[["neutral"]],
                        linetype = "dashed", linewidth = 0.3) +
    ggplot2::geom_line(linewidth = 0.45) +
    ggplot2::geom_point(size = 1.2) +
    ggplot2::scale_colour_manual(values = .plot_survival_values(names(rocs)), name = NULL) +
    ggplot2::labs(
      x = NULL, y = "Time-dependent AUC",
      alt = "Time-dependent AUC of several predictors across evaluation times."
    ) +
    .plot_survival_theme() +
    ggplot2::theme(legend.position = "bottom")
  .plot_with_dimensions(plot, width = 5.5, height = 5)
}

#' Plot survival calibration curves
#'
#' Draws mean predicted survival against Kaplan-Meier observed survival with
#' its confidence interval for each risk group and time, as returned by
#' [surv_calibration()]. Points on the dashed diagonal are well calibrated.
#'
#' @param calibration A data frame from [surv_calibration()].
#' @param time_labels Optional character labels, one per distinct `time`.
#'
#' @return An unprinted standard ggplot object carrying recommended physical
#'   dimensions for [plot_save()].
#' @family plotting
#' @export
plot_surv_calibration <- function(calibration, time_labels = NULL) {
  data <- as.data.frame(calibration, optional = TRUE)
  .plot_require_columns(data, c("time", "predicted", "observed", "lower", "upper"), "calibration")
  .plot_assert_finite_numeric(data$predicted, "`calibration$predicted`")
  .plot_assert_finite_numeric(data$time, "`calibration$time`")
  times <- sort(unique(data$time))
  time_labels <- .plot_time_labels(times, time_labels)
  data$time_label <- factor(time_labels[match(data$time, times)], levels = time_labels)
  data <- data[order(data$time, data$predicted), , drop = FALSE]
  limits <- range(c(data$predicted, data$lower, data$upper, 1), na.rm = TRUE)
  limits[[1L]] <- max(0, floor(limits[[1L]] * 10) / 10)
  plot <- ggplot2::ggplot(data, ggplot2::aes(
    x = .data[["predicted"]], y = .data[["observed"]], colour = .data[["time_label"]]
  )) +
    ggplot2::geom_abline(slope = 1, intercept = 0, colour = .bulkmae_survival[["neutral"]],
                         linetype = "dashed", linewidth = 0.3) +
    ggplot2::geom_errorbar(
      mapping = ggplot2::aes(ymin = .data[["lower"]], ymax = .data[["upper"]]),
      width = 0.015, linewidth = 0.3, na.rm = TRUE
    ) +
    ggplot2::geom_line(linewidth = 0.45, na.rm = TRUE) +
    ggplot2::geom_point(size = 1.2, na.rm = TRUE) +
    ggplot2::scale_colour_manual(values = .plot_survival_values(time_labels), name = NULL) +
    ggplot2::coord_equal(xlim = limits, ylim = limits) +
    ggplot2::labs(
      x = "Predicted survival", y = "Observed survival (KM)",
      alt = paste(
        "Calibration curves of mean predicted against Kaplan-Meier observed",
        "survival by risk group, one line per time, with a dashed diagonal."
      )
    ) +
    .plot_survival_theme() +
    ggplot2::theme(
      legend.position = "inside",
      legend.position.inside = c(0.02, 0.98),
      legend.justification = c(0, 1)
    )
  .plot_with_dimensions(plot, width = 5.5, height = 5.5)
}

#' Plot a decision curve
#'
#' Draws net benefit against threshold probability for each model in a
#' `dcurves` decision-curve result, together with the treat-all and
#' treat-none reference strategies.
#'
#' @param dca A `dca` object, usually from [surv_dca()].
#'
#' @return An unprinted standard ggplot object carrying recommended physical
#'   dimensions for [plot_save()].
#' @family plotting
#' @export
plot_surv_dca <- function(dca) {
  if (!inherits(dca, "dca") || is.null(dca$dca)) {
    stop("`dca` must be a dcurves decision-curve result.", call. = FALSE)
  }
  data <- as.data.frame(dca$dca, optional = TRUE)
  .plot_require_columns(data, c("variable", "label", "threshold", "net_benefit"), "dca$dca")
  data$label <- as.character(data$label)
  references <- unique(data$label[data$variable %in% c("all", "none")])
  models <- unique(data$label[!data$variable %in% c("all", "none")])
  levels <- c(models, references)
  palette <- .plot_survival_values(models)
  palette[data$label[match("all", data$variable)]] <- .bulkmae_survival[["miss"]]
  palette[data$label[match("none", data$variable)]] <- .bulkmae_survival[["ink"]]
  data$label <- factor(data$label, levels = levels)
  model_rows <- !data$variable %in% c("all", "none")
  top <- max(data$net_benefit[is.finite(data$net_benefit)])
  plot <- ggplot2::ggplot(data, ggplot2::aes(
    x = .data[["threshold"]], y = .data[["net_benefit"]], colour = .data[["label"]]
  )) +
    ggplot2::geom_line(linewidth = 0.5, na.rm = TRUE) +
    ggplot2::scale_colour_manual(values = palette[levels], breaks = levels, name = NULL) +
    ggplot2::coord_cartesian(ylim = c(-0.15 * top, 1.05 * top)) +
    ggplot2::labs(
      x = "Threshold probability",
      y = "Net benefit",
      alt = "Decision curves of net benefit against threshold probability for each model and reference strategy."
    ) +
    .plot_survival_theme() +
    ggplot2::theme(
      legend.position = "inside",
      legend.position.inside = c(0.98, 0.98),
      legend.justification = c(1, 1)
    )
  if (!any(model_rows)) stop("`dca` contains no model curves.", call. = FALSE)
  .plot_with_dimensions(plot, width = 6.5, height = 5.5)
}

.bulkmae_survival <- c(
  high = "#B64C3C", low = "#3D6B8C", accent = "#C4A35A", neutral = "#6E6A66",
  ink = "#2C2926", miss = "#D5D1CB", green = "#6E8B74", purple = "#8D6A9F",
  brown = "#5C4A3A"
)

.bulkmae_survival_text_pt <- 6

.plot_survival_values <- function(levels) {
  levels <- as.character(levels)
  palette <- unname(.bulkmae_survival[c("high", "low", "accent", "green", "purple",
                                        "brown", "neutral")])
  if (identical(sort(levels), c("High", "Low"))) {
    return(c(High = .bulkmae_survival[["high"]], Low = .bulkmae_survival[["low"]]))
  }
  if (length(levels) > length(palette)) return(palette)
  stats::setNames(palette[seq_along(levels)], levels)
}

.plot_survival_theme <- function(base_size = 7) {
  ink <- .bulkmae_survival[["ink"]]
  ggplot2::theme_classic(base_size = base_size) +
    ggplot2::theme(
      axis.line = ggplot2::element_line(linewidth = 0.3, colour = ink),
      axis.ticks = ggplot2::element_line(linewidth = 0.3, colour = ink),
      axis.text = ggplot2::element_text(size = base_size, colour = ink),
      axis.title = ggplot2::element_text(size = base_size + 0.5, colour = ink),
      plot.title = ggplot2::element_text(size = base_size + 1, face = "bold", hjust = 0),
      plot.title.position = "plot",
      strip.background = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(size = base_size, face = "bold", colour = ink),
      legend.title = ggplot2::element_text(size = base_size, colour = ink),
      legend.text = ggplot2::element_text(size = base_size - 0.5, colour = ink),
      legend.key.size = grid::unit(3, "mm"),
      legend.background = ggplot2::element_blank(),
      plot.background = ggplot2::element_rect(fill = "white", colour = NA),
      plot.margin = ggplot2::margin(3, 3, 3, 3)
    )
}

.plot_unit_labels <- function(breaks) {
  format(breaks, drop0trailing = TRUE, trim = TRUE)
}

.plot_hr_breaks <- function(values) {
  limits <- range(c(values, 1))
  # Prefer doubling steps; fall back to finer steps only for narrow ranges.
  for (candidates in list(
    c(0.01, 0.1, 1, 10, 100),
    c(0.125, 0.25, 0.5, 1, 2, 4, 8),
    c(0.5, 0.7, 1, 1.5, 2),
    c(0.8, 0.9, 1, 1.1, 1.25)
  )) {
    breaks <- candidates[candidates >= limits[[1L]] & candidates <= limits[[2L]]]
    if (length(breaks) >= 3L) return(breaks)
  }
  breaks
}

.plot_format_p <- function(p) {
  ifelse(
    is.na(p), "P = NA",
    ifelse(p < 0.001, sprintf("P = %.1e", p), sprintf("P = %.3f", p))
  )
}

.plot_km_annotation <- function(test, cox, hr_label, strata) {
  lines <- character()
  if (!is.null(cox)) {
    if (!inherits(cox, "coxph") || !length(stats::coef(cox))) {
      stop("`cox` must be a coxph fit with at least one coefficient.", call. = FALSE)
    }
    estimate <- stats::coef(cox)[[1L]]
    se <- sqrt(as.matrix(stats::vcov(cox))[1L, 1L])
    interval <- exp(estimate + c(-1, 1) * stats::qnorm(0.975) * se)
    lines <- c(lines, sprintf(
      paste0("%s = %.2f (%.2f", intToUtf8(8211L), "%.2f)"),
      hr_label, exp(estimate), interval[[1L]], interval[[2L]]
    ))
  }
  if (!is.null(test)) {
    lines <- c(lines, paste("Log-rank", .plot_format_p(.plot_survdiff_p(test, strata))))
  }
  if (!length(lines)) return(NULL)
  paste(lines, collapse = "\n")
}

.plot_timeroc_parts <- function(roc) {
  if (!inherits(roc, "ipcwsurvivalROC") && !inherits(roc, "ipcwcompetingrisksROC")) {
    stop("`roc` must be a timeROC result.", call. = FALSE)
  }
  fp <- if (!is.null(roc$FP)) roc$FP else roc$FP_1
  auc <- if (!is.null(roc$AUC)) roc$AUC else roc$AUC_1
  if (is.null(fp) || is.null(roc$TP) || is.null(auc)) {
    stop("`roc` does not contain TP, FP, and AUC values.", call. = FALSE)
  }
  list(fp = as.matrix(fp), tp = as.matrix(roc$TP), auc = auc, times = roc$times)
}

.plot_time_labels <- function(times, time_labels) {
  if (is.null(time_labels)) return(paste("t =", format(times, trim = TRUE)))
  if (!is.character(time_labels) || length(time_labels) != length(times) ||
      anyNA(time_labels) || any(!nzchar(time_labels)) || anyDuplicated(time_labels)) {
    stop("`time_labels` must give one unique label per evaluation time.", call. = FALSE)
  }
  time_labels
}

.plot_km_data <- function(fit) {
  if (!inherits(fit, "survfit") || inherits(fit, "survfitms") ||
      !is.null(dim(fit$surv))) {
    stop("`fit` must be a single-state survfit object.", call. = FALSE)
  }
  n <- length(fit$time)
  if (!n) stop("`fit` contains no event times.", call. = FALSE)
  if (is.null(fit$strata)) {
    raw <- rep("All", n)
    strata <- "All"
  } else {
    raw <- rep(names(fit$strata), fit$strata)
    strata <- names(fit$strata)
    single <- !grepl(",", strata, fixed = TRUE) & grepl("=", strata, fixed = TRUE)
    if (all(single) && length(unique(sub("=.*$", "", strata))) == 1L) {
      strata <- sub("^[^=]*=", "", strata)
      raw <- sub("^[^=]*=", "", raw)
    }
  }
  if (anyDuplicated(strata)) stop("Survival strata labels must be unique.", call. = FALSE)
  upper <- if (is.null(fit$upper)) rep(NA_real_, n) else fit$upper
  lower <- if (is.null(fit$lower)) rep(NA_real_, n) else fit$lower
  rows <- data.frame(
    stratum = raw, time = fit$time, surv = fit$surv,
    upper = upper, lower = lower, n_censor = fit$n.censor
  )
  start <- data.frame(stratum = strata, time = 0, surv = 1, upper = 1, lower = 1, n_censor = 0)
  curves <- rbind(start, rows)
  curves$stratum <- factor(curves$stratum, levels = strata)
  curves <- curves[order(curves$stratum, curves$time), , drop = FALSE]
  list(curves = curves, strata = strata)
}

.plot_km_stairs <- function(curves) {
  pieces <- lapply(split(curves, curves$stratum), function(stratum) {
    n <- nrow(stratum)
    if (n < 2L) return(stratum)
    # Each band level holds from its own time until the next step time.
    held <- rep(seq_len(n - 1L), each = 2L)
    data.frame(
      stratum = stratum$stratum[[1L]],
      time = c(as.vector(rbind(stratum$time[-n], stratum$time[-1L])), stratum$time[[n]]),
      upper = c(stratum$upper[held], stratum$upper[[n]]),
      lower = c(stratum$lower[held], stratum$lower[[n]]),
      panel = stratum$panel[[1L]]
    )
  })
  do.call(rbind, pieces)
}

.plot_km_risk_table <- function(fit, extracted, risk_times) {
  max_time <- max(fit$time)
  if (is.null(risk_times)) {
    risk_times <- pretty(c(0, max_time), n = 5L)
    risk_times <- risk_times[risk_times <= max_time]
  }
  .plot_assert_finite_numeric(risk_times, "`risk_times`")
  if (any(risk_times < 0) || anyDuplicated(risk_times)) {
    stop("`risk_times` must be unique and non-negative.", call. = FALSE)
  }
  risk_times <- sort(risk_times)
  summary <- summary(fit, times = risk_times, extend = TRUE)
  strata <- extracted$strata
  stratum <- if (is.null(summary$strata)) {
    rep(strata, length(summary$time))
  } else {
    strata[as.integer(summary$strata)]
  }
  position <- match(stratum, strata)
  data.frame(
    stratum = factor(stratum, levels = strata),
    time = summary$time,
    n_risk = as.character(summary$n.risk),
    y = -0.02 - 0.09 * position
  )
}

.plot_survdiff_p <- function(test, strata) {
  if (!inherits(test, "survdiff")) {
    stop("`test` must be a survival::survdiff() result.", call. = FALSE)
  }
  if (length(test$n) != length(strata)) {
    stop("`test` and `fit` describe different numbers of strata.", call. = FALSE)
  }
  if (!is.null(test$pvalue)) return(test$pvalue)
  stats::pchisq(test$chisq, df = length(test$n) - 1L, lower.tail = FALSE)
}

.plot_cox_table <- function(fit, conf_level) {
  one <- function(model, prefix = NULL) {
    if (!inherits(model, "coxph")) {
      stop("`fit` must be a coxph fit or a list of coxph fits.", call. = FALSE)
    }
    coefficients <- stats::coef(model)
    if (!length(coefficients) || anyNA(coefficients)) {
      stop("Cox fits must have estimable coefficients.", call. = FALSE)
    }
    se <- sqrt(diag(as.matrix(stats::vcov(model))))
    z <- stats::qnorm(1 - (1 - conf_level) / 2)
    ids <- names(coefficients)
    if (!is.null(prefix)) {
      ids <- if (length(ids) == 1L) prefix else paste(prefix, ids, sep = ": ")
    }
    data.frame(
      term_id = ids,
      hazard_ratio = exp(unname(coefficients)),
      lower = exp(unname(coefficients) - z * se),
      upper = exp(unname(coefficients) + z * se),
      p_value = 2 * stats::pnorm(-abs(unname(coefficients) / se))
    )
  }
  if (inherits(fit, "coxph")) {
    table <- one(fit)
  } else if (is.list(fit) && length(fit)) {
    .plot_assert_ids(names(fit), "Names of the `fit` list")
    table <- do.call(rbind, Map(one, fit, names(fit)))
  } else {
    stop("`fit` must be a coxph fit or a named list of coxph fits.", call. = FALSE)
  }
  .plot_assert_ids(table$term_id, "Cox term identifiers")
  rownames(table) <- NULL
  table
}

.plot_assert_cv_glmnet <- function(fit) {
  if (!inherits(fit, "cv.glmnet") || is.null(fit$glmnet.fit)) {
    stop("`fit` must be a cv.glmnet object.", call. = FALSE)
  }
  invisible(fit)
}
