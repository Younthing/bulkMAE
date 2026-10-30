#' Plot Kaplan-Meier survival curves with a number-at-risk table
#'
#' Draws the step curves, pointwise confidence bands, and censoring marks
#' stored in an already-fitted `survfit` object, with a number-at-risk table
#' underneath. The table is a second panel of the same standard ggplot, not a
#' separate figure. The log-rank p-value is printed only when the caller
#' passes the `survdiff` result it was computed from; the constructor does not
#' run a test.
#'
#' @param fit A single-state `survfit` object, usually from [surv_km()].
#' @param test Optional `survival::survdiff()` result for the same strata.
#' @param risk_times Optional non-negative times at which to count subjects
#'   at risk. Defaults to pretty breaks across the observed follow-up.
#' @param conf_int Draw the pointwise confidence band when the fit has one.
#' @param risk_table Draw the number-at-risk panel.
#'
#' @return An unprinted standard ggplot object carrying recommended physical
#'   dimensions for [plot_save()].
#' @family plotting
#' @export
plot_surv_km <- function(
    fit,
    test = NULL,
    risk_times = NULL,
    conf_int = TRUE,
    risk_table = TRUE
) {
  .require_backend("survival", "to read Kaplan-Meier fits")
  .plot_assert_flag(conf_int, "conf_int")
  .plot_assert_flag(risk_table, "risk_table")
  extracted <- .plot_km_data(fit)
  curves <- extracted$curves
  strata <- extracted$strata
  panels <- c("Survival", "Number at risk")
  curves$panel <- factor(panels[[1L]], levels = panels)
  colours <- .plot_discrete_values(factor(strata, levels = strata))
  conf_int <- conf_int && !anyNA(curves$upper[!is.na(curves$surv)]) &&
    !is.null(fit$upper)

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
      inherit.aes = FALSE, alpha = 0.15, colour = NA, na.rm = TRUE
    )
  }
  plot <- plot +
    ggplot2::geom_step(direction = "hv", linewidth = 0.45) +
    ggplot2::geom_point(
      data = curves[curves$n_censor > 0, , drop = FALSE],
      shape = 3, size = 0.8, stroke = 0.35, show.legend = FALSE
    )
  breaks <- seq(0, 1, by = 0.25)
  labels <- format(breaks, nsmall = 2)
  if (risk_table) {
    table <- .plot_km_risk_table(fit, extracted, risk_times)
    table$panel <- factor(panels[[2L]], levels = panels)
    header <- data.frame(
      time = min(table$time), y = -0.035, label = "Number at risk",
      panel = factor(panels[[2L]], levels = panels)
    )
    plot <- plot +
      ggplot2::geom_text(
        data = table,
        mapping = ggplot2::aes(
          x = .data[["time"]], y = .data[["y"]], label = .data[["n_risk"]],
          colour = .data[["stratum"]]
        ),
        size = .bulkmae_text_size_pt, size.unit = "pt", show.legend = FALSE
      ) +
      ggplot2::geom_text(
        data = header,
        mapping = ggplot2::aes(x = .data[["time"]], y = .data[["y"]], label = .data[["label"]]),
        inherit.aes = FALSE, hjust = 0, fontface = "bold",
        size = .bulkmae_text_size_pt, size.unit = "pt"
      ) +
      ggplot2::facet_grid(
        rows = ggplot2::vars(.data[["panel"]]), scales = "free_y", space = "free_y"
      )
    rows <- unique(table[c("stratum", "y")])
    breaks <- c(breaks, rows$y)
    labels <- c(labels, as.character(rows$stratum))
  }
  if (!is.null(test)) {
    annotation <- data.frame(
      time = 0, surv = 0.04,
      label = paste("Log-rank", .plot_format_p(.plot_survdiff_p(test, strata))),
      panel = factor(panels[[1L]], levels = panels)
    )
    plot <- plot + ggplot2::geom_text(
      data = annotation,
      mapping = ggplot2::aes(x = .data[["time"]], y = .data[["surv"]], label = .data[["label"]]),
      inherit.aes = FALSE, hjust = 0, vjust = 0,
      size = .bulkmae_text_size_pt, size.unit = "pt"
    )
  }
  plot <- plot +
    ggplot2::scale_colour_manual(values = colours, breaks = strata, drop = FALSE) +
    ggplot2::scale_fill_manual(values = colours, breaks = strata, drop = FALSE, guide = "none") +
    ggplot2::scale_y_continuous(breaks = breaks, labels = labels,
                                expand = ggplot2::expansion(mult = 0.04, add = c(0.05, 0))) +
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0.04, 0.03))) +
    ggplot2::coord_cartesian(clip = "off") +
    ggplot2::labs(
      x = "Time", y = "Survival probability", colour = extracted$title,
      alt = paste(
        "Kaplan-Meier step curves by stratum with censoring marks",
        if (risk_table) "and a number-at-risk table below." else "."
      )
    ) +
    .plot_survival_theme() +
    ggplot2::theme(
      strip.text = ggplot2::element_blank(),
      legend.position = "top",
      legend.direction = "horizontal",
      panel.spacing.y = grid::unit(0.6, "lines")
    )
  .plot_with_dimensions(
    plot,
    width = 8,
    height = if (risk_table) 7 + 0.35 * length(strata) else 6.5
  )
}

#' Plot hazard ratios from Cox models as a forest plot
#'
#' Draws already-estimated hazard ratios with their confidence intervals on a
#' log axis. The right-hand axis prints each estimate, interval, and Wald
#' p-value, so the figure carries its own table. Pass one multivariable
#' `coxph` fit, or a named list of univariate fits for a screening forest.
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
  table_labels <- sprintf(
    paste0("%.2f (%.2f", intToUtf8(8211L), "%.2f)   %s"),
    data$hazard_ratio, data$lower, data$upper, .plot_format_p(data$p_value)
  )
  names(table_labels) <- as.character(data$term)

  plot <- ggplot2::ggplot(data, ggplot2::aes(
    x = .data[["hazard_ratio"]], y = .data[["term"]]
  )) +
    ggplot2::geom_vline(xintercept = 1, colour = "#999999",
                        linetype = "dashed", linewidth = 0.3) +
    ggplot2::geom_linerange(
      mapping = ggplot2::aes(xmin = .data[["lower"]], xmax = .data[["upper"]]),
      linewidth = 0.4, colour = "#222222"
    ) +
    ggplot2::geom_point(shape = 15, size = 1.6, colour = .bulkmae_colours[["blue"]]) +
    ggplot2::scale_x_log10(breaks = .plot_hr_breaks(c(data$lower, data$upper))) +
    ggplot2::scale_y_discrete(
      sec.axis = ggplot2::dup_axis(
        labels = unname(table_labels[levels(data$term)]),
        name = sprintf("HR (%g%% CI)   p", 100 * conf_level)
      )
    ) +
    ggplot2::labs(
      x = sprintf("Hazard ratio (%g%% CI, log scale)", 100 * conf_level),
      y = NULL,
      alt = paste(
        "Forest plot of Cox hazard ratios with confidence intervals on a",
        "log axis; estimates and p-values are printed on the right."
      )
    ) +
    theme_bulkmae() +
    ggplot2::theme(
      panel.grid.major.y = ggplot2::element_blank(),
      panel.grid.major.x = ggplot2::element_line(colour = "#EBEBEB", linewidth = 0.25),
      axis.line.x = ggplot2::element_line(colour = "#222222", linewidth = 0.3),
      axis.ticks.x = ggplot2::element_line(colour = "#222222", linewidth = 0.3),
      axis.text.y.right = ggplot2::element_text(hjust = 0, colour = "#333333"),
      axis.title.y.right = ggplot2::element_text(angle = 0, vjust = 1, hjust = 0)
    )
  .plot_with_dimensions(
    plot,
    width = 11,
    height = .plot_clamped_dimension(nrow(data), base = 1.6, per_item = 0.42,
                                     minimum = 4, maximum = 18)
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
  top_breaks <- data$log_lambda[index]
  top_labels <- data$nzero[index]
  chosen <- data.frame(
    log_lambda = log(c(fit$lambda.min, fit$lambda.1se)),
    rule = c("lambda.min", "lambda.1se")
  )
  plot <- ggplot2::ggplot(data, ggplot2::aes(x = .data[["log_lambda"]], y = .data[["cvm"]])) +
    ggplot2::geom_vline(
      data = chosen, mapping = ggplot2::aes(xintercept = .data[["log_lambda"]]),
      colour = "#999999", linetype = "dashed", linewidth = 0.3
    ) +
    ggplot2::geom_linerange(
      mapping = ggplot2::aes(ymin = .data[["cvlo"]], ymax = .data[["cvup"]]),
      colour = "#BBBBBB", linewidth = 0.3
    ) +
    ggplot2::geom_point(colour = .bulkmae_colours[["orange"]], size = 0.9) +
    ggplot2::scale_x_continuous(
      sec.axis = ggplot2::sec_axis(
        transform = ~ ., breaks = top_breaks, labels = top_labels,
        name = "Non-zero coefficients"
      )
    ) +
    ggplot2::labs(
      x = "Log(lambda)", y = unname(as.character(fit$name)),
      alt = paste(
        "Cross-validated error with one-standard-error bars along the glmnet",
        "penalty path, with dashed lines at lambda.min and lambda.1se."
      )
    ) +
    .plot_survival_theme()
  .plot_with_dimensions(plot, width = 8, height = 6.5)
}

#' Plot glmnet coefficient paths
#'
#' Draws each feature's coefficient across the penalty path of a `cv.glmnet`
#' fit. Features that are non-zero at the chosen `lambda` are coloured by the
#' sign of that coefficient; all other paths are grey. Text labels are added
#' only for features named in `label_features`.
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
  palette <- c(Negative = .bulkmae_colours[["blue"]], Zero = "#D9D9D9",
               Positive = .bulkmae_colours[["orange"]])
  plot <- ggplot2::ggplot(data, ggplot2::aes(
    x = .data[["log_lambda"]], y = .data[["coefficient"]],
    group = .data[["feature"]], colour = .data[["status"]]
  )) +
    ggplot2::geom_hline(yintercept = 0, colour = "#999999", linewidth = 0.25) +
    ggplot2::geom_vline(xintercept = log(path_lambda[[chosen]]),
                        colour = "#999999", linetype = "dashed", linewidth = 0.3) +
    ggplot2::geom_line(linewidth = 0.35) +
    ggplot2::scale_colour_manual(values = palette, breaks = names(palette), drop = FALSE) +
    ggplot2::labs(
      x = "Log(lambda)", y = "Coefficient",
      colour = paste("At", lambda),
      alt = paste(
        "Coefficient paths along the glmnet penalty, coloured by the sign of",
        "each coefficient at the marked lambda."
      )
    ) +
    .plot_survival_theme() +
    ggplot2::theme(legend.position = "top", legend.direction = "horizontal")
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
        size = .bulkmae_text_size_pt, size.unit = "pt"
      ) +
      ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0.2, 0.03)))
  }
  .plot_with_dimensions(plot, width = 8, height = 6.5)
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
  boundary <- sum(group == "Low") + 0.5
  plot <- ggplot2::ggplot(mapping = ggplot2::aes(x = .data[["rank"]], y = .data[["value"]])) +
    ggplot2::geom_vline(xintercept = boundary, colour = "#999999",
                        linetype = "dashed", linewidth = 0.3) +
    ggplot2::geom_hline(
      data = data.frame(value = cutoff, panel = factor(panels[[1L]], levels = panels)),
      mapping = ggplot2::aes(yintercept = .data[["value"]]),
      colour = "#999999", linetype = "dashed", linewidth = 0.3
    ) +
    ggplot2::geom_point(data = top, mapping = ggplot2::aes(colour = .data[["group"]]),
                        size = 0.7) +
    ggplot2::geom_point(data = bottom, mapping = ggplot2::aes(fill = .data[["status"]]),
                        shape = 21, colour = "transparent", stroke = 0, size = 0.9) +
    ggplot2::facet_grid(rows = ggplot2::vars(.data[["panel"]]), scales = "free_y",
                        switch = "y") +
    ggplot2::scale_colour_manual(
      values = c(Low = .bulkmae_colours[["blue"]], High = .bulkmae_colours[["orange"]]),
      drop = FALSE, name = "Risk group"
    ) +
    ggplot2::scale_fill_manual(
      values = c(Censored = "#999999", Event = "#000000"),
      drop = FALSE, name = "Status"
    ) +
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = 0.01)) +
    ggplot2::labs(
      x = "Samples ranked by risk score", y = NULL,
      alt = paste(
        "Two aligned panels: samples sorted by risk score and split at the",
        "cutoff, and each sample's follow-up time coloured by event status."
      )
    ) +
    .plot_survival_theme() +
    ggplot2::theme(
      strip.placement = "outside",
      strip.text.y.left = ggplot2::element_text(angle = 90, face = "plain", colour = "#222222"),
      legend.position = "top", legend.direction = "horizontal",
      panel.spacing.y = grid::unit(0.8, "lines")
    )
  .plot_with_dimensions(plot, width = 9, height = 8)
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
  if (!inherits(roc, "ipcwsurvivalROC") && !inherits(roc, "ipcwcompetingrisksROC")) {
    stop("`roc` must be a timeROC result.", call. = FALSE)
  }
  fp <- if (!is.null(roc$FP)) roc$FP else roc$FP_1
  tp <- roc$TP
  auc <- if (!is.null(roc$AUC)) roc$AUC else roc$AUC_1
  if (is.null(fp) || is.null(tp) || is.null(auc)) {
    stop("`roc` does not contain TP, FP, and AUC values.", call. = FALSE)
  }
  fp <- as.matrix(fp)
  tp <- as.matrix(tp)
  times <- roc$times
  if (is.null(time_labels)) {
    time_labels <- paste("t =", format(times, trim = TRUE))
  }
  if (!is.character(time_labels) || length(time_labels) != length(times) ||
      anyNA(time_labels) || any(!nzchar(time_labels)) || anyDuplicated(time_labels)) {
    stop("`time_labels` must give one unique label per evaluation time.", call. = FALSE)
  }
  keep <- is.finite(auc)
  if (!any(keep)) stop("No evaluation time has a finite AUC.", call. = FALSE)
  legend <- sprintf("%s (AUC %.3f)", time_labels, auc)
  data <- do.call(rbind, lapply(which(keep), function(index) {
    curve <- data.frame(fp = fp[, index], tp = tp[, index])
    curve <- curve[is.finite(curve$fp) & is.finite(curve$tp), , drop = FALSE]
    curve <- curve[order(curve$fp, curve$tp), , drop = FALSE]
    curve$time <- legend[[index]]
    curve
  }))
  data$time <- factor(data$time, levels = legend[keep])
  plot <- ggplot2::ggplot(data, ggplot2::aes(
    x = .data[["fp"]], y = .data[["tp"]], colour = .data[["time"]]
  )) +
    ggplot2::geom_abline(slope = 1, intercept = 0, colour = "#999999",
                         linetype = "dashed", linewidth = 0.3) +
    ggplot2::geom_path(linewidth = 0.5) +
    ggplot2::scale_colour_manual(values = .plot_discrete_values(data$time)) +
    ggplot2::scale_x_continuous(limits = c(0, 1), expand = ggplot2::expansion(mult = 0.01)) +
    ggplot2::scale_y_continuous(limits = c(0, 1), expand = ggplot2::expansion(mult = 0.01)) +
    ggplot2::coord_equal() +
    ggplot2::labs(
      x = "False positive rate", y = "True positive rate", colour = NULL,
      alt = "Time-dependent ROC curves, one per evaluation time, with AUC in the legend."
    ) +
    .plot_survival_theme() +
    ggplot2::theme(
      legend.position = "inside",
      legend.position.inside = c(0.98, 0.02),
      legend.justification = c(1, 0),
      legend.background = ggplot2::element_rect(fill = "white", colour = NA)
    )
  .plot_with_dimensions(plot, width = 7, height = 7)
}

.plot_survival_theme <- function() {
  theme_bulkmae() +
    ggplot2::theme(
      panel.grid.major = ggplot2::element_blank(),
      axis.line = ggplot2::element_line(colour = "#222222", linewidth = 0.3),
      axis.ticks = ggplot2::element_line(colour = "#222222", linewidth = 0.3),
      plot.background = ggplot2::element_rect(fill = "white", colour = NA)
    )
}

.plot_hr_breaks <- function(values) {
  candidates <- c(0.01, 0.02, 0.05, 0.1, 0.2, 0.3, 0.5, 0.7, 1, 1.5, 2, 3, 5,
                  7, 10, 20, 50, 100)
  limits <- range(c(values, 1))
  breaks <- candidates[candidates >= limits[[1L]] & candidates <= limits[[2L]]]
  while (length(breaks) > 6L) {
    keep <- seq_along(breaks) %% 2L == match(1, breaks) %% 2L
    breaks <- breaks[keep]
  }
  breaks
}

.plot_format_p <- function(p) {
  ifelse(
    is.na(p), "p = NA",
    ifelse(p < 1e-4, "p < 0.0001",
           ifelse(p < 0.001, sprintf("p = %.1e", p), sprintf("p = %.3f", p)))
  )
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
    title <- NULL
  } else {
    raw <- rep(names(fit$strata), fit$strata)
    strata <- names(fit$strata)
    title <- "Stratum"
    single <- !grepl(",", strata, fixed = TRUE) & grepl("=", strata, fixed = TRUE)
    if (all(single)) {
      prefixes <- unique(sub("=.*$", "", strata))
      if (length(prefixes) == 1L) {
        title <- prefixes
        strata <- sub("^[^=]*=", "", strata)
        raw <- sub("^[^=]*=", "", raw)
      }
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
  list(curves = curves, strata = strata, title = title)
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
    y = -0.035 - 0.085 * position
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
