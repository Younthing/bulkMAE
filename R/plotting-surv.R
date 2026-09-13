#' Plot Kaplan-Meier curves from a fitted survfit object
#'
#' Draws already-computed [survival::survfit()] curves as a standard ggplot.
#' When `risk_table = TRUE`, the number at risk is printed under the panel in
#' the same ggplot so the constructor still returns one unprinted plot. The
#' function does not refit Kaplan-Meier estimates or call survminer.
#'
#' @param fit A `survfit` object, usually from [surv_km()].
#' @param conf_int Draw pointwise confidence bands when `fit` contains them.
#' @param censor Mark censoring times.
#' @param risk_table Print the number at risk under the survival panel.
#' @param risk_times Optional non-negative times for the risk table. The
#'   default uses pretty breaks of the observed follow-up.
#'
#' @return An unprinted standard ggplot object carrying recommended physical
#'   dimensions for [plot_save()].
#' @family plotting
#' @export
plot_surv_km <- function(
    fit,
    conf_int = TRUE,
    censor = TRUE,
    risk_table = TRUE,
    risk_times = NULL
) {
  .require_backend("survival", "to plot Kaplan-Meier curves")
  .plot_assert_flag(conf_int, "conf_int")
  .plot_assert_flag(censor, "censor")
  .plot_assert_flag(risk_table, "risk_table")
  curves <- .plot_survfit_frame(fit)
  palette <- .plot_surv_palette(levels(curves$strata))
  stepped <- .plot_surv_stairstep(curves)
  has_interval <- conf_int &&
    any(is.finite(curves$lower) & is.finite(curves$upper))
  censored <- curves[curves$n_censor > 0 & is.finite(curves$surv), , drop = FALSE]

  plot <- ggplot2::ggplot()
  if (has_interval) {
    ribbon <- stepped[is.finite(stepped$lower) & is.finite(stepped$upper), , drop = FALSE]
    if (nrow(ribbon)) {
      plot <- plot +
        ggplot2::geom_ribbon(
          data = ribbon,
          mapping = ggplot2::aes(
            x = .data[["time"]],
            ymin = .data[["lower"]],
            ymax = .data[["upper"]],
            fill = .data[["strata"]]
          ),
          alpha = 0.16,
          colour = NA,
          show.legend = FALSE
        )
    }
  }
  plot <- plot +
    ggplot2::geom_step(
      data = curves,
      mapping = ggplot2::aes(
        x = .data[["time"]],
        y = .data[["surv"]],
        colour = .data[["strata"]]
      ),
      linewidth = 0.45,
      direction = "hv"
    )
  if (censor && nrow(censored)) {
    plot <- plot +
      ggplot2::geom_point(
        data = censored,
        mapping = ggplot2::aes(
          x = .data[["time"]],
          y = .data[["surv"]],
          colour = .data[["strata"]]
        ),
        shape = 3,
        size = 1.1,
        stroke = 0.4,
        show.legend = FALSE
      )
  }

  n_strata <- nlevels(curves$strata)
  y_breaks <- c(0, 0.25, 0.5, 0.75, 1)
  y_limits <- c(0, 1.02)
  if (risk_table) {
    risk <- .plot_surv_risk_table(fit, curves, risk_times)
    row_gap <- 0.11
    table_top <- -0.10
    table_ys <- table_top - row_gap * (seq_len(n_strata) - 1L)
    names(table_ys) <- levels(curves$strata)
    risk$y <- unname(table_ys[as.character(risk$strata)])
    label_x <- min(c(0, curves$time, risk$time), na.rm = TRUE)
    y_limits <- c(min(table_ys) - 0.08, 1.05)
    plot <- plot +
      ggplot2::annotate(
        "segment",
        x = -Inf,
        xend = Inf,
        y = 0,
        yend = 0,
        colour = "#BBBBBB",
        linewidth = 0.2
      ) +
      ggplot2::geom_text(
        data = risk,
        mapping = ggplot2::aes(
          x = .data[["time"]],
          y = .data[["y"]],
          label = .data[["n_risk"]],
          colour = .data[["strata"]]
        ),
        size = .bulkmae_text_size_pt,
        size.unit = "pt",
        vjust = 0.5,
        show.legend = FALSE
      ) +
      ggplot2::annotate(
        "text",
        x = label_x,
        y = unname(table_ys),
        label = names(table_ys),
        hjust = 1.15,
        size = .bulkmae_text_size_pt,
        size.unit = "pt",
        colour = "#333333"
      )
  }

  plot <- plot +
    ggplot2::scale_colour_manual(values = palette, drop = FALSE)
  if (has_interval) {
    plot <- plot +
      ggplot2::scale_fill_manual(values = palette, drop = FALSE)
  }
  plot <- plot +
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0.02, 0.04))) +
    ggplot2::scale_y_continuous(
      breaks = y_breaks,
      labels = function(breaks) {
        ifelse(breaks < 0, "", format(breaks, nsmall = 2L, trim = TRUE))
      }
    ) +
    ggplot2::coord_cartesian(ylim = y_limits, clip = "off") +
    ggplot2::labs(
      x = "Time",
      y = "Survival probability",
      colour = if (n_strata > 1L) "Stratum" else NULL,
      alt = paste0(
        "Kaplan-Meier survival curves",
        if (has_interval) " with pointwise confidence bands" else "",
        if (risk_table) " and a number-at-risk table" else "",
        "."
      )
    ) +
    theme_bulkmae() +
    ggplot2::theme(
      legend.position = if (n_strata > 1L) "right" else "none",
      plot.margin = if (risk_table) {
        ggplot2::margin(4, 8, 8, 28, unit = "pt")
      } else {
        ggplot2::margin(4, 6, 4, 6, unit = "pt")
      }
    )
  .plot_with_dimensions(
    plot,
    width = .plot_clamped_dimension(
      n_strata,
      base = 9,
      per_item = 0.4,
      minimum = 9.5,
      maximum = 14
    ),
    height = if (risk_table) {
      .plot_clamped_dimension(
        n_strata,
        base = 7.2,
        per_item = 0.45,
        minimum = 7.5,
        maximum = 11
      )
    } else {
      6.8
    }
  )
}

#' Plot a Cox hazard-ratio forest from already-computed results
#'
#' Accepts a [survival::coxph()] fit or a data frame from [surv_cox_table()]
#' / [surv_cox_univariable()]. Hazard ratios are drawn on a logarithmic x
#' axis. The constructor does not refit Cox models.
#'
#' @param x A `coxph` object or a data frame with `term`, `hazard_ratio`,
#'   `conf_low`, and `conf_high`.
#' @param term_labels Optional uniquely named labels for a subset or all terms.
#'
#' @return An unprinted standard ggplot object carrying recommended physical
#'   dimensions for [plot_save()].
#' @family plotting
#' @export
plot_surv_forest <- function(x, term_labels = NULL) {
  table <- .plot_surv_forest_table(x)
  table$display <- .plot_surv_term_labels(table$term, term_labels)
  table$display <- factor(table$display, levels = rev(table$display))
  table$hr_label <- sprintf(
    "%.2f (%.2f-%.2f)",
    table$hazard_ratio,
    table$conf_low,
    table$conf_high
  )
  x_max <- max(table$conf_high, 1)
  x_min <- min(table$conf_low, 1)
  if (x_min <= 0) {
    stop("Forest-plot confidence limits must be positive.", call. = FALSE)
  }
  label_x <- exp(log(x_max) + 0.18 * (log(x_max) - log(x_min)))
  plot <- ggplot2::ggplot(table, ggplot2::aes(x = .data[["hazard_ratio"]], y = .data[["display"]])) +
    ggplot2::geom_vline(
      xintercept = 1,
      linetype = 2,
      colour = "#777777",
      linewidth = 0.3
    ) +
    ggplot2::geom_linerange(
      mapping = ggplot2::aes(xmin = .data[["conf_low"]], xmax = .data[["conf_high"]]),
      linewidth = 0.4,
      colour = .bulkmae_colours[["blue"]]
    ) +
    ggplot2::geom_point(
      size = 1.6,
      colour = .bulkmae_colours[["blue"]]
    ) +
    ggplot2::geom_text(
      mapping = ggplot2::aes(x = label_x, label = .data[["hr_label"]]),
      size = .bulkmae_text_size_pt,
      size.unit = "pt",
      hjust = 0,
      colour = "#333333"
    ) +
    ggplot2::scale_x_continuous(
      trans = "log",
      breaks = .plot_surv_hr_breaks(x_min, x_max)
    ) +
    ggplot2::coord_cartesian(clip = "off") +
    ggplot2::labs(
      x = "Hazard ratio",
      y = "Term",
      alt = "Forest plot of Cox hazard ratios with confidence intervals."
    ) +
    theme_bulkmae() +
    ggplot2::theme(
      panel.grid.major.y = ggplot2::element_blank(),
      plot.margin = ggplot2::margin(4, 42, 4, 6, unit = "pt")
    )
  .plot_with_dimensions(
    plot,
    width = .plot_clamped_dimension(
      nrow(table),
      base = 9.5,
      per_item = 0.15,
      minimum = 10,
      maximum = 14
    ),
    height = .plot_clamped_dimension(
      nrow(table),
      base = 3.2,
      per_item = 0.55,
      minimum = 4.5,
      maximum = 12
    )
  )
}

#' Plot a risk-score distribution or outcome-stratified box plot
#'
#' Displays an already-computed, sample-named risk score. Supply `group` to
#' colour a density or to draw boxes by an aligned outcome or risk stratum.
#' The constructor does not compute the score or fit a survival model.
#'
#' @param score A uniquely named finite numeric vector, one value per sample.
#' @param group Optional named vector or factor aligned to `names(score)`.
#' @param type `"density"` or `"box"`. A box plot requires `group`.
#'
#' @return An unprinted standard ggplot object carrying recommended physical
#'   dimensions for [plot_save()].
#' @family plotting
#' @export
plot_surv_risk <- function(score, group = NULL, type = c("density", "box")) {
  type <- match.arg(type)
  .assert_named_numeric(score, "score")
  data <- data.frame(
    sample = names(score),
    score = as.numeric(score),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  if (!is.null(group)) {
    group <- .plot_surv_align_group(group, names(score))
    data$group <- group[data$sample]
  }
  if (identical(type, "box")) {
    if (is.null(group)) {
      stop("`group` is required when `type = \"box\"`.", call. = FALSE)
    }
    palette <- .plot_surv_palette(levels(data$group))
    plot <- ggplot2::ggplot(
      data,
      ggplot2::aes(x = .data[["group"]], y = .data[["score"]], fill = .data[["group"]])
    ) +
      ggplot2::geom_boxplot(
        width = 0.62,
        outlier.size = 0.8,
        colour = "#333333",
        linewidth = 0.3
      ) +
      ggplot2::scale_fill_manual(values = palette, name = "Group", drop = FALSE) +
      ggplot2::labs(
        x = "Group",
        y = "Risk score",
        alt = "Box plot of a risk score stratified by an aligned outcome or risk group."
      ) +
      theme_bulkmae() +
      ggplot2::theme(legend.position = "none")
    return(.plot_with_dimensions(
      plot,
      width = .plot_clamped_dimension(
        nlevels(data$group),
        base = 6.5,
        per_item = 1.1,
        minimum = 7.5,
        maximum = 12
      ),
      height = 6.5
    ))
  }

  if (is.null(group)) {
    plot <- ggplot2::ggplot(data, ggplot2::aes(x = .data[["score"]])) +
      ggplot2::geom_density(
        fill = .bulkmae_colours[["blue"]],
        colour = .bulkmae_colours[["blue"]],
        alpha = 0.28,
        linewidth = 0.4
      ) +
      ggplot2::labs(
        x = "Risk score",
        y = "Density",
        alt = "Density of a sample-named risk score."
      ) +
      theme_bulkmae()
  } else {
    palette <- .plot_surv_palette(levels(data$group))
    plot <- ggplot2::ggplot(
      data,
      ggplot2::aes(x = .data[["score"]], fill = .data[["group"]], colour = .data[["group"]])
    ) +
      ggplot2::geom_density(alpha = 0.22, linewidth = 0.4) +
      ggplot2::scale_fill_manual(values = palette, name = "Group", drop = FALSE) +
      ggplot2::scale_colour_manual(values = palette, name = "Group", drop = FALSE) +
      ggplot2::labs(
        x = "Risk score",
        y = "Density",
        alt = "Overlaid densities of a risk score coloured by an aligned group."
      ) +
      theme_bulkmae()
  }
  .plot_with_dimensions(plot, width = 8.5, height = 6.5)
}

.plot_survfit_frame <- function(fit) {
  if (!inherits(fit, "survfit")) {
    stop("`fit` must be a survfit object.", call. = FALSE)
  }
  fit <- survival::survfit0(fit)
  n <- length(fit$time)
  if (!n) {
    stop("The survfit object contains no time points.", call. = FALSE)
  }
  if (is.null(fit$strata)) {
    strata <- rep("All", n)
    strata_levels <- "All"
  } else {
    raw_levels <- names(fit$strata)
    strata_levels <- .plot_surv_strata_labels(raw_levels)
    strata <- rep(strata_levels, times = unname(fit$strata))
  }
  data <- data.frame(
    time = as.numeric(fit$time),
    surv = as.numeric(fit$surv),
    lower = if (is.null(fit$lower)) rep(NA_real_, n) else as.numeric(fit$lower),
    upper = if (is.null(fit$upper)) rep(NA_real_, n) else as.numeric(fit$upper),
    n_risk = as.numeric(fit$n.risk),
    n_censor = as.numeric(fit$n.censor),
    n_event = as.numeric(fit$n.event),
    strata = factor(strata, levels = unique(strata_levels)),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  if (any(!is.finite(data$time)) || any(!is.finite(data$surv))) {
    stop("survfit times and survival probabilities must be finite.", call. = FALSE)
  }
  data
}

.plot_surv_strata_labels <- function(strata_names) {
  if (!length(strata_names)) {
    return(strata_names)
  }
  prefixes <- sub("=.*$", "", strata_names)
  stripped <- sub("^[^=]+=", "", strata_names)
  if (length(unique(prefixes)) == 1L && !identical(stripped, strata_names)) {
    return(stripped)
  }
  strata_names
}

.plot_surv_stairstep <- function(data) {
  pieces <- lapply(split(data, data$strata, drop = TRUE), function(piece) {
    piece <- piece[order(piece$time), , drop = FALSE]
    n <- nrow(piece)
    if (n == 1L) {
      return(piece)
    }
    idx <- as.vector(rbind(seq_len(n), seq_len(n)))
    out <- piece[idx, , drop = FALSE]
    times <- as.vector(rbind(piece$time, c(piece$time[-1L], piece$time[[n]])))
    out$time <- times
    out[-nrow(out), , drop = FALSE]
  })
  result <- do.call(rbind, pieces)
  rownames(result) <- NULL
  result
}

.plot_surv_risk_table <- function(fit, curves, risk_times) {
  times <- .plot_surv_risk_times(curves$time, risk_times)
  summarised <- summary(fit, times = times, extend = TRUE)
  n <- length(summarised$time)
  if (!n) {
    stop("The number-at-risk table is empty at the requested times.", call. = FALSE)
  }
  if (is.null(summarised$strata)) {
    strata <- factor(rep("All", n), levels = levels(curves$strata))
  } else {
    raw <- as.character(summarised$strata)
    labels <- .plot_surv_strata_labels(raw)
    strata <- factor(labels, levels = levels(curves$strata))
  }
  data.frame(
    time = as.numeric(summarised$time),
    n_risk = as.integer(summarised$n.risk),
    strata = strata,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}

.plot_surv_risk_times <- function(time, times) {
  if (is.null(times)) {
    upper <- max(time)
    times <- pretty(c(0, upper), n = 5L)
    times <- times[times >= 0 & times <= upper + sqrt(.Machine$double.eps)]
    if (!length(times)) {
      times <- c(0, upper)
    }
    if (!0 %in% times) {
      times <- c(0, times)
    }
    times <- unique(times)
  }
  .plot_assert_finite_numeric(times, "`risk_times`")
  if (any(times < 0)) {
    stop("`risk_times` must be non-negative.", call. = FALSE)
  }
  sort(unique(as.numeric(times)))
}

.plot_surv_forest_table <- function(x) {
  table <- if (inherits(x, "coxph")) {
    surv_cox_table(x)
  } else if (is.data.frame(x)) {
    as.data.frame(x, optional = TRUE)
  } else {
    stop("`x` must be a coxph object or a Cox table data frame.", call. = FALSE)
  }
  required <- c("term", "hazard_ratio", "conf_low", "conf_high")
  .plot_require_columns(table, required, "x")
  .plot_assert_ids(as.character(table$term), "`x$term`")
  for (column in c("hazard_ratio", "conf_low", "conf_high")) {
    .plot_assert_finite_numeric(table[[column]], paste0("`x$", column, "`"))
  }
  if (any(table$hazard_ratio <= 0) || any(table$conf_low <= 0) ||
      any(table$conf_high <= 0)) {
    stop("Hazard ratios and confidence limits must be positive.", call. = FALSE)
  }
  if (any(table$conf_low > table$hazard_ratio) ||
      any(table$conf_high < table$hazard_ratio)) {
    stop("Each confidence interval must contain its hazard ratio.", call. = FALSE)
  }
  table$term <- as.character(table$term)
  table
}

.plot_surv_term_labels <- function(terms, labels) {
  if (is.null(labels)) {
    return(terms)
  }
  if (is.null(names(labels)) || anyNA(names(labels)) || any(!nzchar(names(labels))) ||
      anyDuplicated(names(labels))) {
    stop("`term_labels` must be uniquely term-named.", call. = FALSE)
  }
  unknown <- setdiff(names(labels), terms)
  if (length(unknown)) {
    stop(
      "`term_labels` contains unknown terms: ",
      paste(unknown, collapse = ", "),
      ".",
      call. = FALSE
    )
  }
  result <- terms
  names(result) <- terms
  result[names(labels)] <- as.character(labels)
  if (anyNA(result) || any(!nzchar(result))) {
    stop("`term_labels` cannot contain missing or empty labels.", call. = FALSE)
  }
  unname(result)
}

.plot_surv_hr_breaks <- function(x_min, x_max) {
  candidates <- c(0.1, 0.2, 0.25, 0.5, 1, 2, 4, 5, 8, 10, 16, 20)
  selected <- candidates[candidates >= x_min / 1.2 & candidates <= x_max * 1.2]
  if (!length(selected)) {
    selected <- pretty(c(x_min, x_max), n = 4L)
    selected <- selected[selected > 0]
  }
  unique(c(1, selected))
}

.plot_surv_align_group <- function(group, sample_ids) {
  if (is.data.frame(group) || is.matrix(group)) {
    stop("`group` must be a named vector or factor.", call. = FALSE)
  }
  if (is.null(names(group))) {
    stop("`group` must be sample-named.", call. = FALSE)
  }
  .plot_assert_ids(names(group), "Names of `group`")
  if (!setequal(names(group), sample_ids)) {
    stop("`group` names must match `score` names exactly.", call. = FALSE)
  }
  group <- group[sample_ids]
  if (anyNA(group) || any(!nzchar(as.character(group)))) {
    stop("`group` cannot contain missing or empty values.", call. = FALSE)
  }
  if (is.numeric(group) && !is.factor(group) && length(unique(group)) > 8L) {
    stop("`group` must be discrete; convert a continuous variable to a factor first.", call. = FALSE)
  }
  if (is.factor(group)) {
    factor(as.character(group), levels = levels(droplevels(group)))
  } else {
    factor(as.character(group), levels = unique(as.character(group)))
  }
}

.plot_surv_palette <- function(levels) {
  levels <- as.character(levels)
  if (!length(levels)) {
    stop("At least one group is required.", call. = FALSE)
  }
  if (length(levels) > length(.bulkmae_qualitative)) {
    return(unname(.bulkmae_qualitative))
  }
  stats::setNames(unname(.bulkmae_qualitative[seq_along(levels)]), levels)
}
