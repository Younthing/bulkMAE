#' Plot a regulator-by-sample activity heatmap
#'
#' Displays already-inferred activity scores. The fill scale is the backend
#' activity score (for example an ULM *t*-statistic), or a row *z*-score of
#' that score. It is not an enrichment NES or a GO/KEGG bar statistic.
#' Phenotype tracks are drawn as annotation bars; they do not re-run
#' inference.
#'
#' `sources` may be omitted when the result already contains the regulators
#' to display. The function does not select "top" regulators from *p*-values.
#'
#' @param result A data-frame-like result returned by [activity_decouple()],
#'   [activity_progeny()], or [activity_tf()].
#' @param sources Optional regulator/source identifiers to display, in the
#'   requested order before clustering. The default uses every source in
#'   `result`.
#' @param sample_data Optional sample metadata whose row names cover every
#'   plotted sample. Required when `annotation` names metadata columns.
#' @param annotation Optional phenotype track: one or more `sample_data`
#'   column names, or one complete sample-named vector.
#' @param scale Standardize each regulator row or retain activity scores.
#' @param cluster_rows,cluster_columns Reorder regulators or samples by
#'   hierarchical clustering. Dendrograms are not drawn.
#' @inheritParams activity_matrix
#'
#' @return An unprinted standard ggplot object carrying recommended physical
#'   dimensions for [plot_save()].
#' @family plotting
#' @export
plot_activity_heatmap <- function(
    result,
    sources = NULL,
    sample_data = NULL,
    annotation = NULL,
    scale = c("none", "row"),
    cluster_rows = TRUE,
    cluster_columns = TRUE,
    value = "score",
    source = "source",
    sample = "condition",
    statistic = NULL,
    statistic_column = "statistic"
) {
  scale <- match.arg(scale)
  .plot_assert_flag(cluster_rows, "cluster_rows")
  .plot_assert_flag(cluster_columns, "cluster_columns")
  prepared <- .plot_activity_matrix(
    result, sources, statistic, value, source, sample, statistic_column
  )
  matrix <- prepared$matrix
  if (scale == "row") {
    deviations <- apply(matrix, 1L, stats::sd)
    constant <- !is.finite(deviations) | deviations == 0
    if (any(constant)) {
      warning(
        "Removed constant activity rows: ",
        paste(rownames(matrix)[constant], collapse = ", "),
        ".",
        call. = FALSE
      )
      matrix <- matrix[!constant, , drop = FALSE]
    }
    if (!nrow(matrix)) stop("No non-constant activity rows remain.", call. = FALSE)
    matrix <- t(scale(t(matrix)))
  }
  row_order <- rownames(matrix)
  if (cluster_rows && nrow(matrix) > 1L) {
    row_order <- row_order[stats::hclust(stats::dist(matrix))$order]
  }
  column_order <- colnames(matrix)
  if (cluster_columns && ncol(matrix) > 1L) {
    column_order <- column_order[stats::hclust(stats::dist(t(matrix)))$order]
  }
  tracks <- .plot_activity_annotation_tracks(annotation, sample_data, column_order)
  if (!is.null(tracks)) {
    first <- factor(
      tracks[[1L]]$values,
      levels = unique(as.character(tracks[[1L]]$values))
    )
    column_order <- column_order[
      order(first[column_order], match(column_order, column_order))
    ]
  }
  data <- expand.grid(
    source_id = row_order,
    sample = column_order,
    stringsAsFactors = FALSE
  )
  data$value <- matrix[cbind(
    match(data$source_id, rownames(matrix)),
    match(data$sample, colnames(matrix))
  )]
  track_names <- if (is.null(tracks)) character() else {
    vapply(tracks, `[[`, character(1), "name")
  }
  overlap <- intersect(track_names, row_order)
  if (length(overlap)) {
    stop(
      "Annotation track names must not match activity source names: ",
      paste(overlap, collapse = ", "),
      ".",
      call. = FALSE
    )
  }
  y_levels <- c(rev(row_order), track_names)
  data$source_id <- factor(data$source_id, levels = y_levels)
  data$sample <- factor(data$sample, levels = column_order)
  fill_name <- .activity_score_label(
    prepared$result, statistic, statistic_column, scaled = scale == "row"
  )
  plot <- ggplot2::ggplot(data, ggplot2::aes(
    x = .data[["sample"]], y = .data[["source_id"]], fill = .data[["value"]]
  )) +
    ggplot2::geom_tile() +
    ggplot2::scale_fill_gradient2(
      low = "#3B4CC0", mid = "#F7F7F7", high = "#B40426",
      midpoint = 0, name = fill_name
    ) +
    ggplot2::labs(
      x = "Sample", y = "Regulator",
      alt = "A heatmap of inferred regulator activity scores across samples."
    ) +
    theme_bulkmae() +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)
    )
  if (!is.null(tracks)) {
    plot <- .plot_activity_add_annotation(plot, tracks, column_order, y_levels)
  }
  .plot_with_dimensions(
    plot,
    width = .plot_clamped_dimension(
      ncol(matrix), base = 5, per_item = 0.5, minimum = 10, maximum = 18
    ),
    height = .plot_clamped_dimension(
      nrow(matrix) + length(track_names),
      base = 4, per_item = 0.28, minimum = 8, maximum = 20
    )
  )
}

#' Plot ranked regulator activities
#'
#' One value is drawn per regulator. Long-format sample-level results are
#' summarised as the mean activity across samples. Contrast tables from
#' [activity_contrast()] are plotted using `delta` when `value` is omitted
#' or set to `"delta"`. The function does not select regulators from
#' *p*-values; pass `sources` to choose the displayed set and order.
#'
#' @param result A long-format activity table or an [activity_contrast()]
#'   table.
#' @param sources Optional regulator identifiers to display. The default
#'   uses every source in `result`.
#' @inheritParams activity_matrix
#'
#' @return An unprinted standard ggplot object carrying recommended physical
#'   dimensions for [plot_save()].
#' @family plotting
#' @export
plot_activity_rank <- function(
    result,
    sources = NULL,
    value = NULL,
    source = "source",
    sample = "condition",
    statistic = NULL,
    statistic_column = "statistic"
) {
  table <- .plot_activity_rank_table(
    result, sources, value, source, sample, statistic, statistic_column
  )
  table$sign <- factor(
    ifelse(table$value > 0, "Positive", ifelse(table$value < 0, "Negative", "Zero")),
    levels = c("Negative", "Zero", "Positive")
  )
  table$source <- factor(table$source, levels = table$source[order(table$value)])
  x_label <- table$x_label[[1L]]
  plot <- ggplot2::ggplot(table, ggplot2::aes(
    x = .data[["value"]], y = .data[["source"]]
  )) +
    ggplot2::geom_vline(
      xintercept = 0, colour = "#777777", linewidth = 0.3
    ) +
    ggplot2::geom_segment(
      ggplot2::aes(x = 0, xend = .data[["value"]], yend = .data[["source"]]),
      colour = "#777777", linewidth = 0.4
    ) +
    ggplot2::geom_point(
      ggplot2::aes(colour = .data[["sign"]]),
      size = 1.8
    ) +
    ggplot2::scale_colour_manual(
      values = c(
        Negative = .bulkmae_colours[["blue"]],
        Zero = .bulkmae_colours[["grey"]],
        Positive = .bulkmae_colours[["orange"]]
      ),
      breaks = c("Negative", "Zero", "Positive"),
      drop = TRUE,
      name = "Sign"
    ) +
    ggplot2::labs(
      x = x_label, y = "Regulator",
      alt = "A ranked lollipop plot of inferred regulator activity scores."
    ) +
    theme_bulkmae() +
    ggplot2::theme(
      legend.position = "top",
      legend.direction = "horizontal"
    )
  .plot_with_dimensions(
    plot,
    width = 8.5,
    height = .plot_clamped_dimension(
      nrow(table), base = 4, per_item = 0.35, minimum = 6, maximum = 16
    )
  )
}

#' Plot sample-level activity by phenotype
#'
#' Draws already-inferred per-sample activity scores for an explicit set of
#' regulators. Points are the samples; boxes or violins summarise the
#' phenotype groups supplied by the caller. The figure states that only this
#' caller-supplied subset is shown; it is not a *p*-value top-N.
#'
#' @param result A data-frame-like activity result.
#' @param sources Regulator identifiers to display. Required and not chosen
#'   from *p*-values.
#' @param sample_data Optional sample metadata whose row names cover every
#'   plotted sample. Required when `group` is a column name.
#' @param group A metadata column in `sample_data` or a complete sample-named
#'   vector of phenotype labels.
#' @param geom `"boxplot"` or `"violin"`. Both include jittered sample points.
#' @inheritParams activity_matrix
#'
#' @return An unprinted standard ggplot object carrying recommended physical
#'   dimensions for [plot_save()].
#' @family plotting
#' @export
plot_activity_sample <- function(
    result,
    sources,
    sample_data = NULL,
    group,
    geom = c("boxplot", "violin"),
    value = "score",
    source = "source",
    sample = "condition",
    statistic = NULL,
    statistic_column = "statistic"
) {
  geom <- match.arg(geom)
  if (is.null(sources) || !length(sources)) {
    stop("`sources` must explicitly contain at least one regulator.", call. = FALSE)
  }
  prepared <- .plot_activity_matrix(
    result, sources, statistic, value, source, sample, statistic_column
  )
  matrix <- prepared$matrix
  groups <- .plot_mae_annotation(
    group, sample_data, colnames(matrix), "group", "sample"
  )
  data <- expand.grid(
    source = rownames(matrix),
    sample = colnames(matrix),
    stringsAsFactors = FALSE
  )
  data$value <- matrix[cbind(
    match(data$source, rownames(matrix)),
    match(data$sample, colnames(matrix))
  )]
  data$group <- unname(groups[as.character(data$sample)])
  data$source <- factor(data$source, levels = sources)
  if (is.factor(groups)) {
    data$group <- factor(as.character(data$group), levels = levels(droplevels(groups)))
  }
  colours <- .plot_discrete_values(data$group)
  y_label <- .activity_score_label(prepared$result, statistic, statistic_column)
  group_label <- if (is.character(group) && length(group) == 1L && is.null(names(group))) {
    group
  } else {
    "Group"
  }
  plot <- ggplot2::ggplot(data, ggplot2::aes(
    x = .data[["group"]], y = .data[["value"]],
    colour = .data[["group"]], fill = .data[["group"]]
  ))
  plot <- if (geom == "violin") {
    plot + ggplot2::geom_violin(alpha = 0.22, colour = NA, linewidth = 0.3)
  } else {
    plot + ggplot2::geom_boxplot(
      width = 0.45, outlier.shape = NA, alpha = 0.22, linewidth = 0.35
    )
  }
  plot <- plot +
    ggplot2::geom_jitter(
      width = 0.12, height = 0, size = 1.3, alpha = 0.9, stroke = 0
    ) +
    ggplot2::scale_colour_manual(values = colours, name = group_label) +
    ggplot2::scale_fill_manual(values = colours, name = group_label) +
    ggplot2::labs(
      x = group_label, y = y_label,
      caption = paste0(
        "Showing ", length(sources),
        if (length(sources) == 1L) {
          " explicitly selected regulator"
        } else {
          " explicitly selected regulators"
        },
        "; not a p-value top-N."
      ),
      alt = paste(
        "Sample-level inferred activity scores grouped by phenotype",
        "for an explicit regulator subset."
      )
    ) +
    theme_bulkmae() +
    ggplot2::theme(
      legend.position = "top",
      legend.direction = "horizontal"
    )
  if (length(sources) > 1L) {
    plot <- plot + ggplot2::facet_wrap(
      ggplot2::vars(.data[["source"]]),
      scales = "free_y"
    )
  }
  .plot_with_dimensions(
    plot,
    width = if (length(sources) > 1L) {
      .plot_clamped_dimension(
        min(length(sources), 3L), base = 7, per_item = 2.2, minimum = 8.5, maximum = 16
      )
    } else {
      8.5
    },
    height = .plot_clamped_dimension(
      ceiling(length(sources) / 2), base = 5.5, per_item = 3.2, minimum = 6.5, maximum = 16
    )
  )
}

#' Plot a two-group activity contrast volcano
#'
#' Expects the table returned by [activity_contrast()]. The x axis is the
#' mean activity difference; the y axis is `-log10` of the two-group
#' *p*-value from that helper. Colour is the sign of that activity
#' difference (Up / Down / Not significant). A dashed horizontal line marks
#' nominal *p* = 0.05 on the same raw-*p* scale as the y axis. BH FDR is
#' not a y-coordinate and is not used to recolour points.
#'
#' @param result A data frame returned by [activity_contrast()].
#' @param fdr Accepted for interface parity with [plot_de_volcano()]. It is
#'   validated but does not recolour points; the y axis already shows the
#'   two-group raw *p*-value.
#' @param min_abs_effect Minimum absolute activity difference used to
#'   separate Up/Down from Not significant.
#' @param label_sources Regulator identifiers to label; none are selected
#'   automatically.
#' @param source_labels Optional complete source-named labels.
#' @param symmetric_x Centre the x axis on zero with equal negative and
#'   positive limits.
#' @param source,value,p_value,adjusted_p_value Column names in `result`.
#'
#' @return An unprinted standard ggplot object carrying recommended physical
#'   dimensions for [plot_save()].
#' @family plotting
#' @export
plot_activity_volcano <- function(
    result,
    fdr = 0.05,
    min_abs_effect = 0,
    label_sources = NULL,
    source_labels = NULL,
    symmetric_x = TRUE,
    source = "source",
    value = "delta",
    p_value = "p_value",
    adjusted_p_value = "adjusted_p_value"
) {
  .plot_assert_flag(symmetric_x, "symmetric_x")
  .differential_assert_probability(fdr, "fdr")
  if (
    !is.numeric(min_abs_effect) || length(min_abs_effect) != 1L ||
      is.na(min_abs_effect) || !is.finite(min_abs_effect) || min_abs_effect < 0
  ) {
    stop("`min_abs_effect` must be one finite, non-negative number.", call. = FALSE)
  }
  result <- as.data.frame(result, optional = TRUE)
  .de_assert_column(result, source, "source")
  .de_assert_column(result, value, "value")
  .de_assert_column(result, p_value, "p_value")
  .de_assert_column(result, adjusted_p_value, "adjusted_p_value")
  source_ids <- as.character(result[[source]])
  .plot_assert_ids(source_ids, "Activity source identifiers")
  effect <- result[[value]]
  raw_p <- result[[p_value]]
  adj_p <- result[[adjusted_p_value]]
  if (!is.numeric(effect) || any(!is.finite(effect))) {
    stop("Activity differences must contain finite numeric values.", call. = FALSE)
  }
  if (all(is.na(raw_p))) {
    stop("Raw p-values are unavailable in `result`.", call. = FALSE)
  }
  .plot_enrichment_probabilities(raw_p, "Raw p-values", allow_all_missing = TRUE)
  .plot_enrichment_probabilities(adj_p, "Adjusted p-values", allow_all_missing = TRUE)
  positive <- raw_p[is.finite(raw_p) & raw_p > 0]
  floor <- if (length(positive)) min(positive) / 10 else .Machine$double.xmin
  table <- data.frame(
    source = source_ids,
    effect = as.numeric(effect),
    minus_log10_p = -log10(pmax(raw_p, floor, na.rm = FALSE)),
    stringsAsFactors = FALSE
  )
  signed <- is.finite(table$effect) & abs(table$effect) >= min_abs_effect
  direction <- rep("Not significant", nrow(table))
  direction[signed & table$effect < 0] <- "Down"
  direction[signed & table$effect > 0] <- "Up"
  table$direction <- factor(direction, levels = c("Down", "Not significant", "Up"))
  table$label <- .plot_feature_labels(
    table$source,
    label_sources,
    source_labels,
    data.frame(feature_id = table$source, stringsAsFactors = FALSE)
  )
  x_label <- .plot_activity_delta_label(result)
  p_cutoff <- 0.05
  p_cutoff_y <- -log10(p_cutoff)
  plot <- ggplot2::ggplot(table, ggplot2::aes(
    x = .data[["effect"]], y = .data[["minus_log10_p"]]
  )) +
    ggplot2::geom_hline(
      yintercept = p_cutoff_y,
      linetype = 2,
      colour = "#777777"
    ) +
    ggplot2::geom_vline(
      xintercept = c(-min_abs_effect, min_abs_effect),
      linetype = 2, colour = "#777777"
    ) +
    ggplot2::geom_point(
      mapping = ggplot2::aes(colour = .data[["direction"]]),
      alpha = 0.8,
      size = 1.6,
      na.rm = TRUE
    ) +
    .plot_de_scale() +
    ggplot2::labs(
      x = x_label, y = expression(-log[10](italic(p))),
      colour = "Direction",
      caption = "Dashed line: nominal p = 0.05.",
      alt = paste(
        "A volcano plot of two-group activity differences and raw p-values,",
        "coloured by the sign of the activity difference, with a dashed",
        "nominal p = 0.05 line."
      )
    ) +
    .plot_de_theme()
  if (symmetric_x) {
    x_limit <- max(c(abs(table$effect[is.finite(table$effect)]), min_abs_effect))
    if (x_limit == 0) x_limit <- 1
    plot <- plot + ggplot2::coord_cartesian(xlim = c(-x_limit, x_limit))
  }
  labelled <- !is.na(table$label) &
    is.finite(table$effect) & is.finite(table$minus_log10_p)
  if (any(labelled)) {
    y_values <- table$minus_log10_p[is.finite(table$minus_log10_p)]
    y_range <- diff(range(y_values))
    nudge_y <- if (is.finite(y_range) && y_range > 0) y_range * 0.025 else 0
    plot <- plot + ggplot2::geom_text(
      data = table[labelled, , drop = FALSE],
      mapping = ggplot2::aes(
        x = .data[["effect"]],
        y = .data[["minus_log10_p"]],
        label = .data[["label"]]
      ),
      inherit.aes = FALSE,
      nudge_y = nudge_y,
      size = .bulkmae_text_size_pt,
      size.unit = "pt",
      check_overlap = TRUE,
      show.legend = FALSE,
      na.rm = TRUE
    )
  }
  .plot_with_dimensions(plot, width = 8.5, height = 7)
}

.plot_activity_matrix <- function(
    result,
    sources,
    statistic,
    value,
    source,
    sample,
    statistic_column
) {
  result <- as.data.frame(result, optional = TRUE)
  result <- .activity_filter_statistic(result, statistic, statistic_column)
  .activity_warn_enrichment_statistic(result, statistic_column)
  matrix <- activity_matrix(
    result,
    value = value,
    source = source,
    sample = sample,
    statistic = statistic,
    statistic_column = statistic_column
  )
  if (!is.null(sources)) {
    .plot_assert_ids(as.character(sources), "`sources`")
    unknown <- setdiff(sources, rownames(matrix))
    if (length(unknown)) {
      stop("Unknown `sources`: ", paste(unknown, collapse = ", "), ".", call. = FALSE)
    }
    matrix <- matrix[as.character(sources), , drop = FALSE]
  }
  list(matrix = matrix, result = result)
}

.plot_activity_annotation_tracks <- function(annotation, sample_data, sample_ids) {
  if (is.null(annotation)) return(NULL)
  named_vector <- !is.null(names(annotation))
  if (named_vector) {
    values <- .plot_mae_annotation(
      annotation, sample_data, sample_ids, "annotation", "sample"
    )
    return(list(list(name = "annotation", values = values)))
  }
  if (!is.character(annotation) || !length(annotation) || any(!nzchar(annotation))) {
    stop(
      "`annotation` must be metadata column name(s) or a sample-named vector.",
      call. = FALSE
    )
  }
  if (anyDuplicated(annotation)) {
    stop("`annotation` track names must be unique.", call. = FALSE)
  }
  lapply(annotation, function(column) {
    list(
      name = column,
      values = .plot_mae_annotation(
        column, sample_data, sample_ids, "annotation", "sample"
      )
    )
  })
}

.plot_activity_add_annotation <- function(plot, tracks, column_order, y_levels) {
  rows <- lapply(tracks, function(track) {
    values <- as.character(track$values[column_order])
    legend <- if (length(tracks) == 1L) {
      values
    } else {
      paste(track$name, values, sep = ": ")
    }
    pal <- .plot_discrete_values(factor(values, levels = unique(values)))
    fill_hex <- unname(pal[values])
    data.frame(
      sample = factor(column_order, levels = column_order),
      y = factor(track$name, levels = y_levels),
      legend = legend,
      fill_hex = fill_hex,
      stringsAsFactors = FALSE
    )
  })
  ann <- do.call(rbind, rows)
  colours <- stats::setNames(ann$fill_hex, ann$legend)
  colours <- colours[!duplicated(names(colours))]
  legend_name <- if (length(tracks) == 1L) tracks[[1L]]$name else "Phenotype"
  plot +
    ggplot2::geom_tile(
      data = ann,
      mapping = ggplot2::aes(
        x = .data[["sample"]], y = .data[["y"]]
      ),
      fill = ann$fill_hex,
      height = 0.7,
      inherit.aes = FALSE
    ) +
    ggplot2::geom_point(
      data = data.frame(
        sample = ann$sample[[1L]],
        y = ann$y[[1L]],
        legend = names(colours),
        stringsAsFactors = FALSE
      ),
      mapping = ggplot2::aes(
        x = .data[["sample"]], y = .data[["y"]], colour = .data[["legend"]]
      ),
      alpha = 0,
      size = 2,
      inherit.aes = FALSE
    ) +
    ggplot2::scale_colour_manual(values = colours, name = legend_name) +
    ggplot2::guides(
      colour = ggplot2::guide_legend(
        override.aes = list(alpha = 1, shape = 15, size = 3)
      )
    ) +
    ggplot2::scale_y_discrete(limits = y_levels)
}

.plot_activity_rank_table <- function(
    result,
    sources,
    value,
    source,
    sample,
    statistic,
    statistic_column
) {
  result <- as.data.frame(result, optional = TRUE)
  if (is.null(value)) {
    value <- if ("delta" %in% names(result) && !"score" %in% names(result)) {
      "delta"
    } else if ("score" %in% names(result)) {
      "score"
    } else if ("delta" %in% names(result)) {
      "delta"
    } else {
      "score"
    }
  }
  result <- .activity_prepare_result(
    result,
    value = value,
    source = source,
    sample = sample,
    statistic = statistic,
    statistic_column = statistic_column,
    require_sample = value != "delta" && sample %in% names(result)
  )
  if (value != "delta") {
    .activity_warn_enrichment_statistic(result, statistic_column)
  }
  source_ids <- as.character(result[[source]])
  if (anyNA(source_ids) || any(!nzchar(source_ids))) {
    stop("Activity source identifiers cannot be missing or empty.", call. = FALSE)
  }
  scores <- result[[value]]
  if (!is.numeric(scores) || any(!is.finite(scores))) {
    stop("The activity `value` column must contain finite numeric values.", call. = FALSE)
  }
  summarised <- FALSE
  if (sample %in% names(result) && anyDuplicated(source_ids)) {
    sample_ids <- as.character(result[[sample]])
    keys <- paste(source_ids, sample_ids, sep = "\r")
    if (anyDuplicated(keys)) {
      stop("Activity results contain duplicate source/sample rows.", call. = FALSE)
    }
    scores <- tapply(scores, source_ids, mean)
    table <- data.frame(
      source = names(scores),
      value = as.numeric(scores),
      stringsAsFactors = FALSE
    )
    summarised <- TRUE
  } else {
    if (anyDuplicated(source_ids)) {
      stop(
        "Activity rank plots need one value per source or a sample-level grid.",
        call. = FALSE
      )
    }
    table <- data.frame(
      source = source_ids,
      value = as.numeric(scores),
      stringsAsFactors = FALSE
    )
  }
  if (!is.null(sources)) {
    .plot_assert_ids(as.character(sources), "`sources`")
    unknown <- setdiff(sources, table$source)
    if (length(unknown)) {
      stop("Unknown `sources`: ", paste(unknown, collapse = ", "), ".", call. = FALSE)
    }
    table <- table[match(sources, table$source), , drop = FALSE]
  }
  table$x_label <- if (identical(value, "delta")) {
    .plot_activity_delta_label(result)
  } else if (summarised) {
    .activity_score_label(
      result, statistic, statistic_column, prefix = "Mean"
    )
  } else {
    .activity_score_label(result, statistic, statistic_column)
  }
  table
}

.plot_activity_delta_label <- function(result) {
  if (all(c("target", "reference") %in% names(result))) {
    target <- unique(as.character(result$target))
    reference <- unique(as.character(result$reference))
    if (length(target) == 1L && length(reference) == 1L) {
      return(paste0("Activity difference (", target, " \u2212 ", reference, ")"))
    }
  }
  "Activity difference"
}
