#' Plot ranked LINCS connectivity scores
#'
#' Draws already-computed connectivity scores as a lollipop. The x-axis label
#' is the score column actually present in `result` (NCS, Tau, WTCS, or
#' NCSct). Negative scores are labelled `Reverse` and positive scores
#' `Mimic`; that is the sign of the connectivity statistic, not a validated
#' drug effect.
#'
#' @param result A `gessResult` or LINCS table, usually from [drug_lincs()]
#'   or [drug_lincs_table()].
#' @param n Number of already-ranked rows to draw, taking the most extreme
#'   absolute scores when the table is longer.
#' @param score Optional score column forwarded to [drug_lincs_table()].
#' @param compounds Optional compound identifiers. When supplied, rows are
#'   restricted to those names before ranking.
#'
#' @return An unprinted standard ggplot object carrying recommended physical
#'   dimensions for [plot_save()].
#' @family plotting
#' @export
plot_lincs_rank <- function(result, n = 15L, score = NULL, compounds = NULL) {
  table <- .plot_lincs_prepare(result, score, compounds)
  table <- .plot_lincs_top_absolute(table, n)
  score_column <- attr(table, "score_column", exact = TRUE)
  score_label <- attr(table, "score_label", exact = TRUE)
  table$label <- .lincs_row_label(table)
  table$label <- make.unique(table$label)
  table <- table[order(table[[score_column]], decreasing = FALSE), , drop = FALSE]
  table$label <- factor(table$label, levels = table$label)
  palette <- .plot_lincs_direction_scale()

  plot <- ggplot2::ggplot(table, ggplot2::aes(
    x = .data[[score_column]],
    y = .data[["label"]],
    colour = .data[["direction"]]
  )) +
    ggplot2::geom_vline(
      xintercept = 0,
      colour = "#BBBBBB",
      linewidth = 0.3
    ) +
    ggplot2::geom_segment(
      mapping = ggplot2::aes(xend = 0, yend = .data[["label"]]),
      linewidth = 0.4
    ) +
    ggplot2::geom_point(size = 1.6) +
    ggplot2::scale_colour_manual(
      values = palette,
      breaks = names(palette),
      drop = FALSE
    ) +
    ggplot2::labs(
      x = score_label,
      y = "Compound (cell)",
      colour = "Score sign",
      alt = paste(
        "Lollipop ranking of LINCS connectivity scores, with reverse",
        "signatures in blue and mimic signatures in orange."
      )
    ) +
    theme_bulkmae() +
    ggplot2::theme(legend.position = "top", legend.direction = "horizontal")

  .plot_with_dimensions(
    plot,
    width = 9,
    height = .plot_clamped_dimension(
      nrow(table),
      base = 3.2,
      per_item = 0.38,
      minimum = 6,
      maximum = 16
    )
  )
}

#' Plot LINCS connectivity scores as a compound-by-cell heatmap
#'
#' Each tile is one already-computed compound/cell score. The caller supplies
#' the result rows to display; the constructor does not pick "top hits" from a
#' larger unfiltered table unless `compounds` or `cells` is passed.
#'
#' @inheritParams plot_lincs_rank
#' @param cells Optional cell-line identifiers used to subset columns.
#' @param cluster_rows,cluster_columns Reorder compounds or cells by
#'   hierarchical clustering of the displayed scores. Dendrograms are not drawn.
#'
#' @return An unprinted standard ggplot object carrying recommended physical
#'   dimensions for [plot_save()].
#' @family plotting
#' @export
plot_lincs_heatmap <- function(
    result,
    compounds = NULL,
    cells = NULL,
    score = NULL,
    cluster_rows = FALSE,
    cluster_columns = FALSE
) {
  .plot_assert_flag(cluster_rows, "cluster_rows")
  .plot_assert_flag(cluster_columns, "cluster_columns")
  table <- .plot_lincs_prepare(result, score, compounds)
  cell_column <- .lincs_cell_column(table, required = TRUE)
  if (!is.null(cells)) {
    .plot_assert_ids(as.character(cells), "`cells`")
    table <- table[as.character(table[[cell_column]]) %in% cells, , drop = FALSE]
    if (!nrow(table)) {
      stop("None of the requested `cells` are present in `result`.", call. = FALSE)
    }
  }
  score_column <- attr(table, "score_column", exact = TRUE)
  score_label <- attr(table, "score_label", exact = TRUE)
  compound_column <- .lincs_compound_column(table)
  keys <- paste(
    as.character(table[[compound_column]]),
    as.character(table[[cell_column]]),
    sep = "\t"
  )
  if (anyDuplicated(keys)) {
    stop(
      "Heatmap cells must be unique compound and cell pairs.",
      call. = FALSE
    )
  }
  compounds_present <- unique(as.character(table[[compound_column]]))
  cells_present <- unique(as.character(table[[cell_column]]))
  matrix <- matrix(
    NA_real_,
    nrow = length(compounds_present),
    ncol = length(cells_present),
    dimnames = list(compounds_present, cells_present)
  )
  matrix[
    cbind(
      match(as.character(table[[compound_column]]), compounds_present),
      match(as.character(table[[cell_column]]), cells_present)
    )
  ] <- as.numeric(table[[score_column]])
  if (any(!is.finite(matrix))) {
    stop("Heatmap scores must be finite for every displayed tile.", call. = FALSE)
  }
  row_order <- rownames(matrix)
  column_order <- colnames(matrix)
  if (cluster_rows && nrow(matrix) > 1L) {
    row_order <- row_order[stats::hclust(stats::dist(matrix))$order]
  }
  if (cluster_columns && ncol(matrix) > 1L) {
    column_order <- column_order[stats::hclust(stats::dist(t(matrix)))$order]
  }
  data <- expand.grid(
    compound = row_order,
    cell = column_order,
    stringsAsFactors = FALSE
  )
  data$value <- matrix[cbind(
    match(data$compound, rownames(matrix)),
    match(data$cell, colnames(matrix))
  )]
  data$compound <- factor(data$compound, levels = rev(row_order))
  data$cell <- factor(data$cell, levels = column_order)

  plot <- ggplot2::ggplot(data, ggplot2::aes(
    x = .data[["cell"]],
    y = .data[["compound"]],
    fill = .data[["value"]]
  )) +
    ggplot2::geom_tile() +
    ggplot2::scale_fill_gradient2(
      low = "#0072B2",
      mid = "#F7F7F7",
      high = "#D55E00",
      midpoint = 0,
      name = score_label
    ) +
    ggplot2::labs(
      x = "Cell",
      y = "Compound",
      alt = paste(
        "Heatmap of LINCS connectivity scores for selected compounds",
        "across cell lines."
      )
    ) +
    theme_bulkmae() +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)
    )

  .plot_with_dimensions(
    plot,
    width = .plot_clamped_dimension(
      ncol(matrix),
      base = 4.5,
      per_item = 0.7,
      minimum = 7.5,
      maximum = 14
    ),
    height = .plot_clamped_dimension(
      nrow(matrix),
      base = 3.6,
      per_item = 0.42,
      minimum = 6.5,
      maximum = 16
    )
  )
}

#' Plot a waterfall of LINCS connectivity scores
#'
#' Bars run from zero to the already-computed score and are ordered from the
#' strongest reverse signature to the strongest mimic. The axis name follows
#' the resolved score column.
#'
#' @inheritParams plot_lincs_rank
#'
#' @return An unprinted standard ggplot object carrying recommended physical
#'   dimensions for [plot_save()].
#' @family plotting
#' @export
plot_lincs_waterfall <- function(result, n = NULL, score = NULL, compounds = NULL) {
  table <- .plot_lincs_prepare(result, score, compounds)
  if (!is.null(n)) {
    table <- .plot_lincs_top_absolute(table, n)
  }
  score_column <- attr(table, "score_column", exact = TRUE)
  score_label <- attr(table, "score_label", exact = TRUE)
  table$label <- make.unique(.lincs_row_label(table))
  table <- table[order(table[[score_column]], decreasing = FALSE), , drop = FALSE]
  table$label <- factor(table$label, levels = table$label)
  palette <- .plot_lincs_direction_scale()

  plot <- ggplot2::ggplot(table, ggplot2::aes(
    x = .data[["label"]],
    y = .data[[score_column]],
    fill = .data[["direction"]]
  )) +
    ggplot2::geom_hline(
      yintercept = 0,
      colour = "#BBBBBB",
      linewidth = 0.3
    ) +
    ggplot2::geom_col(width = 0.78) +
    ggplot2::scale_fill_manual(
      values = palette,
      breaks = names(palette),
      drop = FALSE
    ) +
    ggplot2::labs(
      x = "Compound (cell)",
      y = score_label,
      fill = "Score sign",
      alt = paste(
        "Waterfall of LINCS connectivity scores ordered from reverse",
        "to mimic signatures."
      )
    ) +
    theme_bulkmae() +
    ggplot2::theme(
      legend.position = "top",
      legend.direction = "horizontal",
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)
    )

  .plot_with_dimensions(
    plot,
    width = .plot_clamped_dimension(
      nrow(table),
      base = 4.8,
      per_item = 0.42,
      minimum = 8,
      maximum = 18
    ),
    height = 7
  )
}

#' Plot query up/down overlap counts for LINCS hits
#'
#' Uses the `N_upset` and `N_downset` columns returned by
#' [signatureSearch::gess_lincs()]. Those counts are the query genes that
#' overlap the extreme ranks of each reference signature. They are not an
#' independent enrichment test.
#'
#' @inheritParams plot_lincs_rank
#'
#' @return An unprinted standard ggplot object carrying recommended physical
#'   dimensions for [plot_save()].
#' @family plotting
#' @export
plot_lincs_overlap <- function(result, n = 12L, score = NULL, compounds = NULL) {
  table <- .plot_lincs_prepare(result, score, compounds)
  .plot_require_columns(table, c("N_upset", "N_downset"), "result")
  .plot_assert_finite_numeric(table$N_upset, "`result$N_upset`")
  .plot_assert_finite_numeric(table$N_downset, "`result$N_downset`")
  if (any(table$N_upset < 0) || any(table$N_downset < 0)) {
    stop("`N_upset` and `N_downset` must be non-negative.", call. = FALSE)
  }
  table <- .plot_lincs_top_absolute(table, n)
  table$label <- make.unique(.lincs_row_label(table))
  score_column <- attr(table, "score_column", exact = TRUE)
  table <- table[order(table[[score_column]], decreasing = FALSE), , drop = FALSE]
  data <- data.frame(
    label = factor(rep(table$label, 2L), levels = table$label),
    set = factor(
      rep(c("Query up", "Query down"), each = nrow(table)),
      levels = c("Query up", "Query down")
    ),
    count = c(table$N_upset, table$N_downset),
    stringsAsFactors = FALSE
  )

  plot <- ggplot2::ggplot(data, ggplot2::aes(
    x = .data[["count"]],
    y = .data[["label"]],
    fill = .data[["set"]]
  )) +
    ggplot2::geom_col(
      position = ggplot2::position_dodge(width = 0.72),
      width = 0.68
    ) +
    ggplot2::scale_fill_manual(
      values = c(
        `Query up` = .bulkmae_colours[["orange"]],
        `Query down` = .bulkmae_colours[["blue"]]
      ),
      breaks = c("Query up", "Query down")
    ) +
    ggplot2::labs(
      x = "Overlapping query genes",
      y = "Compound (cell)",
      fill = "Query set",
      alt = paste(
        "Grouped bars of query up-set and down-set gene counts overlapping",
        "each ranked LINCS signature."
      )
    ) +
    theme_bulkmae() +
    ggplot2::theme(legend.position = "top", legend.direction = "horizontal")

  .plot_with_dimensions(
    plot,
    width = 9,
    height = .plot_clamped_dimension(
      nrow(table),
      base = 3.2,
      per_item = 0.42,
      minimum = 6,
      maximum = 16
    )
  )
}

.plot_lincs_prepare <- function(result, score, compounds) {
  table <- drug_lincs_table(result, score = score)
  if (is.null(compounds)) {
    return(table)
  }
  .plot_assert_ids(as.character(compounds), "`compounds`")
  compound_column <- .lincs_compound_column(table)
  keep <- as.character(table[[compound_column]]) %in% compounds
  if (!any(keep)) {
    stop(
      "None of the requested `compounds` are present in `result`.",
      call. = FALSE
    )
  }
  .lincs_restore_score_attrs(table[keep, , drop = FALSE], table)
}

.plot_lincs_top_absolute <- function(table, n) {
  if (
    !is.numeric(n) || length(n) != 1L || !is.finite(n) || n < 1L ||
      n != as.integer(n)
  ) {
    stop("`n` must be a positive integer.", call. = FALSE)
  }
  n <- min(as.integer(n), nrow(table))
  score_column <- attr(table, "score_column", exact = TRUE)
  values <- abs(as.numeric(table[[score_column]]))
  keep <- utils::head(order(values, decreasing = TRUE), n)
  .lincs_restore_score_attrs(table[sort(keep), , drop = FALSE], table)
}

.lincs_restore_score_attrs <- function(table, source) {
  attr(table, "score_column") <- attr(source, "score_column", exact = TRUE)
  attr(table, "score_label") <- attr(source, "score_label", exact = TRUE)
  attr(table, "bulkmae_lincs_source") <- attr(
    source,
    "bulkmae_lincs_source",
    exact = TRUE
  )
  table
}

.plot_lincs_direction_scale <- function() {
  c(
    Reverse = .bulkmae_colours[["blue"]],
    Unrelated = .bulkmae_colours[["grey"]],
    Mimic = .bulkmae_colours[["orange"]]
  )
}
