#' bulkMAE plot theme
#'
#' A compact, colour-blind-friendly theme used by the package's plotting
#' helpers. Its default text sizes are designed for the recommended physical
#' dimensions used by [plot_save()]. Add ordinary ggplot2 themes or scales to
#' a returned plot to replace any part of it.
#'
#' @param base_size Base font size.
#' @param base_family Base font family.
#'
#' @return A ggplot2 theme object.
#' @family plotting
#' @export
theme_bulkmae <- function(base_size = 6, base_family = "") {
  ggplot2::theme_minimal(base_size = base_size, base_family = base_family) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major.x = ggplot2::element_blank(),
      axis.title = ggplot2::element_text(colour = "#222222"),
      axis.text = ggplot2::element_text(colour = "#333333"),
      strip.text = ggplot2::element_text(face = "bold"),
      legend.position = "right",
      plot.title.position = "plot"
    )
}

#' Save a bulkMAE plot at its recommended physical size
#'
#' Plot themes control text and visual styling, but cannot control the physical
#' size of a graphics device. Each bulkMAE `plot_*()` helper therefore attaches
#' a recommended width and height in centimetres to its otherwise standard
#' ggplot object. `plot_save()` reads that recommendation and delegates to
#' [ggplot2::ggsave()] with `scale = 1` by default. Explicit `width` or `height`
#' values override either recommendation.
#'
#' PDF output uses [grDevices::cairo_pdf()] automatically when Cairo graphics
#' are available. Dimensions remain attached after ordinary ggplot2 `+`
#' operations, so themes, labels, and scales can be changed before saving.
#'
#' The RStudio plot pane and knitr graphics chunks create their graphics
#' devices before a plot is drawn, so they do not consume the attached
#' recommendation. Use `plot_save()` for an exact final output size; set
#' `fig.width` and `fig.height` separately when a rendered document must use a
#' particular preview size.
#'
#' @param filename Output file name, including an extension understood by
#'   [ggplot2::ggsave()].
#' @param plot A ggplot object. Plots not produced by bulkMAE must supply both
#'   `width` and `height`.
#' @param width,height Optional positive physical dimensions. When omitted,
#'   use the recommendation attached by a bulkMAE plotting helper.
#' @param units Physical units for `width` and `height`.
#' @param dpi Raster resolution passed to [ggplot2::ggsave()].
#' @param scale Multiplicative scale passed to [ggplot2::ggsave()]. Keep the
#'   default of one when exact final dimensions matter.
#' @param device Optional graphics device. When `NULL`, infer it from
#'   `filename`, with Cairo PDF preferred when available.
#' @param bg Background colour passed to [ggplot2::ggsave()].
#' @param ... Additional arguments passed to [ggplot2::ggsave()].
#'
#' @return `filename`, invisibly.
#' @family plotting
#' @export
plot_save <- function(
    filename,
    plot,
    width = NULL,
    height = NULL,
    units = c("cm", "in", "mm"),
    dpi = 300,
    scale = 1,
    device = NULL,
    bg = NULL,
    ...
) {
  .assert_scalar_character(filename, "filename")
  if (!inherits(plot, "ggplot")) {
    stop("`plot` must be a ggplot object.", call. = FALSE)
  }
  units <- match.arg(units)
  dimensions <- attr(plot, "bulkmae_dimensions", exact = TRUE)
  if (is.null(width)) {
    if (is.null(dimensions)) {
      stop(
        "`width` and `height` are required for plots not produced by bulkMAE.",
        call. = FALSE
      )
    }
    width <- .plot_convert_cm(dimensions$width, units)
  }
  if (is.null(height)) {
    if (is.null(dimensions)) {
      stop(
        "`width` and `height` are required for plots not produced by bulkMAE.",
        call. = FALSE
      )
    }
    height <- .plot_convert_cm(dimensions$height, units)
  }
  .plot_assert_positive_number(width, "width")
  .plot_assert_positive_number(height, "height")
  .plot_assert_positive_number(scale, "scale")
  if (is.null(device) && identical(tolower(tools::file_ext(filename)), "pdf") &&
      capabilities("cairo")) {
    device <- grDevices::cairo_pdf
  }
  ggplot2::ggsave(
    filename = filename,
    plot = plot,
    device = device,
    scale = scale,
    width = width,
    height = height,
    units = units,
    dpi = dpi,
    bg = bg,
    ...
  )
  invisible(filename)
}

#' Plot sample library QC metrics
#'
#' @param metrics A data frame returned by [qc_library()].
#'
#' @return An unprinted standard ggplot object carrying recommended physical
#'   dimensions for [plot_save()].
#' @family plotting
#' @export
plot_qc_library <- function(metrics) {
  metrics <- as.data.frame(metrics, optional = TRUE)
  required <- c("sample", "library_size", "detected_features", "zero_fraction")
  .plot_require_columns(metrics, required, "metrics")
  .plot_assert_ids(as.character(metrics$sample), "`metrics$sample`")
  for (column in required[-1L]) {
    .plot_assert_finite_numeric(metrics[[column]], paste0("`metrics$", column, "`"))
  }
  if (any(metrics$library_size < 0) || any(metrics$detected_features < 0) ||
      any(metrics$zero_fraction < 0 | metrics$zero_fraction > 1)) {
    stop("Library QC metrics are outside their valid ranges.", call. = FALSE)
  }
  sample_levels <- as.character(metrics$sample)
  data <- data.frame(
    sample = factor(rep(sample_levels, 3L), levels = sample_levels),
    metric = factor(
      rep(c("Library size", "Detected features", "Zero fraction"),
          each = nrow(metrics)),
      levels = c("Library size", "Detected features", "Zero fraction")
    ),
    value = c(metrics$library_size, metrics$detected_features, metrics$zero_fraction)
  )
  plot <- ggplot2::ggplot(data, ggplot2::aes(x = .data[["sample"]], y = .data[["value"]])) +
    ggplot2::geom_col(fill = .bulkmae_colours[["blue"]], width = 0.72) +
    ggplot2::facet_wrap(ggplot2::vars(.data[["metric"]]), scales = "free_y", ncol = 1L) +
    ggplot2::labs(
      x = "Sample", y = "Value",
      alt = "Three aligned bar charts of library size, detected features, and zero fraction by sample."
    ) +
    theme_bulkmae() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
  .plot_with_dimensions(
    plot,
    width = .plot_clamped_dimension(nrow(metrics), base = 4.5, per_item = 0.5,
                                    minimum = 8.5, maximum = 18),
    height = 11.5
  )
}

#' Plot a sample embedding
#'
#' @param embedding A `prcomp` object, limma MDS object, sample-by-dimension
#'   matrix, or `Rtsne` result.
#' @param sample_data Optional sample metadata with unique sample row names
#'   matching the embedding exactly.
#' @param axes Two one-based embedding dimensions.
#' @param colour,shape,label Optional column names in `sample_data` mapped to
#'   colour, shape, and text labels. Shape values must be discrete. The default
#'   discrete colour scale supports up to eight levels; for more levels, replace
#'   it with a standard ggplot2 colour scale before drawing or building.
#'
#' @return An unprinted standard ggplot object carrying recommended physical
#'   dimensions for [plot_save()].
#' @family plotting
#' @export
plot_embedding <- function(
    embedding,
    sample_data = NULL,
    axes = c(1L, 2L),
    colour = NULL,
    shape = NULL,
    label = NULL
) {
  axes <- .plot_axes(axes)
  extracted <- .plot_embedding_data(embedding, axes)
  data <- extracted$data
  metadata_columns <- c(colour, shape, label)
  metadata_columns <- metadata_columns[!vapply(metadata_columns, is.null, logical(1))]
  metadata_columns <- unique(metadata_columns)
  if (length(metadata_columns) && is.null(sample_data)) {
    stop("`sample_data` is required for colour, shape, or label mappings.", call. = FALSE)
  }
  if (!is.null(sample_data)) {
    sample_data <- as.data.frame(sample_data, optional = TRUE)
    .plot_align_frame(sample_data, data$sample, "sample_data")
    .plot_require_columns(sample_data, metadata_columns, "sample_data")
    conflicts <- intersect(metadata_columns, names(data))
    if (length(conflicts)) {
      stop(
        "Mapped `sample_data` columns use reserved embedding names: ",
        paste(conflicts, collapse = ", "),
        ".",
        call. = FALSE
      )
    }
    sample_data <- sample_data[data$sample, metadata_columns, drop = FALSE]
    data <- cbind(data, sample_data, row.names = NULL)
  }
  if (!is.null(colour)) data$.bulkmae_colour <- data[[colour]]
  if (!is.null(shape)) {
    shape_values <- data[[shape]]
    if (!(is.factor(shape_values) || is.character(shape_values) ||
          is.logical(shape_values))) {
      stop(
        "`sample_data$", shape, "` mapped to shape must be discrete; ",
        "convert it to a factor.",
        call. = FALSE
      )
    }
    data$.bulkmae_shape <- shape_values
  }
  mapping <- if (!is.null(colour) && !is.null(shape)) {
    ggplot2::aes(
      x = .data[["x"]], y = .data[["y"]],
      colour = .data[[".bulkmae_colour"]], shape = .data[[".bulkmae_shape"]]
    )
  } else if (!is.null(colour)) {
    ggplot2::aes(x = .data[["x"]], y = .data[["y"]], colour = .data[[".bulkmae_colour"]])
  } else if (!is.null(shape)) {
    ggplot2::aes(x = .data[["x"]], y = .data[["y"]], shape = .data[[".bulkmae_shape"]])
  } else {
    ggplot2::aes(x = .data[["x"]], y = .data[["y"]])
  }
  plot <- ggplot2::ggplot(data, mapping) +
    ggplot2::geom_point(size = 2.8, alpha = 0.9) +
    ggplot2::labs(
      x = extracted$x_label, y = extracted$y_label,
      colour = colour, shape = shape,
      alt = "A two-dimensional sample embedding with samples joined to metadata by identifier."
    ) +
    theme_bulkmae()
  if (!is.null(colour) && (is.factor(data[[colour]]) || is.character(data[[colour]]) ||
      is.logical(data[[colour]]))) {
    plot <- plot + ggplot2::scale_colour_manual(values = .plot_discrete_values(data[[colour]]))
  } else if (!is.null(colour)) {
    .plot_assert_finite_numeric(data[[colour]], paste0("`sample_data$", colour, "`"))
    plot <- plot + ggplot2::scale_colour_viridis_c()
  }
  if (!is.null(label)) {
    labels <- as.character(data[[label]])
    if (anyNA(labels)) stop("Embedding labels cannot be missing.", call. = FALSE)
    data$.bulkmae_label <- labels
    plot <- plot + ggplot2::geom_text(
      data = data,
      mapping = ggplot2::aes(label = .data[[".bulkmae_label"]]),
      nudge_y = diff(range(data$y)) * 0.025,
      size = .bulkmae_text_size_pt,
      size.unit = "pt",
      check_overlap = TRUE,
      show.legend = FALSE
    )
  }
  .plot_with_dimensions(plot, width = 8.5, height = 7)
}

#' Plot a sample correlation matrix
#'
#' @param correlation A square, sample-named correlation matrix.
#' @param cluster Reorder rows and columns using hierarchical clustering.
#' @param show_names Show sample names on both axes.
#'
#' @return An unprinted standard ggplot object carrying recommended physical
#'   dimensions for [plot_save()].
#' @family plotting
#' @export
plot_qc_correlation <- function(correlation, cluster = FALSE, show_names = TRUE) {
  correlation <- as.matrix(correlation)
  .plot_assert_flag(cluster, "cluster")
  .plot_assert_flag(show_names, "show_names")
  if (!is.numeric(correlation) || nrow(correlation) != ncol(correlation) || !nrow(correlation)) {
    stop("`correlation` must be a non-empty numeric square matrix.", call. = FALSE)
  }
  .plot_assert_ids(rownames(correlation), "Correlation row names")
  .plot_assert_ids(colnames(correlation), "Correlation column names")
  if (!setequal(rownames(correlation), colnames(correlation))) {
    stop("Correlation row and column names must contain the same samples.", call. = FALSE)
  }
  correlation <- correlation[rownames(correlation), rownames(correlation), drop = FALSE]
  if (any(!is.finite(correlation)) || any(correlation < -1 | correlation > 1)) {
    stop("`correlation` must contain finite values between -1 and 1.", call. = FALSE)
  }
  order <- rownames(correlation)
  if (cluster && length(order) > 1L) {
    order <- order[stats::hclust(stats::as.dist(1 - correlation))$order]
  }
  data <- expand.grid(row = order, column = order, stringsAsFactors = FALSE)
  data$value <- correlation[cbind(match(data$row, rownames(correlation)),
                                  match(data$column, colnames(correlation)))]
  data$row <- factor(data$row, levels = rev(order))
  data$column <- factor(data$column, levels = order)
  plot <- ggplot2::ggplot(data, ggplot2::aes(
    x = .data[["column"]], y = .data[["row"]], fill = .data[["value"]]
  )) +
    ggplot2::geom_tile() +
    ggplot2::scale_fill_gradient2(
      low = "#3B4CC0", mid = "#F7F7F7", high = "#B40426",
      midpoint = 0, limits = c(-1, 1), name = "Correlation"
    ) +
    ggplot2::coord_equal() +
    ggplot2::labs(
      x = "Sample", y = "Sample",
      alt = "A heatmap of pairwise sample correlations."
    ) +
    theme_bulkmae() +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, vjust = 1)
    )
  if (!show_names) {
    plot <- plot + ggplot2::theme(
      axis.text = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank()
    )
  }
  panel_size <- .plot_clamped_dimension(
    nrow(correlation), base = 5, per_item = 0.25,
    minimum = 8.5, maximum = 17
  )
  .plot_with_dimensions(plot, width = panel_size + 1, height = panel_size)
}

#' Plot robust sample outlier distances
#'
#' @param outliers A data frame returned by [qc_outliers()].
#'
#' @return An unprinted standard ggplot object carrying recommended physical
#'   dimensions for [plot_save()].
#' @family plotting
#' @export
plot_qc_outliers <- function(outliers) {
  outliers <- as.data.frame(outliers, optional = TRUE)
  required <- c("sample", "robust_distance", "cutoff", "flagged")
  .plot_require_columns(outliers, required, "outliers")
  .plot_assert_ids(as.character(outliers$sample), "`outliers$sample`")
  .plot_assert_finite_numeric(outliers$robust_distance, "`outliers$robust_distance`")
  .plot_assert_finite_numeric(outliers$cutoff, "`outliers$cutoff`")
  if (!is.logical(outliers$flagged) || anyNA(outliers$flagged)) {
    stop("`outliers$flagged` must be a non-missing logical vector.", call. = FALSE)
  }
  data <- outliers
  data$sample <- factor(data$sample, levels = data$sample)
  plot <- ggplot2::ggplot(data, ggplot2::aes(
    x = .data[["sample"]], y = .data[["robust_distance"]],
    colour = .data[["flagged"]]
  )) +
    ggplot2::geom_hline(yintercept = unique(data$cutoff), linetype = 2, colour = "#666666") +
    ggplot2::geom_point(size = 3) +
    ggplot2::scale_colour_manual(
      values = c(`FALSE` = .bulkmae_colours[["blue"]], `TRUE` = .bulkmae_colours[["orange"]]),
      breaks = c(FALSE, TRUE), labels = c("No", "Yes"), name = "Flagged"
    ) +
    ggplot2::labs(
      x = "Sample", y = "Robust distance",
      alt = "Robust sample distances with the screening cutoff shown as a dashed line."
    ) +
    theme_bulkmae() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
  .plot_with_dimensions(
    plot,
    width = .plot_clamped_dimension(nrow(data), base = 4.5, per_item = 0.5,
                                    minimum = 8.5, maximum = 18),
    height = 6.5
  )
}

#' Plot a differential-expression volcano plot
#'
#' @param result Input accepted by [de_table()].
#' @param coef Optional coefficient passed to [de_table()].
#' @param fdr Maximum adjusted p-value.
#' @param min_abs_effect Minimum absolute effect estimate.
#' @param label_features Feature identifiers to label; no features are selected
#'   for labels automatically.
#' @param feature_labels Optional complete feature-named labels or a named
#'   feature-label column in a data frame result.
#' @param symmetric_x Centre the x axis on zero with equal negative and
#'   positive limits. Set to `FALSE` to use ggplot2's automatic x range.
#'
#' @return An unprinted standard ggplot object carrying recommended physical
#'   dimensions for [plot_save()].
#' @family plotting
#' @export
plot_de_volcano <- function(
    result,
    coef = NULL,
    fdr = 0.05,
    min_abs_effect = 1,
    label_features = NULL,
    feature_labels = NULL,
    symmetric_x = TRUE
) {
  .plot_assert_flag(symmetric_x, "symmetric_x")
  table <- .plot_de_data(result, coef, fdr, min_abs_effect)
  if (all(is.na(table$p_value))) stop("Raw p-values are unavailable in `result`.", call. = FALSE)
  positive <- table$p_value[is.finite(table$p_value) & table$p_value > 0]
  floor <- if (length(positive)) min(positive) / 10 else .Machine$double.xmin
  table$minus_log10_p <- -log10(pmax(table$p_value, floor, na.rm = FALSE))
  table$label <- .plot_feature_labels(table$feature_id, label_features, feature_labels, result)
  plot <- ggplot2::ggplot(table, ggplot2::aes(
    x = .data[["effect"]], y = .data[["minus_log10_p"]],
    colour = .data[["direction"]]
  )) +
    ggplot2::geom_point(alpha = 0.72, size = 1.6, na.rm = TRUE) +
    ggplot2::geom_vline(xintercept = c(-min_abs_effect, min_abs_effect), linetype = 2, colour = "#777777") +
    .plot_de_scale() +
    ggplot2::labs(
      x = "Effect estimate", y = expression(-log[10](italic(p))),
      colour = "Direction",
      alt = "A volcano plot of effect estimates and raw p-values, coloured by FDR and effect-size thresholds."
    ) +
    .plot_de_theme()
  if (symmetric_x) {
    x_limit <- max(c(abs(table$effect[is.finite(table$effect)]), min_abs_effect))
    if (x_limit == 0) x_limit <- 1
    plot <- plot + ggplot2::coord_cartesian(xlim = c(-x_limit, x_limit))
  }
  plot <- .plot_add_feature_text(plot, table, "minus_log10_p")
  .plot_with_dimensions(plot, width = 8.5, height = 7)
}

#' Plot a differential-expression MA plot
#'
#' @inheritParams plot_de_volcano
#' @param abundance Optional complete feature-named abundance vector.
#' @param abundance_column Optional abundance column in a data-frame result.
#' @param abundance_transform Display abundance unchanged, on a base-10 log
#'   axis, or choose the backend convention automatically. Zero abundances are
#'   omitted from a log10 plot; negative abundances are rejected.
#'
#' @return An unprinted standard ggplot object carrying recommended physical
#'   dimensions for [plot_save()].
#' @family plotting
#' @export
plot_de_ma <- function(
    result,
    coef = NULL,
    fdr = 0.05,
    min_abs_effect = 1,
    label_features = NULL,
    feature_labels = NULL,
    abundance = NULL,
    abundance_column = NULL,
    abundance_transform = c("auto", "identity", "log10")
) {
  abundance_transform <- match.arg(abundance_transform)
  if (!is.null(abundance) && !is.null(abundance_column)) {
    stop("Supply only one of `abundance` and `abundance_column`.", call. = FALSE)
  }
  table <- .plot_de_data(result, coef, fdr, min_abs_effect)
  abundance_info <- .plot_de_abundance(result, table$feature_id, abundance, abundance_column)
  table$abundance <- abundance_info$value
  table$label <- .plot_feature_labels(table$feature_id, label_features, feature_labels, result)
  transform <- if (abundance_transform == "auto") abundance_info$transform else abundance_transform
  if (transform == "log10") {
    if (any(table$abundance < 0, na.rm = TRUE)) {
      stop("A log10 abundance axis requires non-negative values.", call. = FALSE)
    }
    zero_abundance <- !is.na(table$abundance) & table$abundance == 0
    if (any(zero_abundance)) {
      warning(
        "Removed ", sum(zero_abundance),
        " zero abundance value(s) from the log10 MA plot.",
        call. = FALSE
      )
      table$abundance[zero_abundance] <- NA_real_
    }
  }
  plot <- ggplot2::ggplot(table, ggplot2::aes(
    x = .data[["abundance"]], y = .data[["effect"]],
    colour = .data[["direction"]]
  )) +
    ggplot2::geom_point(alpha = 0.72, size = 1.6, na.rm = TRUE) +
    ggplot2::geom_hline(yintercept = c(-min_abs_effect, 0, min_abs_effect),
                        linetype = c(2, 1, 2), colour = "#777777") +
    .plot_de_scale() +
    ggplot2::labs(
      x = "Mean abundance", y = "Effect estimate", colour = "Direction",
      alt = "An MA plot of mean abundance and effect estimates, coloured by FDR and effect-size thresholds."
    ) +
    .plot_de_theme()
  if (transform == "log10") plot <- plot + ggplot2::scale_x_log10()
  plot <- .plot_add_feature_text(plot, table, "effect")
  .plot_with_dimensions(plot, width = 8.5, height = 7)
}

#' Plot selected assay features as a heatmap
#'
#' @param x A `MultiAssayExperiment`.
#' @param experiment Experiment name.
#' @param assay Assay name or index.
#' @param features Explicit feature identifiers to plot.
#' @param scale Standardize each feature row or retain assay values.
#' @param cluster_rows,cluster_columns Reorder features or samples by
#'   hierarchical clustering. Dendrograms are not drawn.
#' @param column_split Optional sample metadata column or complete sample-named
#'   vector used to split the x axis.
#' @param feature_label Optional feature metadata column or complete
#'   feature-named labels. Duplicate display labels are made unique.
#'
#' @return An unprinted standard ggplot object carrying recommended physical
#'   dimensions for [plot_save()].
#' @family plotting
#' @export
plot_assay_heatmap <- function(
    x,
    experiment,
    assay,
    features,
    scale = c("row", "none"),
    cluster_rows = TRUE,
    cluster_columns = TRUE,
    column_split = NULL,
    feature_label = NULL
) {
  scale <- match.arg(scale)
  .plot_assert_flag(cluster_rows, "cluster_rows")
  .plot_assert_flag(cluster_columns, "cluster_columns")
  if (is.null(features) || !length(features)) {
    stop("`features` must explicitly contain at least one feature.", call. = FALSE)
  }
  matrix <- .select_features(.pull_matrix(x, experiment, assay), features)
  if (any(!is.finite(matrix))) stop("The heatmap assay must contain finite values.", call. = FALSE)
  labels <- .plot_mae_annotation(
    feature_label, mae_feature_data(x, experiment), rownames(matrix),
    "feature_label", "feature"
  )
  if (is.null(labels)) labels <- stats::setNames(rownames(matrix), rownames(matrix))
  if (scale == "row") {
    deviations <- apply(matrix, 1L, stats::sd)
    constant <- !is.finite(deviations) | deviations == 0
    if (any(constant)) {
      warning(
        "Removed constant heatmap rows: ",
        paste(rownames(matrix)[constant], collapse = ", "),
        ".", call. = FALSE
      )
      matrix <- matrix[!constant, , drop = FALSE]
      labels <- labels[rownames(matrix)]
    }
    if (!nrow(matrix)) stop("No non-constant heatmap rows remain.", call. = FALSE)
    matrix <- t(scale(t(matrix)))
  }
  row_order <- rownames(matrix)
  if (cluster_rows && nrow(matrix) > 1L) {
    row_order <- row_order[stats::hclust(stats::dist(matrix))$order]
  }
  split <- .plot_mae_annotation(
    column_split, mae_samples(x, experiment), colnames(matrix),
    "column_split", "sample"
  )
  column_order <- colnames(matrix)
  if (cluster_columns && ncol(matrix) > 1L) {
    column_order <- column_order[stats::hclust(stats::dist(t(matrix)))$order]
  }
  if (!is.null(split)) {
    split <- factor(split, levels = unique(as.character(split)))
    column_order <- column_order[order(split[column_order], match(column_order, column_order))]
  }
  data <- expand.grid(feature_id = row_order, sample = column_order, stringsAsFactors = FALSE)
  data$value <- matrix[cbind(match(data$feature_id, rownames(matrix)),
                             match(data$sample, colnames(matrix)))]
  display_labels <- stats::setNames(
    make.unique(as.character(labels[row_order])),
    row_order
  )
  data$feature <- factor(
    unname(display_labels[data$feature_id]),
    levels = rev(unname(display_labels[row_order]))
  )
  data$sample <- factor(data$sample, levels = column_order)
  if (!is.null(split)) data$split <- split[as.character(data$sample)]
  plot <- ggplot2::ggplot(data, ggplot2::aes(
    x = .data[["sample"]], y = .data[["feature"]], fill = .data[["value"]]
  )) +
    ggplot2::geom_tile() +
    ggplot2::labs(
      x = "Sample", y = "Feature",
      alt = "A heatmap of explicitly selected assay features across samples."
    ) +
    theme_bulkmae() +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)
    )
  plot <- if (scale == "row") {
    plot + ggplot2::scale_fill_gradient2(
      low = "#3B4CC0", mid = "#F7F7F7", high = "#B40426",
      midpoint = 0, name = "Row z-score"
    )
  } else {
    plot + ggplot2::scale_fill_viridis_c(name = "Assay value")
  }
  if (!is.null(split)) {
    plot <- plot + ggplot2::facet_grid(cols = ggplot2::vars(.data[["split"]]),
                                      scales = "free_x", space = "free_x")
  }
  .plot_with_dimensions(
    plot,
    width = .plot_clamped_dimension(ncol(matrix), base = 5, per_item = 0.5,
                                    minimum = 10, maximum = 18),
    height = .plot_clamped_dimension(nrow(matrix), base = 4, per_item = 0.28,
                                     minimum = 8, maximum = 20)
  )
}

.bulkmae_colours <- c(
  orange = "#D55E00", blue = "#0072B2", green = "#009E73",
  yellow = "#F0E442", sky = "#56B4E9", vermillion = "#E69F00",
  purple = "#CC79A7", black = "#000000", grey = "#999999"
)

.bulkmae_qualitative <- c(
  vermillion = "#D55E00", blue = "#0072B2", green = "#009E73",
  yellow = "#F0E442", sky = "#56B4E9", orange = "#E69F00",
  purple = "#CC79A7", black = "#000000"
)

.bulkmae_text_size_pt <- 6

.plot_require_columns <- function(x, columns, argument) {
  missing <- setdiff(columns, names(x))
  if (length(missing)) stop("`", argument, "` is missing columns: ", paste(missing, collapse = ", "), ".", call. = FALSE)
}

.plot_assert_ids <- function(ids, label) {
  if (is.null(ids) || !length(ids) || anyNA(ids) || any(!nzchar(ids)) || anyDuplicated(ids)) {
    stop(label, " must be unique and non-missing.", call. = FALSE)
  }
  invisible(ids)
}

.plot_assert_finite_numeric <- function(x, label) {
  if (!is.numeric(x) || any(!is.finite(x))) stop(label, " must contain finite numeric values.", call. = FALSE)
}

.plot_assert_positive_number <- function(x, argument) {
  if (!is.numeric(x) || length(x) != 1L || is.na(x) ||
      !is.finite(x) || x <= 0) {
    stop("`", argument, "` must be one finite positive number.", call. = FALSE)
  }
  invisible(x)
}

.plot_with_dimensions <- function(plot, width, height) {
  .plot_assert_positive_number(width, "width")
  .plot_assert_positive_number(height, "height")
  attr(plot, "bulkmae_dimensions") <- list(
    width = as.numeric(width),
    height = as.numeric(height),
    units = "cm"
  )
  plot
}

.plot_clamped_dimension <- function(
    item_count,
    base,
    per_item,
    minimum,
    maximum
) {
  min(maximum, max(minimum, base + per_item * item_count))
}

.plot_convert_cm <- function(value, units) {
  switch(
    units,
    cm = value,
    `in` = value / 2.54,
    mm = value * 10
  )
}

.plot_assert_flag <- function(x, argument) {
  if (!is.logical(x) || length(x) != 1L || is.na(x)) stop("`", argument, "` must be TRUE or FALSE.", call. = FALSE)
}

.plot_axes <- function(axes) {
  if (!is.numeric(axes) || length(axes) != 2L || anyNA(axes) ||
      any(!is.finite(axes)) || any(axes != trunc(axes)) || any(axes < 1L) ||
      anyDuplicated(axes)) stop("`axes` must contain two distinct positive integers.", call. = FALSE)
  as.integer(axes)
}

.plot_embedding_data <- function(embedding, axes) {
  labels <- paste("Dimension", axes)
  if (inherits(embedding, "prcomp")) {
    coordinates <- embedding$x
    variance <- embedding$sdev^2 / sum(embedding$sdev^2) * 100
    labels <- sprintf("PC%d (%.1f%%)", axes, variance[axes])
  } else if (methods::is(embedding, "MDS")) {
    coordinates <- sweep(
      as.matrix(embedding$eigen.vectors), 2L,
      sqrt(pmax(as.numeric(embedding$eigen.values), 0)), "*"
    )
    rownames(coordinates) <- rownames(embedding$distance.matrix.squared)
    variance <- as.numeric(embedding$var.explained) * 100
    labels <- sprintf("MDS%d (%.1f%%)", axes, variance[axes])
  } else if (inherits(embedding, "Rtsne") && is.matrix(embedding$Y)) {
    coordinates <- embedding$Y
    labels <- paste("t-SNE", axes)
  } else if (is.matrix(embedding) || is.data.frame(embedding)) {
    coordinates <- as.matrix(embedding)
  } else {
    stop("Unsupported `embedding`; use prcomp, limma MDS, Rtsne, or a matrix.", call. = FALSE)
  }
  if (!is.numeric(coordinates) || max(axes) > ncol(coordinates)) {
    stop("Requested `axes` are unavailable in `embedding`.", call. = FALSE)
  }
  samples <- rownames(coordinates)
  .plot_assert_ids(samples, "Embedding sample names")
  values <- coordinates[, axes, drop = FALSE]
  if (any(!is.finite(values))) stop("Selected embedding axes contain non-finite values.", call. = FALSE)
  list(
    data = data.frame(sample = samples, x = values[, 1L], y = values[, 2L], check.names = FALSE),
    x_label = labels[[1L]], y_label = labels[[2L]]
  )
}

.plot_align_frame <- function(
    frame,
    ids,
    argument,
    ids_label = "embedding samples"
) {
  .plot_assert_ids(rownames(frame), paste0("`", argument, "` row names"))
  if (!setequal(rownames(frame), ids)) {
    stop(
      "`", argument, "` row names must match ", ids_label, " exactly.",
      call. = FALSE
    )
  }
  invisible(frame)
}

.plot_discrete_values <- function(x) {
  levels <- if (is.factor(x)) levels(droplevels(x)) else unique(as.character(x))
  if (length(levels) > length(.bulkmae_qualitative)) {
    # An unnamed, deliberately short scale lets the plot constructor return a
    # standard ggplot. Its build fails clearly unless the caller first replaces
    # the scale with `+ scale_colour_*()`, avoiding silent colour recycling.
    return(unname(.bulkmae_qualitative))
  }
  stats::setNames(
    unname(.bulkmae_qualitative[seq_along(levels)]),
    levels
  )
}

.plot_de_data <- function(result, coef, fdr, min_abs_effect) {
  .differential_assert_probability(fdr, "fdr")
  if (!is.numeric(min_abs_effect) || length(min_abs_effect) != 1L ||
      is.na(min_abs_effect) || !is.finite(min_abs_effect) || min_abs_effect < 0) {
    stop("`min_abs_effect` must be one finite, non-negative number.", call. = FALSE)
  }
  table <- de_table(result, coef = coef)
  statistic_type <- attr(table, "statistic_type", exact = TRUE)
  if (!is.null(statistic_type) &&
      statistic_type %in% c("nondirectional_lrt", "nondirectional_f_or_lr")) {
    stop(
      "This result has no single effect direction. Use a DESeq2 Wald result or a one-degree-of-freedom test.",
      call. = FALSE
    )
  }
  significant <- !is.na(table$adjusted_p_value) & table$adjusted_p_value <= fdr &
    is.finite(table$effect) & abs(table$effect) >= min_abs_effect
  direction <- rep("Not significant", nrow(table))
  direction[significant & table$effect < 0] <- "Down"
  direction[significant & table$effect > 0] <- "Up"
  table$direction <- factor(direction, levels = c("Down", "Not significant", "Up"))
  table
}

.plot_de_scale <- function() {
  ggplot2::scale_colour_manual(
    values = c(Down = .bulkmae_colours[["blue"]],
               `Not significant` = .bulkmae_colours[["grey"]],
               Up = .bulkmae_colours[["orange"]]),
    breaks = c("Down", "Not significant", "Up"), drop = FALSE,
    guide = ggplot2::guide_legend(nrow = 1, byrow = TRUE)
  )
}

.plot_de_theme <- function() {
  theme_bulkmae() +
    ggplot2::theme(
      legend.position = "top",
      legend.direction = "horizontal"
    )
}

.plot_feature_labels <- function(ids, selected, labels, result) {
  output <- rep(NA_character_, length(ids))
  if (is.null(selected)) return(output)
  .plot_assert_ids(as.character(selected), "`label_features`")
  unknown <- setdiff(selected, ids)
  if (length(unknown)) stop("Unknown `label_features`: ", paste(unknown, collapse = ", "), ".", call. = FALSE)
  resolved <- stats::setNames(ids, ids)
  if (!is.null(labels)) {
    if (is.character(labels) && length(labels) == 1L && is.data.frame(result) && labels %in% names(result)) {
      result_ids <- if ("feature_id" %in% names(result)) as.character(result$feature_id) else rownames(result)
      .plot_assert_ids(result_ids, "Result feature identifiers")
      resolved <- stats::setNames(as.character(result[[labels]]), result_ids)
    } else {
      if (is.null(names(labels))) stop("`feature_labels` must be feature-named.", call. = FALSE)
      .plot_assert_ids(names(labels), "Names of `feature_labels`")
      if (!setequal(names(labels), ids)) stop("`feature_labels` must name every result feature exactly once.", call. = FALSE)
      resolved <- stats::setNames(as.character(labels[ids]), ids)
    }
    if (anyNA(resolved) || any(!nzchar(resolved))) stop("`feature_labels` cannot contain missing or empty labels.", call. = FALSE)
  }
  output[match(selected, ids)] <- unname(resolved[selected])
  output
}

.plot_add_feature_text <- function(plot, data, y_column) {
  x_column <- if (y_column == "effect") "abundance" else "effect"
  labelled <- !is.na(data$label) &
    is.finite(data[[x_column]]) & is.finite(data[[y_column]])
  if (!any(labelled)) return(plot)
  y_values <- data[[y_column]][is.finite(data[[y_column]])]
  y_range <- diff(range(y_values))
  nudge_y <- if (is.finite(y_range) && y_range > 0) y_range * 0.025 else 0
  plot + ggplot2::geom_text(
    data = data[labelled, , drop = FALSE],
    mapping = ggplot2::aes(
      x = .data[[x_column]],
      y = .data[[y_column]], label = .data[["label"]]
    ),
    nudge_y = nudge_y,
    size = .bulkmae_text_size_pt,
    size.unit = "pt",
    check_overlap = TRUE, show.legend = FALSE, na.rm = TRUE
  )
}

.plot_de_abundance <- function(result, ids, abundance, abundance_column) {
  transform <- "identity"
  if (!is.null(abundance)) {
    if (!is.numeric(abundance) || is.null(names(abundance))) stop("`abundance` must be a feature-named numeric vector.", call. = FALSE)
    .plot_assert_ids(names(abundance), "Names of `abundance`")
    if (!setequal(names(abundance), ids)) stop("`abundance` must name every result feature exactly once.", call. = FALSE)
    value <- abundance[ids]
  } else {
    source <- NULL
    column <- abundance_column
    if (methods::is(result, "DESeqResults")) {
      source <- as.data.frame(result, optional = TRUE)
      if (is.null(column)) column <- "baseMean"
      transform <- "log10"
    } else if (inherits(result, c("DGELRT", "DGEExact"))) {
      source <- as.data.frame(result$table, optional = TRUE)
      if (is.null(column)) column <- "logCPM"
    } else if (inherits(result, "TopTags") && !is.null(result$table)) {
      source <- as.data.frame(result$table, optional = TRUE)
      if (is.null(column)) column <- "logCPM"
    } else if (inherits(result, "MArrayLM")) {
      if (is.null(result$Amean)) stop("The limma result does not contain `Amean`; supply `abundance`.", call. = FALSE)
      value <- stats::setNames(as.numeric(result$Amean), rownames(result$coefficients))[ids]
    } else if (is.data.frame(result) || is.matrix(result)) {
      source <- as.data.frame(result, optional = TRUE)
      candidates <- intersect(c("baseMean", "logCPM", "Amean", "AveExpr", "mean_abundance"), names(source))
      if (is.null(column)) {
        if (length(candidates) != 1L) stop("Supply `abundance` or `abundance_column`; automatic abundance detection is ambiguous or unavailable.", call. = FALSE)
        column <- candidates
      }
      if (identical(column, "baseMean")) transform <- "log10"
    } else {
      stop("Supply `abundance` for this result class.", call. = FALSE)
    }
    if (exists("source", inherits = FALSE) && !is.null(source)) {
      .plot_require_columns(source, column, "result")
      source_ids <- rownames(source)
      if ("feature_id" %in% names(source)) source_ids <- as.character(source$feature_id)
      .plot_assert_ids(source_ids, "Abundance feature identifiers")
      if (!setequal(source_ids, ids)) stop("Abundance rows must match differential-result features exactly.", call. = FALSE)
      value <- stats::setNames(source[[column]], source_ids)[ids]
    }
  }
  .plot_assert_finite_numeric(value, "Abundance values")
  list(value = as.numeric(value), transform = transform)
}

.plot_mae_annotation <- function(value, data, ids, argument, axis) {
  if (is.null(value)) return(NULL)
  metadata_column <- is.character(value) && length(value) == 1L && is.null(names(value))
  if (metadata_column) {
    .plot_require_columns(data, value, if (axis == "sample") "sample metadata" else "feature metadata")
    if (!all(ids %in% rownames(data))) stop("Metadata do not contain all selected ", axis, "s.", call. = FALSE)
    result <- data[ids, value, drop = TRUE]
    names(result) <- ids
  } else {
    if (is.null(names(value))) stop("`", argument, "` must be a metadata column or a named vector.", call. = FALSE)
    result <- value
  }
  .plot_assert_ids(names(result), paste0("Names of `", argument, "`"))
  if (!setequal(names(result), ids)) stop("`", argument, "` must describe the selected ", axis, "s exactly.", call. = FALSE)
  result <- result[ids]
  if (anyNA(result) || any(!nzchar(as.character(result)))) stop("`", argument, "` cannot contain missing or empty values.", call. = FALSE)
  result
}

.plot_enrichment_probabilities <- function(
    values,
    label,
    allow_all_missing = FALSE
) {
  if (!is.numeric(values) || any(is.infinite(values)) ||
      any(values < 0 | values > 1, na.rm = TRUE)) {
    stop(label, " must contain probabilities between zero and one.", call. = FALSE)
  }
  if (!allow_all_missing && !any(is.finite(values))) {
    stop(label, " are unavailable.", call. = FALSE)
  }
  invisible(values)
}

.plot_enrichment_labels <- function(data, labels) {
  if (is.null(labels)) return(make.unique(data$term_label))
  if (is.null(names(labels))) {
    stop("`term_labels` must be term-named.", call. = FALSE)
  }
  .plot_assert_ids(names(labels), "Names of `term_labels`")
  if (!setequal(names(labels), data$term_id)) {
    stop("`term_labels` must describe the selected terms exactly.", call. = FALSE)
  }
  labels <- as.character(labels[data$term_id])
  if (anyNA(labels) || any(!nzchar(labels))) {
    stop("`term_labels` cannot contain missing or empty labels.", call. = FALSE)
  }
  make.unique(labels)
}

utils::globalVariables(".data")
