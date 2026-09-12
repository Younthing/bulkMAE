#' Plot WGCNA soft-threshold diagnostics
#'
#' @param fit A native list returned by [coexpr_pick_power()].
#'
#' @return An unprinted standard ggplot object carrying recommended physical
#'   dimensions for [plot_save()].
#' @family plotting
#' @export
plot_coexpr_power <- function(fit) {
  if (!is.list(fit) || is.null(fit$fitIndices)) {
    stop("`fit` must be a pickSoftThreshold list containing `fitIndices`.", call. = FALSE)
  }
  indices <- as.data.frame(fit$fitIndices, optional = TRUE)
  required <- c("Power", "SFT.R.sq", "mean.k.")
  .plot_require_columns(indices, required, "fit$fitIndices")
  .plot_assert_finite_numeric(indices$Power, "`fit$fitIndices$Power`")
  .plot_assert_finite_numeric(indices$SFT.R.sq, "`fit$fitIndices$SFT.R.sq`")
  .plot_assert_finite_numeric(indices$mean.k., "`fit$fitIndices$mean.k.`")
  data <- rbind(
    data.frame(
      power = indices$Power,
      value = indices$SFT.R.sq,
      metric = "Scale-free R\u00b2",
      stringsAsFactors = FALSE
    ),
    data.frame(
      power = indices$Power,
      value = indices$mean.k.,
      metric = "Mean connectivity",
      stringsAsFactors = FALSE
    )
  )
  data$metric <- factor(
    data$metric,
    levels = c("Scale-free R\u00b2", "Mean connectivity")
  )
  plot <- ggplot2::ggplot(data, ggplot2::aes(
    x = .data[["power"]], y = .data[["value"]]
  )) +
    ggplot2::geom_line(colour = .bulkmae_qualitative[["blue"]], linewidth = 0.4) +
    ggplot2::geom_point(colour = .bulkmae_qualitative[["blue"]], size = 1.2) +
    ggplot2::facet_wrap(ggplot2::vars(.data[["metric"]]), scales = "free_y") +
    ggplot2::labs(
      x = "Soft-threshold power",
      y = "Diagnostic value",
      alt = paste(
        "Scale-free topology R-squared and mean connectivity",
        "across candidate WGCNA soft-threshold powers."
      )
    ) +
    theme_bulkmae()
  .plot_with_dimensions(plot, width = 12, height = 5.5)
}

#' Plot WGCNA module sizes
#'
#' @param modules A feature-named module vector from [coexpr_modules()].
#'
#' @return An unprinted standard ggplot object carrying recommended physical
#'   dimensions for [plot_save()].
#' @family plotting
#' @export
plot_coexpr_modules <- function(modules) {
  modules <- .plot_named_labels(modules, "modules")
  counts <- as.data.frame(table(module = unname(modules)), stringsAsFactors = FALSE)
  names(counts) <- c("module", "n")
  counts$n <- as.integer(counts$n)
  counts$module <- factor(counts$module, levels = counts$module[order(-counts$n)])
  fills <- .plot_module_colours(levels(counts$module))
  plot <- ggplot2::ggplot(counts, ggplot2::aes(
    x = .data[["module"]], y = .data[["n"]], fill = .data[["module"]]
  )) +
    ggplot2::geom_col(width = 0.72, colour = "#222222", linewidth = 0.2) +
    ggplot2::scale_fill_manual(values = fills, guide = "none") +
    ggplot2::labs(
      x = "Module",
      y = "Features",
      alt = "A bar chart of feature counts in each WGCNA module."
    ) +
    theme_bulkmae() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
  .plot_with_dimensions(
    plot,
    width = .plot_clamped_dimension(
      nrow(counts), base = 6, per_item = 0.7, minimum = 8, maximum = 16
    ),
    height = 6.5
  )
}

#' Plot module–trait correlations
#'
#' @param trait A list returned by [coexpr_module_trait()] with `correlation`
#'   and `p_value` matrices.
#'
#' @return An unprinted standard ggplot object carrying recommended physical
#'   dimensions for [plot_save()].
#' @family plotting
#' @export
plot_coexpr_trait <- function(trait) {
  if (
    !is.list(trait) || is.null(trait$correlation) || is.null(trait$p_value)
  ) {
    stop("`trait` must contain `correlation` and `p_value` matrices.", call. = FALSE)
  }
  correlation <- as.matrix(trait$correlation)
  p_value <- as.matrix(trait$p_value)
  if (
    !is.numeric(correlation) || !is.numeric(p_value) ||
      !identical(dim(correlation), dim(p_value)) ||
      !identical(dimnames(correlation), dimnames(p_value))
  ) {
    stop("`correlation` and `p_value` must be numeric matrices with the same dimnames.", call. = FALSE)
  }
  .plot_assert_ids(rownames(correlation), "Module names")
  .plot_assert_ids(colnames(correlation), "Trait names")
  if (any(!is.finite(correlation) & !is.na(correlation))) {
    stop("`correlation` must contain finite or missing values.", call. = FALSE)
  }
  data <- expand.grid(
    module = rownames(correlation),
    trait = colnames(correlation),
    stringsAsFactors = FALSE
  )
  data$correlation <- correlation[cbind(
    match(data$module, rownames(correlation)),
    match(data$trait, colnames(correlation))
  )]
  data$p_value <- p_value[cbind(
    match(data$module, rownames(p_value)),
    match(data$trait, colnames(p_value))
  )]
  data$label <- ifelse(
    is.finite(data$correlation),
    sprintf("%.2f", data$correlation),
    ""
  )
  data$module <- factor(data$module, levels = rev(rownames(correlation)))
  data$trait <- factor(data$trait, levels = colnames(correlation))
  plot <- ggplot2::ggplot(data, ggplot2::aes(
    x = .data[["trait"]], y = .data[["module"]], fill = .data[["correlation"]]
  )) +
    ggplot2::geom_tile(colour = "#222222", linewidth = 0.15) +
    ggplot2::geom_text(
      ggplot2::aes(label = .data[["label"]]),
      size = .bulkmae_text_size_pt,
      size.unit = "pt",
      colour = "#111111"
    ) +
    ggplot2::scale_fill_gradient2(
      low = "#3B4CC0", mid = "#F7F7F7", high = "#B40426",
      midpoint = 0, limits = c(-1, 1), name = "Correlation",
      na.value = "#EEEEEE"
    ) +
    ggplot2::labs(
      x = "Trait",
      y = "Module",
      alt = "A heatmap of Pearson correlations between module eigengenes and sample traits."
    ) +
    theme_bulkmae() +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)
    )
  .plot_with_dimensions(
    plot,
    width = .plot_clamped_dimension(
      ncol(correlation), base = 4, per_item = 1.2, minimum = 8, maximum = 16
    ),
    height = .plot_clamped_dimension(
      nrow(correlation), base = 3.5, per_item = 0.7, minimum = 6, maximum = 16
    )
  )
}

#' Plot module membership against gene significance
#'
#' @param membership A feature-by-module matrix from [coexpr_membership()].
#' @param gene_significance A feature-named numeric vector on the same IDs.
#' @param module One module label or membership column name.
#'
#' @return An unprinted standard ggplot object carrying recommended physical
#'   dimensions for [plot_save()].
#' @family plotting
#' @export
plot_coexpr_membership <- function(membership, gene_significance, module) {
  membership <- .coexpr_membership_matrix(membership)
  if (
    !is.numeric(gene_significance) || is.null(names(gene_significance))
  ) {
    stop("`gene_significance` must be a named numeric vector.", call. = FALSE)
  }
  .plot_assert_ids(names(gene_significance), "Gene-significance names")
  if (!setequal(names(gene_significance), rownames(membership))) {
    stop(
      "`gene_significance` names must match membership features exactly.",
      call. = FALSE
    )
  }
  if (any(!is.finite(gene_significance))) {
    stop("`gene_significance` must contain finite values.", call. = FALSE)
  }
  column <- .coexpr_match_module_column(module, colnames(membership))
  features <- rownames(membership)
  data <- data.frame(
    feature = features,
    membership = unname(membership[features, column]),
    gene_significance = unname(gene_significance[features]),
    check.names = FALSE
  )
  plot <- ggplot2::ggplot(data, ggplot2::aes(
    x = .data[["membership"]], y = .data[["gene_significance"]]
  )) +
    ggplot2::geom_point(
      colour = .bulkmae_qualitative[["blue"]],
      size = 1.1,
      alpha = 0.85
    ) +
    ggplot2::geom_smooth(
      method = "lm",
      se = FALSE,
      colour = .bulkmae_qualitative[["vermillion"]],
      linewidth = 0.4,
      formula = y ~ x
    ) +
    ggplot2::labs(
      x = "Module membership",
      y = "Gene significance",
      alt = paste(
        "A scatter plot of feature-module membership against",
        "feature-trait gene significance."
      )
    ) +
    theme_bulkmae()
  .plot_with_dimensions(plot, width = 8, height = 8)
}

#' Plot a consensus-clustering consensus matrix
#'
#' @param results A native [cluster_consensus()] list, or a square consensus
#'   matrix with sample dimnames.
#' @param k Cluster count used when `results` is a ConsensusClusterPlus list.
#' @param show_names Draw sample axis labels.
#'
#' @return An unprinted standard ggplot object carrying recommended physical
#'   dimensions for [plot_save()].
#' @family plotting
#' @export
plot_cluster_consensus <- function(results, k = NULL, show_names = FALSE) {
  .plot_assert_flag(show_names, "show_names")
  matrix <- .plot_consensus_matrix(results, k)
  order <- rownames(matrix)
  data <- expand.grid(row = order, column = order, stringsAsFactors = FALSE)
  data$value <- matrix[cbind(
    match(data$row, rownames(matrix)),
    match(data$column, colnames(matrix))
  )]
  data$row <- factor(data$row, levels = rev(order))
  data$column <- factor(data$column, levels = order)
  plot <- ggplot2::ggplot(data, ggplot2::aes(
    x = .data[["column"]], y = .data[["row"]], fill = .data[["value"]]
  )) +
    ggplot2::geom_tile() +
    ggplot2::scale_fill_gradient(
      low = "#FFFFFF", high = "#08306B",
      limits = c(0, 1), name = "Consensus"
    ) +
    ggplot2::coord_equal() +
    ggplot2::labs(
      x = "Sample",
      y = "Sample",
      alt = "A heatmap of pairwise consensus-clustering co-assignment values."
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
    nrow(matrix), base = 5, per_item = 0.25, minimum = 8.5, maximum = 17
  )
  .plot_with_dimensions(plot, width = panel_size + 1, height = panel_size)
}

#' Plot consensus CDFs across candidate k
#'
#' @param results A native [cluster_consensus()] list.
#'
#' @return An unprinted standard ggplot object carrying recommended physical
#'   dimensions for [plot_save()].
#' @family plotting
#' @export
plot_cluster_cdf <- function(results) {
  if (!is.list(results) || length(results) < 2L) {
    stop("`results` must be a ConsensusClusterPlus result list.", call. = FALSE)
  }
  rows <- list()
  for (index in seq_along(results)) {
    if (index < 2L) {
      next
    }
    solution <- results[[index]]
    if (!is.list(solution) || is.null(solution$consensusMatrix)) {
      next
    }
    matrix <- as.matrix(solution$consensusMatrix)
    values <- matrix[upper.tri(matrix, diag = FALSE)]
    values <- values[is.finite(values)]
    if (!length(values)) {
      next
    }
    if (any(values < 0 | values > 1)) {
      stop("Consensus values must lie in [0, 1].", call. = FALSE)
    }
    ordered <- sort(values)
    rows[[length(rows) + 1L]] <- data.frame(
      k = as.character(index),
      consensus = ordered,
      cdf = seq_along(ordered) / length(ordered),
      stringsAsFactors = FALSE
    )
  }
  if (!length(rows)) {
    stop("No consensus matrices were available for k >= 2.", call. = FALSE)
  }
  data <- do.call(rbind, rows)
  data$k <- factor(data$k, levels = unique(data$k))
  fills <- .plot_discrete_values(data$k)
  plot <- ggplot2::ggplot(data, ggplot2::aes(
    x = .data[["consensus"]], y = .data[["cdf"]], colour = .data[["k"]]
  )) +
    ggplot2::geom_step(linewidth = 0.45) +
    ggplot2::scale_colour_manual(values = fills, name = "k") +
    ggplot2::labs(
      x = "Consensus",
      y = "CDF",
      alt = "Empirical cumulative distributions of pairwise consensus values for each k."
    ) +
    theme_bulkmae()
  .plot_with_dimensions(plot, width = 8.5, height = 6.5)
}

#' Plot subtype sample counts
#'
#' @param classes A sample-named vector from [cluster_consensus_classes()] or
#'   [cluster_nmf_classes()].
#' @param group Optional sample-named vector or factor used to split counts.
#'   Names must match `classes` exactly.
#'
#' @return An unprinted standard ggplot object carrying recommended physical
#'   dimensions for [plot_save()].
#' @family plotting
#' @export
plot_cluster_sizes <- function(classes, group = NULL) {
  classes <- .plot_named_labels(classes, "classes")
  data <- data.frame(
    sample = names(classes),
    subtype = as.character(unname(classes)),
    stringsAsFactors = FALSE
  )
  if (is.null(group)) {
    data$subtype <- factor(data$subtype, levels = unique(data$subtype))
    plot <- ggplot2::ggplot(data, ggplot2::aes(x = .data[["subtype"]])) +
      ggplot2::geom_bar(
        fill = unname(.bulkmae_qualitative[["blue"]]),
        colour = "#222222",
        linewidth = 0.2,
        width = 0.72
      ) +
      ggplot2::labs(
        x = "Subtype",
        y = "Samples",
        alt = "A bar chart of sample counts in each molecular subtype."
      )
  } else {
    group <- .plot_align_named_labels(group, names(classes), "group")
    data$group <- as.character(unname(group[data$sample]))
    data$subtype <- factor(data$subtype, levels = unique(data$subtype))
    data$group <- factor(data$group, levels = unique(data$group))
    fills <- .plot_discrete_values(data$group)
    plot <- ggplot2::ggplot(data, ggplot2::aes(
      x = .data[["subtype"]], fill = .data[["group"]]
    )) +
      ggplot2::geom_bar(colour = "#222222", linewidth = 0.2, width = 0.72) +
      ggplot2::scale_fill_manual(values = fills, name = "Group") +
      ggplot2::labs(
        x = "Subtype",
        y = "Samples",
        alt = "A stacked bar chart of subtype sample counts split by a sample group."
      )
  }
  plot <- plot + theme_bulkmae()
  .plot_with_dimensions(
    plot,
    width = .plot_clamped_dimension(
      length(unique(data$subtype)), base = 5, per_item = 0.9,
      minimum = 7, maximum = 14
    ),
    height = 6.5
  )
}

#' Plot stacked cell-type fractions
#'
#' @param fractions A cell-type-by-sample matrix from [deconv_fractions()].
#'
#' @return An unprinted standard ggplot object carrying recommended physical
#'   dimensions for [plot_save()].
#' @family plotting
#' @export
plot_deconv_stacked <- function(fractions) {
  data <- .plot_fraction_frame(fractions)
  totals <- tapply(data$fraction, data$sample, sum)
  y_name <- if (all(abs(totals - 1) <= 0.05)) "Fraction" else "Estimated value"
  cells <- unique(as.character(data$cell_type))
  data$cell_type <- factor(data$cell_type, levels = cells)
  fills <- .plot_discrete_values(data$cell_type)
  plot <- ggplot2::ggplot(data, ggplot2::aes(
    x = .data[["sample"]], y = .data[["fraction"]], fill = .data[["cell_type"]]
  )) +
    ggplot2::geom_col(width = 0.85, colour = "#222222", linewidth = 0.15) +
    ggplot2::scale_fill_manual(values = fills, name = "Cell type") +
    ggplot2::labs(
      x = "Sample",
      y = y_name,
      alt = "Stacked bars of estimated cell-type values for each sample."
    ) +
    theme_bulkmae() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
  .plot_with_dimensions(
    plot,
    width = .plot_clamped_dimension(
      length(unique(data$sample)), base = 5, per_item = 0.45,
      minimum = 8, maximum = 16
    ),
    height = 7
  )
}

#' Plot cell-type fractions by sample group
#'
#' @param fractions A cell-type-by-sample matrix from [deconv_fractions()].
#' @param group A sample-named vector or factor. Names must match the
#'   fraction columns exactly.
#'
#' @return An unprinted standard ggplot object carrying recommended physical
#'   dimensions for [plot_save()].
#' @family plotting
#' @export
plot_deconv_box <- function(fractions, group) {
  data <- .plot_fraction_frame(fractions)
  group <- .plot_align_named_labels(group, levels(data$sample), "group")
  data$group <- as.character(unname(group[as.character(data$sample)]))
  data$group <- factor(data$group, levels = unique(data$group))
  fills <- .plot_discrete_values(data$group)
  plot <- ggplot2::ggplot(data, ggplot2::aes(
    x = .data[["cell_type"]], y = .data[["fraction"]], fill = .data[["group"]]
  )) +
    ggplot2::geom_boxplot(
      outlier.size = 0.6,
      linewidth = 0.3,
      width = 0.7
    ) +
    ggplot2::scale_fill_manual(values = fills, name = "Group") +
    ggplot2::labs(
      x = "Cell type",
      y = "Estimated value",
      alt = "Boxplots of estimated cell-type values grouped by a sample annotation."
    ) +
    theme_bulkmae() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
  .plot_with_dimensions(
    plot,
    width = .plot_clamped_dimension(
      length(unique(data$cell_type)), base = 5, per_item = 1.1,
      minimum = 8, maximum = 16
    ),
    height = 7
  )
}

#' Plot cell-type fractions as a heatmap
#'
#' @param fractions A cell-type-by-sample matrix from [deconv_fractions()].
#' @param column_split Optional sample-named vector used to facet columns,
#'   for example subtype labels. Names must cover every sample.
#'
#' @return An unprinted standard ggplot object carrying recommended physical
#'   dimensions for [plot_save()].
#' @family plotting
#' @export
plot_deconv_heatmap <- function(fractions, column_split = NULL) {
  data <- .plot_fraction_frame(fractions)
  data$cell_type <- factor(
    data$cell_type,
    levels = rev(levels(data$cell_type))
  )
  if (!is.null(column_split)) {
    split <- .plot_align_named_labels(
      column_split,
      levels(data$sample),
      "column_split"
    )
    data$split <- factor(
      as.character(unname(split[as.character(data$sample)])),
      levels = unique(as.character(unname(split)))
    )
  }
  plot <- ggplot2::ggplot(data, ggplot2::aes(
    x = .data[["sample"]], y = .data[["cell_type"]], fill = .data[["fraction"]]
  )) +
    ggplot2::geom_tile() +
    ggplot2::scale_fill_viridis_c(name = "Estimated value") +
    ggplot2::labs(
      x = "Sample",
      y = "Cell type",
      alt = "A heatmap of estimated cell-type values across samples."
    ) +
    theme_bulkmae() +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)
    )
  if (!is.null(column_split)) {
    plot <- plot + ggplot2::facet_grid(
      cols = ggplot2::vars(.data[["split"]]),
      scales = "free_x",
      space = "free_x"
    )
  }
  .plot_with_dimensions(
    plot,
    width = .plot_clamped_dimension(
      length(levels(data$sample)), base = 5, per_item = 0.45,
      minimum = 8, maximum = 16
    ),
    height = .plot_clamped_dimension(
      length(levels(data$cell_type)), base = 4, per_item = 0.45,
      minimum = 6, maximum = 14
    )
  )
}

.plot_named_labels <- function(value, argument) {
  if (
    (!is.atomic(value) && !is.factor(value)) || is.null(names(value))
  ) {
    stop("`", argument, "` must be a named atomic vector.", call. = FALSE)
  }
  names <- names(value)
  .plot_assert_ids(names, paste0("`", argument, "` names"))
  if (anyNA(value) || any(!nzchar(as.character(value)))) {
    stop("`", argument, "` values cannot be missing or empty.", call. = FALSE)
  }
  stats::setNames(value, names)
}

.plot_align_named_labels <- function(value, ids, argument) {
  value <- .plot_named_labels(value, argument)
  if (!setequal(names(value), ids)) {
    stop(
      "`", argument, "` names must match the plotted samples exactly.",
      call. = FALSE
    )
  }
  value[ids]
}

.plot_fraction_frame <- function(fractions) {
  fractions <- as.matrix(fractions)
  if (!is.numeric(fractions) || !nrow(fractions) || !ncol(fractions)) {
    stop("`fractions` must be a non-empty numeric matrix.", call. = FALSE)
  }
  .plot_assert_ids(rownames(fractions), "Cell-type names")
  .plot_assert_ids(colnames(fractions), "Sample names")
  if (any(!is.finite(fractions)) || any(fractions < 0)) {
    stop("`fractions` must contain finite, non-negative values.", call. = FALSE)
  }
  data <- expand.grid(
    cell_type = rownames(fractions),
    sample = colnames(fractions),
    stringsAsFactors = FALSE
  )
  data$fraction <- fractions[cbind(
    match(data$cell_type, rownames(fractions)),
    match(data$sample, colnames(fractions))
  )]
  data$cell_type <- factor(data$cell_type, levels = rownames(fractions))
  data$sample <- factor(data$sample, levels = colnames(fractions))
  data
}

.plot_consensus_matrix <- function(results, k) {
  if (is.matrix(results)) {
    matrix <- results
  } else {
    if (is.null(k)) {
      stop("`k` is required when `results` is a consensus-clustering list.", call. = FALSE)
    }
    classes <- cluster_consensus_classes(results, k)
    solution <- results[[as.integer(k)]]
    if (is.null(solution$consensusMatrix)) {
      stop("The selected result has no `consensusMatrix` element.", call. = FALSE)
    }
    matrix <- as.matrix(solution$consensusMatrix)
    if (is.null(rownames(matrix)) || is.null(colnames(matrix))) {
      samples <- names(classes)
      if (length(samples) != nrow(matrix) || nrow(matrix) != ncol(matrix)) {
        stop("Consensus matrix dimensions do not match `consensusClass`.", call. = FALSE)
      }
      rownames(matrix) <- samples
      colnames(matrix) <- samples
    }
    tree <- solution$consensusTree
    if (inherits(tree, "hclust") && length(tree$order) == nrow(matrix)) {
      order <- rownames(matrix)[tree$order]
      matrix <- matrix[order, order, drop = FALSE]
    } else {
      order <- names(sort(classes))
      matrix <- matrix[order, order, drop = FALSE]
    }
  }
  if (!is.numeric(matrix) || nrow(matrix) != ncol(matrix) || !nrow(matrix)) {
    stop("Consensus values must be a non-empty numeric square matrix.", call. = FALSE)
  }
  .plot_assert_ids(rownames(matrix), "Consensus row names")
  .plot_assert_ids(colnames(matrix), "Consensus column names")
  if (!setequal(rownames(matrix), colnames(matrix))) {
    stop("Consensus row and column names must contain the same samples.", call. = FALSE)
  }
  matrix <- matrix[rownames(matrix), rownames(matrix), drop = FALSE]
  if (any(!is.finite(matrix)) || any(matrix < 0 | matrix > 1)) {
    stop("Consensus values must be finite and lie in [0, 1].", call. = FALSE)
  }
  matrix
}

.plot_module_colours <- function(modules) {
  standard <- c(
    grey = "#A9A9A9", turquoise = "#40E0D0", blue = "#0000FF",
    brown = "#A52A2A", yellow = "#FFFF00", green = "#00FF00",
    red = "#FF0000", black = "#000000", pink = "#FFC0CB",
    magenta = "#FF00FF", purple = "#800080", greenyellow = "#ADFF2F",
    tan = "#D2B48C", salmon = "#FA8072", cyan = "#00FFFF",
    midnightblue = "#191970", lightcyan = "#E0FFFF", grey60 = "#999999",
    lightgreen = "#90EE90", lightyellow = "#FFFFE0", royalblue = "#4169E1",
    darkred = "#8B0000", darkgreen = "#006400", darkturquoise = "#00CED1",
    darkgrey = "#505050", orange = "#FFA500", darkorange = "#FF8C00",
    white = "#F7F7F7", skyblue = "#87CEEB"
  )
  labels <- c(
    "0" = "grey", "1" = "turquoise", "2" = "blue", "3" = "brown",
    "4" = "yellow", "5" = "green", "6" = "red", "7" = "black",
    "8" = "pink", "9" = "magenta", "10" = "purple"
  )
  mapped <- as.character(modules)
  mapped[mapped %in% names(labels)] <- unname(labels[mapped[mapped %in% names(labels)]])
  fills <- standard[mapped]
  missing <- is.na(fills)
  if (any(missing)) {
    fallback <- .plot_discrete_values(factor(modules[missing], levels = unique(modules[missing])))
    fills[missing] <- unname(fallback[as.character(modules[missing])])
  }
  stats::setNames(unname(fills), modules)
}
