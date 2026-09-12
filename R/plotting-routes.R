#' Plot WGCNA soft-threshold diagnostics
#'
#' Two panels share the candidate-power x axis: (a) scale-free topology
#' model fit and (b) mean connectivity. A dashed reference line marks the
#' conventional R-squared cut. The helper does not choose a power.
#'
#' @param fit A native list returned by [coexpr_pick_power()].
#' @param r_squared Horizontal reference drawn only on panel (a).
#'
#' @return An unprinted standard ggplot object carrying recommended physical
#'   dimensions for [plot_save()].
#' @family plotting
#' @export
plot_coexpr_power <- function(fit, r_squared = 0.8) {
  if (!is.list(fit) || is.null(fit$fitIndices)) {
    stop("`fit` must be a pickSoftThreshold list containing `fitIndices`.", call. = FALSE)
  }
  if (
    !is.numeric(r_squared) || length(r_squared) != 1L || is.na(r_squared) ||
      !is.finite(r_squared) || r_squared < 0 || r_squared > 1
  ) {
    stop("`r_squared` must be one finite number in [0, 1].", call. = FALSE)
  }
  indices <- as.data.frame(fit$fitIndices, optional = TRUE)
  required <- c("Power", "SFT.R.sq", "mean.k.")
  .plot_require_columns(indices, required, "fit$fitIndices")
  .plot_assert_finite_numeric(indices$Power, "`fit$fitIndices$Power`")
  .plot_assert_finite_numeric(indices$SFT.R.sq, "`fit$fitIndices$SFT.R.sq`")
  .plot_assert_finite_numeric(indices$mean.k., "`fit$fitIndices$mean.k.`")
  panel_levels <- c("(a) Scale-free R\u00b2", "(b) Mean connectivity")
  data <- rbind(
    data.frame(
      power = indices$Power,
      value = indices$SFT.R.sq,
      panel = panel_levels[[1L]],
      stringsAsFactors = FALSE
    ),
    data.frame(
      power = indices$Power,
      value = indices$mean.k.,
      panel = panel_levels[[2L]],
      stringsAsFactors = FALSE
    )
  )
  data$panel <- factor(data$panel, levels = panel_levels)
  reference <- data.frame(
    panel = factor(panel_levels[[1L]], levels = panel_levels),
    y = r_squared
  )
  plot <- ggplot2::ggplot(data, ggplot2::aes(
    x = .data[["power"]], y = .data[["value"]]
  )) +
    ggplot2::geom_hline(
      data = reference,
      mapping = ggplot2::aes(yintercept = .data[["y"]]),
      colour = "#666666",
      linetype = "dashed",
      linewidth = 0.35
    ) +
    ggplot2::geom_line(colour = .bulkmae_qualitative[["blue"]], linewidth = 0.4) +
    ggplot2::geom_point(colour = .bulkmae_qualitative[["blue"]], size = 1.2) +
    ggplot2::facet_wrap(
      ggplot2::vars(.data[["panel"]]),
      scales = "free_y",
      nrow = 1,
      strip.position = "left"
    ) +
    ggplot2::labs(
      x = "Soft-threshold power",
      y = NULL,
      alt = paste(
        "Dual-panel WGCNA soft-threshold diagnostics: scale-free R-squared",
        "with a reference line, and mean connectivity versus power."
      )
    ) +
    theme_bulkmae() +
    ggplot2::theme(
      strip.placement = "outside",
      strip.background = ggplot2::element_blank()
    )
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
  data$stars <- .plot_significance_stars(data$p_value)
  data$label <- ifelse(
    is.finite(data$correlation),
    ifelse(
      nzchar(data$stars),
      paste0(sprintf("%.2f", data$correlation), "\n", data$stars),
      sprintf("%.2f", data$correlation)
    ),
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
      alt = paste(
        "A heatmap of Pearson correlations between module eigengenes",
        "and sample traits, with significance stars in each cell."
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

#' Plot a WGCNA gene dendrogram with a module colour bar
#'
#' @param dendrogram An `hclust` or `dendrogram` object, or a
#'   [coexpr_wgcna()] list containing `dendrograms`.
#' @param modules Optional feature-named module vector used as a colour bar.
#'   When `dendrogram` is a WGCNA list, [coexpr_modules()] is used by default.
#'
#' @return An unprinted standard ggplot object carrying recommended physical
#'   dimensions for [plot_save()].
#' @family plotting
#' @export
plot_coexpr_dendrogram <- function(dendrogram, modules = NULL) {
  block_labels <- NULL
  if (is.list(dendrogram) && !is.null(dendrogram$dendrograms)) {
    if (is.null(modules) && !is.null(dendrogram$colors)) {
      modules <- coexpr_modules(dendrogram)
    }
    if (!is.null(dendrogram$colors)) {
      feature_names <- names(coexpr_modules(dendrogram))
      block <- dendrogram$blockGenes
      if (is.list(block) && length(block) && is.numeric(block[[1L]])) {
        block_labels <- feature_names[as.integer(block[[1L]])]
      } else {
        block_labels <- feature_names
      }
    }
    trees <- dendrogram$dendrograms
    if (!length(trees)) {
      stop("The WGCNA fit has no `dendrograms`.", call. = FALSE)
    }
    if (length(trees) > 1L) {
      warning("Using the first WGCNA block dendrogram.", call. = FALSE)
    }
    dendrogram <- trees[[1L]]
  }
  tree <- .plot_as_hclust(dendrogram)
  if (is.null(tree$labels) || anyNA(tree$labels) || any(!nzchar(as.character(tree$labels)))) {
    tree$labels <- block_labels
  }
  labels <- tree$labels
  if (is.null(labels) || anyNA(labels) || any(!nzchar(as.character(labels)))) {
    stop("The dendrogram must have unique, non-missing leaf labels.", call. = FALSE)
  }
  .plot_assert_ids(labels, "Dendrogram leaf labels")
  segments <- .plot_hclust_segments(tree)
  leaf_order <- labels[tree$order]
  segments$component <- factor("Dendrogram", levels = c("Dendrogram", "Module"))
  plot <- ggplot2::ggplot() +
    ggplot2::geom_segment(
      data = segments,
      mapping = ggplot2::aes(
        x = .data[["x"]], xend = .data[["xend"]],
        y = .data[["y"]], yend = .data[["yend"]]
      ),
      colour = "#222222",
      linewidth = 0.25
    )
  if (!is.null(modules)) {
    modules <- .plot_named_labels(modules, "modules")
    if (!setequal(names(modules), labels)) {
      stop("`modules` names must match dendrogram leaves exactly.", call. = FALSE)
    }
    bar <- data.frame(
      feature = leaf_order,
      x = seq_along(leaf_order),
      module = as.character(unname(modules[leaf_order])),
      component = factor("Module", levels = c("Dendrogram", "Module")),
      stringsAsFactors = FALSE
    )
    fills <- .plot_module_colours(unique(bar$module))
    plot <- plot +
      ggplot2::geom_tile(
        data = bar,
        mapping = ggplot2::aes(
          x = .data[["x"]], y = 0.5, fill = .data[["module"]]
        ),
        height = 1,
        colour = "#222222",
        linewidth = 0.05
      ) +
      ggplot2::scale_fill_manual(values = fills, name = "Module")
  }
  plot <- plot +
    ggplot2::facet_grid(
      rows = ggplot2::vars(.data[["component"]]),
      scales = "free_y",
      space = "free_y"
    ) +
    ggplot2::scale_x_continuous(expand = c(0.01, 0)) +
    ggplot2::labs(
      x = "Feature",
      y = "Height",
      alt = "A gene dendrogram with an optional WGCNA module colour bar."
    ) +
    theme_bulkmae() +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_blank()
    )
  .plot_with_dimensions(
    plot,
    width = .plot_clamped_dimension(
      length(leaf_order), base = 8, per_item = 0.04, minimum = 10, maximum = 18
    ),
    height = 7
  )
}

#' Plot a topological-overlap matrix
#'
#' @param tom A feature-named square TOM or adjacency matrix.
#' @param modules Optional feature-named module vector used to order features
#'   and draw annotation bars.
#' @param show_names Draw feature axis labels.
#'
#' @return An unprinted standard ggplot object carrying recommended physical
#'   dimensions for [plot_save()].
#' @family plotting
#' @export
plot_coexpr_tom <- function(tom, modules = NULL, show_names = FALSE) {
  .plot_assert_flag(show_names, "show_names")
  tom <- as.matrix(tom)
  if (!is.numeric(tom) || nrow(tom) != ncol(tom) || !nrow(tom)) {
    stop("`tom` must be a non-empty numeric square matrix.", call. = FALSE)
  }
  .plot_assert_ids(rownames(tom), "TOM row names")
  .plot_assert_ids(colnames(tom), "TOM column names")
  if (!setequal(rownames(tom), colnames(tom))) {
    stop("TOM row and column names must contain the same features.", call. = FALSE)
  }
  tom <- tom[rownames(tom), rownames(tom), drop = FALSE]
  if (any(!is.finite(tom)) || any(tom < 0 | tom > 1)) {
    stop("TOM values must be finite and lie in [0, 1].", call. = FALSE)
  }
  order <- rownames(tom)
  if (!is.null(modules)) {
    modules <- .plot_named_labels(modules, "modules")
    if (!setequal(names(modules), rownames(tom))) {
      stop("`modules` names must match TOM features exactly.", call. = FALSE)
    }
    order <- names(modules)[order(as.character(unname(modules)), method = "radix")]
    tom <- tom[order, order, drop = FALSE]
  } else if (nrow(tom) > 2L) {
    order <- rownames(tom)[stats::hclust(stats::dist(tom))$order]
    tom <- tom[order, order, drop = FALSE]
  }
  data <- expand.grid(row = order, column = order, stringsAsFactors = FALSE)
  data$value <- tom[cbind(
    match(data$row, rownames(tom)),
    match(data$column, colnames(tom))
  )]
  data$row <- factor(data$row, levels = rev(order))
  data$column <- factor(data$column, levels = order)
  data$component <- factor("TOM", levels = c("Module", "TOM"))
  plot <- ggplot2::ggplot(data, ggplot2::aes(
    x = .data[["column"]], y = .data[["row"]], fill = .data[["value"]]
  )) +
    ggplot2::geom_tile() +
    ggplot2::scale_fill_gradient(
      low = "#FFFFFF", high = "#B40426",
      limits = c(0, 1), name = "TOM"
    )
  if (!is.null(modules)) {
    bar <- data.frame(
      column = factor(order, levels = order),
      row = factor("Module", levels = "Module"),
      module = as.character(unname(modules[order])),
      component = factor("Module", levels = c("Module", "TOM")),
      stringsAsFactors = FALSE
    )
    fills <- .plot_module_colours(unique(bar$module))
    plot <- plot +
      ggplot2::geom_tile(
        data = bar,
        mapping = ggplot2::aes(
          x = .data[["column"]], y = .data[["row"]], colour = .data[["module"]]
        ),
        fill = unname(fills[bar$module]),
        linewidth = 0,
        inherit.aes = FALSE
      ) +
      ggplot2::scale_colour_manual(values = fills, name = "Module")
  }
  plot <- plot +
    ggplot2::facet_grid(
      rows = ggplot2::vars(.data[["component"]]),
      scales = "free_y",
      space = "free_y"
    ) +
    ggplot2::labs(
      x = "Feature",
      y = "Feature",
      alt = "A heatmap of topological overlap with an optional module colour bar."
    ) +
    theme_bulkmae() +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)
    )
  if (!show_names) {
    plot <- plot + ggplot2::theme(
      axis.text = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank()
    )
  }
  panel_size <- .plot_clamped_dimension(
    nrow(tom), base = 6, per_item = 0.12, minimum = 8.5, maximum = 16
  )
  .plot_with_dimensions(
    plot,
    width = panel_size + 1.2,
    height = panel_size + if (is.null(modules)) 0 else 1
  )
}

#' Plot a consensus-clustering consensus matrix
#'
#' @param results A native [cluster_consensus()] list, or a square consensus
#'   matrix with sample dimnames.
#' @param k Cluster count used when `results` is a ConsensusClusterPlus list.
#' @param show_names Draw sample axis labels.
#' @param annotation Optional sample-named vector, sample-row data frame, or
#'   named list of sample-named vectors drawn as annotation bars.
#'
#' @return An unprinted standard ggplot object carrying recommended physical
#'   dimensions for [plot_save()].
#' @family plotting
#' @export
plot_cluster_consensus <- function(
    results,
    k = NULL,
    show_names = FALSE,
    annotation = NULL
) {
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
  component_levels <- if (is.null(annotation)) "Consensus" else c("Annotation", "Consensus")
  data$component <- factor("Consensus", levels = component_levels)
  plot <- ggplot2::ggplot(data, ggplot2::aes(
    x = .data[["column"]], y = .data[["row"]], fill = .data[["value"]]
  )) +
    ggplot2::geom_tile() +
    ggplot2::scale_fill_gradient(
      low = "#FFFFFF", high = "#08306B",
      limits = c(0, 1), name = "Consensus"
    )
  if (!is.null(annotation)) {
    tracks <- .plot_annotation_frame(annotation, order, "annotation")
    styled <- .plot_annotation_style(tracks, order)
    styled$data$component <- factor("Annotation", levels = component_levels)
    plot <- plot +
      ggplot2::geom_tile(
        data = styled$data,
        mapping = ggplot2::aes(
          x = .data[["sample"]],
          y = .data[["track"]],
          colour = .data[["legend"]]
        ),
        fill = styled$data$fill_colour,
        linewidth = 0.15,
        inherit.aes = FALSE
      ) +
      ggplot2::scale_colour_manual(values = styled$legend_fills, name = "Annotation") +
      ggplot2::facet_grid(
        rows = ggplot2::vars(.data[["component"]]),
        scales = "free_y",
        space = "free_y"
      )
  } else {
    plot <- plot + ggplot2::coord_equal()
  }
  plot <- plot +
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
  .plot_with_dimensions(
    plot,
    width = panel_size + 1.4,
    height = panel_size + if (is.null(annotation)) 0 else 1.2
  )
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
  pairs <- .cluster_consensus_pair_values(results)
  rows <- lapply(names(pairs), function(k) {
    values <- sort(pairs[[k]])
    unique_x <- c(0, unique(values), 1)
    cdf <- vapply(unique_x, function(x) mean(values <= x), numeric(1))
    data.frame(
      k = k,
      consensus = unique_x,
      cdf = cdf,
      stringsAsFactors = FALSE
    )
  })
  data <- do.call(rbind, rows)
  data$k <- factor(data$k, levels = unique(data$k))
  fills <- .plot_discrete_values(data$k)
  plot <- ggplot2::ggplot(data, ggplot2::aes(
    x = .data[["consensus"]], y = .data[["cdf"]], colour = .data[["k"]]
  )) +
    ggplot2::geom_step(linewidth = 0.45, direction = "hv") +
    ggplot2::scale_colour_manual(values = fills, name = "k") +
    ggplot2::scale_x_continuous(limits = c(0, 1), expand = c(0.01, 0)) +
    ggplot2::scale_y_continuous(limits = c(0, 1), expand = c(0.02, 0)) +
    ggplot2::labs(
      x = "Consensus",
      y = "CDF",
      alt = "Empirical cumulative distributions of pairwise consensus values for each k."
    ) +
    theme_bulkmae()
  .plot_with_dimensions(plot, width = 8.5, height = 6.5)
}

#' Plot relative delta area of consensus CDFs
#'
#' @param results A native [cluster_consensus()] list.
#'
#' @return An unprinted standard ggplot object carrying recommended physical
#'   dimensions for [plot_save()].
#' @family plotting
#' @export
plot_cluster_delta <- function(results) {
  table <- cluster_consensus_delta_area(results)
  table$k <- factor(table$k, levels = table$k)
  plot <- ggplot2::ggplot(table, ggplot2::aes(
    x = .data[["k"]], y = .data[["delta"]]
  )) +
    ggplot2::geom_line(
      ggplot2::aes(group = 1),
      colour = .bulkmae_qualitative[["blue"]],
      linewidth = 0.4
    ) +
    ggplot2::geom_point(colour = .bulkmae_qualitative[["blue"]], size = 1.4) +
    ggplot2::labs(
      x = "k",
      y = "Relative delta area",
      alt = "Relative change in consensus CDF area across candidate k."
    ) +
    theme_bulkmae()
  .plot_with_dimensions(plot, width = 8, height = 6.5)
}

#' Plot the proportion of ambiguous clustering
#'
#' @param results A native [cluster_consensus()] list.
#' @inheritParams cluster_consensus_pac
#'
#' @return An unprinted standard ggplot object carrying recommended physical
#'   dimensions for [plot_save()].
#' @family plotting
#' @export
plot_cluster_pac <- function(results, lower = 0.1, upper = 0.9) {
  pac <- cluster_consensus_pac(results, lower = lower, upper = upper)
  data <- data.frame(
    k = factor(names(pac), levels = names(pac)),
    pac = unname(pac),
    stringsAsFactors = FALSE
  )
  plot <- ggplot2::ggplot(data, ggplot2::aes(x = .data[["k"]], y = .data[["pac"]])) +
    ggplot2::geom_col(
      fill = unname(.bulkmae_qualitative[["blue"]]),
      colour = "#222222",
      linewidth = 0.2,
      width = 0.72
    ) +
    ggplot2::labs(
      x = "k",
      y = "PAC",
      alt = "Proportion of ambiguous pairwise consensus values at each k."
    ) +
    theme_bulkmae()
  .plot_with_dimensions(plot, width = 8, height = 6.5)
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
      y = .plot_fraction_axis_name(data),
      alt = "Stacked bars of estimated cell-type values for each sample."
    ) +
    .plot_deconv_theme()
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
      outlier.shape = NA,
      linewidth = 0.3,
      width = 0.7
    ) +
    ggplot2::geom_jitter(
      ggplot2::aes(colour = .data[["group"]]),
      position = ggplot2::position_jitterdodge(
        jitter.width = 0.18,
        dodge.width = 0.75
      ),
      size = 0.7,
      alpha = 0.8,
      show.legend = FALSE
    ) +
    ggplot2::scale_fill_manual(values = fills, name = "Group") +
    ggplot2::scale_colour_manual(values = fills, guide = "none") +
    ggplot2::labs(
      x = "Cell type",
      y = "Estimated fraction",
      alt = "Boxplots with jittered points of estimated cell-type fractions grouped by a sample annotation."
    ) +
    .plot_deconv_theme()
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
    ggplot2::scale_fill_viridis_c(name = .plot_fraction_axis_name(data)) +
    ggplot2::labs(
      x = "Sample",
      y = "Cell type",
      alt = "A heatmap of estimated cell-type values across samples."
    ) +
    .plot_deconv_theme()
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

.plot_significance_stars <- function(p_value) {
  ifelse(
    !is.finite(p_value),
    "",
    ifelse(
      p_value < 0.001,
      "***",
      ifelse(p_value < 0.01, "**", ifelse(p_value < 0.05, "*", ""))
    )
  )
}

.plot_deconv_theme <- function() {
  theme_bulkmae() +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)
    )
}

.plot_fraction_axis_name <- function(data) {
  totals <- tapply(data$fraction, data$sample, sum)
  if (all(abs(totals - 1) <= 0.05)) "Estimated fraction" else "Estimated value"
}

.plot_annotation_frame <- function(annotation, ids, argument) {
  if (is.atomic(annotation) && is.null(dim(annotation))) {
    annotation <- .plot_align_named_labels(annotation, ids, argument)
    return(data.frame(
      sample = ids,
      track = "Annotation",
      label = as.character(unname(annotation[ids])),
      stringsAsFactors = FALSE
    ))
  }
  if (is.data.frame(annotation) || is.matrix(annotation)) {
    annotation <- as.data.frame(annotation, optional = TRUE)
    if (!ncol(annotation)) {
      stop("`", argument, "` must contain at least one column.", call. = FALSE)
    }
    if (is.null(rownames(annotation))) {
      stop("`", argument, "` must have sample row names.", call. = FALSE)
    }
    .plot_assert_ids(rownames(annotation), paste0("`", argument, "` row names"))
    .plot_assert_ids(names(annotation), paste0("`", argument, "` column names"))
    if (!setequal(rownames(annotation), ids)) {
      stop(
        "`", argument, "` row names must match the plotted samples exactly.",
        call. = FALSE
      )
    }
    annotation <- annotation[ids, , drop = FALSE]
    rows <- lapply(names(annotation), function(track) {
      values <- as.character(annotation[[track]])
      if (anyNA(values) || any(!nzchar(values))) {
        stop("`", argument, "` cannot contain missing or empty values.", call. = FALSE)
      }
      data.frame(
        sample = ids,
        track = track,
        label = values,
        stringsAsFactors = FALSE
      )
    })
    return(do.call(rbind, rows))
  }
  if (is.list(annotation)) {
    if (
      is.null(names(annotation)) || anyNA(names(annotation)) ||
        any(!nzchar(names(annotation))) || anyDuplicated(names(annotation))
    ) {
      stop("`", argument, "` list names must be unique and non-empty.", call. = FALSE)
    }
    rows <- lapply(names(annotation), function(track) {
      value <- .plot_align_named_labels(
        annotation[[track]],
        ids,
        paste0(argument, "$", track)
      )
      data.frame(
        sample = ids,
        track = track,
        label = as.character(unname(value[ids])),
        stringsAsFactors = FALSE
      )
    })
    return(do.call(rbind, rows))
  }
  stop(
    "`", argument, "` must be a named vector, data frame, or named list.",
    call. = FALSE
  )
}

.plot_annotation_colours <- function(labels) {
  labels <- unique(as.character(labels))
  n <- length(labels)
  if (n <= length(.bulkmae_qualitative)) {
    return(stats::setNames(unname(.bulkmae_qualitative[seq_len(n)]), labels))
  }
  hues <- grDevices::hcl(
    h = seq(15, 375, length.out = n + 1L)[seq_len(n)],
    c = 65,
    l = 50
  )
  stats::setNames(hues, labels)
}

.plot_annotation_style <- function(tracks, ids) {
  tracks$sample <- factor(tracks$sample, levels = ids)
  tracks$track <- factor(tracks$track, levels = rev(unique(as.character(tracks$track))))
  tracks$legend <- ifelse(
    length(unique(tracks$track)) == 1L,
    as.character(tracks$label),
    paste0(tracks$track, ": ", tracks$label)
  )
  fills <- .plot_annotation_colours(tracks$legend)
  tracks$fill_colour <- unname(fills[tracks$legend])
  list(data = tracks, legend_fills = fills)
}

.plot_as_hclust <- function(tree) {
  if (inherits(tree, "hclust")) {
    return(tree)
  }
  if (inherits(tree, "dendrogram")) {
    return(stats::as.hclust(tree))
  }
  stop("`dendrogram` must be an hclust, dendrogram, or WGCNA result list.", call. = FALSE)
}

.plot_hclust_segments <- function(tree) {
  n <- length(tree$order)
  leaf_x <- integer(n)
  leaf_x[tree$order] <- seq_len(n)
  merge <- tree$merge
  height <- tree$height
  node_x <- numeric(nrow(merge))
  rows <- vector("list", nrow(merge))
  for (i in seq_len(nrow(merge))) {
    child <- function(index) {
      if (index < 0L) {
        c(x = leaf_x[[-index]], y = 0)
      } else {
        c(x = node_x[[index]], y = height[[index]])
      }
    }
    left <- child(merge[i, 1L])
    right <- child(merge[i, 2L])
    node_x[[i]] <- mean(c(left[["x"]], right[["x"]]))
    rows[[i]] <- data.frame(
      x = c(left[["x"]], right[["x"]], left[["x"]]),
      xend = c(left[["x"]], right[["x"]], right[["x"]]),
      y = c(left[["y"]], right[["y"]], height[[i]]),
      yend = c(height[[i]], height[[i]], height[[i]])
    )
  }
  do.call(rbind, rows)
}
