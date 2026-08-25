#' Plot one classic GSEA running-score profile
#'
#' Displays the running enrichment score, gene-set hit positions with a rank
#' heat band, and the ranked metric as three intrinsic parts of one standard
#' ggplot. The parts share a rank-position axis and use fixed relative panel
#' heights of 0.50, 0.20, and 0.30. This function does not call private
#' enrichplot helpers or rerun GSEA.
#'
#' A native clusterProfiler gseaResult supplies its ranked vector, gene sets,
#' exponent, and leading-edge membership automatically. For an fgsea or
#' generic GSEA table, ranks, gene_sets, and exponent are all required. Those
#' inputs must reproduce the reported enrichment score.
#'
#' @param result A native clusterProfiler gseaResult or a GSEA result table.
#' @param term One explicit term identifier to display.
#' @param ranks Optional decreasing, gene-named, finite numeric rank vector.
#' @param gene_sets Optional uniquely named list mapping terms to genes.
#' @param exponent Optional finite non-negative GSEA weighting exponent.
#' @param term_label Optional single non-empty display label for term.
#'
#' @return An unprinted standard ggplot object carrying recommended physical
#'   dimensions for [plot_save()].
#' @family plotting
#' @export
plot_gsea_classic <- function(
    result,
    term,
    ranks = NULL,
    gene_sets = NULL,
    exponent = NULL,
    term_label = NULL
) {
  term <- .plot_gsea_terms(term, singular = TRUE)
  view <- .plot_gsea_result_view(result)
  row <- .plot_gsea_select(view$data, term)
  native <- identical(view$backend, "clusterProfiler")

  if (!native && is.null(ranks)) {
    stop("'ranks' is required for a tabular GSEA result.", call. = FALSE)
  }
  if (!native && is.null(gene_sets)) {
    stop("'gene_sets' is required for a tabular GSEA result.", call. = FALSE)
  }
  if (!native && is.null(exponent)) {
    stop("'exponent' is required for a tabular GSEA result.", call. = FALSE)
  }
  if (is.null(ranks)) ranks <- view$ranks
  if (is.null(gene_sets)) gene_sets <- view$gene_sets
  if (is.null(exponent)) exponent <- view$exponent

  ranks <- .plot_gsea_ranks(ranks)
  gene_sets <- .plot_gsea_gene_sets(gene_sets, term)
  exponent <- .plot_gsea_exponent(exponent)
  running <- .plot_gsea_running_score(ranks, gene_sets[[term]], exponent)
  .plot_gsea_assert_es(row$enrichment_score, running$enrichment_score, term)

  raw_label <- if (is.null(term_label)) row$term_label else term_label
  raw_label <- .plot_gsea_single_label(raw_label, "term_label")
  display_label <- if (is.null(term_label)) {
    .plot_gsea_classic_default_label(raw_label)
  } else {
    .plot_gsea_wrap_label(raw_label, width = 35L)
  }
  panels <- c("running", "hits", "ranked")
  bands <- list(
    score = c(0, 0.50),
    hits = c(1, 1.20),
    metric = c(2, 2.30)
  )
  panel_factor <- function(x) factor(x, levels = panels)

  score_data_range <- .plot_gsea_clean_range(
    range(c(running$data$running_score, 0))
  )
  score_span <- diff(score_data_range)
  score_range <- score_data_range + c(-0.05, 0.05) * score_span
  metric_range <- range(c(ranks, 0))
  score_data <- data.frame(
    position = running$data$position,
    value = .plot_gsea_map_band(
      running$data$running_score,
      score_range,
      bands$score
    ),
    running_score = running$data$running_score,
    panel = panel_factor(panels[[1L]])
  )
  score_zero <- .plot_gsea_map_band(0, score_range, bands$score)
  annotation_y <- score_data_range[[2L]] - 0.12 * score_span
  annotation_data <- data.frame(
    position = length(ranks) * 0.96,
    value = .plot_gsea_map_band(
      annotation_y,
      score_range,
      bands$score
    ),
    label = .plot_gsea_classic_statistics(row),
    panel = panel_factor(panels[[1L]])
  )
  hit_data <- data.frame(
    position = running$data$position[running$data$hit],
    y = .plot_gsea_map_band(0, c(0, 1), bands$hits),
    yend = .plot_gsea_map_band(1, c(0, 1), bands$hits),
    panel = panel_factor(panels[[2L]])
  )
  rank_band <- .plot_gsea_rank_band(ranks, strip_bins = 200L)
  heat_data <- data.frame(
    xmin = rank_band$xmin,
    xmax = rank_band$xmax,
    ymin = .plot_gsea_map_band(0, c(0, 1), bands$hits),
    ymax = .plot_gsea_map_band(0.30, c(0, 1), bands$hits),
    rank_metric = rank_band$rank_metric,
    panel = panel_factor(panels[[2L]])
  )
  metric_polygon <- data.frame(
    position = c(1, seq_along(ranks), length(ranks)),
    value = .plot_gsea_map_band(
      c(0, as.numeric(ranks), 0),
      metric_range,
      bands$metric
    ),
    panel = panel_factor(panels[[3L]])
  )
  metric_zero <- .plot_gsea_map_band(0, metric_range, bands$metric)
  hit_data$panel <- panel_factor(panels[[2L]])
  blank_data <- do.call(rbind, lapply(seq_along(panels), function(index) {
    band <- bands[[index]]
    data.frame(
      position = c(1, length(ranks)),
      value = band,
      panel = panel_factor(panels[[index]])
    )
  }))
  ticks <- .plot_gsea_classic_ticks(
    score_range,
    bands$score,
    metric_range,
    bands$metric
  )

  plot <- ggplot2::ggplot(
    blank_data,
    ggplot2::aes(x = .data[["position"]], y = .data[["value"]])
  ) +
    ggplot2::geom_blank() +
    ggplot2::geom_line(
      data = score_data,
      ggplot2::aes(
        x = .data[["position"]],
        y = .data[["value"]],
        colour = .data[["running_score"]]
      ),
      linewidth = 0.6,
      show.legend = FALSE
    ) +
    ggplot2::geom_hline(
      data = data.frame(
        panel = panel_factor(panels[[1L]]),
        zero = score_zero
      ),
      ggplot2::aes(yintercept = .data[["zero"]]),
      inherit.aes = FALSE,
      colour = "black",
      linetype = "dashed",
      linewidth = 0.3
    ) +
    ggplot2::geom_text(
      data = annotation_data,
      ggplot2::aes(
        x = .data[["position"]],
        y = .data[["value"]],
        label = .data[["label"]]
      ),
      inherit.aes = FALSE,
      hjust = 1,
      vjust = 1,
      colour = "grey20",
      fontface = "italic",
      lineheight = 0.95,
      size = .bulkmae_text_size_pt,
      size.unit = "pt"
    ) +
    ggplot2::geom_rect(
      data = heat_data,
      ggplot2::aes(
        xmin = .data[["xmin"]],
        xmax = .data[["xmax"]],
        ymin = .data[["ymin"]],
        ymax = .data[["ymax"]],
        fill = .data[["rank_metric"]]
      ),
      inherit.aes = FALSE,
      colour = NA,
      alpha = 0.8,
      show.legend = FALSE
    ) +
    ggplot2::geom_segment(
      data = hit_data,
      ggplot2::aes(
        x = .data[["position"]],
        xend = .data[["position"]],
        y = .data[["y"]],
        yend = .data[["yend"]]
      ),
      inherit.aes = FALSE,
      colour = "black",
      linewidth = 0.2
    ) +
    ggplot2::geom_polygon(
      data = metric_polygon,
      ggplot2::aes(
        x = .data[["position"]],
        y = .data[["value"]]
      ),
      inherit.aes = FALSE,
      fill = "grey70"
    ) +
    ggplot2::geom_hline(
      data = data.frame(
        panel = panel_factor(panels[[3L]]),
        zero = metric_zero
      ),
      ggplot2::aes(yintercept = .data[["zero"]]),
      inherit.aes = FALSE,
      colour = "black",
      linetype = "dashed",
      linewidth = 0.35
    ) +
    ggplot2::scale_colour_gradient(
      low = "#76BA99",
      high = "#EB4747",
      guide = "none"
    ) +
    ggplot2::scale_fill_gradient2(
      low = "#08519C",
      mid = "white",
      high = "#A50F15",
      midpoint = 0,
      guide = "none"
    ) +
    ggplot2::scale_x_continuous(
      breaks = pretty(c(1, length(ranks)), n = 4),
      expand = ggplot2::expansion(mult = 0)
    ) +
    ggplot2::scale_y_continuous(
      breaks = ticks$breaks,
      labels = ticks$labels,
      expand = ggplot2::expansion(mult = 0)
    ) +
    ggplot2::facet_grid(
      rows = ggplot2::vars(.data[["panel"]]),
      scales = "free_y",
      space = "free_y",
      switch = "y",
      labeller = ggplot2::labeller(panel = c(
        running = "Enrichment Score",
        hits = "",
        ranked = "Ranked List"
      ))
    ) +
    ggplot2::labs(
      x = "Rank in Ordered Dataset",
      y = NULL,
      title = display_label,
      alt = paste0(
        "A classic GSEA plot for ", raw_label,
        ", showing the running enrichment score, gene-set hit positions, ",
        "a rank heat band, and the ranked metric."
      )
    ) +
    ggplot2::theme_bw(base_size = 6) +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      panel.border = ggplot2::element_rect(
        colour = "black",
        fill = NA,
        linewidth = 0.25
      ),
      panel.spacing.y = grid::unit(0, "cm"),
      strip.placement = "outside",
      strip.background = ggplot2::element_blank(),
      strip.text.y.left = ggplot2::element_text(
        angle = 90,
        face = "plain",
        size = 6
      ),
      axis.text = ggplot2::element_text(colour = "black"),
      plot.title = ggplot2::element_text(
        hjust = 0.5,
        face = "bold",
        size = 6
      ),
      plot.margin = ggplot2::margin(0, 0, 0, 0),
      legend.position = "none"
    )
  .plot_with_dimensions(plot, width = 5.9, height = 5.3)
}

#' Plot GSEA rank-metric distributions as ridges
#'
#' Draws density ridges and gene-position barcodes for explicitly selected
#' terms. This visualizes where member genes lie in the ranked metric; it does
#' not refit GSEA or select terms.
#'
#' @param result A native clusterProfiler gseaResult or a GSEA result table.
#' @param terms Unique explicit term identifiers, in top-to-bottom order.
#' @param ranks Optional decreasing, gene-named, finite numeric rank vector.
#' @param gene_sets Optional uniquely named list mapping terms to genes.
#' @param membership Genes used for each density: leading_edge or gene_set.
#' @param term_labels Optional complete term-named labels for selected terms.
#' @param show_statistics Show available NES and p-value text.
#'
#' @return An unprinted standard ggplot object carrying recommended physical
#'   dimensions for [plot_save()].
#' @family plotting
#' @export
plot_gsea_ridge <- function(
    result,
    terms,
    ranks = NULL,
    gene_sets = NULL,
    membership = c("leading_edge", "gene_set"),
    term_labels = NULL,
    show_statistics = FALSE
) {
  .require_backend("ggridges", "to draw GSEA ridge plots")
  terms <- .plot_gsea_terms(terms)
  membership <- match.arg(membership)
  .plot_assert_flag(show_statistics, "show_statistics")
  view <- .plot_gsea_result_view(result)
  selected <- .plot_gsea_select(view$data, terms)
  native <- identical(view$backend, "clusterProfiler")
  if (!native && is.null(ranks)) {
    stop("'ranks' is required for a tabular GSEA result.", call. = FALSE)
  }
  if (is.null(ranks)) ranks <- view$ranks
  ranks <- .plot_gsea_ranks(ranks)

  if (membership == "gene_set") {
    if (!native && is.null(gene_sets)) {
      stop(
        "'gene_sets' is required with membership = 'gene_set' for a tabular GSEA result.",
        call. = FALSE
      )
    }
    if (is.null(gene_sets)) gene_sets <- view$gene_sets
    gene_sets <- .plot_gsea_gene_sets(gene_sets, terms)
    memberships <- gene_sets[terms]
  } else {
    memberships <- stats::setNames(selected$leading_edge, selected$term_id)
    missing_membership <- !lengths(memberships)
    if (any(missing_membership)) {
      stop(
        "Leading-edge membership is unavailable for terms: ",
        paste(names(memberships)[missing_membership], collapse = ", "),
        ".",
        call. = FALSE
      )
    }
  }
  memberships <- lapply(terms, function(term) {
    members <- as.character(memberships[[term]])
    .plot_assert_ids(members, paste0("Membership genes for '", term, "'"))
    if (membership == "leading_edge") {
      unknown <- setdiff(members, names(ranks))
      if (length(unknown)) {
        stop(
          "Leading-edge genes are absent from 'ranks' for term '", term,
          "': ", paste(unknown, collapse = ", "), ".",
          call. = FALSE
        )
      }
    }
    intersect(members, names(ranks))
  })
  names(memberships) <- terms
  too_small <- lengths(memberships) < 2L
  if (any(too_small)) {
    stop(
      "At least two ranked membership genes are required for terms: ",
      paste(terms[too_small], collapse = ", "),
      ".",
      call. = FALSE
    )
  }
  constant_membership <- vapply(memberships, function(members) {
    length(unique(unname(ranks[members]))) < 2L
  }, logical(1))
  if (any(constant_membership)) {
    stop(
      "At least two distinct rank metric values are required for terms: ",
      paste(terms[constant_membership], collapse = ", "),
      ".",
      call. = FALSE
    )
  }

  labels <- if (is.null(term_labels)) {
    vapply(seq_along(terms), function(index) {
      .plot_gsea_ridge_default_label(
        selected$term_id[[index]],
        selected$term_label[[index]]
      )
    }, character(1))
  } else {
    .plot_enrichment_labels(selected, term_labels)
  }
  labels <- vapply(
    labels,
    .plot_gsea_wrap_label,
    character(1),
    width = 40L
  )
  palette <- stats::setNames(.plot_gsea_palette(length(terms)), terms)
  baselines <- rev(seq_along(terms))
  names(baselines) <- terms
  density_parts <- vector("list", length(terms))
  hit_parts <- vector("list", length(terms))
  statistic_parts <- vector("list", length(terms))
  ridge_height <- 0.36
  for (index in seq_along(terms)) {
    term <- terms[[index]]
    values <- unname(ranks[memberships[[term]]])
    density <- .plot_gsea_density(values)
    height <- density$y / max(density$y) * ridge_height
    density_parts[[index]] <- data.frame(
      pathway = labels[[index]],
      term_id = term,
      x = density$x,
      y = baselines[[term]],
      height = height,
      colour = palette[[term]],
      stringsAsFactors = FALSE
    )
    hit_parts[[index]] <- data.frame(
      pathway = labels[[index]],
      term_id = term,
      x = values,
      y = baselines[[term]],
      colour = palette[[term]],
      stringsAsFactors = FALSE
    )
    statistic_parts[[index]] <- data.frame(
      pathway = labels[[index]],
      term_id = term,
      y = baselines[[term]] + ridge_height * 0.92,
      label = .plot_gsea_ridge_statistics(selected[index, , drop = FALSE]),
      colour = palette[[term]],
      stringsAsFactors = FALSE
    )
  }
  density_data <- do.call(rbind, density_parts)
  hit_data <- do.call(rbind, hit_parts)
  statistic_data <- do.call(rbind, statistic_parts)
  x_min <- min(density_data$x, na.rm = TRUE)
  x_max <- max(density_data$x, na.rm = TRUE)
  if (all(hit_data$x >= 0)) x_min <- min(0, x_min)
  x_pad <- max((x_max - x_min) * 0.04, 0.15)
  x_limits <- c(x_min - x_pad * 0.2, x_max + x_pad)
  statistic_data$x <- x_limits[[2L]] - diff(x_limits) * 0.02
  statistic_data <- statistic_data[nzchar(statistic_data$label), , drop = FALSE]

  plot <- ggplot2::ggplot() +
    ggridges::geom_ridgeline(
      data = density_data,
      ggplot2::aes(
        x = .data[["x"]],
        y = .data[["y"]],
        height = .data[["height"]],
        group = .data[["pathway"]],
        colour = .data[["colour"]],
        fill = .data[["colour"]]
      ),
      linewidth = 0.35,
      alpha = 0.18,
      scale = 1,
      min_height = 0
    ) +
    ggplot2::geom_segment(
      data = hit_data,
      ggplot2::aes(
        x = .data[["x"]],
        xend = .data[["x"]],
        y = .data[["y"]] - 0.30,
        yend = .data[["y"]] - 0.09,
        colour = .data[["colour"]]
      ),
      inherit.aes = FALSE,
      linewidth = 0.20,
      alpha = 0.52,
      lineend = "round"
    ) +
    ggplot2::scale_colour_identity() +
    ggplot2::scale_fill_identity() +
    ggplot2::scale_x_continuous(
      breaks = pretty(x_limits, n = 3),
      expand = ggplot2::expansion(mult = c(0.01, 0.02))
    ) +
    ggplot2::scale_y_continuous(
      breaks = baselines,
      labels = labels,
      limits = c(min(baselines) - 0.45, max(baselines) + 0.58),
      expand = ggplot2::expansion(mult = 0)
    ) +
    ggplot2::coord_cartesian(xlim = x_limits, clip = "off") +
    ggplot2::labs(
      x = NULL,
      y = NULL,
      alt = paste0(
        "Density ridges and gene-position barcodes for ", length(terms),
        " explicitly selected GSEA terms using ",
        if (membership == "leading_edge") "leading-edge" else "gene-set",
        " membership."
      )
    ) +
    ggplot2::theme_classic(base_size = 6) +
    ggplot2::theme(
      panel.border = ggplot2::element_rect(
        colour = "black",
        fill = NA,
        linewidth = 0.25
      ),
      axis.line.y = ggplot2::element_blank(),
      axis.ticks.y = ggplot2::element_blank(),
      axis.ticks.length = grid::unit(1.2, "mm"),
      axis.text.y = ggplot2::element_text(
        size = 6,
        colour = "black",
        hjust = 1,
        margin = ggplot2::margin(r = 1.3, unit = "mm")
      ),
      axis.text.x = ggplot2::element_text(size = 6, colour = "black"),
      plot.margin = ggplot2::margin(
        t = 1.5,
        r = 2.2,
        b = 1.5,
        l = 1.2,
        unit = "mm"
      ),
      legend.position = "none"
    )
  if (show_statistics && nrow(statistic_data)) {
    plot <- plot + ggplot2::geom_text(
      data = statistic_data,
      ggplot2::aes(
        x = .data[["x"]],
        y = .data[["y"]],
        label = .data[["label"]]
      ),
      inherit.aes = FALSE,
      hjust = 1,
      vjust = 1,
      size = .bulkmae_text_size_pt,
      size.unit = "pt",
      lineheight = 0.95
    )
  }
  .plot_with_dimensions(
    plot,
    width = 12,
    height = 5
  )
}

.plot_gsea_result_view <- function(result) {
  if (methods::is(result, "gseaResult")) {
    table <- as.data.frame(result, optional = TRUE)
    if (!nrow(table)) stop("'result' contains no GSEA terms.", call. = FALSE)
    .plot_require_columns(
      table,
      c(
        "ID", "Description", "enrichmentScore", "NES", "pvalue",
        "p.adjust", "core_enrichment"
      ),
      "result"
    )
    data <- .plot_gsea_make_view(
      term_id = table$ID,
      term_label = table$Description,
      enrichment_score = table$enrichmentScore,
      normalized_score = table$NES,
      p_value = table$pvalue,
      adjusted_p_value = table$p.adjust,
      q_value = if ("qvalue" %in% names(table)) {
        table$qvalue
      } else {
        rep(NA_real_, nrow(table))
      },
      leading_edge = table$core_enrichment
    )
    return(list(
      backend = "clusterProfiler",
      data = data,
      ranks = methods::slot(result, "geneList"),
      gene_sets = methods::slot(result, "geneSets"),
      exponent = methods::slot(result, "params")[["exponent"]]
    ))
  }
  if (!(is.data.frame(result) || is.matrix(result))) {
    stop(
      "Unsupported 'result'; use a clusterProfiler gseaResult or GSEA result table.",
      call. = FALSE
    )
  }
  table <- as.data.frame(result, optional = TRUE)
  if (!nrow(table)) stop("'result' contains no GSEA terms.", call. = FALSE)
  id_column <- .plot_gsea_find_column(
    table,
    c("pathway", "ID", "term_id", "term"),
    "term identifier"
  )
  score_column <- .plot_gsea_find_column(
    table,
    c("ES", "enrichmentScore", "enrichment_score"),
    "enrichment score"
  )
  label_column <- .plot_gsea_find_column(
    table,
    c("Description", "description", "term_label"),
    "term label",
    required = FALSE
  )
  normalized_column <- .plot_gsea_find_column(
    table,
    c("NES", "normalized_score"),
    "normalized enrichment score",
    required = FALSE
  )
  p_column <- .plot_gsea_find_column(
    table,
    c("pval", "pvalue", "p_value"),
    "p-value",
    required = FALSE
  )
  adjusted_column <- .plot_gsea_find_column(
    table,
    c("padj", "p.adjust", "adjusted_p_value"),
    "adjusted p-value",
    required = FALSE
  )
  q_column <- .plot_gsea_find_column(
    table,
    c("qvalue", "qvalues", "q_value"),
    "q-value",
    required = FALSE
  )
  leading_column <- .plot_gsea_find_column(
    table,
    c("leadingEdge", "leading_edge", "core_enrichment"),
    "leading-edge membership",
    required = FALSE
  )
  row_count <- nrow(table)
  data <- .plot_gsea_make_view(
    term_id = table[[id_column]],
    term_label = if (is.null(label_column)) table[[id_column]] else
      table[[label_column]],
    enrichment_score = table[[score_column]],
    normalized_score = if (is.null(normalized_column)) {
      rep(NA_real_, row_count)
    } else {
      table[[normalized_column]]
    },
    p_value = if (is.null(p_column)) rep(NA_real_, row_count) else
      table[[p_column]],
    adjusted_p_value = if (is.null(adjusted_column)) {
      rep(NA_real_, row_count)
    } else {
      table[[adjusted_column]]
    },
    q_value = if (is.null(q_column)) rep(NA_real_, row_count) else
      table[[q_column]],
    leading_edge = if (is.null(leading_column)) {
      rep(list(character()), row_count)
    } else {
      table[[leading_column]]
    }
  )
  list(
    backend = "table",
    data = data,
    ranks = NULL,
    gene_sets = NULL,
    exponent = NULL
  )
}

.plot_gsea_make_view <- function(
    term_id,
    term_label,
    enrichment_score,
    normalized_score,
    p_value,
    adjusted_p_value,
    q_value,
    leading_edge
) {
  term_id <- as.character(term_id)
  .plot_assert_ids(term_id, "GSEA term identifiers")
  term_label <- as.character(term_label)
  if (length(term_label) != length(term_id) || anyNA(term_label) ||
      any(!nzchar(term_label))) {
    stop("GSEA term labels must be non-missing and non-empty.", call. = FALSE)
  }
  .plot_gsea_numeric(enrichment_score, length(term_id), "Enrichment scores")
  .plot_gsea_optional_numeric(
    normalized_score,
    length(term_id),
    "Normalized enrichment scores"
  )
  if (is.logical(q_value) && all(is.na(q_value))) {
    q_value <- rep(NA_real_, length(q_value))
  }
  .plot_enrichment_probabilities(p_value, "Raw p-values", allow_all_missing = TRUE)
  .plot_enrichment_probabilities(
    adjusted_p_value,
    "Adjusted p-values",
    allow_all_missing = TRUE
  )
  .plot_enrichment_probabilities(
    q_value,
    "Q-values",
    allow_all_missing = TRUE
  )
  if (length(p_value) != length(term_id) ||
      length(adjusted_p_value) != length(term_id) ||
      length(q_value) != length(term_id)) {
    stop("GSEA evidence columns must match the number of terms.", call. = FALSE)
  }
  leading_edge <- .plot_gsea_memberships(leading_edge, length(term_id))
  data.frame(
    term_id = term_id,
    term_label = term_label,
    enrichment_score = as.numeric(enrichment_score),
    normalized_score = as.numeric(normalized_score),
    p_value = as.numeric(p_value),
    adjusted_p_value = as.numeric(adjusted_p_value),
    q_value = as.numeric(q_value),
    leading_edge = I(leading_edge),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}

.plot_gsea_find_column <- function(table, candidates, label, required = TRUE) {
  found <- intersect(candidates, names(table))
  if (length(found) > 1L) {
    stop(
      "'result' has ambiguous ", label, " columns: ",
      paste(found, collapse = ", "), ".",
      call. = FALSE
    )
  }
  if (!length(found)) {
    if (!required) return(NULL)
    stop(
      "'result' is missing a recognized ", label, " column (",
      paste(candidates, collapse = ", "), ").",
      call. = FALSE
    )
  }
  found[[1L]]
}

.plot_gsea_numeric <- function(values, length_expected, label) {
  if (!is.numeric(values) || length(values) != length_expected ||
      any(!is.finite(values))) {
    stop(label, " must contain one finite numeric value per term.", call. = FALSE)
  }
  invisible(values)
}

.plot_gsea_optional_numeric <- function(values, length_expected, label) {
  if (!is.numeric(values) || length(values) != length_expected ||
      any(is.infinite(values))) {
    stop(label, " must contain finite numeric values or missing values.",
         call. = FALSE)
  }
  invisible(values)
}

.plot_gsea_memberships <- function(values, length_expected) {
  if (is.list(values)) {
    if (length(values) != length_expected) {
      stop("Leading-edge memberships must match the number of terms.",
           call. = FALSE)
    }
    output <- lapply(values, as.character)
  } else {
    values <- as.character(values)
    if (length(values) != length_expected) {
      stop("Leading-edge memberships must match the number of terms.",
           call. = FALSE)
    }
    output <- lapply(values, function(value) {
      if (is.na(value) || !nzchar(value)) character() else
        strsplit(value, "/", fixed = TRUE)[[1L]]
    })
  }
  lapply(output, function(members) {
    if (anyNA(members) || any(!nzchar(members)) || anyDuplicated(members)) {
      stop(
        "Leading-edge gene identifiers must be unique, non-missing, and non-empty.",
        call. = FALSE
      )
    }
    members
  })
}

.plot_gsea_terms <- function(terms, singular = FALSE) {
  if (is.null(terms) || !length(terms)) {
    stop(
      if (singular) "'term' must explicitly identify one term." else
        "'terms' must explicitly contain at least one term identifier.",
      call. = FALSE
    )
  }
  terms <- as.character(terms)
  .plot_assert_ids(terms, if (singular) "'term'" else "'terms'")
  if (singular && length(terms) != 1L) {
    stop("'term' must explicitly identify one term.", call. = FALSE)
  }
  terms
}

.plot_gsea_select <- function(data, terms) {
  unknown <- setdiff(terms, data$term_id)
  if (length(unknown)) {
    stop("Unknown GSEA terms: ", paste(unknown, collapse = ", "), ".",
         call. = FALSE)
  }
  data[match(terms, data$term_id), , drop = FALSE]
}

.plot_gsea_ranks <- function(ranks) {
  if (!is.numeric(ranks) || is.null(names(ranks))) {
    stop("'ranks' must be a gene-named numeric vector.", call. = FALSE)
  }
  if (length(ranks) < 2L) {
    stop("'ranks' must contain at least two genes.", call. = FALSE)
  }
  .plot_assert_ids(names(ranks), "Names of 'ranks'")
  .plot_assert_finite_numeric(ranks, "'ranks'")
  if (any(diff(unname(ranks)) > 0)) {
    stop("'ranks' must be sorted in decreasing order.", call. = FALSE)
  }
  ranks
}

.plot_gsea_gene_sets <- function(gene_sets, terms) {
  if (!is.list(gene_sets) || is.null(names(gene_sets))) {
    stop("'gene_sets' must be a term-named list.", call. = FALSE)
  }
  .plot_assert_ids(names(gene_sets), "Names of 'gene_sets'")
  unknown <- setdiff(terms, names(gene_sets))
  if (length(unknown)) {
    stop(
      "'gene_sets' is missing selected terms: ",
      paste(unknown, collapse = ", "), ".",
      call. = FALSE
    )
  }
  for (term in terms) {
    members <- as.character(gene_sets[[term]])
    .plot_assert_ids(members, paste0("Genes in gene set '", term, "'"))
    gene_sets[[term]] <- members
  }
  gene_sets
}

.plot_gsea_exponent <- function(exponent) {
  if (!is.numeric(exponent) || length(exponent) != 1L || is.na(exponent) ||
      !is.finite(exponent) || exponent < 0) {
    stop("'exponent' must be one finite non-negative number.", call. = FALSE)
  }
  as.numeric(exponent)
}

.plot_gsea_running_score <- function(ranks, members, exponent) {
  ranks <- .plot_gsea_ranks(ranks)
  exponent <- .plot_gsea_exponent(exponent)
  members <- as.character(members)
  .plot_assert_ids(members, "Gene-set members")
  members <- intersect(members, names(ranks))
  hit <- names(ranks) %in% members
  hit_count <- sum(hit)
  if (!hit_count) {
    stop("The gene set has no members represented in 'ranks'.", call. = FALSE)
  }
  if (hit_count == length(ranks)) {
    stop("The gene set cannot contain every gene in 'ranks'.", call. = FALSE)
  }
  hit_weight <- numeric(length(ranks))
  hit_weight[hit] <- abs(ranks[hit])^exponent
  normalizer <- sum(hit_weight)
  if (!is.finite(normalizer) || normalizer <= 0) {
    stop(
      "The weighted gene-set hits sum to zero; the running score is undefined.",
      call. = FALSE
    )
  }
  cumulative_hit <- cumsum(hit_weight / normalizer)
  cumulative_miss <- cumsum((!hit) / (length(ranks) - hit_count))
  running_score <- cumulative_hit - cumulative_miss
  maximum <- max(running_score)
  minimum <- min(running_score)
  enrichment_score <- if (abs(maximum) > abs(minimum)) maximum else minimum
  list(
    enrichment_score = enrichment_score,
    data = data.frame(
      position = seq_along(ranks),
      gene_id = names(ranks),
      rank_metric = as.numeric(ranks),
      hit = hit,
      running_score = running_score,
      stringsAsFactors = FALSE
    )
  )
}

.plot_gsea_assert_es <- function(reported, computed, term) {
  tolerance <- 1e-6 * max(1, abs(reported))
  if (abs(reported - computed) > tolerance) {
    stop(
      "'ranks', 'gene_sets', and 'exponent' do not reproduce the reported ES ",
      "for term '", term, "' (reported ", signif(reported, 7),
      ", computed ", signif(computed, 7), ").",
      call. = FALSE
    )
  }
  invisible(computed)
}

.plot_gsea_single_label <- function(label, argument) {
  if (!is.character(label) || length(label) != 1L || is.na(label) ||
      !nzchar(label)) {
    stop("'", argument, "' must be one non-empty character value.",
         call. = FALSE)
  }
  label
}

.plot_gsea_map_band <- function(values, source_range, target_range) {
  source_span <- diff(source_range)
  if (!is.finite(source_span) || source_span == 0) {
    return(rep(mean(target_range), length(values)))
  }
  target_range[[1L]] +
    (values - source_range[[1L]]) / source_span * diff(target_range)
}

.plot_gsea_classic_ticks <- function(
    score_range,
    score_band,
    metric_range,
    metric_band
) {
  score_values <- .plot_gsea_pretty_values(score_range)
  metric_values <- .plot_gsea_pretty_values(metric_range)
  list(
    breaks = c(
      .plot_gsea_map_band(score_values, score_range, score_band),
      .plot_gsea_map_band(metric_values, metric_range, metric_band)
    ),
    labels = c(
      .plot_gsea_number_labels(score_values),
      .plot_gsea_number_labels(metric_values)
    )
  )
}

.plot_gsea_pretty_values <- function(value_range) {
  if (diff(value_range) == 0) return(value_range[[1L]])
  values <- pretty(value_range, n = 5)
  values <- values[values >= value_range[[1L]] & values <= value_range[[2L]]]
  zero <- if (value_range[[1L]] <= 0 && value_range[[2L]] >= 0) 0 else
    numeric()
  sort(unique(c(values, zero)))
}

.plot_gsea_clean_range <- function(value_range) {
  tolerance <- sqrt(.Machine$double.eps) * max(1, abs(value_range))
  value_range[abs(value_range) <= tolerance] <- 0
  value_range
}

.plot_gsea_wrap_label <- function(label, width) {
  stringr::str_wrap(label, width = width)
}

.plot_gsea_classic_default_label <- function(label) {
  pieces <- strsplit(label, "_", fixed = TRUE)[[1L]]
  if (length(pieces) > 1L) pieces <- pieces[-1L]
  label <- paste(stringr::str_to_title(pieces), collapse = " ")
  .plot_gsea_wrap_label(label, width = 35L)
}

.plot_gsea_ridge_default_label <- function(term_id, label) {
  if (is.na(label) || !nzchar(label) || identical(label, term_id)) {
    label <- gsub("^REACTOME_", "", term_id)
    label <- tools::toTitleCase(tolower(gsub("_", " ", label, fixed = TRUE)))
  }
  label <- gsub("[[:space:]]+", " ", trimws(label))
  reactome <- grepl("REACTOME", term_id, ignore.case = TRUE) ||
    grepl("^\\[Reactome\\]", label) ||
    grepl("reactome", label, ignore.case = TRUE)
  if (reactome && !grepl("^\\[Reactome\\]", label)) {
    label <- paste0("[Reactome] ", label)
  }
  label
}

.plot_gsea_number_labels <- function(values) {
  vapply(values, function(value) {
    if (value == 0) return("0")
    format(signif(value, 3), trim = TRUE, scientific = abs(value) >= 1e4)
  }, character(1))
}

.plot_gsea_classic_statistics <- function(row) {
  paste0(
    "NES: ", .plot_gsea_round_or_na(row$normalized_score[[1L]], 2L), "\n",
    "P.adj: ", .plot_gsea_threshold_probability(
      row$adjusted_p_value[[1L]]
    ), "\n",
    "FDR: ", .plot_gsea_threshold_probability(row$q_value[[1L]])
  )
}

.plot_gsea_ridge_statistics <- function(row) {
  paste0(
    "NES = ", if (is.finite(row$normalized_score[[1L]])) {
      sprintf("%.3f", row$normalized_score[[1L]])
    } else {
      "NA"
    }, "\n",
    "P.adj ", .plot_gsea_ridge_probability(
      row$adjusted_p_value[[1L]]
    ), "\n",
    "FDR ", .plot_gsea_ridge_probability(row$q_value[[1L]])
  )
}

.plot_gsea_round_or_na <- function(value, digits) {
  if (!is.finite(value)) return("NA")
  as.character(round(value, digits))
}

.plot_gsea_threshold_probability <- function(value) {
  if (!is.finite(value)) return("NA")
  if (value < 0.001) return("< 0.001")
  as.character(round(value, 3L))
}

.plot_gsea_ridge_probability <- function(value) {
  if (!is.finite(value)) return("NA")
  if (value < 0.001) return("< 0.001")
  paste0("= ", formatC(value, format = "f", digits = 3L))
}

.plot_gsea_rank_band <- function(ranks, strip_bins = 200L) {
  limits <- as.numeric(stats::quantile(ranks, c(0.1, 0.9), names = FALSE))
  trimmed <- pmin(pmax(as.numeric(ranks), limits[[1L]]), limits[[2L]])
  if (limits[[1L]] == limits[[2L]]) {
    return(data.frame(
      xmin = 1,
      xmax = length(ranks),
      rank_metric = limits[[1L]]
    ))
  }
  strip_bins <- max(as.integer(strip_bins), 2L)
  breaks <- seq(
    min(trimmed),
    max(trimmed),
    length.out = strip_bins + 1L
  )
  bins <- cut(
    trimmed,
    breaks = breaks,
    include.lowest = TRUE,
    labels = FALSE
  )
  positions <- split(seq_along(ranks), bins)
  data.frame(
    xmin = vapply(positions, min, numeric(1)),
    xmax = vapply(positions, max, numeric(1)),
    rank_metric = vapply(positions, function(index) {
      mean(trimmed[index])
    }, numeric(1))
  )
}

.plot_gsea_density <- function(values) {
  if (length(unique(values)) < 2L) {
    stop(
      "A ridge density requires at least two distinct rank metric values.",
      call. = FALSE
    )
  }
  stats::density(
    values,
    n = 512L,
    adjust = 1
  )
}

.plot_gsea_palette <- function(count) {
  base <- c(
    "#F6A282", "#4B63A7", "#25C6C2", "#FF6B5C",
    "#63D0EC", "#A78BCA", "#7EC87E"
  )
  if (count > length(base)) {
    stop(
      "The default GSEA ridge palette supports at most seven terms.",
      call. = FALSE
    )
  }
  base[seq_len(count)]
}
