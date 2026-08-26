#' Plot explicitly selected ORA terms as bubbles
#'
#' Displays one of three distinct over-representation metrics without
#' treating them as interchangeable: rich factor is overlap divided by the
#' background term size, gene ratio is overlap divided by the input list size,
#' and fold enrichment is gene ratio divided by background ratio. Point area
#' represents overlap count and colour represents finite `-log10()` evidence;
#' adjusted evidence is labelled compactly as `adj p`.
#'
#' @param result A clusterProfiler `enrichResult` or an ORA result data frame
#'   containing `ID`, `Description`, `GeneRatio`, `BgRatio`, `pvalue`,
#'   `p.adjust`, and `Count`.
#' @param terms Unique term identifiers to display, in the desired order.
#'   Terms are never selected automatically.
#' @param x ORA metric displayed on the horizontal axis.
#' @param p_value Evidence column to display. `"auto"` uses adjusted p-values
#'   only when every selected term has one; otherwise the whole plot uses raw
#'   p-values.
#' @param term_labels Optional complete term-named character vector describing
#'   the selected `terms`.
#'
#' @return An unprinted standard ggplot object carrying recommended physical
#'   dimensions for [plot_save()].
#' @family plotting
#' @export
plot_ora_bubble <- function(
    result,
    terms,
    x = c("rich_factor", "gene_ratio", "fold_enrichment"),
    p_value = c("auto", "adjusted", "raw"),
    term_labels = NULL
) {
  terms <- .plot_explicit_terms(terms)
  x <- match.arg(x)
  p_value <- match.arg(p_value)
  view <- .plot_ora_data(result)
  selected <- .plot_ora_select(view, terms, term_labels, p_value)
  data <- selected$data
  data$x_value <- data[[x]]
  data$display_term <- .plot_wrap_labels(data$term_label, width = 40L)
  data$display_term <- make.unique(data$display_term)
  data$term <- factor(data$display_term, levels = rev(data$display_term))

  x_label <- switch(
    x,
    rich_factor = "Rich factor (%)",
    gene_ratio = "Gene ratio (%)",
    fold_enrichment = "Fold enrichment"
  )
  plot <- ggplot2::ggplot(
    data,
    ggplot2::aes(
      x = .data[["x_value"]],
      y = .data[["term"]],
      size = .data[["count"]],
      colour = .data[["minus_log10_p"]]
    )
  ) +
    ggplot2::geom_point(alpha = 0.92) +
    ggplot2::scale_size_continuous(
      range = c(0.8, 2.8),
      name = "Gene count"
    ) +
    ggplot2::scale_colour_gradient(
      low = "#D9EAF7",
      high = "#B2182B",
      name = .plot_ora_evidence_legend(selected$evidence_label)
    ) +
    ggplot2::guides(
      size = ggplot2::guide_legend(order = 1),
      colour = ggplot2::guide_colourbar(order = 2)
    ) +
    ggplot2::labs(
      x = x_label,
      y = NULL,
      size = "Gene count",
      colour = .plot_ora_evidence_legend(selected$evidence_label),
      alt = paste(
        "A bubble plot of explicitly selected over-representation terms;",
        x_label,
        "is on the horizontal axis, overlap is point area, and evidence is colour."
      )
    ) +
    ggplot2::theme_bw(
      base_size = .bulkmae_text_size_pt,
      base_line_size = 0.2,
      base_rect_size = 0.2
    ) +
    ggplot2::theme(
      panel.grid.major = ggplot2::element_line(
        colour = "grey90", linewidth = 0.2
      ),
      panel.grid.minor = ggplot2::element_blank(),
      panel.border = ggplot2::element_rect(
        colour = "black", fill = NA, linewidth = 0.2
      ),
      axis.text = ggplot2::element_text(
        colour = "black", size = .bulkmae_text_size_pt
      ),
      axis.title = ggplot2::element_text(
        colour = "black", size = .bulkmae_text_size_pt
      ),
      axis.ticks = ggplot2::element_line(colour = "black", linewidth = 0.2),
      axis.ticks.length = grid::unit(0.8, "mm"),
      plot.title = ggplot2::element_text(
        colour = "black", face = "bold", hjust = 0.5,
        size = .bulkmae_text_size_pt
      ),
      legend.position = "right",
      legend.title = ggplot2::element_text(size = .bulkmae_text_size_pt),
      legend.text = ggplot2::element_text(size = .bulkmae_text_size_pt),
      legend.key.height = grid::unit(2, "mm"),
      legend.key.width = grid::unit(2, "mm"),
      plot.background = ggplot2::element_rect(
        fill = "white", colour = NA
      ),
      plot.margin = ggplot2::margin(0, 0, 0, 0, unit = "mm")
    )
  if (x %in% c("rich_factor", "gene_ratio")) {
    plot <- plot + ggplot2::scale_x_continuous(
      labels = .plot_percent_labels,
      n.breaks = 4,
      expand = ggplot2::expansion(mult = c(0.01, 0.08))
    )
  } else {
    plot <- plot + ggplot2::scale_x_continuous(
      n.breaks = 4,
      expand = ggplot2::expansion(mult = c(0.01, 0.08))
    )
  }
  .plot_with_dimensions(plot, width = 10, height = 8)
}

#' Plot an ORA term-feature community network
#'
#' Builds a deterministic community layout from explicitly selected ORA terms
#' and their enriched-feature membership. Term-to-term edges represent the
#' Jaccard coefficient of the complete enriched-feature sets, not ontology or
#' full-pathway similarity. A shared feature is drawn once inside every
#' selected term community that contains it, matching the reference community
#' grammar. Optional feature values control every corresponding visual node's
#' size and alpha by absolute magnitude. Values do not affect term selection.
#'
#' @param result A clusterProfiler `enrichResult` or compatible ORA result data
#'   frame. See [plot_ora_bubble()].
#' @param terms Unique term identifiers to display, in the desired order.
#' @param membership Optional term-named list of enriched feature identifiers.
#'   When omitted, membership is read from the result's `geneID` column. An
#'   explicit list must describe the selected terms exactly.
#' @param feature_values Optional feature-named finite numeric vector. It must
#'   cover every displayed feature; additional values are ignored.
#' @param features Optional explicit feature identifiers to display. Supply a
#'   character vector for one global feature set, or a term-named list to
#'   select displayed membership independently within each selected term.
#'   `NULL` displays all enriched features. Selection never changes overlap
#'   statistics computed from the complete membership.
#' @param term_labels Optional complete term-named labels.
#' @param feature_labels Optional feature-named labels covering every displayed
#'   feature. Additional labels are ignored.
#' @param label_features Explicit displayed feature identifiers to label.
#'   `NULL` labels no feature nodes.
#' @param p_value Evidence column used for term-node size.
#'
#' @return An unprinted standard ggplot object carrying recommended physical
#'   dimensions for [plot_save()].
#' @family plotting
#' @export
plot_ora_network <- function(
    result,
    terms,
    membership = NULL,
    feature_values = NULL,
    features = NULL,
    term_labels = NULL,
    feature_labels = NULL,
    label_features = NULL,
    p_value = c("auto", "adjusted", "raw")
) {
  graph <- .plot_ora_graph(
    result = result,
    terms = terms,
    membership = membership,
    feature_values = feature_values,
    features = features,
    term_labels = term_labels,
    feature_labels = feature_labels,
    label_features = label_features,
    p_value = p_value
  )
  layout <- .plot_ora_community_layout(graph)
  plot <- .plot_ora_community_layers(graph, layout)
  .plot_with_dimensions(plot, width = 8, height = 8)
}

#' Plot an ORA radial term-feature network
#'
#' Places terms on an inner ring and unique enriched features on an outer ring.
#' Membership curves retain every selected term-feature relation and inner
#' term-edge width and alpha show the number of shared enriched features. For
#' layout only, each outer feature is assigned to the most significant adjacent
#' term, with selected-term order breaking ties; shared features remain
#' connected to every adjacent term.
#'
#' @inheritParams plot_ora_network
#'
#' @return An unprinted standard ggplot object carrying recommended physical
#'   dimensions for [plot_save()].
#' @family plotting
#' @export
plot_ora_radial <- function(
    result,
    terms,
    membership = NULL,
    feature_values = NULL,
    features = NULL,
    term_labels = NULL,
    feature_labels = NULL,
    label_features = NULL,
    p_value = c("auto", "adjusted", "raw")
) {
  graph <- .plot_ora_graph(
    result = result,
    terms = terms,
    membership = membership,
    feature_values = feature_values,
    features = features,
    term_labels = term_labels,
    feature_labels = feature_labels,
    label_features = label_features,
    p_value = p_value
  )
  layout <- .plot_ora_radial_layout(graph)
  plot <- .plot_ora_radial_layers(graph, layout)
  .plot_with_dimensions(plot, width = 8, height = 8)
}

.plot_explicit_terms <- function(terms) {
  if (is.null(terms) || !length(terms)) {
    stop("`terms` must explicitly contain at least one term identifier.",
         call. = FALSE)
  }
  terms <- as.character(terms)
  .plot_assert_ids(terms, "`terms`")
  terms
}

.plot_ora_data <- function(result) {
  if (methods::is(result, "gseaResult")) {
    stop("`result` is a GSEA result; use a GSEA plotting function.", call. = FALSE)
  }
  if (!(methods::is(result, "enrichResult") || is.data.frame(result) ||
        is.matrix(result))) {
    stop(
      "Unsupported `result`; use a clusterProfiler ORA result or compatible data frame.",
      call. = FALSE
    )
  }
  table <- as.data.frame(result, optional = TRUE)
  if (!nrow(table)) stop("`result` contains no ORA terms.", call. = FALSE)
  required <- c(
    "ID", "Description", "GeneRatio", "BgRatio", "pvalue", "p.adjust", "Count"
  )
  .plot_require_columns(table, required, "result")

  term_id <- as.character(table$ID)
  .plot_assert_ids(term_id, "ORA term identifiers")
  term_label <- as.character(table$Description)
  if (length(term_label) != length(term_id) || anyNA(term_label) ||
      any(!nzchar(term_label))) {
    stop("ORA term descriptions must be non-missing and non-empty.", call. = FALSE)
  }
  .plot_enrichment_probabilities(table$pvalue, "Raw p-values", allow_all_missing = TRUE)
  .plot_enrichment_probabilities(
    table$p.adjust, "Adjusted p-values", allow_all_missing = TRUE
  )
  count <- table$Count
  if (!is.numeric(count) || any(!is.finite(count)) || any(count <= 0) ||
      any(count != trunc(count))) {
    stop("ORA overlap counts must be finite positive integers.", call. = FALSE)
  }

  gene_ratio <- .plot_parse_ora_ratio(table$GeneRatio, "`GeneRatio`")
  background_ratio <- .plot_parse_ora_ratio(table$BgRatio, "`BgRatio`")
  if (any(gene_ratio$numerator != count)) {
    stop("Each `GeneRatio` numerator must equal its ORA `Count`.", call. = FALSE)
  }
  if (any(background_ratio$numerator < count)) {
    stop("Each background term size must be at least its ORA `Count`.",
         call. = FALSE)
  }

  parsed_membership <- NULL
  if ("geneID" %in% names(table)) {
    parsed_membership <- .plot_parse_ora_membership(table$geneID, term_id)
  }
  list(
    data = data.frame(
      term_id = term_id,
      term_label = term_label,
      p_value = as.numeric(table$pvalue),
      adjusted_p_value = as.numeric(table$p.adjust),
      count = as.numeric(count),
      gene_ratio = gene_ratio$value,
      rich_factor = count / background_ratio$numerator,
      fold_enrichment = gene_ratio$value / background_ratio$value,
      stringsAsFactors = FALSE,
      check.names = FALSE
    ),
    membership = parsed_membership
  )
}

.plot_parse_ora_ratio <- function(x, label) {
  value <- trimws(as.character(x))
  match <- regexec("^([0-9]+)[[:space:]]*/[[:space:]]*([1-9][0-9]*)$", value)
  pieces <- regmatches(value, match)
  valid <- lengths(pieces) == 3L
  if (anyNA(value) || any(!valid)) {
    stop(label, " must contain ratios written as non-negative-integer/positive-integer.",
         call. = FALSE)
  }
  numerator <- as.numeric(vapply(pieces, `[[`, character(1), 2L))
  denominator <- as.numeric(vapply(pieces, `[[`, character(1), 3L))
  if (any(!is.finite(numerator)) || any(!is.finite(denominator)) ||
      any(numerator > denominator)) {
    stop(label, " must contain finite ratios between zero and one.", call. = FALSE)
  }
  list(
    numerator = numerator,
    denominator = denominator,
    value = numerator / denominator
  )
}

.plot_parse_ora_membership <- function(x, terms) {
  if (is.list(x) && !is.data.frame(x)) {
    values <- lapply(x, as.character)
  } else {
    values <- strsplit(as.character(x), "/", fixed = TRUE)
  }
  if (length(values) != length(terms)) {
    stop("The ORA `geneID` column is inconsistent with its terms.", call. = FALSE)
  }
  values <- lapply(values, function(ids) unique(trimws(ids)))
  invalid <- vapply(
    values,
    function(ids) !length(ids) || anyNA(ids) || any(!nzchar(ids)),
    logical(1)
  )
  if (any(invalid)) {
    stop("ORA enriched-feature membership cannot contain missing or empty IDs.",
         call. = FALSE)
  }
  names(values) <- terms
  values
}

.plot_ora_select <- function(view, terms, term_labels, p_value) {
  unknown <- setdiff(terms, view$data$term_id)
  if (length(unknown)) {
    stop("Unknown `terms`: ", paste(unknown, collapse = ", "), ".", call. = FALSE)
  }
  data <- view$data[match(terms, view$data$term_id), , drop = FALSE]
  labels <- .plot_enrichment_labels(data, term_labels)
  data$term_label <- labels
  evidence <- .plot_ora_evidence(data, p_value)
  data$minus_log10_p <- evidence$value
  list(data = data, evidence_label = evidence$label)
}

.plot_ora_evidence <- function(data, p_value) {
  p_value <- match.arg(p_value, c("auto", "adjusted", "raw"))
  adjusted_available <- all(is.finite(data$adjusted_p_value))
  evidence <- if (p_value == "auto") {
    if (adjusted_available) "adjusted" else "raw"
  } else {
    p_value
  }
  if (evidence == "adjusted" && !adjusted_available) {
    stop("Adjusted p-values are unavailable for the selected terms.", call. = FALSE)
  }
  values <- if (evidence == "adjusted") {
    data$adjusted_p_value
  } else {
    data$p_value
  }
  if (!any(is.finite(values))) {
    stop("No usable ", evidence, " p-values are available for the selected terms.",
         call. = FALSE)
  }
  unavailable <- !is.finite(values)
  if (any(unavailable)) {
    stop(
      "Selected terms have unavailable ", evidence, " p-values: ",
      paste(data$term_id[unavailable], collapse = ", "), ".",
      call. = FALSE
    )
  }
  positive <- values[values > 0]
  floor <- if (length(positive)) {
    max(min(positive) / 10, .Machine$double.xmin)
  } else {
    .Machine$double.xmin
  }
  list(
    value = -log10(pmax(values, floor)),
    label = if (evidence == "adjusted") {
      "-log10 adjusted p-value"
    } else {
      "-log10 p-value"
    }
  )
}

.plot_percent_labels <- function(x) {
  labels <- rep(NA_character_, length(x))
  finite <- is.finite(x)
  labels[finite] <- sprintf("%.1f%%", 100 * x[finite])
  labels
}

.plot_ora_evidence_legend <- function(label) {
  if (identical(label, "-log10 adjusted p-value")) {
    return(expression(-log[10]("adj p")))
  }
  expression(-log[10]("p-value"))
}

.plot_ora_graph <- function(
    result,
    terms,
    membership,
    feature_values,
    features,
    term_labels,
    feature_labels,
    label_features,
  p_value
) {
  terms <- .plot_explicit_terms(terms)
  p_value <- match.arg(p_value, c("auto", "adjusted", "raw"))
  view <- .plot_ora_data(result)
  selected <- .plot_ora_select(view, terms, term_labels, p_value)
  complete_membership <- .plot_ora_membership(
    view = view,
    selected = selected$data,
    terms = terms,
    membership = membership
  )
  complete_features <- unique(unlist(complete_membership, use.names = FALSE))
  if (is.null(features)) {
    features <- complete_features
    display_membership <- complete_membership
  } else if (is.list(features) && !is.data.frame(features)) {
    if (is.null(names(features))) {
      stop("A list supplied as `features` must be term-named.", call. = FALSE)
    }
    .plot_assert_ids(names(features), "Names of `features`")
    if (!setequal(names(features), terms)) {
      stop("A list supplied as `features` must describe selected terms exactly.",
           call. = FALSE)
    }
    display_membership <- lapply(features[terms], function(ids) {
      ids <- as.character(ids)
      .plot_assert_ids(ids, "Each element of `features`")
      ids
    })
    for (term in terms) {
      unknown <- setdiff(display_membership[[term]], complete_membership[[term]])
      if (length(unknown)) {
        stop(
          "Unknown displayed features for term ", term, ": ",
          paste(unknown, collapse = ", "), ".",
          call. = FALSE
        )
      }
    }
    features <- unique(unlist(display_membership, use.names = FALSE))
  } else {
    features <- as.character(features)
    .plot_assert_ids(features, "`features`")
    unknown <- setdiff(features, complete_features)
    if (length(unknown)) {
      stop(
        "Unknown `features`: ", paste(unknown, collapse = ", "), ".",
        call. = FALSE
      )
    }
    display_membership <- lapply(
      complete_membership,
      function(ids) features[features %in% ids]
    )
  }
  names(display_membership) <- terms
  empty <- names(display_membership)[!lengths(display_membership)]
  if (length(empty)) {
    stop(
      "`features` leaves selected terms without displayed nodes: ",
      paste(empty, collapse = ", "), ".",
      call. = FALSE
    )
  }

  values <- NULL
  if (!is.null(feature_values)) {
    .assert_named_numeric(feature_values, "feature_values")
    missing_values <- setdiff(features, names(feature_values))
    if (length(missing_values)) {
      stop(
        "`feature_values` is missing displayed features: ",
        paste(missing_values, collapse = ", "), ".",
        call. = FALSE
      )
    }
    values <- as.numeric(feature_values[features])
    names(values) <- features
  }

  labels <- features
  names(labels) <- features
  if (!is.null(feature_labels)) {
    if (is.null(names(feature_labels))) {
      stop("`feature_labels` must be feature-named.", call. = FALSE)
    }
    .plot_assert_ids(names(feature_labels), "Names of `feature_labels`")
    missing_labels <- setdiff(features, names(feature_labels))
    if (length(missing_labels)) {
      stop(
        "`feature_labels` is missing displayed features: ",
        paste(missing_labels, collapse = ", "), ".",
        call. = FALSE
      )
    }
    labels <- as.character(feature_labels[features])
    names(labels) <- features
    if (anyNA(labels) || any(!nzchar(labels))) {
      stop("`feature_labels` cannot contain missing or empty labels.",
           call. = FALSE)
    }
  }

  if (is.null(label_features)) {
    label_features <- character()
  } else {
    label_features <- as.character(label_features)
    .plot_assert_ids(label_features, "`label_features`")
    unknown_labels <- setdiff(label_features, features)
    if (length(unknown_labels)) {
      stop(
        "Unknown `label_features`: ", paste(unknown_labels, collapse = ", "), ".",
        call. = FALSE
      )
    }
  }

  membership_edges <- do.call(
    rbind,
    lapply(terms, function(term) {
      data.frame(
        term_id = term,
        feature_id = display_membership[[term]],
        stringsAsFactors = FALSE
      )
    })
  )
  rownames(membership_edges) <- NULL
  overlap <- .plot_ora_overlap(complete_membership, terms)
  colours <- .plot_ora_palette(length(terms))
  names(colours) <- terms

  list(
    terms = terms,
    term_data = selected$data,
    features = features,
    membership = display_membership,
    complete_membership = complete_membership,
    membership_edges = membership_edges,
    overlap = overlap,
    feature_values = values,
    feature_labels = labels,
    label_features = label_features,
    colours = colours,
    evidence_label = selected$evidence_label
  )
}

.plot_ora_membership <- function(view, selected, terms, membership) {
  if (is.null(membership)) {
    if (is.null(view$membership)) {
      stop(
        "`membership` is required because `result` has no `geneID` column.",
        call. = FALSE
      )
    }
    values <- view$membership[terms]
  } else {
    if (!is.list(membership) || is.data.frame(membership) ||
        is.null(names(membership))) {
      stop("`membership` must be a term-named list.", call. = FALSE)
    }
    .plot_assert_ids(names(membership), "Names of `membership`")
    if (!setequal(names(membership), terms)) {
      stop("`membership` must describe the selected terms exactly.", call. = FALSE)
    }
    values <- membership[terms]
  }
  values <- lapply(values, function(ids) unique(as.character(ids)))
  invalid <- vapply(
    values,
    function(ids) !length(ids) || anyNA(ids) || any(!nzchar(ids)),
    logical(1)
  )
  if (any(invalid)) {
    stop("Every selected term must have non-missing enriched-feature membership.",
         call. = FALSE)
  }
  observed_count <- lengths(values)
  expected_count <- selected$count[match(terms, selected$term_id)]
  if (any(observed_count != expected_count)) {
    bad <- terms[observed_count != expected_count]
    stop(
      "Enriched-feature membership does not match ORA `Count` for terms: ",
      paste(bad, collapse = ", "), ".",
      call. = FALSE
    )
  }
  names(values) <- terms
  values
}

.plot_ora_overlap <- function(membership, terms) {
  if (length(terms) < 2L) {
    return(data.frame(
      from = character(), to = character(), shared_n = integer(),
      jaccard = numeric(), stringsAsFactors = FALSE
    ))
  }
  pairs <- utils::combn(terms, 2L, simplify = FALSE)
  edges <- lapply(pairs, function(pair) {
    shared <- intersect(membership[[pair[[1L]]]], membership[[pair[[2L]]]])
    union <- union(membership[[pair[[1L]]]], membership[[pair[[2L]]]])
    data.frame(
      from = pair[[1L]],
      to = pair[[2L]],
      shared_n = length(shared),
      jaccard = length(shared) / length(union),
      stringsAsFactors = FALSE
    )
  })
  result <- do.call(rbind, edges)
  rownames(result) <- NULL
  result[result$shared_n > 0L, , drop = FALSE]
}

.plot_ora_palette_base <- function() {
  c("#6C63B8", "#2BB3A3", "#F2A65A", "#D1495B", "#8DAA22")
}

.plot_ora_palette <- function(n) {
  grDevices::colorRampPalette(.plot_ora_palette_base())(n)
}

.plot_ora_radial_palette <- function(n) {
  base <- .plot_ora_palette_base()
  if (n <= length(base)) return(base[seq_len(n)])
  grDevices::colorRampPalette(base)(n)
}

.plot_wrap_labels <- function(x, width = 24L) {
  vapply(
    x,
    function(label) paste(strwrap(label, width = width), collapse = "\n"),
    character(1),
    USE.NAMES = FALSE
  )
}

.plot_wrap_words <- function(x, words = 3L) {
  vapply(
    as.character(x),
    function(label) {
      pieces <- strsplit(trimws(label), "[[:space:]]+")[[1L]]
      groups <- split(pieces, ceiling(seq_along(pieces) / words))
      paste(vapply(groups, paste, character(1), collapse = " "), collapse = "\n")
    },
    character(1),
    USE.NAMES = FALSE
  )
}

.plot_rescale <- function(x, to = c(0, 1)) {
  if (!length(x)) return(numeric())
  range <- range(x, finite = TRUE)
  if (!all(is.finite(range))) stop("Cannot rescale non-finite values.", call. = FALSE)
  if (diff(range) == 0) return(rep(mean(to), length(x)))
  to[[1L]] + (x - range[[1L]]) / diff(range) * diff(to)
}

.plot_pt_to_mm <- function(x) {
  as.numeric(x) / (72.27 / 25.4)
}

.plot_ora_shape_map <- function(terms) {
  shapes <- rep(c(21, 22, 23, 24, 25), length.out = length(terms))
  stats::setNames(shapes, terms)
}

.plot_ora_primary_terms <- function(graph) {
  term_order <- stats::setNames(seq_along(graph$terms), graph$terms)
  evidence <- stats::setNames(
    graph$term_data$minus_log10_p,
    graph$term_data$term_id
  )
  result <- vapply(graph$features, function(feature) {
    adjacent <- graph$membership_edges$term_id[
      graph$membership_edges$feature_id == feature
    ]
    best <- max(evidence[adjacent])
    candidates <- adjacent[evidence[adjacent] == best]
    candidates[[which.min(term_order[candidates])]]
  }, character(1))
  stats::setNames(result, graph$features)
}

.plot_ora_network_style <- function() {
  list(
    layout_iterations = 800L,
    layout_step = 0.04,
    layout_min_dist = 7,
    layout_max_dist = 10,
    layout_attraction = 0.008,
    layout_repulsion = 0.35,
    layout_hard_repulsion = 0.2,
    layout_overlap_power = 0.3,
    layout_max_radius = 3.3,
    module_edge_width = c(0.8, 1.3),
    module_edge_alpha = c(0.35, 0.55),
    module_edge_colour = "#c7d0e4",
    module_edge_gap = 0.12,
    hub_gene_min_dist = 0.45,
    inner_compression = 0.5,
    shared_gene_push = 1.5,
    shared_align_jitter = 1.2,
    canvas_limit = 3.3,
    bg_n_rings = 15L,
    bg_alpha = c(0.001, 0.018),
    bg_radius_scale = 1.5,
    feature_size = c(2, 3.5),
    hub_size = c(3.5, 5),
    feature_alpha = c(0.35, 0.58),
    feature_min_dist = 0.20
  )
}

.plot_ora_network_centres <- function(graph, style) {
  n <- length(graph$terms)
  if (n == 1L) {
    return(data.frame(
      term_id = graph$terms, center_x = 0, center_y = 0,
      stringsAsFactors = FALSE
    ))
  }
  golden <- pi * (3 - sqrt(5))
  index <- seq_len(n)
  radius <- sqrt(index / n) * sqrt(n)
  angle <- index * golden
  x <- radius * cos(angle) + 0.08 * sin(index * sqrt(2))
  y <- radius * sin(angle) + 0.08 * cos(index * sqrt(3))
  edges <- graph$overlap
  if (nrow(edges)) {
    max_overlap <- max(edges$jaccard)
    edges$i <- match(edges$from, graph$terms)
    edges$j <- match(edges$to, graph$terms)
    edges$weight <- (edges$jaccard / max_overlap)^style$layout_overlap_power
    edges$target <- style$layout_max_dist -
      (style$layout_max_dist - style$layout_min_dist) * edges$weight
  }
  for (iteration in seq_len(style$layout_iterations)) {
    dx_total <- numeric(n)
    dy_total <- numeric(n)
    for (i in seq_len(n - 1L)) {
      for (j in seq.int(i + 1L, n)) {
        dx <- x[[j]] - x[[i]]
        dy <- y[[j]] - y[[i]]
        distance <- sqrt(dx^2 + dy^2) + 1e-6
        ux <- dx / distance
        uy <- dy / distance
        force <- min(style$layout_repulsion / distance^2, 1)
        dx_total[[i]] <- dx_total[[i]] - ux * force
        dy_total[[i]] <- dy_total[[i]] - uy * force
        dx_total[[j]] <- dx_total[[j]] + ux * force
        dy_total[[j]] <- dy_total[[j]] + uy * force
        if (distance < style$layout_min_dist) {
          hard <- (style$layout_min_dist - distance) *
            style$layout_hard_repulsion
          dx_total[[i]] <- dx_total[[i]] - ux * hard
          dy_total[[i]] <- dy_total[[i]] - uy * hard
          dx_total[[j]] <- dx_total[[j]] + ux * hard
          dy_total[[j]] <- dy_total[[j]] + uy * hard
        }
      }
    }
    if (nrow(edges)) {
      for (edge in seq_len(nrow(edges))) {
        i <- edges$i[[edge]]
        j <- edges$j[[edge]]
        dx <- x[[j]] - x[[i]]
        dy <- y[[j]] - y[[i]]
        distance <- sqrt(dx^2 + dy^2) + 1e-6
        force <- style$layout_attraction * edges$weight[[edge]] *
          (distance - edges$target[[edge]])
        ux <- dx / distance
        uy <- dy / distance
        dx_total[[i]] <- dx_total[[i]] + ux * force
        dy_total[[i]] <- dy_total[[i]] + uy * force
        dx_total[[j]] <- dx_total[[j]] - ux * force
        dy_total[[j]] <- dy_total[[j]] - uy * force
      }
    }
    step <- style$layout_step *
      (1 - 0.65 * iteration / style$layout_iterations)
    x <- x + step * dx_total
    y <- y + step * dy_total
    x <- x - mean(x)
    y <- y - mean(y)
  }
  maximum <- max(sqrt(x^2 + y^2))
  if (maximum > style$layout_max_radius) {
    x <- x * style$layout_max_radius / maximum
    y <- y * style$layout_max_radius / maximum
  }
  data.frame(
    term_id = graph$terms, center_x = x, center_y = y,
    stringsAsFactors = FALSE
  )
}

.plot_ora_repel_features <- function(features, minimum, iterations = 90L) {
  if (nrow(features) < 2L || minimum <= 0) return(features)
  groups <- if ("term_id" %in% names(features)) {
    split(
      seq_len(nrow(features)),
      factor(features$term_id, levels = unique(features$term_id))
    )
  } else {
    list(seq_len(nrow(features)))
  }
  result <- lapply(groups, function(rows) {
    data <- features[rows, , drop = FALSE]
    if (nrow(data) < 2L) return(data)
    x <- data$x
    y <- data$y
    anchor_x <- x
    anchor_y <- y
    golden <- pi * (3 - sqrt(5))
    for (iteration in seq_len(iterations)) {
      dx_total <- numeric(length(x))
      dy_total <- numeric(length(y))
      for (i in seq_len(length(x) - 1L)) {
        for (j in seq.int(i + 1L, length(x))) {
          dx <- x[[j]] - x[[i]]
          dy <- y[[j]] - y[[i]]
          distance <- sqrt(dx^2 + dy^2)
          if (distance < 1e-8) {
            ux <- cos((i + j) * golden)
            uy <- sin((i + j) * golden)
            distance <- 0
          } else {
            ux <- dx / distance
            uy <- dy / distance
          }
          if (distance < minimum) {
            push <- (minimum - distance) * 0.5
            dx_total[[i]] <- dx_total[[i]] - ux * push
            dy_total[[i]] <- dy_total[[i]] - uy * push
            dx_total[[j]] <- dx_total[[j]] + ux * push
            dy_total[[j]] <- dy_total[[j]] + uy * push
          }
        }
      }
      x <- x + dx_total + (anchor_x - x) * 0.02
      y <- y + dy_total + (anchor_y - y) * 0.02
    }
    data$x <- x
    data$y <- y
    data
  })
  result <- do.call(rbind, result)
  rownames(result) <- NULL
  result
}

.plot_ora_community_layout <- function(graph) {
  style <- .plot_ora_network_style()
  centres <- .plot_ora_network_centres(graph, style)
  terms <- graph$term_data
  terms$center_x <- centres$center_x[match(terms$term_id, centres$term_id)]
  terms$center_y <- centres$center_y[match(terms$term_id, centres$term_id)]
  terms$x <- terms$center_x
  terms$y <- terms$center_y
  terms$colour <- unname(graph$colours[terms$term_id])
  terms$shape_id <- terms$term_id
  terms$module_radius <- .plot_rescale(sqrt(terms$count), c(0.7, 1.1))
  terms$hub_size <- .plot_rescale(terms$minus_log10_p, style$hub_size)
  terms$label <- .plot_wrap_words(terms$term_label, words = 3L)

  degree <- table(graph$membership_edges$feature_id)
  magnitudes <- if (is.null(graph$feature_values)) {
    rep(1, length(graph$features))
  } else {
    abs(unname(graph$feature_values[graph$features]))
  }
  feature_size <- .plot_rescale(magnitudes, style$feature_size)
  feature_alpha <- .plot_rescale(magnitudes, style$feature_alpha)
  names(feature_size) <- graph$features
  names(feature_alpha) <- graph$features
  golden <- pi * (3 - sqrt(5))
  rows <- vector("list", length(graph$terms))
  for (term_index in seq_along(graph$terms)) {
    term <- graph$terms[[term_index]]
    term_row <- match(term, terms$term_id)
    term_features <- graph$membership_edges$feature_id[
      graph$membership_edges$term_id == term
    ]
    shared <- unname(degree[term_features]) > 1L
    term_order <- order(shared, seq_along(term_features))
    term_features <- term_features[term_order]
    shared <- shared[term_order]
    term_rows <- vector("list", length(term_features))
    for (feature_index in seq_along(term_features)) {
      feature <- term_features[[feature_index]]
      adjacent <- graph$membership_edges$term_id[
        graph$membership_edges$feature_id == feature
      ]
      other_rows <- match(setdiff(adjacent, term), terms$term_id)
      if (length(other_rows)) {
        target_angle <- atan2(
          mean(terms$center_y[other_rows]) - terms$center_y[[term_row]],
          mean(terms$center_x[other_rows]) - terms$center_x[[term_row]]
        )
        jitter <- style$shared_align_jitter *
          sin((term_index + feature_index) * sqrt(2))
        local_angle <- target_angle + jitter
      } else {
        local_angle <- term_index * golden + feature_index * golden
      }
      raw_radius <- style$hub_gene_min_dist +
        (terms$module_radius[[term_row]] - style$hub_gene_min_dist) *
        sqrt(feature_index / length(term_features))
      push <- if (shared[[feature_index]]) style$shared_gene_push else 1
      local_radius <- raw_radius * style$inner_compression * push
      term_rows[[feature_index]] <- data.frame(
        term_id = term,
        feature_id = feature,
        primary_term = term,
        x = terms$center_x[[term_row]] + cos(local_angle) * local_radius,
        y = terms$center_y[[term_row]] + sin(local_angle) * local_radius,
        degree = unname(degree[[feature]]),
        label = unname(graph$feature_labels[[feature]]),
        value = if (is.null(graph$feature_values)) NA_real_ else
          unname(graph$feature_values[[feature]]),
        plot_size = unname(feature_size[[feature]]),
        feature_alpha = unname(feature_alpha[[feature]]),
        colour = unname(graph$colours[[term]]),
        shape_id = term,
        stringsAsFactors = FALSE
      )
    }
    rows[[term_index]] <- do.call(rbind, term_rows)
  }
  features <- do.call(rbind, rows)
  rownames(features) <- NULL
  features <- .plot_ora_repel_features(
    features, style$feature_min_dist, iterations = 90L
  )

  membership <- graph$membership_edges
  membership$x_term <- terms$center_x[match(membership$term_id, terms$term_id)]
  membership$y_term <- terms$center_y[match(membership$term_id, terms$term_id)]
  membership_key <- paste(membership$term_id, membership$feature_id, sep = "\r")
  feature_key <- paste(features$term_id, features$feature_id, sep = "\r")
  feature_rows <- match(membership_key, feature_key)
  membership$x_feature <- features$x[feature_rows]
  membership$y_feature <- features$y[feature_rows]

  overlap <- graph$overlap
  if (nrow(overlap)) {
    overlap$x <- terms$center_x[match(overlap$from, terms$term_id)]
    overlap$y <- terms$center_y[match(overlap$from, terms$term_id)]
    overlap$xend <- terms$center_x[match(overlap$to, terms$term_id)]
    overlap$yend <- terms$center_y[match(overlap$to, terms$term_id)]
    strength <- overlap$jaccard
    overlap$edge_width <- .plot_rescale(strength, style$module_edge_width)
    overlap$edge_alpha <- .plot_rescale(strength, style$module_edge_alpha)
    dx <- overlap$xend - overlap$x
    dy <- overlap$yend - overlap$y
    distance <- sqrt(dx^2 + dy^2) + 1e-6
    gap <- pmin(style$module_edge_gap, distance * 0.45)
    overlap$x_plot <- overlap$x + dx / distance * gap
    overlap$y_plot <- overlap$y + dy / distance * gap
    overlap$xend_plot <- overlap$xend - dx / distance * gap
    overlap$yend_plot <- overlap$yend - dy / distance * gap
  }

  feature_bottom <- vapply(graph$terms, function(term) {
    ids <- features$term_id == term
    if (any(ids)) min(features$y[ids]) else
      terms$center_y[match(term, terms$term_id)]
  }, numeric(1))
  terms$label_x <- terms$center_x
  terms$label_y <- pmin(
    feature_bottom,
    terms$center_y - terms$module_radius *
      (style$bg_radius_scale - 0.1) * 1.1 * style$inner_compression
  ) - 0.18

  rings <- seq_len(style$bg_n_rings)
  halos <- do.call(rbind, lapply(seq_len(nrow(terms)), function(i) {
    fraction <- seq(1.1, 0.3, length.out = style$bg_n_rings) *
      style$inner_compression
    data.frame(
      term_id = terms$term_id[[i]],
      x0 = terms$center_x[[i]],
      y0 = terms$center_y[[i]],
      a = terms$module_radius[[i]] * style$bg_radius_scale * fraction,
      b = terms$module_radius[[i]] * style$bg_radius_scale * fraction,
      angle = 0,
      colour = terms$colour[[i]],
      halo_alpha = seq(
        style$bg_alpha[[1L]], style$bg_alpha[[2L]],
        length.out = style$bg_n_rings
      ),
      ring = rings,
      stringsAsFactors = FALSE
    )
  }))
  rownames(halos) <- NULL
  list(
    terms = terms,
    features = features,
    membership = membership,
    overlap = overlap,
    halos = halos,
    shape_map = .plot_ora_shape_map(graph$terms),
    style = style
  )
}

.plot_ora_community_layers <- function(graph, layout) {
  plot <- ggplot2::ggplot()
  if (nrow(layout$overlap)) {
    plot <- plot + ggplot2::geom_curve(
      data = layout$overlap,
      ggplot2::aes(
        x = .data[["x_plot"]], y = .data[["y_plot"]],
        xend = .data[["xend_plot"]], yend = .data[["yend_plot"]],
        linewidth = .data[["edge_width"]], alpha = .data[["edge_alpha"]]
      ),
      colour = layout$style$module_edge_colour,
      curvature = 0,
      lineend = "square",
      inherit.aes = FALSE
    )
  }
  plot <- plot +
    ggforce::geom_ellipse(
      data = layout$halos,
      ggplot2::aes(
        x0 = .data[["x0"]], y0 = .data[["y0"]],
        a = .data[["a"]], b = .data[["b"]], angle = .data[["angle"]],
        fill = .data[["colour"]], alpha = .data[["halo_alpha"]],
        group = interaction(.data[["term_id"]], .data[["ring"]])
      ),
      colour = NA,
      inherit.aes = FALSE
    ) +
    ggplot2::geom_segment(
      data = layout$membership,
      ggplot2::aes(
        x = .data[["x_term"]], y = .data[["y_term"]],
        xend = .data[["x_feature"]], yend = .data[["y_feature"]]
      ),
      colour = "#beb0b2",
      linewidth = 0.45,
      alpha = 0.35,
      inherit.aes = FALSE
    ) +
    ggplot2::geom_point(
      data = layout$features,
      ggplot2::aes(
        x = .data[["x"]], y = .data[["y"]],
        size = .data[["plot_size"]], alpha = .data[["feature_alpha"]],
        fill = .data[["colour"]], shape = .data[["shape_id"]]
      ),
      colour = "white",
      stroke = 0,
      inherit.aes = FALSE
    ) +
    ggplot2::geom_point(
      data = layout$terms,
      ggplot2::aes(
        x = .data[["center_x"]], y = .data[["center_y"]],
        size = .data[["hub_size"]], fill = .data[["colour"]],
        shape = .data[["shape_id"]]
      ),
      colour = "white",
      stroke = 0.4,
      inherit.aes = FALSE
    ) +
    ggrepel::geom_text_repel(
      data = layout$terms,
      ggplot2::aes(
        x = .data[["label_x"]], y = .data[["label_y"]],
        label = .data[["label"]]
      ),
      vjust = 1,
      size = .plot_pt_to_mm(.bulkmae_text_size_pt),
      lineheight = 0.9,
      colour = "black",
      fontface = "plain",
      bg.color = NA,
      bg.r = 0.1,
      box.padding = 0.1,
      point.padding = 0,
      direction = "y",
      segment.size = 0.2,
      segment.alpha = 0.3,
      max.overlaps = Inf,
      seed = 2026,
      inherit.aes = FALSE
    )
  if (length(graph$label_features)) {
    labels <- layout$features[
      layout$features$feature_id %in% graph$label_features, , drop = FALSE
    ]
    plot <- plot + ggrepel::geom_text_repel(
      data = labels,
      ggplot2::aes(
        x = .data[["x"]], y = .data[["y"]], label = .data[["label"]]
      ),
      size = .plot_pt_to_mm(.bulkmae_text_size_pt),
      colour = "#474545",
      box.padding = 0.08,
      point.padding = 0.05,
      segment.size = 0.15,
      segment.alpha = 0.3,
      max.overlaps = Inf,
      seed = 2027,
      inherit.aes = FALSE
    )
  }
  plot +
    ggplot2::scale_fill_identity() +
    ggplot2::scale_alpha_identity() +
    ggplot2::scale_size_identity() +
    ggplot2::scale_linewidth_identity() +
    ggplot2::scale_shape_manual(values = layout$shape_map, guide = "none") +
    ggplot2::coord_fixed(
      xlim = c(-layout$style$canvas_limit, layout$style$canvas_limit),
      ylim = c(-layout$style$canvas_limit, layout$style$canvas_limit),
      clip = "off"
    ) +
    ggplot2::labs(
      x = NULL,
      y = NULL,
      alt = paste(
        "A community network of selected ORA terms and term-specific enriched",
        "feature nodes. Term-feature lines show membership, pale term-term lines",
        "show Jaccard overlap, feature size and alpha show absolute feature values",
        "when supplied, and term size shows enrichment evidence."
      )
    ) +
    ggplot2::theme_void(base_size = .bulkmae_text_size_pt) +
    ggplot2::theme(
      legend.position = "none",
      plot.background = ggplot2::element_rect(fill = "white", colour = NA),
      plot.margin = ggplot2::margin(5, 0, 0, 0, unit = "pt"),
      plot.title = ggplot2::element_text(
        face = "plain", hjust = 0.5, size = .bulkmae_text_size_pt,
        margin = ggplot2::margin(b = 2)
      )
    )
}

.plot_ora_radial_style <- function() {
  list(
    outer_radius = 6.3,
    inner_radius = 4.1,
    label_radius = 6.65,
    group_gap = 0.075,
    angle_offset = pi / 2,
    canvas_limit = 7.2,
    term_size = c(6, 12),
    term_alpha = c(0.35, 0.55),
    feature_size = c(1.5, 3),
    feature_alpha = c(0.45, 0.77),
    repel_min_dist = 2.4,
    repel_iterations = 500L,
    repel_step = 0.1,
    repel_anchor = 0.06,
    repel_radial = 0.06,
    overlap_width = c(0.8, 1.4),
    overlap_alpha = c(0.25, 0.40),
    overlap_power = 0.85
  )
}

.plot_ora_repel_terms <- function(terms, style) {
  if (nrow(terms) < 2L) return(terms)
  x <- terms$x
  y <- terms$y
  golden <- pi * (3 - sqrt(5))
  for (iteration in seq_len(style$repel_iterations)) {
    dx_total <- numeric(length(x))
    dy_total <- numeric(length(y))
    for (i in seq_len(length(x) - 1L)) {
      for (j in seq.int(i + 1L, length(x))) {
        dx <- x[[j]] - x[[i]]
        dy <- y[[j]] - y[[i]]
        distance <- sqrt(dx^2 + dy^2)
        if (distance < 1e-8) {
          ux <- cos((i + j) * golden)
          uy <- sin((i + j) * golden)
          distance <- 0
        } else {
          ux <- dx / distance
          uy <- dy / distance
        }
        if (distance < style$repel_min_dist) {
          push <- (style$repel_min_dist - distance) * 0.5
          dx_total[[i]] <- dx_total[[i]] - ux * push
          dy_total[[i]] <- dy_total[[i]] - uy * push
          dx_total[[j]] <- dx_total[[j]] + ux * push
          dy_total[[j]] <- dy_total[[j]] + uy * push
        }
      }
    }
    dx_total <- dx_total + (terms$x_anchor - x) * style$repel_anchor
    dy_total <- dy_total + (terms$y_anchor - y) * style$repel_anchor
    x <- x + style$repel_step * dx_total
    y <- y + style$repel_step * dy_total
    radius <- sqrt(x^2 + y^2) + 1e-8
    x <- x + style$repel_radial * (style$inner_radius - radius) * x / radius
    y <- y + style$repel_radial * (style$inner_radius - radius) * y / radius
  }
  terms$x <- x
  terms$y <- y
  terms
}

.plot_ora_radial_layout <- function(graph) {
  style <- .plot_ora_radial_style()
  colours <- .plot_ora_radial_palette(length(graph$terms))
  names(colours) <- graph$terms
  primary <- .plot_ora_primary_terms(graph)
  group_count <- vapply(
    graph$terms,
    function(term) sum(primary == term),
    integer(1)
  )
  span_count <- pmax(group_count, 1L)
  available <- 2 * pi - style$group_gap * length(graph$terms)
  span <- available * span_count / sum(span_count)
  start <- style$angle_offset + cumsum(c(
    0,
    utils::head(span + style$group_gap, -1L)
  ))
  group_arcs <- data.frame(
    term_id = graph$terms,
    feature_n = group_count,
    arc_start = start,
    arc_end = start + span,
    arc_mid = start + span / 2,
    stringsAsFactors = FALSE
  )

  terms <- graph$term_data
  terms$angle_anchor <- group_arcs$arc_mid[
    match(terms$term_id, group_arcs$term_id)
  ]
  terms$x_anchor <- style$inner_radius * cos(terms$angle_anchor)
  terms$y_anchor <- style$inner_radius * sin(terms$angle_anchor)
  terms$x <- terms$x_anchor
  terms$y <- terms$y_anchor
  terms <- .plot_ora_repel_terms(terms, style)
  terms$angle <- atan2(terms$y, terms$x)
  terms$colour <- unname(colours[terms$term_id])
  terms$pathway_size <- .plot_rescale(terms$minus_log10_p, style$term_size)
  terms$pathway_alpha <- .plot_rescale(terms$minus_log10_p, style$term_alpha)
  terms$label <- vapply(
    terms$term_label,
    stringr::str_wrap,
    character(1),
    width = 20L,
    USE.NAMES = FALSE
  )

  features <- data.frame(
    feature_id = graph$features,
    primary_term = unname(primary[graph$features]),
    label = unname(graph$feature_labels[graph$features]),
    value = if (is.null(graph$feature_values)) NA_real_ else
      unname(graph$feature_values[graph$features]),
    stringsAsFactors = FALSE
  )
  features$abs_value <- if (is.null(graph$feature_values)) 1 else abs(features$value)
  features$angle <- NA_real_
  for (term in graph$terms) {
    rows <- which(features$primary_term == term)
    if (!length(rows)) next
    rows <- rows[order(-features$abs_value[rows], features$label[rows])]
    arc <- group_arcs[group_arcs$term_id == term, , drop = FALSE]
    features$angle[rows] <- arc$arc_start +
      (arc$arc_end - arc$arc_start) *
      (seq_along(rows) - 0.5) / length(rows)
  }
  features$x <- style$outer_radius * cos(features$angle)
  features$y <- style$outer_radius * sin(features$angle)
  features$label_x <- style$label_radius * cos(features$angle)
  features$label_y <- style$label_radius * sin(features$angle)
  angle_degrees <- (features$angle * 180 / pi) %% 360
  flip <- angle_degrees > 90 & angle_degrees < 270
  features$label_angle <- ifelse(flip, angle_degrees + 180, angle_degrees)
  features$label_hjust <- ifelse(flip, 1, 0)
  features$feature_size <- .plot_rescale(features$abs_value, style$feature_size)
  features$feature_alpha <- .plot_rescale(features$abs_value, style$feature_alpha)
  features$colour <- unname(colours[features$primary_term])

  membership <- graph$membership_edges
  membership$x_term <- terms$x[match(membership$term_id, terms$term_id)]
  membership$y_term <- terms$y[match(membership$term_id, terms$term_id)]
  membership$x_feature <- features$x[match(membership$feature_id, features$feature_id)]
  membership$y_feature <- features$y[match(membership$feature_id, features$feature_id)]
  membership$term_colour <- unname(colours[membership$term_id])

  overlap <- graph$overlap
  if (nrow(overlap)) {
    overlap$x <- terms$x[match(overlap$from, terms$term_id)]
    overlap$y <- terms$y[match(overlap$from, terms$term_id)]
    overlap$xend <- terms$x[match(overlap$to, terms$term_id)]
    overlap$yend <- terms$y[match(overlap$to, terms$term_id)]
    strength <- overlap$shared_n^style$overlap_power
    overlap$edge_width <- .plot_rescale(strength, style$overlap_width)
    overlap$edge_alpha <- .plot_rescale(strength, style$overlap_alpha)
  }
  list(
    terms = terms,
    features = features,
    membership = membership,
    overlap = overlap,
    group_arcs = group_arcs,
    style = style
  )
}

.plot_ora_radial_layers <- function(graph, layout) {
  plot <- ggplot2::ggplot()
  if (nrow(layout$overlap)) {
    plot <- plot + ggplot2::geom_curve(
      data = layout$overlap,
      ggplot2::aes(
        x = .data[["x"]], y = .data[["y"]],
        xend = .data[["xend"]], yend = .data[["yend"]],
        linewidth = .data[["edge_width"]], alpha = .data[["edge_alpha"]]
      ),
      colour = "#c0c7d5",
      curvature = 0,
      lineend = "round",
      inherit.aes = FALSE
    )
  }
  plot <- plot +
    ggplot2::geom_curve(
      data = layout$membership,
      ggplot2::aes(
        x = .data[["x_term"]], y = .data[["y_term"]],
        xend = .data[["x_feature"]], yend = .data[["y_feature"]],
        colour = .data[["term_colour"]]
      ),
      linewidth = 0.25,
      alpha = 0.17,
      curvature = 0.1,
      lineend = "round",
      inherit.aes = FALSE
    ) +
    ggplot2::geom_point(
      data = layout$terms,
      ggplot2::aes(
        x = .data[["x"]], y = .data[["y"]],
        size = .data[["pathway_size"]], alpha = .data[["pathway_alpha"]],
        fill = .data[["colour"]]
      ),
      shape = 21,
      colour = "white",
      stroke = 0.45,
      inherit.aes = FALSE
    ) +
    ggplot2::geom_text(
      data = layout$terms,
      ggplot2::aes(
        x = .data[["x"]], y = .data[["y"]], label = .data[["label"]]
      ),
      size = .bulkmae_text_size_pt,
      size.unit = "pt",
      fontface = "plain",
      lineheight = 0.9,
      colour = "grey22",
      hjust = 0.5,
      vjust = 0.5,
      inherit.aes = FALSE
    ) +
    ggplot2::geom_point(
      data = layout$features,
      ggplot2::aes(
        x = .data[["x"]], y = .data[["y"]],
        size = .data[["feature_size"]], alpha = .data[["feature_alpha"]],
        fill = .data[["colour"]]
      ),
      shape = 21,
      colour = "white",
      stroke = 0.28,
      inherit.aes = FALSE
    )
  if (length(graph$label_features)) {
    labels <- layout$features[
      match(graph$label_features, layout$features$feature_id), , drop = FALSE
    ]
    plot <- plot + ggplot2::geom_text(
      data = labels,
      ggplot2::aes(
        x = .data[["label_x"]], y = .data[["label_y"]],
        label = .data[["label"]], angle = .data[["label_angle"]],
        hjust = .data[["label_hjust"]]
      ),
      size = .bulkmae_text_size_pt,
      size.unit = "pt",
      colour = "#474545",
      vjust = 0.5,
      inherit.aes = FALSE
    )
  }
  plot +
    ggplot2::scale_fill_identity() +
    ggplot2::scale_colour_identity() +
    ggplot2::scale_alpha_identity() +
    ggplot2::scale_size_identity() +
    ggplot2::scale_linewidth_identity() +
    ggplot2::coord_fixed(
      xlim = c(-layout$style$canvas_limit, layout$style$canvas_limit),
      ylim = c(-layout$style$canvas_limit, layout$style$canvas_limit),
      clip = "off"
    ) +
    ggplot2::labs(
      x = NULL,
      y = NULL,
      alt = paste(
        "A radial network of selected ORA terms on an inner ring and unique",
        "enriched features on an outer ring. Curves show term membership, pale",
        "inner-line width and alpha show shared-feature counts, feature size and",
        "alpha show absolute feature values when supplied, and term size shows",
        "enrichment evidence."
      )
    ) +
    ggplot2::theme_void(base_size = .bulkmae_text_size_pt) +
    ggplot2::theme(
      legend.position = "none",
      plot.background = ggplot2::element_rect(fill = "white", colour = NA),
      plot.margin = ggplot2::margin(5, 5, 5, 5, unit = "pt"),
      plot.title = ggplot2::element_text(
        face = "plain", hjust = 0.1, size = .bulkmae_text_size_pt,
        margin = ggplot2::margin(b = 0)
      )
    )
}
