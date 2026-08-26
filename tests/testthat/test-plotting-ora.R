.toy_ora_plot_result <- function() {
  data.frame(
    ID = c("t1", "t2", "t3"),
    Description = c("Immune response", "Cell cycle", "Stress response"),
    GeneRatio = c("3/10", "3/10", "2/10"),
    BgRatio = c("10/100", "12/100", "8/100"),
    pvalue = c(0, 0.01, 0.03),
    p.adjust = c(0.001, 0.02, 0.04),
    Count = c(3, 3, 2),
    geneID = c("g1/g2/g3", "g3/g4/g5", "g2/g6"),
    stringsAsFactors = FALSE
  )
}

test_that("ORA bubble keeps ratio metrics statistically distinct", {
  result <- .toy_ora_plot_result()
  rich <- plot_ora_bubble(result, terms = c("t3", "t1", "t2"))
  ratio <- plot_ora_bubble(result, terms = c("t1", "t2", "t3"), x = "gene_ratio")
  fold <- plot_ora_bubble(
    result, terms = c("t1", "t2", "t3"), x = "fold_enrichment"
  )

  expect_s3_class(rich, "ggplot")
  expect_no_error(ggplot2::ggplot_build(rich))
  expect_true(nzchar(ggplot2::get_alt_text(rich)))
  expect_identical(rich$data$term_id, c("t3", "t1", "t2"))
  expect_equal(rich$data$x_value, c(2 / 8, 3 / 10, 3 / 12))
  expect_equal(ratio$data$x_value, c(3 / 10, 3 / 10, 2 / 10))
  expect_equal(fold$data$x_value, c(3, 2.5, 2.5))
  expect_identical(rich$labels$x, "Rich factor (%)")
  expect_identical(ratio$labels$x, "Gene ratio (%)")
  expect_identical(fold$labels$x, "Fold enrichment")
  expect_identical(rich$labels$size, "Gene count")
  expect_true(all(is.finite(rich$data$minus_log10_p)))
  expect_identical(rich$theme$legend.position, "right")
  expect_equal(rich$guides$guides$size$params$order, 1)
  expect_equal(rich$guides$guides$colour$params$order, 2)
  expect_equal(rich$layers[[1L]]$aes_params$alpha, 0.92)
  expect_identical(
    rich$scales$get_scales("colour")$palette(c(0, 1)),
    c("#D9EAF7", "#B2182B")
  )
  expect_equal(
    rich$scales$get_scales("size")$palette(c(0, 1)),
    c(0.8, 2.8)
  )
  expect_true(ggplot2::is_waiver(rich$scales$get_scales("size")$breaks))
  expect_identical(
    rich$scales$get_scales("x")$labels(c(0.05, 0.055)),
    c("5.0%", "5.5%")
  )
  expect_identical(
    deparse1(rich$scales$get_scales("colour")$name),
    "expression(-log[10](\"adj p\"))"
  )
  expect_equal(
    rich$scales$get_scales("x")$expand,
    c(0.01, 0, 0.08, 0)
  )
  expect_equal(rich$scales$get_scales("x")$n.breaks, 4)
  expect_identical(rich$theme$panel.grid.major$colour, "grey90")
  expect_equal(rich$theme$panel.grid.major$linewidth, 0.2)
  expect_identical(rich$theme$panel.border$colour, "black")
  expect_equal(rich$theme$panel.border$linewidth, 0.2)
  expect_identical(rich$theme$plot.background$fill, "white")
  expect_equal(as.numeric(rich$theme$axis.ticks.length), 0.8)
  expect_equal(rich$theme$plot.title$size, 6)
  expect_equal(rich$theme$plot.title$hjust, 0.5)
  expect_equal(
    attr(rich, "bulkmae_dimensions", exact = TRUE),
    list(width = 10, height = 8, units = "cm")
  )
  expect_identical(
    attr(rich + ggplot2::theme_classic(), "bulkmae_dimensions", exact = TRUE),
    attr(rich, "bulkmae_dimensions", exact = TRUE)
  )
})

test_that("ORA bubble validates term joins, evidence, and ratio contracts", {
  result <- .toy_ora_plot_result()
  expect_error(plot_ora_bubble(result, NULL), "explicitly")
  expect_error(plot_ora_bubble(result, c("t1", "t1")), "unique")
  expect_error(plot_ora_bubble(result, "missing"), "Unknown `terms`")
  expect_error(
    plot_ora_bubble(result, c("t1", "t2"), term_labels = c(t1 = "A")),
    "selected terms exactly"
  )

  partial <- result
  partial$p.adjust[[2L]] <- NA_real_
  fallback <- plot_ora_bubble(partial, c("t1", "t2"))
  expect_identical(
    deparse1(fallback$scales$get_scales("colour")$name),
    "expression(-log[10](\"p-value\"))"
  )
  expect_error(
    plot_ora_bubble(partial, c("t1", "t2"), p_value = "adjusted"),
    "unavailable"
  )

  malformed <- result
  malformed$GeneRatio[[1L]] <- "3 of 10"
  expect_error(plot_ora_bubble(malformed, "t1"), "non-negative-integer")
  mismatch <- result
  mismatch$GeneRatio[[1L]] <- "2/10"
  expect_error(plot_ora_bubble(mismatch, "t1"), "numerator.*Count")
  invalid_background <- result
  invalid_background$BgRatio[[1L]] <- "2/100"
  expect_error(plot_ora_bubble(invalid_background, "t1"), "at least.*Count")
})

test_that("ORA bubble uses the reference 40-character term wrapping", {
  result <- .toy_ora_plot_result()[1L, , drop = FALSE]
  result$Description <- "1234567890123456789012345678901234567890"

  plot <- plot_ora_bubble(result, terms = "t1")

  expect_identical(as.character(plot$data$term), result$Description)
})

test_that("ORA graph has unique features and one complete Jaccard contract", {
  result <- .toy_ora_plot_result()
  values <- c(g6 = -0.4, g4 = -2, g2 = -1, g1 = 2, g5 = 1.2, g3 = 0.5)
  graph <- bulkMAE:::.plot_ora_graph(
    result = result,
    terms = c("t1", "t2", "t3"),
    membership = NULL,
    feature_values = values,
    features = NULL,
    term_labels = NULL,
    feature_labels = NULL,
    label_features = c("g2", "g3"),
    p_value = "auto"
  )

  expect_identical(graph$features, paste0("g", 1:6))
  expect_equal(unname(graph$feature_values), c(2, -1, 0.5, -2, 1.2, -0.4))
  expect_equal(nrow(graph$membership_edges), 8L)
  expect_equal(
    graph$overlap$jaccard[graph$overlap$from == "t1" & graph$overlap$to == "t2"],
    1 / 5
  )
  expect_equal(
    graph$overlap$jaccard[graph$overlap$from == "t1" & graph$overlap$to == "t3"],
    1 / 4
  )

  subset_graph <- bulkMAE:::.plot_ora_graph(
    result = result,
    terms = c("t1", "t2", "t3"),
    membership = NULL,
    feature_values = values,
    features = c("g1", "g3", "g4", "g6"),
    term_labels = NULL,
    feature_labels = NULL,
    label_features = NULL,
    p_value = "auto"
  )
  expect_equal(subset_graph$overlap$jaccard, graph$overlap$jaccard)
  expect_equal(nrow(subset_graph$membership_edges), 5L)

  per_term_graph <- bulkMAE:::.plot_ora_graph(
    result = result,
    terms = c("t1", "t2", "t3"),
    membership = NULL,
    feature_values = values,
    features = list(
      t1 = c("g1", "g2"),
      t2 = c("g3", "g4"),
      t3 = "g6"
    ),
    term_labels = NULL,
    feature_labels = NULL,
    label_features = NULL,
    p_value = "auto"
  )
  expect_identical(lengths(per_term_graph$membership), c(t1 = 2L, t2 = 2L, t3 = 1L))
  expect_equal(nrow(per_term_graph$membership_edges), 5L)
  expect_equal(per_term_graph$overlap$jaccard, graph$overlap$jaccard)
})

test_that("ORA community layout keeps one visual feature node per term", {
  result <- .toy_ora_plot_result()
  values <- c(g1 = 2, g2 = -1, g3 = 0.5, g4 = -2, g5 = 1.2, g6 = -0.4)
  graph <- bulkMAE:::.plot_ora_graph(
    result, c("t1", "t2", "t3"), NULL, values, NULL,
    NULL, NULL, NULL, "auto"
  )

  layout <- bulkMAE:::.plot_ora_community_layout(graph)
  node_keys <- paste(layout$features$term_id, layout$features$feature_id, sep = "\r")
  edge_keys <- paste(
    graph$membership_edges$term_id,
    graph$membership_edges$feature_id,
    sep = "\r"
  )

  expect_equal(nrow(layout$features), nrow(graph$membership_edges))
  expect_setequal(node_keys, edge_keys)
  expect_equal(anyDuplicated(node_keys), 0L)
  expect_equal(sum(layout$features$feature_id == "g2"), 2L)
  expect_equal(sum(layout$features$feature_id == "g3"), 2L)
  expect_equal(nrow(layout$membership), nrow(layout$features))
})

test_that("ORA community plot has the reference white background and no gene labels", {
  result <- .toy_ora_plot_result()
  values <- c(g1 = 2, g2 = -1, g3 = 0.5, g4 = -2, g5 = 1.2, g6 = -0.4)
  plot <- plot_ora_network(
    result,
    terms = c("t1", "t2", "t3"),
    feature_values = values
  )

  expect_identical(plot$theme$plot.background$fill, "white")
  expect_identical(
    unname(vapply(
      plot$layers, function(layer) class(layer$geom)[[1L]], character(1)
    )),
    c(
      "GeomCurve", "GeomCircle", "GeomSegment", "GeomPoint", "GeomPoint",
      "GeomTextRepel"
    )
  )
})

test_that("ORA radial plot uses shared counts and a white background", {
  result <- .toy_ora_plot_result()
  values <- c(g1 = 2, g2 = -1, g3 = 0.5, g4 = -2, g5 = 1.2, g6 = -0.4)
  graph <- bulkMAE:::.plot_ora_graph(
    result, c("t1", "t2", "t3"), NULL, values, NULL,
    NULL, NULL, names(values), "auto"
  )
  graph$overlap <- data.frame(
    from = c("t1", "t1"),
    to = c("t2", "t3"),
    shared_n = c(10L, 2L),
    jaccard = c(0.1, 0.8),
    stringsAsFactors = FALSE
  )

  layout <- bulkMAE:::.plot_ora_radial_layout(graph)
  plot <- bulkMAE:::.plot_ora_radial_layers(graph, layout)

  expect_gt(layout$overlap$edge_width[[1L]], layout$overlap$edge_width[[2L]])
  expect_gt(layout$overlap$edge_alpha[[1L]], layout$overlap$edge_alpha[[2L]])
  background <- ggplot2::calc_element("plot.background", plot$theme)
  expect_identical(background$fill, "white")
  expect_true(is.na(background$colour))
  expect_equal(nrow(plot$layers[[6L]]$data), length(values))
})

test_that("ORA community and radial plots reproduce the reference visual grammar", {
  result <- .toy_ora_plot_result()
  values <- c(g1 = 2, g2 = -1, g3 = 0.5, g4 = -2, g5 = 1.2, g6 = -0.4)
  network <- plot_ora_network(
    result,
    terms = c("t1", "t2", "t3"),
    feature_values = values,
    label_features = c("g2", "g3")
  )
  radial <- plot_ora_radial(
    result,
    terms = c("t1", "t2", "t3"),
    feature_values = values,
    label_features = names(values)
  )

  for (plot in list(network, radial)) {
    expect_s3_class(plot, "ggplot")
    expect_no_error(ggplot2::ggplot_build(plot))
    expect_true(nzchar(ggplot2::get_alt_text(plot)))
    expect_identical(plot$theme$legend.position, "none")
    expect_equal(plot$theme$text$size, 6)
    expect_equal(
      attr(plot, "bulkmae_dimensions", exact = TRUE),
      list(width = 8, height = 8, units = "cm")
    )
  }

  graph <- bulkMAE:::.plot_ora_graph(
    result, c("t1", "t2", "t3"), NULL, values, NULL,
    NULL, NULL, names(values), "auto"
  )
  first <- bulkMAE:::.plot_ora_community_layout(graph)
  second <- bulkMAE:::.plot_ora_community_layout(graph)
  expect_identical(first$terms, second$terms)
  expect_identical(first$features, second$features)
  expect_identical(
    unname(graph$colours),
    c("#6C63B8", "#F2A65A", "#8DAA22")
  )
  expect_equal(nrow(first$halos), length(graph$terms) * 15L)
  expect_equal(first$halos$a, first$halos$b)
  expect_equal(range(first$halos$halo_alpha), c(0.001, 0.018))
  expect_equal(range(first$features$plot_size), c(2, 3.5))
  expect_equal(range(first$features$feature_alpha), c(0.35, 0.58))
  expect_equal(range(first$terms$hub_size), c(3.5, 5))
  expect_equal(range(first$overlap$edge_width), c(0.8, 1.3))
  expect_equal(range(first$overlap$edge_alpha), c(0.35, 0.55))
  expect_identical(
    unname(vapply(
      network$layers, function(layer) class(layer$geom)[[1L]], character(1)
    )),
    c(
      "GeomCurve", "GeomCircle", "GeomSegment", "GeomPoint", "GeomPoint",
      "GeomTextRepel", "GeomTextRepel"
    )
  )
  expect_identical(network$layers[[1L]]$aes_params$colour, "#c7d0e4")
  expect_equal(network$layers[[1L]]$geom_params$curvature, 0)
  expect_identical(network$layers[[1L]]$geom_params$lineend, "square")
  expect_identical(network$layers[[3L]]$aes_params$colour, "#beb0b2")
  expect_equal(network$layers[[3L]]$aes_params$linewidth, 0.45)
  expect_equal(network$layers[[3L]]$aes_params$alpha, 0.35)
  expect_equal(network$layers[[6L]]$aes_params$size, 6 / (72.27 / 25.4))
  expect_identical(network$layers[[6L]]$aes_params$fontface, "plain")
  expect_true(is.na(network$layers[[6L]]$aes_params$bg.colour))
  expect_equal(network$layers[[6L]]$aes_params$bg.r, 0.1)
  expect_equal(network$theme$plot.title$size, 6)
  expect_equal(network$theme$plot.title$hjust, 0.5)

  radial_layout <- bulkMAE:::.plot_ora_radial_layout(graph)
  expect_identical(
    radial_layout$terms$colour,
    c("#6C63B8", "#2BB3A3", "#F2A65A")
  )
  expect_identical(sort(radial_layout$features$feature_id), sort(graph$features))
  expect_equal(anyDuplicated(radial_layout$features$feature_id), 0L)
  expect_equal(
    sum(radial_layout$membership$feature_id == "g2"),
    2L
  )
  expect_equal(nrow(radial_layout$membership), nrow(graph$membership_edges))
  expect_equal(
    sqrt(radial_layout$features$x^2 + radial_layout$features$y^2),
    rep(6.3, length(graph$features))
  )
  expect_equal(
    sqrt(radial_layout$features$label_x^2 + radial_layout$features$label_y^2),
    rep(6.65, length(graph$features))
  )
  expect_true(all(abs(
    sqrt(radial_layout$terms$x^2 + radial_layout$terms$y^2) - 4.1
  ) < 0.05))
  expect_equal(
    sum(radial_layout$group_arcs$arc_end - radial_layout$group_arcs$arc_start) +
      nrow(radial_layout$group_arcs) * 0.075,
    2 * pi
  )
  expect_equal(range(radial_layout$features$feature_size), c(1.5, 3))
  expect_equal(range(radial_layout$features$feature_alpha), c(0.45, 0.77))
  expect_equal(range(radial_layout$terms$pathway_size), c(6, 12))
  expect_equal(range(radial_layout$terms$pathway_alpha), c(0.35, 0.55))
  expect_equal(radial_layout$overlap$edge_width, rep(1.1, 2L))
  expect_equal(radial_layout$overlap$edge_alpha, rep(0.325, 2L))
  left <- (radial_layout$features$angle * 180 / pi) %% 360 > 90 &
    (radial_layout$features$angle * 180 / pi) %% 360 < 270
  expect_identical(radial_layout$features$label_hjust, ifelse(left, 1, 0))
  expect_identical(
    unname(vapply(
      radial$layers, function(layer) class(layer$geom)[[1L]], character(1)
    )),
    c("GeomCurve", "GeomCurve", "GeomPoint", "GeomText", "GeomPoint", "GeomText")
  )
  expect_identical(radial$layers[[1L]]$aes_params$colour, "#c0c7d5")
  expect_equal(radial$layers[[2L]]$aes_params$linewidth, 0.25)
  expect_equal(radial$layers[[2L]]$aes_params$alpha, 0.17)
  expect_equal(radial$layers[[2L]]$geom_params$curvature, 0.1)
  expect_equal(radial$layers[[3L]]$aes_params$stroke, 0.45)
  expect_equal(radial$layers[[4L]]$aes_params$size, 6)
  expect_identical(radial$layers[[4L]]$geom_params$size.unit, "pt")
  expect_identical(radial$layers[[4L]]$aes_params$fontface, "plain")
  expect_equal(radial$layers[[4L]]$aes_params$hjust, 0.5)
  expect_equal(radial$layers[[4L]]$aes_params$vjust, 0.5)
  expect_equal(radial$layers[[5L]]$aes_params$stroke, 0.28)
  expect_equal(radial$layers[[6L]]$aes_params$size, 6)
  expect_identical(radial$layers[[6L]]$geom_params$size.unit, "pt")
  expect_equal(radial$theme$plot.title$size, 6)
  expect_equal(radial$theme$plot.title$hjust, 0.1)
})

test_that("ORA graph inputs fail before IDs can be guessed", {
  result <- .toy_ora_plot_result()
  values <- c(g1 = 2, g2 = -1, g3 = 0.5, g4 = -2, g5 = 1.2, g6 = -0.4)
  expect_error(
    plot_ora_network(result, c("t1", "t2"), membership = list(t1 = c("g1"))),
    "selected terms exactly"
  )
  wrong_count <- list(
    t1 = c("g1", "g2"),
    t2 = c("g3", "g4", "g5")
  )
  expect_error(
    plot_ora_network(result, c("t1", "t2"), membership = wrong_count),
    "does not match"
  )
  expect_error(
    plot_ora_network(
      result, c("t1", "t2", "t3"),
      feature_values = values[names(values) != "g6"]
    ),
    "missing displayed features"
  )
  expect_error(
    plot_ora_network(
      result, c("t1", "t2", "t3"), features = c("g1", "g3", "g4")
    ),
    "without displayed nodes.*t3"
  )
  expect_error(
    plot_ora_radial(
      result, c("t1", "t2"), label_features = "missing"
    ),
    "Unknown `label_features`"
  )
  expect_error(
    plot_ora_network(
      result, c("t1", "t2"),
      feature_labels = c(g1 = "A", g2 = "B")
    ),
    "missing displayed features"
  )
})

test_that("ORA graph palettes scale to reference-sized term selections", {
  term_id <- paste0("t", seq_len(10L))
  feature_id <- paste0("g", seq_len(10L))
  result <- data.frame(
    ID = term_id,
    Description = paste("Reference term", seq_len(10L)),
    GeneRatio = rep("1/10", 10L),
    BgRatio = rep("10/100", 10L),
    pvalue = seq(0.001, 0.01, length.out = 10L),
    p.adjust = seq(0.002, 0.02, length.out = 10L),
    Count = rep(1, 10L),
    geneID = feature_id,
    stringsAsFactors = FALSE
  )
  values <- stats::setNames(seq_along(feature_id), feature_id)
  graph <- bulkMAE:::.plot_ora_graph(
    result, term_id, NULL, values, NULL,
    NULL, NULL, NULL, "auto"
  )

  expect_length(graph$colours, 10L)
  expect_equal(anyDuplicated(graph$colours), 0L)
  expect_no_error(ggplot2::ggplot_build(
    plot_ora_network(result, term_id, feature_values = values)
  ))
  expect_no_error(ggplot2::ggplot_build(
    plot_ora_radial(result, term_id, feature_values = values)
  ))
})

test_that("native clusterProfiler ORA results feed all ORA plots", {
  suppressMessages(suppressWarnings(skip_if_not_installed("clusterProfiler")))
  sets <- list(
    signal = paste0("g", 1:20),
    shared = paste0("g", 11:30),
    background = paste0("g", 31:50)
  )
  result <- suppressMessages(suppressWarnings(enrich_ora(
    genes = paste0("g", 1:18),
    gene_sets = sets,
    universe = paste0("g", 1:80),
    pvalueCutoff = 1,
    qvalueCutoff = 1,
    minGSSize = 5,
    maxGSSize = 30
  )))
  table <- as.data.frame(result)
  terms <- utils::head(table$ID, 2L)
  features <- unique(unlist(strsplit(table$geneID[match(terms, table$ID)], "/")))
  values <- stats::setNames(seq_along(features), features)

  expect_no_error(ggplot2::ggplot_build(plot_ora_bubble(result, terms)))
  expect_no_error(ggplot2::ggplot_build(
    plot_ora_network(result, terms, feature_values = values)
  ))
  expect_no_error(ggplot2::ggplot_build(
    plot_ora_radial(result, terms, feature_values = values)
  ))
})
