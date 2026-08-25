.gsea_plot_fixture <- function() {
  ranks <- c(g1 = 3, g2 = 2, g3 = 1, g4 = -1, g5 = -2, g6 = -3)
  gene_sets <- list(
    positive = c("g1", "g2", "g5"),
    negative = c("g2", "g4", "g6")
  )
  scores <- vapply(gene_sets, function(members) {
    bulkMAE:::.plot_gsea_running_score(
      ranks,
      members,
      exponent = 1
    )$enrichment_score
  }, numeric(1))
  result <- data.frame(
    pathway = names(gene_sets),
    ES = unname(scores),
    NES = c(1.8, -1.5),
    pval = c(0.001, 0.02),
    padj = c(0.002, 0.03),
    qvalue = c(0.003, 0.04),
    size = lengths(gene_sets),
    stringsAsFactors = FALSE
  )
  result$leadingEdge <- I(list(c("g1", "g2"), c("g4", "g6")))
  list(result = result, ranks = ranks, gene_sets = gene_sets)
}

test_that("the pure running score matches the weighted GSEA definition", {
  fixture <- .gsea_plot_fixture()
  running <- bulkMAE:::.plot_gsea_running_score(
    fixture$ranks,
    fixture$gene_sets$positive,
    exponent = 1
  )

  expect_equal(
    running$data$running_score,
    c(3 / 7, 5 / 7, 5 / 7 - 1 / 3, 5 / 7 - 2 / 3,
      1 - 2 / 3, 0)
  )
  expect_equal(running$enrichment_score, 5 / 7)
  expect_identical(
    running$data$gene_id[running$data$hit],
    fixture$gene_sets$positive
  )

  expect_error(
    bulkMAE:::.plot_gsea_running_score(
      fixture$ranks,
      names(fixture$ranks),
      exponent = 1
    ),
    "cannot contain every gene"
  )
  expect_error(
    bulkMAE:::.plot_gsea_running_score(
      fixture$ranks,
      "absent",
      exponent = 1
    ),
    "no members"
  )
})

test_that("classic GSEA is one ggplot with intrinsic 50-20-30 panels", {
  fixture <- .gsea_plot_fixture()
  plot <- plot_gsea_classic(
    fixture$result,
    term = "positive",
    ranks = fixture$ranks,
    gene_sets = fixture$gene_sets,
    exponent = 1,
    term_label = "Positive response"
  )
  built <- ggplot2::ggplot_build(plot)

  expect_s3_class(plot, "ggplot")
  expect_no_error(ggplot2::ggplot_build(plot))
  expect_true(nzchar(ggplot2::get_alt_text(plot)))
  expect_identical(plot$labels$x, "Rank in Ordered Dataset")
  expect_identical(plot$labels$title, "Positive response")
  expect_null(plot$labels$subtitle)
  expect_s3_class(plot$facet, "FacetGrid")
  expect_equal(
    vapply(
      built$layout$panel_scales_y,
      function(scale) diff(scale$range$range),
      numeric(1)
    ),
    c(0.5, 0.2, 0.3)
  )
  expect_true(any(vapply(
    plot$layers,
    function(layer) inherits(layer$geom, "GeomRect"),
    logical(1)
  )))
  expect_true(any(vapply(
    plot$layers,
    function(layer) inherits(layer$geom, "GeomSegment"),
    logical(1)
  )))
  expect_true(any(vapply(
    plot$layers,
    function(layer) inherits(layer$geom, "GeomPolygon"),
    logical(1)
  )))
  text_layers <- plot$layers[vapply(
    plot$layers,
    function(layer) inherits(layer$geom, "GeomText"),
    logical(1)
  )]
  expect_true(all(vapply(
    text_layers,
    function(layer) identical(layer$geom_params$size.unit, "pt"),
    logical(1)
  )))
  expect_true(all(vapply(
    text_layers,
    function(layer) identical(layer$aes_params$size, 6),
    logical(1)
  )))
  expect_length(text_layers, 1L)
  expect_identical(
    text_layers[[1L]]$data$label,
    "NES: 1.8\nP.adj: 0.002\nFDR: 0.003"
  )
  expect_identical(text_layers[[1L]]$aes_params$fontface, "italic")

  curve_layer <- plot$layers[[
    which(vapply(
      plot$layers,
      function(layer) inherits(layer$geom, "GeomLine"),
      logical(1)
    ))
  ]]
  expect_identical(curve_layer$aes_params$linewidth, 0.6)
  expect_identical(
    plot$scales$get_scales("colour")$palette(c(0, 1)),
    c("#76BA99", "#EB4747")
  )
  expect_identical(
    plot$scales$get_scales("fill")$palette(c(0, 0.5, 1)),
    c("#08519C", "#FFFFFF", "#A50F15")
  )
  heat_layer <- plot$layers[[
    which(vapply(
      plot$layers,
      function(layer) inherits(layer$geom, "GeomRect"),
      logical(1)
    ))
  ]]
  expect_lte(nrow(heat_layer$data), 200L)
  expect_equal(
    range(heat_layer$data$rank_metric),
    as.numeric(stats::quantile(fixture$ranks, c(0.1, 0.9)))
  )
  expect_identical(heat_layer$aes_params$alpha, 0.8)
  polygon_layer <- plot$layers[[
    which(vapply(
      plot$layers,
      function(layer) inherits(layer$geom, "GeomPolygon"),
      logical(1)
    ))
  ]]
  expect_identical(polygon_layer$aes_params$fill, "grey70")

  y_scale <- plot$scales$get_scales("y")
  expect_true(all(c("-2", "0", "2") %in% y_scale$labels))
  dimensions <- attr(plot, "bulkmae_dimensions", exact = TRUE)
  expect_equal(
    dimensions,
    list(width = 5.9, height = 5.3, units = "cm")
  )
  expect_identical(
    attr(plot + ggplot2::theme_classic(), "bulkmae_dimensions", exact = TRUE),
    dimensions
  )
})

test_that("classic GSEA rejects incomplete or inconsistent tabular inputs", {
  fixture <- .gsea_plot_fixture()

  expect_error(
    plot_gsea_classic(fixture$result, "positive"),
    "'ranks' is required"
  )
  expect_error(
    plot_gsea_classic(
      fixture$result,
      "positive",
      ranks = fixture$ranks
    ),
    "'gene_sets' is required"
  )
  expect_error(
    plot_gsea_classic(
      fixture$result,
      "positive",
      ranks = fixture$ranks,
      gene_sets = fixture$gene_sets
    ),
    "'exponent' is required"
  )
  mismatched <- fixture$result
  mismatched$ES[[1L]] <- mismatched$ES[[1L]] + 0.1
  expect_error(
    plot_gsea_classic(
      mismatched,
      "positive",
      ranks = fixture$ranks,
      gene_sets = fixture$gene_sets,
      exponent = 1
    ),
    "do not reproduce the reported ES"
  )
  expect_error(
    plot_gsea_classic(
      fixture$result,
      "positive",
      ranks = rev(fixture$ranks),
      gene_sets = fixture$gene_sets,
      exponent = 1
    ),
    "sorted in decreasing order"
  )
  duplicate_ranks <- fixture$ranks
  names(duplicate_ranks)[[2L]] <- names(duplicate_ranks)[[1L]]
  expect_error(
    plot_gsea_classic(
      fixture$result,
      "positive",
      ranks = duplicate_ranks,
      gene_sets = fixture$gene_sets,
      exponent = 1
    ),
    "unique and non-missing"
  )
  expect_error(
    plot_gsea_classic(
      fixture$result,
      "positive",
      ranks = fixture$ranks,
      gene_sets = fixture$gene_sets["negative"],
      exponent = 1
    ),
    "missing selected terms"
  )
  expect_error(
    plot_gsea_classic(
      fixture$result,
      "positive",
      ranks = fixture$ranks,
      gene_sets = fixture$gene_sets,
      exponent = -1
    ),
    "non-negative"
  )
})

test_that("GSEA term identifiers and labels are never joined by position", {
  fixture <- .gsea_plot_fixture()

  expect_error(plot_gsea_classic(fixture$result, NULL), "explicitly")
  expect_error(
    plot_gsea_classic(fixture$result, c("positive", "negative")),
    "identify one term"
  )
  expect_error(
    plot_gsea_classic(fixture$result, "absent"),
    "Unknown GSEA terms"
  )

  duplicate <- fixture$result
  duplicate$pathway[[2L]] <- duplicate$pathway[[1L]]
  expect_error(
    plot_gsea_ridge(
      duplicate,
      "positive",
      ranks = fixture$ranks
    ),
    "term identifiers.*unique"
  )
  expect_error(
    plot_gsea_ridge(
      fixture$result,
      c("positive", "negative"),
      ranks = fixture$ranks,
      term_labels = c(positive = "Only one")
    ),
    "selected terms exactly"
  )
})

test_that("GSEA ridges preserve term order, memberships, and 6 pt labels", {
  fixture <- .gsea_plot_fixture()
  plot <- plot_gsea_ridge(
    fixture$result,
    terms = c("negative", "positive"),
    ranks = fixture$ranks,
    term_labels = c(
      positive = "Positive response",
      negative = "Negative response"
    ),
    show_statistics = TRUE
  )
  built <- ggplot2::ggplot_build(plot)

  expect_s3_class(plot, "ggplot")
  expect_no_error(ggplot2::ggplot_build(plot))
  expect_true(nzchar(ggplot2::get_alt_text(plot)))
  expect_identical(
    unique(plot$layers[[1L]]$data$term_id),
    c("negative", "positive")
  )
  expect_identical(
    unname(plot$scales$get_scales("y")$labels),
    c("Negative response", "Positive response")
  )
  expect_null(plot$labels$x)
  expect_null(plot$labels$y)
  expect_identical(
    unique(built$data[[1L]]$colour),
    c("#F6A282", "#4B63A7")
  )
  expect_s3_class(plot$layers[[1L]]$geom, "GeomRidgeline")
  expect_identical(plot$layers[[1L]]$aes_params$linewidth, 0.35)
  expect_identical(plot$layers[[1L]]$aes_params$alpha, 0.18)
  expect_identical(plot$layers[[1L]]$aes_params$scale, 1)
  expect_true(any(vapply(
    plot$layers,
    function(layer) inherits(layer$geom, "GeomSegment"),
    logical(1)
  )))
  text_layer <- plot$layers[[
    which(vapply(
      plot$layers,
      function(layer) inherits(layer$geom, "GeomText"),
      logical(1)
    ))
  ]]
  expect_identical(text_layer$geom_params$size.unit, "pt")
  expect_identical(text_layer$aes_params$size, 6)
  expect_identical(
    text_layer$data$label,
    c(
      "NES = -1.500\nP.adj = 0.030\nFDR = 0.040",
      "NES = 1.800\nP.adj = 0.002\nFDR = 0.003"
    )
  )
  expect_equal(
    text_layer$data$y,
    c(2, 1) + 0.36 * 0.92
  )
  barcode_layer <- plot$layers[[
    which(vapply(
      plot$layers,
      function(layer) inherits(layer$geom, "GeomSegment"),
      logical(1)
    ))
  ]]
  expect_identical(barcode_layer$aes_params$linewidth, 0.2)
  expect_identical(barcode_layer$aes_params$alpha, 0.52)
  expected_density <- stats::density(
    unname(fixture$ranks[c("g4", "g6")]),
    n = 512L,
    adjust = 1
  )
  first_density <- plot$layers[[1L]]$data
  first_density <- first_density[first_density$term_id == "negative", ]
  expect_equal(first_density$x, expected_density$x)
  expect_equal(max(first_density$height), 0.36)
  expect_equal(
    attr(plot, "bulkmae_dimensions", exact = TRUE),
    list(width = 12, height = 5, units = "cm")
  )

  gene_set_plot <- plot_gsea_ridge(
    fixture$result,
    terms = c("positive", "negative"),
    ranks = fixture$ranks,
    gene_sets = fixture$gene_sets,
    membership = "gene_set",
    show_statistics = FALSE
  )
  expect_no_error(ggplot2::ggplot_build(gene_set_plot))
  expect_false(any(vapply(
    gene_set_plot$layers,
    function(layer) inherits(layer$geom, "GeomText"),
    logical(1)
  )))
})

test_that("GSEA ridge recognizes the qvalues spelling used by reference tables", {
  fixture <- .gsea_plot_fixture()
  names(fixture$result)[names(fixture$result) == "qvalue"] <- "qvalues"

  plot <- plot_gsea_ridge(
    fixture$result,
    terms = c("negative", "positive"),
    ranks = fixture$ranks,
    show_statistics = TRUE
  )
  text_layer <- plot$layers[[
    which(vapply(
      plot$layers,
      function(layer) inherits(layer$geom, "GeomText"),
      logical(1)
    ))
  ]]

  expect_identical(
    text_layer$data$label,
    c(
      "NES = -1.500\nP.adj = 0.030\nFDR = 0.040",
      "NES = 1.800\nP.adj = 0.002\nFDR = 0.003"
    )
  )
})

test_that("GSEA ridges reject unavailable or mismatched memberships", {
  fixture <- .gsea_plot_fixture()

  expect_error(
    plot_gsea_ridge(fixture$result, terms = "positive"),
    "'ranks' is required"
  )
  expect_error(
    plot_gsea_ridge(
      fixture$result,
      terms = "positive",
      ranks = fixture$ranks,
      membership = "gene_set"
    ),
    "'gene_sets' is required"
  )
  no_leading <- fixture$result
  no_leading$leadingEdge <- NULL
  expect_error(
    plot_gsea_ridge(
      no_leading,
      terms = "positive",
      ranks = fixture$ranks
    ),
    "Leading-edge membership is unavailable"
  )
  unknown_leading <- fixture$result
  unknown_leading$leadingEdge[[1L]] <- c("g1", "absent")
  expect_error(
    plot_gsea_ridge(
      unknown_leading,
      terms = "positive",
      ranks = fixture$ranks
    ),
    "absent from 'ranks'"
  )
  one_gene <- fixture$result
  one_gene$leadingEdge[[1L]] <- "g1"
  expect_error(
    plot_gsea_ridge(
      one_gene,
      terms = "positive",
      ranks = fixture$ranks
    ),
    "At least two"
  )
  tied_ranks <- fixture$ranks
  tied_ranks[c("g1", "g2")] <- 3
  expect_error(
    plot_gsea_ridge(
      fixture$result,
      terms = "positive",
      ranks = tied_ranks
    ),
    "two distinct rank metric values"
  )
  expect_error(
    plot_gsea_ridge(
      fixture$result,
      terms = "positive",
      ranks = fixture$ranks,
      show_statistics = NA
    ),
    "TRUE or FALSE"
  )
  expect_error(
    bulkMAE:::.plot_gsea_palette(8L),
    "at most seven terms"
  )
})

test_that("GSEA source palettes and balanced label wrapping are preserved", {
  expect_identical(
    bulkMAE:::.plot_gsea_palette(7L),
    c(
      "#F6A282", "#4B63A7", "#25C6C2", "#FF6B5C",
      "#63D0EC", "#A78BCA", "#7EC87E"
    )
  )
  expect_identical(
    bulkMAE:::.plot_gsea_wrap_label(
      paste(
        "Genes encoding proteins involved in processing of drugs and",
        "other xenobiotics."
      ),
      width = 40L
    ),
    paste(
      "Genes encoding proteins involved",
      "in processing of drugs and other",
      "xenobiotics.",
      sep = "\n"
    )
  )
  expect_identical(
    bulkMAE:::.plot_gsea_classic_default_label(
      "HALLMARK_SIGNALING_OF_THE_CELL"
    ),
    "Signaling Of The Cell"
  )
})

test_that("native clusterProfiler GSEA supplies classic and ridge inputs", {
  suppressMessages(suppressWarnings(skip_if_not_installed("clusterProfiler")))
  suppressMessages(suppressWarnings(skip_if_not_installed("fgsea")))
  ranks <- seq(4, -4, length.out = 80L)
  names(ranks) <- paste0("g", seq_len(80L))
  gene_sets <- list(
    signal = paste0("g", 1:20),
    background = paste0("g", 21:40)
  )
  set.seed(11L)
  result <- suppressMessages(suppressWarnings(enrich_gsea(
    ranks = ranks,
    gene_sets = gene_sets,
    pvalueCutoff = 1,
    minGSSize = 5,
    maxGSSize = 30,
    eps = 0,
    verbose = FALSE,
    seed = TRUE
  )))

  classic <- plot_gsea_classic(result, term = "signal")
  ridge <- plot_gsea_ridge(
    result,
    terms = c("signal", "background")
  )

  expect_s3_class(classic, "ggplot")
  expect_s3_class(ridge, "ggplot")
  expect_no_error(ggplot2::ggplot_build(classic))
  expect_no_error(ggplot2::ggplot_build(ridge))
})
