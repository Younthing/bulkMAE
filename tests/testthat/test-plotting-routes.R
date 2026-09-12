.route_power_fit <- function() {
  list(
    fitIndices = data.frame(
      Power = 1:5,
      SFT.R.sq = c(0.21, 0.48, 0.73, 0.86, 0.88),
      mean.k. = c(18, 11, 7, 4.2, 2.8)
    ),
    powerEstimate = 4
  )
}

.route_trait <- function() {
  list(
    correlation = matrix(
      c(0.62, -0.41, 0.08, 0.21, 0.54, -0.33),
      nrow = 3L,
      dimnames = list(
        c("MEturquoise", "MEblue", "MEgrey"),
        c("condition", "batch")
      )
    ),
    p_value = matrix(
      c(0.01, 0.04, 0.41, 0.22, 0.03, 0.11),
      nrow = 3L,
      dimnames = list(
        c("MEturquoise", "MEblue", "MEgrey"),
        c("condition", "batch")
      )
    )
  )
}

.route_consensus <- function() {
  samples <- paste0("sample", 1:6)
  make_matrix <- function(seed) {
    set.seed(seed)
    value <- matrix(runif(36L, 0.2, 0.95), nrow = 6L)
    value <- (value + t(value)) / 2
    diag(value) <- 1
    dimnames(value) <- list(samples, samples)
    value
  }
  results <- vector("list", 3L)
  results[[2]] <- list(
    consensusMatrix = make_matrix(11L),
    consensusClass = stats::setNames(rep(1:2, each = 3L), samples)
  )
  results[[3]] <- list(
    consensusMatrix = make_matrix(12L),
    consensusClass = stats::setNames(rep(1:3, each = 2L), samples)
  )
  results
}

.route_fractions <- function() {
  matrix(
    c(
      0.4, 0.35, 0.2, 0.25,
      0.35, 0.4, 0.3, 0.25,
      0.25, 0.25, 0.5, 0.5
    ),
    nrow = 3L,
    dimnames = list(
      c("B_cell", "T_cell", "Myeloid"),
      paste0("sample", 1:4)
    )
  )
}

expect_route_plot <- function(plot) {
  expect_s3_class(plot, "ggplot")
  expect_no_error(ggplot2::ggplot_build(plot))
  expect_true(nzchar(ggplot2::get_alt_text(plot)))
  dimensions <- attr(plot, "bulkmae_dimensions", exact = TRUE)
  expect_identical(dimensions$units, "cm")
  expect_true(dimensions$width > 0)
  expect_true(dimensions$height > 0)
}

test_that("P0 co-expression plots build with stable labels", {
  modules <- stats::setNames(
    rep(c("turquoise", "blue", "grey"), length.out = 12L),
    paste0("gene", 1:12)
  )
  membership <- matrix(
    seq(-0.8, 0.8, length.out = 36L),
    nrow = 12L,
    dimnames = list(names(modules), c("MEturquoise", "MEblue", "MEgrey"))
  )
  gene_significance <- stats::setNames(seq(-0.5, 0.5, length.out = 12L), names(modules))
  power_plot <- plot_coexpr_power(.route_power_fit())
  module_plot <- plot_coexpr_modules(modules)
  trait_plot <- plot_coexpr_trait(.route_trait())
  membership_plot <- plot_coexpr_membership(
    membership,
    gene_significance,
    module = "blue"
  )

  for (plot in list(power_plot, module_plot, trait_plot, membership_plot)) {
    expect_route_plot(plot)
  }
  expect_identical(power_plot$labels$x, "Soft-threshold power")
  expect_identical(module_plot$labels$y, "Features")
  expect_identical(trait_plot$scales$get_scales("fill")$name, "Correlation")
  expect_identical(membership_plot$labels$x, "Module membership")
  expect_identical(trait_plot$layers[[2L]]$aes_params$size, 6)
  expect_identical(trait_plot$layers[[2L]]$geom_params$size.unit, "pt")
  expect_error(
    plot_coexpr_membership(membership, gene_significance[-1L], "blue"),
    "match membership features exactly"
  )
})

test_that("P0 subtype plots join consensus and class labels by name", {
  results <- .route_consensus()
  classes <- results[[2]]$consensusClass
  group <- stats::setNames(
    rep(c("control", "treated"), each = 3L),
    names(classes)
  )
  consensus_plot <- plot_cluster_consensus(results, k = 2L, show_names = TRUE)
  cdf_plot <- plot_cluster_cdf(results)
  size_plot <- plot_cluster_sizes(classes, group = group)
  matrix_plot <- plot_cluster_consensus(results[[2]]$consensusMatrix)

  for (plot in list(consensus_plot, cdf_plot, size_plot, matrix_plot)) {
    expect_route_plot(plot)
  }
  expect_identical(consensus_plot$scales$get_scales("fill")$name, "Consensus")
  expect_identical(cdf_plot$labels$y, "CDF")
  expect_identical(size_plot$labels$x, "Subtype")
  expect_error(plot_cluster_consensus(results), "`k` is required")
  expect_error(
    plot_cluster_sizes(classes, group = group[-1L]),
    "match the plotted samples exactly"
  )
})

test_that("P0 deconvolution plots use cell-type by sample fractions", {
  fractions <- .route_fractions()
  group <- stats::setNames(
    c("control", "control", "treated", "treated"),
    colnames(fractions)
  )
  stacked <- plot_deconv_stacked(fractions)
  boxed <- plot_deconv_box(fractions, group)
  heat <- plot_deconv_heatmap(fractions, column_split = group)

  for (plot in list(stacked, boxed, heat)) {
    expect_route_plot(plot)
  }
  expect_identical(stacked$labels$y, "Fraction")
  expect_identical(boxed$labels$x, "Cell type")
  expect_identical(heat$labels$y, "Cell type")
  expect_error(
    plot_deconv_box(fractions, group[-1L]),
    "match the plotted samples exactly"
  )
  expect_error(
    plot_deconv_stacked(fractions * -1),
    "non-negative"
  )
})
