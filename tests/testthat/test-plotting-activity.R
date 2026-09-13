make_activity_plot_result <- function() {
  samples <- paste0("s", 1:6)
  sources <- c("TF_A", "TF_B", "TF_C")
  result <- expand.grid(
    source = sources,
    condition = samples,
    stringsAsFactors = FALSE
  )
  result$statistic <- "ulm"
  scores <- list(
    TF_A = c(-1.2, -0.8, -1.0, 1.1, 0.9, 1.3),
    TF_B = c(0.4, 0.2, 0.3, -0.5, -0.1, -0.4),
    TF_C = c(-0.2, 0.1, 0.0, 0.2, -0.1, 0.3)
  )
  result$score <- NA_real_
  for (src in sources) {
    result$score[result$source == src] <- scores[[src]]
  }
  result$p_value <- 0.2
  result
}

make_activity_plot_samples <- function() {
  data.frame(
    dex = factor(rep(c("untrt", "trt"), each = 3), levels = c("untrt", "trt")),
    cell = factor(rep(c("N61311", "N052611", "N080611"), 2)),
    row.names = paste0("s", 1:6),
    check.names = FALSE
  )
}

test_that("activity plots return buildable ggplots with activity axis labels", {
  result <- make_activity_plot_result()
  samples <- make_activity_plot_samples()
  heatmap <- plot_activity_heatmap(
    result,
    annotation = c("dex", "cell"),
    sample_data = samples,
    cluster_rows = FALSE,
    cluster_columns = FALSE
  )
  rank <- plot_activity_rank(result, sources = c("TF_C", "TF_A", "TF_B"))
  sample_plot <- plot_activity_sample(
    result,
    sources = c("TF_A", "TF_B"),
    sample_data = samples,
    group = "dex"
  )
  contrast <- activity_contrast(
    result, samples, group = "dex", reference = "untrt", target = "trt"
  )
  volcano <- plot_activity_volcano(
    contrast, label_sources = "TF_A"
  )

  for (plot in list(heatmap, rank, sample_plot, volcano)) {
    expect_s3_class(plot, "ggplot")
    expect_no_error(ggplot2::ggplot_build(plot))
    expect_true(nzchar(ggplot2::get_alt_text(plot)))
    dimensions <- attr(plot, "bulkmae_dimensions", exact = TRUE)
    expect_identical(dimensions$units, "cm")
    expect_true(dimensions$width > 0)
    expect_true(dimensions$height > 0)
  }
  expect_identical(heatmap$labels$x, "Sample")
  expect_identical(heatmap$labels$y, "Regulator")
  expect_identical(heatmap$scales$get_scales("fill")$name, "ULM activity score")
  expect_identical(rank$labels$x, "Mean ULM activity score")
  expect_identical(rank$labels$y, "Regulator")
  expect_identical(sample_plot$labels$y, "ULM activity score")
  expect_identical(sample_plot$labels$x, "dex")
  expect_match(sample_plot$labels$caption, "explicitly selected regulators")
  expect_match(sample_plot$labels$caption, "not a p-value top-N")
  expect_match(volcano$labels$x, "Activity difference \\(trt")
  expect_identical(as.character(rank$data$source), c("TF_C", "TF_A", "TF_B"))
})

test_that("activity heatmap annotation joins shuffled sample metadata by name", {
  result <- make_activity_plot_result()
  samples <- make_activity_plot_samples()
  samples <- samples[rev(rownames(samples)), , drop = FALSE]
  plot <- plot_activity_heatmap(
    result,
    sources = c("TF_A", "TF_B"),
    annotation = "dex",
    sample_data = samples,
    cluster_rows = FALSE,
    cluster_columns = FALSE
  )
  built <- ggplot2::ggplot_build(plot)
  expect_s3_class(plot, "ggplot")
  expect_true("untrt" %in% unique(as.character(samples$dex)))
  expect_gt(length(built$data), 1L)
  expect_error(
    plot_activity_heatmap(
      result,
      annotation = "dex",
      sample_data = samples[-1L, , drop = FALSE]
    ),
    "selected samples"
  )
  expect_error(
    plot_activity_heatmap(result, sources = "missing"),
    "Unknown `sources`"
  )
})

test_that("activity rank and sample plots keep explicit source order contracts", {
  result <- make_activity_plot_result()
  samples <- make_activity_plot_samples()
  rank <- plot_activity_rank(result, sources = c("TF_B", "TF_A"))
  expect_setequal(as.character(rank$data$source), c("TF_B", "TF_A"))
  expect_error(plot_activity_sample(result, sources = character(), group = "dex"),
               "explicitly contain")
  expect_error(
    plot_activity_sample(
      result, sources = "TF_A", sample_data = samples[-1L, , drop = FALSE],
      group = "dex"
    ),
    "selected samples"
  )
  named_group <- stats::setNames(as.character(samples$dex), rownames(samples))
  named_plot <- plot_activity_sample(
    result, sources = "TF_A", group = named_group[rev(names(named_group))]
  )
  expect_no_error(ggplot2::ggplot_build(named_plot))
  expect_identical(
    as.character(named_plot$data$group),
    as.character(named_group[as.character(named_plot$data$sample)])
  )
})

test_that("activity contrast volcano uses two-group activity p-values", {
  result <- make_activity_plot_result()
  samples <- make_activity_plot_samples()
  contrast <- activity_contrast(
    result, samples, group = "dex", reference = "untrt", target = "trt"
  )
  expect_identical(contrast$source, c("TF_A", "TF_B", "TF_C"))
  expect_true(contrast$delta[contrast$source == "TF_A"] > 0)
  expect_true(all(c("p_value", "adjusted_p_value", "delta") %in% names(contrast)))
  named <- stats::setNames(as.character(samples$dex), rownames(samples))
  named_contrast <- activity_contrast(
    result, group = named, reference = "untrt", target = "trt"
  )
  expect_equal(named_contrast$delta, contrast$delta)
  none <- activity_contrast(
    result, samples, group = "dex", test = "none"
  )
  expect_true(all(is.na(none$p_value)))
  expect_error(plot_activity_volcano(none), "Raw p-values are unavailable")
  expect_error(
    activity_contrast(result, samples, group = "cell"),
    "reference.*target"
  )

  volcano <- plot_activity_volcano(contrast, fdr = 1, min_abs_effect = 0)
  expect_true(all(c("Down", "Not significant", "Up") %in% levels(volcano$data$direction)))
  expect_setequal(as.character(volcano$data$direction), c("Up", "Down"))
  point_layer <- Filter(
    function(layer) inherits(layer$geom, "GeomPoint"),
    volcano$layers
  )[[1L]]
  expect_false(is.null(point_layer$mapping$colour))
  built <- ggplot2::ggplot_build(volcano)
  point_idx <- which(vapply(
    volcano$layers,
    function(layer) inherits(layer$geom, "GeomPoint"),
    logical(1)
  ))
  point_colours <- unique(built$data[[point_idx]]$colour)
  expect_setequal(point_colours, c("#0072B2", "#D55E00"))
  hline <- Filter(
    function(layer) inherits(layer$geom, "GeomHline"),
    volcano$layers
  )[[1L]]
  yintercept <- hline$aes_params$yintercept
  if (is.null(yintercept)) yintercept <- hline$data$yintercept
  expect_equal(yintercept, -log10(0.05))
  expect_match(volcano$labels$caption, "p = 0.05")
  expect_error(
    plot_activity_volcano(contrast, label_sources = "missing"),
    "Unknown `label_features`"
  )
})

test_that("activity plots warn when given enrichment-style statistics", {
  result <- make_activity_plot_result()
  result$statistic <- "ora"
  expect_warning(
    plot_activity_rank(result),
    "enrichment-style"
  )
})

test_that("activity plot dimensions survive ggplot additions", {
  result <- make_activity_plot_result()
  plot <- plot_activity_rank(result, sources = c("TF_A", "TF_B"))
  dimensions <- attr(plot, "bulkmae_dimensions", exact = TRUE)
  modified <- plot + ggplot2::labs(title = "Diagnostic activity")
  expect_identical(
    attr(modified, "bulkmae_dimensions", exact = TRUE),
    dimensions
  )
  expect_identical(plot$layers[[3L]]$aes_params$size, 1.8)
})
