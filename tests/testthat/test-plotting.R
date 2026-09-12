test_that("QC plots return buildable ggplot objects with stable labels", {
  mae <- make_toy_mae(n_features = 8L, n_samples = 4L)
  library_plot <- plot_qc_library(qc_library(mae, "rna"))
  correlation_plot <- plot_qc_correlation(
    qc_correlation(mae, "rna", "log_expression"), cluster = TRUE
  )
  pca <- reduce_pca(mae, "rna", "log_expression")
  outlier_plot <- plot_qc_outliers(qc_outliers(pca, components = 1:2))

  for (plot in list(library_plot, correlation_plot, outlier_plot)) {
    expect_s3_class(plot, "ggplot")
    expect_no_error(ggplot2::ggplot_build(plot))
    expect_true(nzchar(ggplot2::get_alt_text(plot)))
    dimensions <- attr(plot, "bulkmae_dimensions", exact = TRUE)
    expect_identical(dimensions$units, "cm")
    expect_true(dimensions$width > 0)
    expect_true(dimensions$height > 0)
  }
  expect_identical(library_plot$labels$x, "Sample")
  expect_identical(correlation_plot$scales$get_scales("fill")$name, "Correlation")
  expect_identical(correlation_plot$theme$axis.text.x$angle, 45)
  expect_identical(outlier_plot$labels$y, "Robust distance")
})

test_that("embedding metadata are joined by sample identifier", {
  mae <- make_toy_mae(n_features = 12L, n_samples = 4L)
  pca <- reduce_pca(mae, "rna", "log_expression")
  samples <- mae_samples(mae, "rna")
  samples$display <- paste("Sample", seq_len(nrow(samples)))
  samples <- samples[rev(rownames(samples)), , drop = FALSE]
  plot <- plot_embedding(
    pca, samples, colour = "condition", shape = "batch", label = "display"
  )

  expect_s3_class(plot, "ggplot")
  expect_no_error(ggplot2::ggplot_build(plot))
  expect_match(plot$labels$x, "^PC1 \\(")
  expect_identical(as.character(plot$data$condition),
                   as.character(samples[plot$data$sample, "condition"]))
  expect_error(
    plot_embedding(pca, samples[-1L, , drop = FALSE], colour = "condition"),
    "match embedding samples exactly"
  )
  duplicated <- matrix(seq_len(10), nrow = 5L,
                       dimnames = list(c(rownames(samples), rownames(samples)[1L]), c("a", "b")))
  expect_error(plot_embedding(pca, duplicated), "match embedding samples exactly|unique and non-missing")
  reserved <- samples
  reserved$x <- seq_len(nrow(reserved))
  expect_no_error(ggplot2::ggplot_build(
    plot_embedding(pca, reserved, colour = "condition")
  ))
  expect_error(plot_embedding(pca, reserved, label = "x"), "reserved embedding names")
  reserved$numeric_shape <- seq_len(nrow(reserved))
  expect_error(
    plot_embedding(pca, reserved, shape = "numeric_shape"),
    "shape must be discrete"
  )
  expect_identical(plot$layers[[2L]]$aes_params$size, 6)
  expect_identical(plot$layers[[2L]]$geom_params$size.unit, "pt")

  many_coordinates <- matrix(
    seq_len(18L),
    ncol = 2L,
    dimnames = list(paste0("s", seq_len(9L)), c("x", "y"))
  )
  many_samples <- data.frame(
    group = factor(letters[seq_len(9L)]),
    row.names = rownames(many_coordinates)
  )
  many_plot <- plot_embedding(
    many_coordinates,
    many_samples,
    colour = "group"
  )
  expect_s3_class(many_plot, "ggplot")
  expect_error(
    ggplot2::ggplot_build(many_plot),
    "Insufficient values in manual scale"
  )
  expect_no_error(suppressMessages(ggplot2::ggplot_build(
    many_plot + ggplot2::scale_colour_viridis_d()
  )))
})

test_that("recommended physical dimensions survive ggplot additions and drive export", {
  mae <- make_toy_mae(n_features = 8L, n_samples = 4L)
  plot <- plot_qc_library(qc_library(mae, "rna"))
  dimensions <- attr(plot, "bulkmae_dimensions", exact = TRUE)
  modified <- plot +
    ggplot2::theme_classic() +
    ggplot2::labs(title = "QC") +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = 0.02))

  expect_identical(
    attr(modified, "bulkmae_dimensions", exact = TRUE),
    dimensions
  )

  output <- tempfile(fileext = ".png")
  on.exit(unlink(output), add = TRUE)
  expect_identical(plot_save(output, modified, dpi = 100), output)
  expect_true(file.exists(output))
  header <- readBin(output, what = "raw", n = 24L)
  png_width <- sum(as.integer(header[17:20]) * 256^(3:0))
  png_height <- sum(as.integer(header[21:24]) * 256^(3:0))
  expect_equal(png_width, round(dimensions$width / 2.54 * 100), tolerance = 1)
  expect_equal(png_height, round(dimensions$height / 2.54 * 100), tolerance = 1)

  ordinary <- ggplot2::ggplot(data.frame(x = 1, y = 1), ggplot2::aes(x, y)) +
    ggplot2::geom_point()
  expect_error(plot_save(tempfile(fileext = ".png"), ordinary),
               "width.*height.*required")
})

test_that("plot families carry publication-oriented dimension presets", {
  mae <- make_toy_mae(n_features = 8L, n_samples = 4L)
  pca <- reduce_pca(mae, "rna", "log_expression")
  dimensions <- function(plot) attr(plot, "bulkmae_dimensions", exact = TRUE)

  expect_equal(dimensions(plot_qc_library(qc_library(mae, "rna")))[c("width", "height")],
               list(width = 8.5, height = 11.5))
  expect_equal(dimensions(plot_embedding(pca))[c("width", "height")],
               list(width = 8.5, height = 7))
  expect_equal(dimensions(plot_qc_correlation(
    qc_correlation(mae, "rna", "log_expression")
  ))[c("width", "height")], list(width = 9.5, height = 8.5))
  expect_equal(dimensions(plot_qc_outliers(
    qc_outliers(pca, components = 1:2)
  ))[c("width", "height")], list(width = 8.5, height = 6.5))
})

test_that("matrix and limma MDS embeddings are supported", {
  coordinates <- matrix(1:12, nrow = 4L,
                        dimnames = list(paste0("s", 1:4), paste0("D", 1:3)))
  expect_no_error(ggplot2::ggplot_build(plot_embedding(coordinates, axes = c(3, 1))))

  skip_if_not_installed("limma")
  set.seed(10)
  mds <- limma::plotMDS(matrix(rnorm(80), nrow = 20,
                              dimnames = list(NULL, paste0("s", 1:4))), plot = FALSE)
  plot <- plot_embedding(mds)
  expect_match(plot$labels$x, "^MDS1")
  expect_no_error(ggplot2::ggplot_build(plot))
})

test_that("DE plots share threshold classification and explicit labels", {
  result <- data.frame(
    feature_id = paste0("g", 1:5),
    effect = c(-2, -0.5, 0, 1.5, 3),
    standard_error = rep(0.2, 5),
    statistic = c(-10, -2.5, 0, 7.5, 15),
    p_value = c(0, 0.01, NA, 0.04, 1e-8),
    adjusted_p_value = c(0.001, 0.02, NA, 0.04, 0.001),
    abundance = c(10, 20, 30, 40, 50),
    row.names = paste0("g", 1:5),
    check.names = FALSE
  )
  volcano <- plot_de_volcano(
    result, fdr = 0.05, min_abs_effect = 1,
    label_features = "g5",
    feature_labels = stats::setNames(paste0("Gene ", 1:5), paste0("g", 1:5))
  )
  ma <- plot_de_ma(
    result, fdr = 0.05, min_abs_effect = 1,
    abundance_column = "abundance", label_features = "g5"
  )

  expect_identical(as.character(volcano$data$direction),
                   c("Down", "Not significant", "Not significant", "Up", "Up"))
  expect_identical(as.character(ma$data$direction), as.character(volcano$data$direction))
  expect_true(is.finite(volcano$data$minus_log10_p[1L]))
  expect_identical(volcano$data$label[5L], "Gene 5")
  expect_false(any(vapply(
    volcano$layers,
    function(layer) inherits(layer$geom, "GeomHline"),
    logical(1)
  )))
  expect_equal(volcano$coordinates$limits$x, c(-3, 3))
  expect_identical(volcano$theme$legend.position, "top")
  expect_identical(volcano$theme$legend.direction, "horizontal")
  count_layer <- Filter(
    function(layer) inherits(layer$geom, "GeomText"),
    volcano$layers
  )[[1L]]
  expect_identical(count_layer$aes_params$size, 6)
  expect_identical(count_layer$geom_params$size.unit, "pt")
  expect_setequal(count_layer$data$label, c("Down 1", "Up 2"))
  expect_true(any(vapply(
    volcano$layers,
    function(layer) inherits(layer$geom, "GeomTextRepel"),
    logical(1)
  )))
  expect_true(any(vapply(
    ma$layers,
    function(layer) inherits(layer$geom, "GeomTextRepel"),
    logical(1)
  )))
  expect_no_error(ggplot2::ggplot_build(volcano))
  expect_no_error(ggplot2::ggplot_build(ma))
  asymmetric <- plot_de_volcano(
    result, fdr = 0.05, min_abs_effect = 1, symmetric_x = FALSE
  )
  expect_null(asymmetric$coordinates$limits$x)
  expect_error(plot_de_volcano(result, symmetric_x = NA), "TRUE or FALSE")
  zero_significant <- result
  zero_significant$p_value[3L] <- 1e-6
  zero_significant$adjusted_p_value[3L] <- 1e-5
  expect_false(de_selected(zero_significant, min_abs_effect = 0)[["g3"]])
  expect_identical(
    as.character(plot_de_volcano(
      zero_significant, min_abs_effect = 0
    )$data$direction[3L]),
    "Not significant"
  )
  expect_identical(names(de_table(result)), c(
    "feature_id", "effect", "standard_error", "statistic", "p_value",
    "adjusted_p_value"
  ))
})

test_that("DE plots reject non-directional tests and ambiguous abundance", {
  result <- data.frame(
    feature_id = c("g1", "g2"), effect = c(-1, 1), F = c(4, 9),
    p_value = c(0.02, 0.01), adjusted_p_value = c(0.03, 0.02),
    baseMean = c(10, 20), Amean = c(2, 3), row.names = c("g1", "g2")
  )
  expect_error(plot_de_volcano(result), "no single effect direction")

  directional <- result[c("feature_id", "effect", "p_value", "adjusted_p_value",
                          "baseMean", "Amean")]
  expect_error(plot_de_ma(directional), "ambiguous")
  expect_error(plot_de_volcano(directional, label_features = "missing"),
               "Unknown `label_features`")
})

test_that("DE plots handle missing FDR and reject unusable result inputs", {
  missing_fdr <- data.frame(
    feature_id = c("g1", "g2"), effect = c(-2, 2),
    p_value = c(0.01, 0.02), adjusted_p_value = c(NA_real_, 0.02),
    abundance = c(10, 20), row.names = c("g1", "g2")
  )
  plot <- plot_de_volcano(missing_fdr)
  expect_identical(as.character(plot$data$direction),
                   c("Not significant", "Up"))
  expect_no_error(ggplot2::ggplot_build(plot))

  all_na_p <- missing_fdr
  all_na_p$p_value <- NA_real_
  all_na_p$adjusted_p_value <- NA_real_
  expect_error(plot_de_volcano(all_na_p), "Raw p-values are unavailable")
  expect_error(plot_de_volcano(list(effect = 1)), "must be a DESeqResults")
  zero_abundance <- missing_fdr
  zero_abundance$abundance[1L] <- 0
  expect_warning(
    zero_plot <- plot_de_ma(
      zero_abundance, abundance_column = "abundance",
      abundance_transform = "log10"
    ),
    "Removed 1 zero abundance"
  )
  expect_true(is.na(zero_plot$data$abundance[1L]))
  expect_no_error(ggplot2::ggplot_build(zero_plot))
  negative_abundance <- missing_fdr
  negative_abundance$abundance[1L] <- -1
  expect_error(
    plot_de_ma(
      negative_abundance, abundance_column = "abundance",
      abundance_transform = "log10"
    ),
    "non-negative"
  )
  factor_abundance <- missing_fdr
  factor_abundance$abundance <- factor(c("10", "20"))
  expect_error(
    plot_de_ma(
      factor_abundance, abundance_column = "abundance",
      abundance_transform = "identity"
    ),
    "numeric values"
  )
})

test_that("assay heatmap aligns annotations and removes constant scaled rows", {
  mae <- make_toy_mae(n_features = 5L, n_samples = 4L)
  values <- mae_pull_assay(mae, "rna", "log_expression")
  values[1L, ] <- 2
  mae <- mae_add_assay(mae, "rna", values, name = "heat")
  symbols <- stats::setNames(paste0("Symbol", 1:5), rownames(values))
  feature_data <- data.frame(
    symbol = unname(symbols[rev(names(symbols))]),
    row.names = rev(names(symbols))
  )
  mae <- mae_add_feature_data(mae, "rna", feature_data)

  expect_warning(
    plot <- plot_assay_heatmap(
      mae, "rna", "heat", features = rownames(values)[1:4],
      column_split = "condition", feature_label = "symbol"
    ),
    "Removed constant heatmap rows"
  )
  expect_s3_class(plot, "ggplot")
  expect_no_error(ggplot2::ggplot_build(plot))
  expect_false("Symbol1" %in% as.character(plot$data$feature))
  expect_error(plot_assay_heatmap(mae, "rna", "heat", features = character()),
               "explicitly contain")
  expect_error(
    plot_assay_heatmap(
      mae, "rna", "heat", features = rownames(values)[1:2],
      feature_label = stats::setNames(c("a", "b"), c("gene1", "missing"))
    ),
    "selected features exactly"
  )
  duplicate_labels <- stats::setNames(
    rep("Repeated", 2L),
    rownames(values)[2:3]
  )
  duplicate_plot <- plot_assay_heatmap(
    mae, "rna", "heat", features = rownames(values)[2:3],
    scale = "none", cluster_rows = FALSE, cluster_columns = FALSE,
    feature_label = duplicate_labels
  )
  expect_setequal(levels(duplicate_plot$data$feature), c("Repeated", "Repeated.1"))
  expect_no_error(ggplot2::ggplot_build(duplicate_plot))
})

test_that("assay expression joins groups by sample identifier", {
  mae <- make_toy_mae(n_features = 5L, n_samples = 4L)
  features <- c("gene1", "gene3")
  plot <- plot_assay_expression(
    mae, "rna", "log_expression",
    features = features,
    colour = "condition",
    feature_label = stats::setNames(c("Gene 1", "Gene 3"), features)
  )
  expect_s3_class(plot, "ggplot")
  expect_no_error(ggplot2::ggplot_build(plot))
  expect_identical(levels(plot$data$feature), c("Gene 1", "Gene 3"))
  expect_identical(
    as.character(plot$data$group),
    as.character(mae_samples(mae, "rna")[as.character(plot$data$sample), "condition"])
  )
  expect_error(
    plot_assay_expression(mae, "rna", "log_expression", features = character()),
    "explicitly contain"
  )
  expect_error(
    plot_assay_expression(mae, "rna", "log_expression", features = "missing"),
    "Unknown features"
  )
  expect_error(
    plot_assay_expression(
      mae, "rna", "log_expression", features = "gene2", geom = "violin"
    ),
    "at least two samples"
  )
  expect_no_error(ggplot2::ggplot_build(
    plot_assay_expression(
      mae, "rna", "log_expression",
      features = "gene2", colour = "condition", geom = "violin"
    )
  ))
})

test_that("backend MA abundance conventions do not change de_table", {
  mae <- make_model_mae(n_features = 80L, n_samples = 8L)

  skip_if_not_installed("edgeR")
  edge <- de_edger(mae, "rna", ~ batch + condition, filter = FALSE,
                   robust = FALSE, coef = "conditiontreated")
  expect_no_error(ggplot2::ggplot_build(plot_de_ma(edge)))
  expect_identical(names(de_table(edge)), c(
    "feature_id", "effect", "standard_error", "statistic", "p_value",
    "adjusted_p_value"
  ))

  skip_if_not_installed("limma")
  limma_fit <- de_limma(mae, "rna", ~ batch + condition, assay = "log_expression")
  expect_no_error(ggplot2::ggplot_build(plot_de_ma(limma_fit, coef = "conditiontreated")))

  skip_if_not_installed("DESeq2")
  deseq_fit <- de_deseq2(mae, "rna", ~ batch + condition,
                         fitType = "mean", quiet = TRUE)
  deseq_result <- de_deseq2_results(
    deseq_fit, contrast = c("condition", "treated", "control")
  )
  deseq_plot <- plot_de_ma(deseq_result)
  expect_s3_class(deseq_plot$scales$get_scales("x"), "ScaleContinuousPosition")
  expect_no_error(ggplot2::ggplot_build(deseq_plot))
})

test_that("an unnamed DESeqTransform assay can be added", {
  skip_if_not_installed("DESeq2")
  mae <- make_model_mae(n_features = 1200L, n_samples = 8L)
  transformed <- transform_vst(
    mae, "rna", blind = FALSE, design = ~ batch + condition, fit_type = "mean"
  )
  expect_null(SummarizedExperiment::assayNames(transformed))
  updated <- mae_add_assay(mae, "rna", transformed, name = "vst")
  expect_true("vst" %in% mae_assays(updated, "rna"))
  expect_equal(
    mae_pull_assay(updated, "rna", "vst"),
    SummarizedExperiment::assay(transformed, 1L)
  )
})
