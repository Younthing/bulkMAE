.lincs_example_table <- function() {
  drug_lincs_table(drug_lincs_example())
}

test_that("LINCS plots return buildable ggplots with score-aware labels", {
  table <- .lincs_example_table()
  rank <- plot_lincs_rank(table, n = 10L)
  heatmap <- plot_lincs_heatmap(
    table,
    compounds = c("toy_cp_01", "toy_cp_02", "toy_cp_07")
  )
  waterfall <- plot_lincs_waterfall(table, n = 12L)
  overlap <- plot_lincs_overlap(table, n = 8L)

  for (plot in list(rank, heatmap, waterfall, overlap)) {
    expect_s3_class(plot, "ggplot")
    expect_no_error(ggplot2::ggplot_build(plot))
    expect_true(nzchar(ggplot2::get_alt_text(plot)))
    dimensions <- attr(plot, "bulkmae_dimensions", exact = TRUE)
    expect_identical(dimensions$units, "cm")
    expect_true(dimensions$width > 0)
    expect_true(dimensions$height > 0)
  }

  expect_identical(rank$labels$x, "Normalized connectivity score (NCS)")
  expect_identical(rank$labels$y, "Compound (cell)")
  expect_identical(heatmap$labels$x, "Cell")
  expect_identical(heatmap$labels$y, "Compound")
  expect_identical(heatmap$scales$get_scales("fill")$name,
                   "Normalized connectivity score (NCS)")
  expect_identical(waterfall$labels$y, "Normalized connectivity score (NCS)")
  expect_identical(overlap$labels$x, "Overlapping query genes")
})

test_that("LINCS plot axes follow the requested score column", {
  table <- drug_lincs_example()
  tau <- plot_lincs_rank(table, n = 8L, score = "Tau")
  wtcs <- plot_lincs_waterfall(table, n = 8L, score = "WTCS")

  expect_identical(tau$labels$x, "Tau")
  expect_identical(wtcs$labels$y, "Weighted connectivity score (WTCS)")
  expect_no_error(ggplot2::ggplot_build(tau))
  expect_no_error(ggplot2::ggplot_build(wtcs))
})

test_that("LINCS direction colours stay colour-blind and complete", {
  plot <- plot_lincs_rank(drug_lincs_example(), n = 12L)
  built <- ggplot2::ggplot_build(plot)
  scale <- plot$scales$get_scales("colour")
  expect_identical(
    scale$palette(3),
    c(Reverse = "#0072B2", Unrelated = "#999999", Mimic = "#D55E00")
  )
  expect_setequal(
    as.character(unique(built$data[[3]]$colour)),
    c("#0072B2", "#D55E00")
  )
})

test_that("LINCS recommended sizes survive ggplot additions and plot_save", {
  plot <- plot_lincs_rank(drug_lincs_example(), n = 8L)
  dimensions <- attr(plot, "bulkmae_dimensions", exact = TRUE)
  modified <- plot +
    ggplot2::theme_classic() +
    ggplot2::labs(title = "Diagnostic LINCS ranking")

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
})

test_that("LINCS heatmap joins compound and cell by name", {
  table <- drug_lincs_example()
  compounds <- c("toy_cp_08", "toy_cp_01")
  cells <- c("PC3", "MCF7")
  plot <- plot_lincs_heatmap(table, compounds = compounds, cells = cells)
  built <- ggplot2::ggplot_build(plot)
  tile <- built$data[[1]]
  expect_identical(levels(plot$data$compound), rev(compounds))
  expect_identical(levels(plot$data$cell), cells)
  expect_equal(nrow(tile), 4L)

  expect_error(
    plot_lincs_heatmap(table, compounds = "missing_cp"),
    "requested `compounds`"
  )
  expect_error(
    plot_lincs_heatmap(table, cells = "HEPG2"),
    "requested `cells`"
  )
})

test_that("LINCS plot constructors reject broken score tables", {
  table <- drug_lincs_example()
  table$NCS[1L] <- Inf
  expect_error(plot_lincs_rank(table), "finite values")

  duplicated <- drug_lincs_example()[c(1L, 1L, 2L), , drop = FALSE]
  expect_error(plot_lincs_heatmap(duplicated), "unique compound and cell")

  no_overlap <- drug_lincs_example()
  no_overlap$N_upset <- NULL
  expect_error(plot_lincs_overlap(no_overlap), "N_upset")
  expect_error(plot_lincs_rank(drug_lincs_example(), n = 0L), "positive integer")
})
