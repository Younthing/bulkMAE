test_that("representative plotting families remain visually stable", {
  skip_if_not_installed("vdiffr")

  mae <- make_toy_mae(n_features = 8L, n_samples = 4L)
  metrics <- qc_library(mae, "rna")
  pca <- reduce_pca(mae, "rna", "log_expression")
  de <- data.frame(
    feature_id = paste0("g", 1:5),
    effect = c(-2, -0.5, 0, 1.5, 3),
    standard_error = rep(0.2, 5),
    statistic = c(-10, -2.5, 0, 7.5, 15),
    p_value = c(0.001, 0.1, 1, 0.04, 1e-8),
    adjusted_p_value = c(0.005, 0.2, 1, 0.05, 0.001),
    abundance = c(10, 20, 30, 40, 50),
    row.names = paste0("g", 1:5)
  )

  vdiffr::expect_doppelganger("library QC", plot_qc_library(metrics))
  vdiffr::expect_doppelganger("PCA embedding", plot_embedding(pca))
  vdiffr::expect_doppelganger(
    "DE volcano",
    plot_de_volcano(de, fdr = 0.05, min_abs_effect = 1)
  )
  vdiffr::expect_doppelganger(
    "assay heatmap",
    plot_assay_heatmap(
      mae, "rna", "log_expression", features = paste0("gene", 1:5),
      cluster_rows = FALSE, cluster_columns = FALSE
    )
  )
})
