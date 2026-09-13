test_that("survival plots return buildable ggplots with stable labels", {
  skip_if_not_installed("survival")
  mae <- mae_simulate(n_features = 40L, n_samples = 24L, seed = 12L)
  samples <- mae_samples(mae, "rna")
  score <- stats::setNames(samples$risk_score, rownames(samples))
  groups <- surv_risk_groups(score)
  derived <- data.frame(
    risk_group = groups,
    row.names = names(groups)
  )
  mae <- mae_add_sample_data(mae, "rna", derived)
  km <- surv_km(
    mae,
    "rna",
    surv_formula("survival_time", "event", predictors = "risk_group")
  )
  cox <- surv_cox(
    mae,
    "rna",
    surv_formula("survival_time", "event", predictors = c("risk_score", "age"))
  )
  table <- surv_cox_table(cox)
  km_plot <- plot_surv_km(km)
  forest_plot <- plot_surv_forest(cox)
  box_plot <- plot_surv_risk(score, group = groups, type = "box")
  density_plot <- plot_surv_risk(score, group = groups)

  for (plot in list(km_plot, forest_plot, box_plot, density_plot)) {
    expect_s3_class(plot, "ggplot")
    expect_no_error(ggplot2::ggplot_build(plot))
    expect_true(nzchar(ggplot2::get_alt_text(plot)))
    dimensions <- attr(plot, "bulkmae_dimensions", exact = TRUE)
    expect_identical(dimensions$units, "cm")
    expect_true(dimensions$width > 0)
    expect_true(dimensions$height > 0)
  }

  expect_identical(km_plot$labels$x, "Time (days)")
  expect_identical(km_plot$labels$y, "Survival probability")
  expect_identical(km_plot$labels$colour, "Stratum")
  expect_identical(forest_plot$labels$x, "Hazard ratio")
  expect_identical(forest_plot$labels$y, "Variable")
  expect_identical(box_plot$labels$x, "Group")
  expect_identical(box_plot$labels$y, "Risk score")
  expect_identical(density_plot$labels$x, "Risk score")
  expect_identical(density_plot$labels$y, "Density")

  km_text <- vapply(km_plot$layers, function(layer) {
    inherits(layer$geom, "GeomText") && isTRUE(grepl("Log-rank", layer$aes_params$label))
  }, logical(1))
  expect_true(any(km_text))
  box_text <- vapply(box_plot$layers, function(layer) {
    inherits(layer$geom, "GeomText")
  }, logical(1))
  expect_true(any(box_text))
  expect_true(any(vapply(box_plot$layers, function(layer) {
    inherits(layer$geom, "GeomPoint")
  }, logical(1))))
  expect_true(any(grepl("p =|p <", forest_plot$data$hr_label)))
  expect_identical(
    as.character(forest_plot$data$display),
    c("Age", "Risk score")[match(as.character(forest_plot$data$term), c("age", "risk_score"))]
  )

  forest_from_table <- plot_surv_forest(table)
  expect_no_error(ggplot2::ggplot_build(forest_from_table))
  expect_identical(
    attr(forest_plot + ggplot2::labs(title = "Cox"), "bulkmae_dimensions", exact = TRUE),
    attr(forest_plot, "bulkmae_dimensions", exact = TRUE)
  )
})

test_that("Kaplan-Meier risk tables and strata labels stay name-aligned", {
  skip_if_not_installed("survival")
  mae <- mae_simulate(n_features = 40L, n_samples = 16L, seed = 13L)
  km <- surv_km(
    mae,
    "rna",
    surv_formula("survival_time", "event", predictors = "condition")
  )
  plot <- plot_surv_km(km, risk_times = c(0, 50, 100))
  built <- ggplot2::ggplot_build(plot)
  expect_identical(levels(plot$layers[[2L]]$data$strata), c("control", "treated"))
  risk_layer <- vapply(plot$layers, function(layer) {
    inherits(layer$geom, "GeomText")
  }, logical(1))
  expect_true(any(risk_layer))
  risk_data <- plot$layers[[which(risk_layer)[[1L]]]]$data
  expect_true(all(risk_data$time %in% c(0, 50, 100)))
  expect_true(all(risk_data$n_risk >= 0L))

  no_table <- plot_surv_km(km, risk_table = FALSE, conf_int = FALSE, censor = FALSE)
  expect_no_error(ggplot2::ggplot_build(no_table))
  no_table_text <- vapply(no_table$layers, function(layer) {
    inherits(layer$geom, "GeomText")
  }, logical(1))
  expect_true(any(no_table_text))
  expect_true(grepl(
    "Log-rank",
    no_table$layers[[which(no_table_text)[[1L]]]]$aes_params$label
  ))
  hidden_p <- plot_surv_km(
    km,
    risk_table = FALSE,
    conf_int = FALSE,
    censor = FALSE,
    show_pvalue = FALSE
  )
  expect_false(any(vapply(hidden_p$layers, function(layer) {
    inherits(layer$geom, "GeomText")
  }, logical(1))))
})

test_that("forest and risk plots reject misaligned or invalid inputs", {
  skip_if_not_installed("survival")
  score <- c(a = -0.4, b = 0.1, c = 1.2, d = 0.7)
  group <- c(b = "event", a = "censor", d = "event", c = "censor")
  shuffled <- group[rev(names(group))]
  aligned_plot <- plot_surv_risk(score, shuffled, type = "box")
  expect_no_error(ggplot2::ggplot_build(aligned_plot))
  expect_identical(
    as.character(aligned_plot$data$group),
    as.character(group[aligned_plot$data$sample])
  )
  expect_error(plot_surv_risk(score, type = "box"), "group")
  expect_error(plot_surv_risk(score, group = group[-1L]), "match `score` names")
  expect_error(plot_surv_risk(unname(score)), "named numeric")

  bad_table <- data.frame(
    term = c("age", "age"),
    hazard_ratio = c(1.1, 0.8),
    conf_low = c(0.9, 0.5),
    conf_high = c(1.4, 1.2)
  )
  expect_error(plot_surv_forest(bad_table), "unique and non-missing")
  expect_error(
    plot_surv_forest(data.frame(
      term = "age",
      hazard_ratio = 1.2,
      conf_low = 1.5,
      conf_high = 2
    )),
    "contain its hazard ratio"
  )
  expect_error(plot_surv_km(list(time = 1)), "survfit")
  expect_error(
    plot_surv_forest(
      data.frame(
        term = "age",
        hazard_ratio = 1.2,
        conf_low = 0.8,
        conf_high = 1.6
      ),
      term_labels = c(missing = "Age")
    ),
    "unknown terms"
  )
})

test_that("survival plot dimensions are publication-oriented", {
  skip_if_not_installed("survival")
  mae <- mae_simulate(n_features = 40L, n_samples = 16L, seed = 14L)
  samples <- mae_samples(mae, "rna")
  km <- surv_km(mae, "rna", surv_formula("survival_time", "event"))
  cox <- surv_cox(
    mae,
    "rna",
    surv_formula("survival_time", "event", predictors = "age")
  )
  dimensions <- function(plot) attr(plot, "bulkmae_dimensions", exact = TRUE)
  expect_equal(
    dimensions(plot_surv_km(km, risk_table = FALSE))[c("width", "height")],
    list(width = 9.5, height = 6.8)
  )
  expect_equal(
    dimensions(plot_surv_risk(stats::setNames(samples$risk_score, rownames(samples))))[c("width", "height")],
    list(width = 8.5, height = 6.5)
  )
  expect_equal(
    dimensions(plot_surv_forest(cox))$units,
    "cm"
  )

  output <- tempfile(fileext = ".png")
  on.exit(unlink(output), add = TRUE)
  plot <- plot_surv_forest(cox)
  plot_dims <- dimensions(plot)
  expect_identical(plot_save(output, plot, dpi = 100), output)
  header <- readBin(output, what = "raw", n = 24L)
  png_width <- sum(as.integer(header[17:20]) * 256^(3:0))
  png_height <- sum(as.integer(header[21:24]) * 256^(3:0))
  expect_equal(png_width, round(plot_dims$width / 2.54 * 100), tolerance = 1)
  expect_equal(png_height, round(plot_dims$height / 2.54 * 100), tolerance = 1)
})

test_that("log-rank on the KM plot matches survdiff", {
  skip_if_not_installed("survival")
  mae <- mae_simulate(n_features = 40L, n_samples = 24L, seed = 15L)
  samples <- mae_samples(mae, "rna")
  groups <- surv_risk_groups(stats::setNames(samples$risk_score, rownames(samples)))
  mae <- mae_add_sample_data(
    mae,
    "rna",
    data.frame(risk_group = groups, row.names = names(groups))
  )
  samples <- mae_samples(mae, "rna")
  km <- surv_km(
    mae,
    "rna",
    surv_formula("survival_time", "event", predictors = "risk_group")
  )
  expected <- survival::survdiff(
    survival::Surv(survival_time, event) ~ risk_group,
    data = samples
  )
  expected_p <- stats::pchisq(expected$chisq, df = length(expected$n) - 1L, lower.tail = FALSE)
  expect_equal(bulkMAE:::.plot_surv_logrank_p(km), expected_p, tolerance = 1e-8)
})

test_that("AUC plot hook accepts a diagnostic table", {
  table <- data.frame(
    time = c(90, 180, 360),
    auc = c(0.62, 0.70, 0.66),
    n_case = c(8L, 14L, 19L),
    n_control = c(30L, 22L, 11L),
    estimator = "cumulative_dynamic"
  )
  plot <- plot_surv_auc(table)
  expect_s3_class(plot, "ggplot")
  expect_no_error(ggplot2::ggplot_build(plot))
  expect_identical(plot$labels$x, "Time (days)")
  expect_identical(plot$labels$y, "AUC(t)")
  expect_match(ggplot2::get_alt_text(plot), "cumulative_dynamic")
  expect_error(plot_surv_auc(list(auc = 0.6)), "timeROC result or an AUC table")
})
