.survival_lung <- function() {
  data <- survival::lung
  data$sex <- factor(data$sex, levels = c(1, 2), labels = c("Male", "Female"))
  data$status <- data$status - 1L
  data
}

.survival_glmnet <- function() {
  .with_seed(20261030, {
    x <- matrix(stats::rnorm(120 * 12), 120, dimnames = list(NULL, paste0("g", 1:12)))
    y <- survival::Surv(stats::rexp(120, exp(x[, 1] - x[, 2])), stats::rbinom(120, 1, 0.7))
    glmnet::cv.glmnet(x, y, family = "cox", nfolds = 5)
  })
}

.expect_bulkmae_plot <- function(plot) {
  expect_s3_class(plot, "ggplot")
  expect_no_error(ggplot2::ggplot_build(plot))
  expect_true(nzchar(ggplot2::get_alt_text(plot)))
  dimensions <- attr(plot, "bulkmae_dimensions", exact = TRUE)
  expect_identical(dimensions$units, "cm")
  expect_true(dimensions$width > 0 && dimensions$height > 0)
}

test_that("plot_surv_km draws strata, censor marks, risk table, and log-rank p", {
  skip_if_not_installed("survival")
  data <- .survival_lung()
  fit <- survival::survfit(survival::Surv(time, status) ~ sex, data = data)
  test <- survival::survdiff(survival::Surv(time, status) ~ sex, data = data)
  cox <- survival::coxph(survival::Surv(time, status) ~ sex, data = data)
  plot <- plot_surv_km(fit, test = test, cox = cox, hr_label = "HR",
                       risk_table = TRUE, conf_int = TRUE,
                       risk_times = c(0, 250, 500, 750))
  .expect_bulkmae_plot(plot)
  expect_identical(plot$labels$y, "Survival probability")
  built <- ggplot2::ggplot_build(plot)
  layers <- vapply(plot$layers, function(layer) class(layer$geom)[[1L]], character(1))
  expect_true(all(c("GeomRibbon", "GeomStep", "GeomPoint", "GeomText") %in% layers))
  risk_layer <- which(vapply(plot$layers, function(layer) {
    is.data.frame(layer$data) && "n_risk" %in% names(layer$data)
  }, logical(1)))
  risk <- built$data[[risk_layer]]
  expected <- summary(fit, times = c(0, 250, 500, 750), extend = TRUE)$n.risk
  expect_setequal(as.integer(risk$label), as.integer(expected))
  labels <- vapply(plot$layers, function(layer) {
    data <- layer$data
    if (is.data.frame(data) && "label" %in% names(data)) paste(data$label, collapse = "|") else ""
  }, character(1))
  hr <- sprintf("%.2f", exp(stats::coef(cox)[[1L]]))
  expect_true(any(grepl(paste0("^HR = ", hr, " .*
Log-rank P = 0.001$"), labels)))
  expect_identical(
    unname(plot$scales$get_scales("colour")$palette(2)),
    c("#B64C3C", "#3D6B8C")
  )
  # Confidence ribbons stay in the survival panel.
  ribbon <- built$data[[which(layers == "GeomRibbon")]]
  expect_identical(unique(as.integer(ribbon$PANEL)), 1L)
})

test_that("plot_surv_km handles unstratified fits and rejects bad inputs", {
  skip_if_not_installed("survival")
  data <- .survival_lung()
  fit <- survival::survfit(survival::Surv(time, status) ~ 1, data = data)
  plot <- plot_surv_km(fit, risk_table = FALSE)
  .expect_bulkmae_plot(plot)
  expect_error(plot_surv_km(list()), "single-state survfit")
  other <- survival::survdiff(survival::Surv(time, status) ~ ph.ecog, data = data)
  stratified <- survival::survfit(survival::Surv(time, status) ~ sex, data = data)
  expect_error(plot_surv_km(stratified, test = other), "different numbers of strata")
  expect_error(plot_surv_km(stratified, risk_table = TRUE, risk_times = c(-1, 10)), "non-negative")
  expect_error(plot_surv_km(stratified, cox = list()), "coxph fit")
})

test_that("plot_surv_forest keeps right-hand estimates aligned with rows", {
  skip_if_not_installed("survival")
  data <- .survival_lung()
  fit <- survival::coxph(survival::Surv(time, status) ~ age + sex + ph.ecog, data = data)
  plot <- plot_surv_forest(fit, term_labels = c(age = "Age", sexFemale = "Female", ph.ecog = "ECOG"))
  .expect_bulkmae_plot(plot)
  expect_identical(levels(plot$data$term), c("ECOG", "Female", "Age"))
  hr <- exp(stats::coef(fit))
  expect_equal(plot$data$hazard_ratio, unname(hr))
  secondary <- plot$scales$get_scales("y")$secondary.axis
  expect_match(secondary$labels[[3L]], sprintf("^%.2f ", hr[["age"]]))
  expect_match(secondary$labels[[1L]], sprintf("^%.2f ", hr[["ph.ecog"]]))
  expect_identical(plot$labels$x, "Hazard ratio (95% CI)")
  # Rows run top to bottom, so ECOG sits at y = 1 and age at y = 3.
  built <- ggplot2::ggplot_build(plot)$data[[3L]]
  expect_identical(built$colour[built$y == 1], "#B64C3C")
  expect_identical(built$colour[built$y == 3], "#6E6A66")
})

test_that("plot_surv_forest accepts a named list of univariate fits", {
  skip_if_not_installed("survival")
  data <- .survival_lung()
  fits <- list(
    Age = survival::coxph(survival::Surv(time, status) ~ age, data = data),
    Sex = survival::coxph(survival::Surv(time, status) ~ sex, data = data)
  )
  plot <- plot_surv_forest(fits, terms = c("Sex", "Age"))
  .expect_bulkmae_plot(plot)
  expect_identical(as.character(plot$data$term_id), c("Sex", "Age"))
  expect_error(plot_surv_forest(fits, terms = "Grade"), "Unknown `terms`")
  expect_error(plot_surv_forest(unname(fits)), "unique and non-missing")
})

test_that("glmnet plots mark the chosen penalty and colour by coefficient sign", {
  skip_if_not_installed("glmnet")
  skip_if_not_installed("survival")
  fit <- .survival_glmnet()
  cv <- plot_glmnet_cv(fit)
  .expect_bulkmae_plot(cv)
  expect_identical(cv$labels$y, "Partial Likelihood Deviance")
  vlines <- ggplot2::ggplot_build(cv)$data[[1L]]$xintercept
  expect_equal(sort(vlines), sort(log(c(fit$lambda.min, fit$lambda.1se))))

  path <- plot_glmnet_path(fit, lambda = "lambda.1se", label_features = c("g1", "g2"))
  .expect_bulkmae_plot(path)
  chosen <- as.matrix(stats::coef(fit, s = "lambda.1se"))[, 1]
  observed <- unique(path$data[c("feature", "status")])
  expected <- ifelse(chosen > 0, "Positive", ifelse(chosen < 0, "Negative", "Zero"))
  expect_identical(
    as.character(observed$status[match(names(chosen), observed$feature)]),
    unname(expected)
  )
  expect_error(plot_glmnet_path(fit, label_features = "nope"), "Unknown `label_features`")
  expect_error(plot_glmnet_cv(list()), "cv.glmnet")
})

test_that("plot_surv_risk joins metadata by sample name and splits at the cutoff", {
  score <- c(b = 2, a = -1, c = 0.5, d = 3)
  samples <- data.frame(
    time = c(10, 20, 5, 1), event = c(0, 1, 1, 0),
    row.names = c("a", "b", "c", "d")
  )
  plot <- plot_surv_risk(score, samples, "time", "event", cutoff = 1)
  .expect_bulkmae_plot(plot)
  built <- ggplot2::ggplot_build(plot)
  top <- built$data[[3L]]
  bottom <- built$data[[4L]]
  expect_equal(top$x[order(top$y)], 1:4)
  # Sample b has rank 3 and was followed for 20 time units with an event.
  expect_equal(bottom$y[bottom$x == 3], 20)
  expect_identical(bottom$fill[bottom$x == 3], "#2C2926")
  expect_error(plot_surv_risk(score, samples[1:3, ], "time", "event", 1), "exactly")
  samples$event[1] <- 2
  expect_error(plot_surv_risk(score, samples, "time", "event", 1), "coded 0")
  expect_error(plot_surv_risk(unname(score), samples, "time", "event", 1), "sample-named")
})

test_that("plot_surv_roc draws one curve per time with AUC labels", {
  skip_if_not_installed("timeROC")
  skip_if_not_installed("survival")
  data <- .survival_lung()
  data <- data[stats::complete.cases(data[c("time", "status", "age")]), ]
  rownames(data) <- paste0("patient", seq_len(nrow(data)))
  expression <- matrix(
    0, nrow = 1L, ncol = nrow(data),
    dimnames = list("feature", rownames(data))
  )
  mae <- mae_from_matrix(expression, data[c("time", "status", "age")])
  roc <- surv_roc(mae, "rna", "time", "status", "age", times = c(180, 365))
  plot <- plot_surv_roc(roc, time_labels = c("6 months", "1 year"))
  .expect_bulkmae_plot(plot)
  expect_identical(
    levels(plot$data$time),
    sprintf("%s  AUC = %.3f", c("6 months", "1 year"), roc$AUC)
  )
  expect_error(plot_surv_roc(roc, time_labels = "one"), "one unique label")
  expect_error(plot_surv_roc(list()), "timeROC result")
})

.survival_mae <- function() {
  data <- .survival_lung()
  data <- data[stats::complete.cases(data[c("time", "status", "age", "sex", "ph.ecog")]), ]
  rownames(data) <- paste0("patient", seq_len(nrow(data)))
  expression <- matrix(
    0, nrow = 1L, ncol = nrow(data),
    dimnames = list("feature", rownames(data))
  )
  mae_from_matrix(expression, data[c("time", "status", "age", "sex", "ph.ecog")])
}

test_that("plot_surv_auc compares predictors across shared times", {
  skip_if_not_installed("timeROC")
  skip_if_not_installed("survival")
  mae <- .survival_mae()
  rocs <- list(
    Age = surv_roc(mae, "rna", "time", "status", "age", times = c(180, 365)),
    ECOG = surv_roc(mae, "rna", "time", "status", "ph.ecog", times = c(180, 365))
  )
  plot <- plot_surv_auc(rocs, time_labels = c("6 months", "1 year"))
  .expect_bulkmae_plot(plot)
  expect_equal(plot$data$auc[plot$data$predictor == "ECOG"], unname(rocs$ECOG$AUC))
  expect_identical(levels(plot$data$predictor), c("Age", "ECOG"))
  other <- surv_roc(mae, "rna", "time", "status", "age", times = c(100, 365))
  expect_error(plot_surv_auc(list(A = rocs$Age, B = other)), "same evaluation times")
  expect_error(plot_surv_auc(rocs$Age), "named list")
})

test_that("surv_calibration groups predicted survival and plots against KM", {
  skip_if_not_installed("survival")
  mae <- .survival_mae()
  fit <- surv_cox(mae, "rna", surv_formula("time", "status", c("age", "sex", "ph.ecog")),
                  model = TRUE)
  calibration <- surv_calibration(mae, "rna", fit, "time", "status",
                                  times = c(180, 365), groups = 3L)
  expect_identical(
    names(calibration),
    c("time", "group", "n", "predicted", "observed", "lower", "upper")
  )
  expect_identical(sum(calibration$n[calibration$time == 180]), nrow(mae_samples(mae, "rna")))
  expect_true(all(calibration$predicted > 0 & calibration$predicted < 1))
  plot <- plot_surv_calibration(calibration, time_labels = c("180 d", "365 d"))
  .expect_bulkmae_plot(plot)
  expect_identical(levels(plot$data$time_label), c("180 d", "365 d"))
  expect_error(surv_calibration(mae, "rna", fit, "time", "status", times = 0), "positive")
  expect_error(plot_surv_calibration(calibration[, -4]), "missing columns")
})

test_that("surv_dca returns a native dca object that plots with references", {
  skip_if_not_installed("survival")
  skip_if_not_installed("dcurves")
  mae <- .survival_mae()
  models <- list(
    Full = surv_cox(mae, "rna", surv_formula("time", "status", c("age", "sex", "ph.ecog")),
                    model = TRUE),
    Age = surv_cox(mae, "rna", surv_formula("time", "status", "age"), model = TRUE)
  )
  decision <- surv_dca(mae, "rna", models, "time", "status", horizon = 365,
                       thresholds = seq(0.1, 0.8, by = 0.1))
  expect_s3_class(decision, "dca")
  expect_setequal(as.character(unique(decision$dca$label)),
                  c("Full", "Age", "Treat All", "Treat None"))
  plot <- plot_surv_dca(decision)
  .expect_bulkmae_plot(plot)
  palette <- plot$scales$get_scales("colour")$palette(4)
  expect_identical(unname(palette[c("Full", "Treat None")]), c("#B64C3C", "#2C2926"))
  expect_error(surv_dca(mae, "rna", models[[1L]], "time", "status", 365), "named list")
  expect_error(plot_surv_dca(list()), "decision-curve")
})
