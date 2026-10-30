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
  plot <- plot_surv_km(fit, test = test, risk_times = c(0, 250, 500, 750))
  .expect_bulkmae_plot(plot)
  expect_identical(plot$labels$y, "Survival probability")
  expect_identical(plot$labels$colour, "sex")
  built <- ggplot2::ggplot_build(plot)
  layers <- vapply(plot$layers, function(layer) class(layer$geom)[[1L]], character(1))
  expect_true(all(c("GeomRibbon", "GeomStep", "GeomPoint", "GeomText") %in% layers))
  risk <- built$data[[which(layers == "GeomText")[[1L]]]]
  expected <- summary(fit, times = c(0, 250, 500, 750), extend = TRUE)$n.risk
  expect_setequal(as.integer(risk$label), as.integer(expected))
  labels <- vapply(plot$layers, function(layer) {
    data <- layer$data
    if (is.data.frame(data) && "label" %in% names(data)) paste(data$label, collapse = "|") else ""
  }, character(1))
  expect_true(any(grepl("^Log-rank p = 0.001$", labels)))
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
  expect_error(plot_surv_km(stratified, risk_times = c(-1, 10)), "non-negative")
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
  expect_identical(plot$labels$x, "Hazard ratio (95% CI, log scale)")
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
  expect_identical(bottom$fill[bottom$x == 3], "#000000")
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
    sprintf("%s (AUC %.3f)", c("6 months", "1 year"), roc$AUC)
  )
  expect_error(plot_surv_roc(roc, time_labels = "one"), "one unique label")
  expect_error(plot_surv_roc(list()), "timeROC result")
})
