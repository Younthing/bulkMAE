test_that("decoupleR method discovery is a real offline backend call", {
  skip_if_not_installed("decoupleR", minimum_version = "2.0.0")

  methods <- suppressMessages(suppressWarnings(activity_methods()))

  expect_s3_class(methods, "data.frame")
  expect_gt(nrow(methods), 0L)
  expect_true(all(c("Function", "Name") %in% names(methods)))
  expect_true("run_ulm" %in% methods$Function)
})

test_that("decoupleR infers local weighted activities for every sample", {
  skip_if_not_installed("decoupleR", minimum_version = "2.0.0")
  mae <- make_simulated_mae(seed = 1701L)
  samples <- rownames(mae_samples(mae, "rna"))
  network <- data.frame(
    source = rep(c("regulator_a", "regulator_b"), each = 30L),
    target = sprintf("gene%04d", seq_len(60L)),
    weight = rep(c(-1, 1), 30L),
    stringsAsFactors = FALSE
  )

  result <- suppressMessages(suppressWarnings(activity_decouple(
    mae,
    "rna",
    network = network,
    assay = "log_expression",
    statistics = "ulm",
    method_args = list(NULL),
    consensus = FALSE,
    mor = "weight",
    min_size = 10L
  )))

  expect_s3_class(result, "data.frame")
  expect_identical(nrow(result), 2L * length(samples))
  expect_true(all(c(
    "statistic", "source", "condition", "score", "p_value"
  ) %in% names(result)))
  expect_identical(unique(result$statistic), "ulm")
  expect_setequal(result$source, c("regulator_a", "regulator_b"))
  expect_setequal(result$condition, samples)
  expect_true(all(is.finite(result$score)))
  expect_true(all(is.finite(result$p_value)))
})

test_that("resource-backed activity functions are opt-in online tests", {
  run_online <- identical(
    tolower(Sys.getenv("BULKMAE_RUN_ONLINE_TESTS", "false")),
    "true"
  )
  skip_if(
    !run_online,
    paste(
      "decoupleR resources may use remote services; set",
      "BULKMAE_RUN_ONLINE_TESTS=true to opt in."
    )
  )
  skip_if_not_installed("decoupleR", minimum_version = "2.0.0")

  resource <- activity_resource("PROGENy", organism = "human")
  progeny_network <- decoupleR::get_progeny(organism = "human", top = 50L)
  tf_network <- decoupleR::get_collectri(
    organism = "human",
    split_complexes = FALSE
  )
  network_features <- unique(c(
    progeny_network$target,
    tf_network$target
  ))
  network_features <- utils::head(network_features, 3000L)

  mae <- make_simulated_mae(
    n_features = max(1200L, length(network_features)),
    seed = 1701L
  )
  experiment <- mae_pull_experiment(mae, "rna")
  remaining <- setdiff(rownames(experiment), network_features)
  rownames(experiment) <- c(network_features, remaining)[seq_len(nrow(experiment))]
  online_mae <- mae_create(
    list(rna = experiment),
    mae_samples(mae, "rna")
  )

  progeny <- activity_progeny(
    online_mae,
    "rna",
    assay = "log_expression",
    organism = "human",
    top = 50L
  )
  tf <- activity_tf(
    online_mae,
    "rna",
    assay = "log_expression",
    organism = "human",
    resource = "collectri"
  )

  expect_s3_class(resource, "data.frame")
  expect_s3_class(progeny, "data.frame")
  expect_s3_class(tf, "data.frame")
  expect_gt(nrow(resource), 0L)
  expect_gt(nrow(progeny), 0L)
  expect_gt(nrow(tf), 0L)
})
