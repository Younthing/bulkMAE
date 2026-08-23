test_that("mae_create_experiment validates every assay", {
  samples <- data.frame(group = c("a", "b"), row.names = c("s1", "s2"))
  first <- matrix(
    1:4,
    nrow = 2,
    dimnames = list(c("g1", "g2"), c("s1", "s2"))
  )
  reordered <- first[, c("s2", "s1"), drop = FALSE]

  expect_error(mae_create_experiment(list(), samples), "at least one assay")
  expect_error(
    mae_create_experiment(
      stats::setNames(list(first, first), c("counts", "counts")),
      samples
    ),
    "unique, non-empty names"
  )
  expect_error(
    mae_create_experiment(list(first = first, second = reordered), samples),
    "identical dimensions and dimnames"
  )

  duplicated <- first
  rownames(duplicated) <- c("g1", "g1")
  expect_error(
    mae_create_experiment(duplicated, samples),
    "feature names must be non-empty and unique"
  )
})

test_that("named sample vectors and designs are aligned by identifier", {
  samples <- data.frame(
    group = c("a", "b", "a", "b"),
    row.names = paste0("s", 1:4)
  )
  value <- c(s4 = 4, s2 = 2, s1 = 1, s3 = 3)
  aligned <- bulkMAE:::.column_or_vector(value, samples, "value")

  expect_identical(unname(aligned), c(1, 2, 3, 4))
  expect_identical(names(aligned), rownames(samples))

  design <- stats::model.matrix(~ group, samples[c(4, 2, 1, 3), , drop = FALSE])
  checked <- bulkMAE:::.validate_design_matrix(design, samples)
  expect_identical(rownames(checked), rownames(samples))
})

test_that("model construction rejects missing and aliased covariates", {
  missing <- data.frame(
    group = c("a", NA, "b"),
    row.names = paste0("s", 1:3)
  )
  expect_error(
    bulkMAE:::.model_matrix(~ group, missing),
    "Could not construct the model frame"
  )

  aliased <- data.frame(
    first = c(0, 0, 1, 1),
    second = c(0, 0, 1, 1),
    row.names = paste0("s", 1:4)
  )
  expect_error(
    bulkMAE:::.model_matrix(~ first + second, aliased),
    "not full rank"
  )
})

test_that("local random seeds leave caller state unchanged", {
  set.seed(912)
  state <- .Random.seed
  first <- bulkMAE:::.with_seed(17, stats::runif(4))
  second <- bulkMAE:::.with_seed(17, stats::runif(4))

  expect_identical(first, second)
  expect_identical(.Random.seed, state)
})

test_that("local random seeds do not create caller state", {
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had_seed) {
    old_seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  }
  on.exit({
    if (had_seed) {
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)
  if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
    rm(".Random.seed", envir = .GlobalEnv)
  }

  bulkMAE:::.with_seed(19, stats::runif(2))
  expect_false(exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
})

test_that("PCA checks finite values and removes constant features", {
  mae <- make_toy_mae(n_features = 4L, n_samples = 6L)
  fit <- reduce_pca(mae, "rna", "log_expression")
  expect_s3_class(fit, "prcomp")
  expect_identical(rownames(fit$x), paste0("sample", 1:6))

  matrix <- mae_pull_assay(mae, "rna", "log_expression")
  matrix[,] <- 1
  samples <- mae_samples(mae, "rna")
  constant <- mae_create(
    list(rna = mae_create_experiment(matrix, samples, assay_name = "log")),
    samples
  )
  expect_error(reduce_pca(constant, "rna", "log"), "No variable features")
})

test_that("embedding wrappers expose sample-safe defaults", {
  expect_null(formals(reduce_umap)$neighbors)
  expect_null(formals(reduce_tsne)$perplexity)
})

test_that("decouple networks make signed weights explicit", {
  network <- data.frame(
    regulator = rep(c("r1", "r2"), each = 3),
    gene = paste0("g", 1:6),
    weight = c(-1, 1, 0.5, -0.5, 1, 1),
    stringsAsFactors = FALSE
  )
  normalized <- bulkMAE:::.prepare_decouple_network(
    network,
    features = paste0("g", 1:6),
    source = "regulator",
    target = "gene",
    mor = "weight",
    min_size = 2L
  )

  expect_identical(normalized$mor, network$weight)
  expect_error(
    bulkMAE:::.assert_decouple_mor(network, "mlm"),
    "require a normalized `mor`"
  )
  expect_error(
    bulkMAE:::.validate_decouple_consensus(
      c("ora", "ulm"),
      consensus = TRUE,
      consensus_statistics = NULL
    ),
    "Do not form one consensus"
  )
  expect_silent(
    bulkMAE:::.validate_decouple_consensus(
      "ulm",
      consensus = TRUE,
      consensus_statistics = "norm_ulm"
    )
  )
})

test_that("reference labels align by cell identifier", {
  labels <- c(cell3 = "B", cell1 = "A", cell2 = "A")
  aligned <- bulkMAE:::.align_reference_labels(
    labels,
    c("cell1", "cell2", "cell3"),
    "labels"
  )

  expect_identical(aligned, c("A", "A", "B"))
  expect_error(
    bulkMAE:::.align_reference_labels(
      labels[-1],
      c("cell1", "cell2", "cell3"),
      "labels"
    ),
    "one value per reference cell"
  )
})

test_that("clinical metadata arguments name exactly one column", {
  data <- data.frame(time = 1:3, event = c(0, 1, 1))
  expect_silent(bulkMAE:::.assert_metadata_column(data, "time", "time"))
  expect_error(
    bulkMAE:::.assert_metadata_column(data, c("time", "event"), "time"),
    "exactly one metadata column"
  )
})
