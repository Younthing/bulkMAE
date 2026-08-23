make_simulated_mae <- function(
    n_features = 1200L,
    n_samples = 24L,
    seed = 20260824L
) {
  if (n_features < 1200L) {
    stop("`n_features` must be at least 1200 for the shared fixture.")
  }
  if (n_samples != 24L) {
    stop("The shared fixture currently requires exactly 24 samples.")
  }

  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  old_kind <- RNGkind()
  if (had_seed) {
    old_seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  }
  on.exit({
    do.call(RNGkind, as.list(old_kind))
    if (had_seed) {
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)

  set.seed(seed)
  feature_names <- sprintf("gene%04d", seq_len(n_features))
  sample_names <- sprintf("sample%02d", seq_len(n_samples))

  subject_id <- rep(sprintf("subject%02d", seq_len(12L)), each = 2L)
  condition <- rep(rep(c("control", "treated"), each = 6L), each = 2L)
  batch <- rep(c("batch_a", "batch_b"), times = 12L)
  outcome <- as.integer(condition == "treated")
  risk_score <- as.numeric(scale(0.65 * outcome + stats::rnorm(n_samples)))
  age <- as.integer(round(stats::runif(n_samples, min = 42, max = 78)))
  event_template <- c(rep(1L, 8L), rep(0L, 4L))
  event <- c(sample(event_template), sample(event_template))
  follow_up <- stats::rexp(n_samples, rate = 1 / 650)
  follow_up <- follow_up * exp(-0.25 * risk_score)
  follow_up <- round(pmax(follow_up, 14), digits = 1L)

  sample_data <- data.frame(
    condition = condition,
    batch = batch,
    subject_id = subject_id,
    time = follow_up,
    event = event,
    risk_score = risk_score,
    outcome = outcome,
    age = age,
    row.names = sample_names,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  baseline <- exp(stats::rnorm(n_features, mean = log(80), sd = 0.65))
  condition_effect <- numeric(n_features)
  condition_effect[seq_len(80L)] <- log(2)
  batch_effect <- numeric(n_features)
  batch_effect[81:140] <- log(1.35)
  library_factor <- exp(stats::rnorm(n_samples, mean = 0, sd = 0.18))

  log_mean <- outer(log(baseline), rep(1, n_samples))
  log_mean <- log_mean + outer(condition_effect, outcome)
  log_mean <- log_mean + outer(batch_effect, as.integer(batch == "batch_b"))
  log_mean <- sweep(log_mean, 2L, log(library_factor), FUN = "+")
  mean_matrix <- exp(log_mean)
  counts <- matrix(
    stats::rnbinom(length(mean_matrix), mu = as.vector(mean_matrix), size = 8),
    nrow = n_features,
    dimnames = list(feature_names, sample_names)
  )
  storage.mode(counts) <- "integer"

  counts_per_million <- sweep(counts, 2L, colSums(counts), FUN = "/") * 1e6
  log_expression <- log2(counts_per_million + 0.5)
  gene_length <- stats::runif(n_features, min = 500, max = 3000)
  reads_per_kb <- counts / (gene_length / 1000)
  tpm <- sweep(reads_per_kb, 2L, colSums(reads_per_kb), FUN = "/") * 1e6
  weights <- 1 / sqrt(log2(counts + 2))

  row_data <- data.frame(
    gene_length = gene_length,
    row.names = feature_names,
    check.names = FALSE
  )
  assay_col_data <- data.frame(row.names = sample_names)
  experiment <- mae_create_experiment(
    assays = list(
      counts = counts,
      log_expression = log_expression,
      tpm = tpm,
      weights = weights
    ),
    col_data = assay_col_data,
    row_data = row_data
  )
  mae_create(list(rna = experiment), sample_data)
}
