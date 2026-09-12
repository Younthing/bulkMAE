#' Simulate a self-contained bulk-transcriptomics MAE
#'
#' This function is intended for examples, teaching, backend smoke tests, and
#' reproducible bug reports. Its synthetic annotations and outcomes have no
#' biological meaning and must not be used as reference data for real studies.
#'
#' @param n_features Number of simulated expression features.
#' @param n_samples Even number of simulated samples.
#' @param seed Local random seed; the caller's random-number state is restored.
#' @param experiment Experiment name in the returned MAE.
#'
#' @return A [MultiAssayExperiment::MultiAssayExperiment] containing integer
#'   `counts`, `log_expression`, and `tpm` assays; feature lengths and synthetic
#'   identifiers in row metadata; and common design, prediction, and survival
#'   variables in sample metadata.
#'
#' @examples
#' example_mae <- mae_simulate(n_features = 100, n_samples = 8, seed = 1)
#' mae_assays(example_mae, "rna")
#' @export
mae_simulate <- function(
    n_features = 300L,
    n_samples = 12L,
    seed = 1L,
    experiment = "rna"
) {
  if (
    !is.numeric(n_features) || length(n_features) != 1L ||
      !is.finite(n_features) || n_features != trunc(n_features) ||
      n_features < 20L || n_features > .Machine$integer.max
  ) {
    stop("`n_features` must be an integer of at least 20.", call. = FALSE)
  }
  if (
    !is.numeric(n_samples) || length(n_samples) != 1L ||
      !is.finite(n_samples) || n_samples != trunc(n_samples) ||
      n_samples < 8L || n_samples %% 2L != 0L ||
      n_samples > .Machine$integer.max
  ) {
    stop("`n_samples` must be an even integer of at least 8.", call. = FALSE)
  }
  if (
    !is.numeric(seed) || length(seed) != 1L || is.na(seed) ||
      !is.finite(seed) || seed != trunc(seed) || seed < 0 ||
      seed > .Machine$integer.max
  ) {
    stop("`seed` must be one non-negative integer.", call. = FALSE)
  }
  .assert_scalar_character(experiment, "experiment")
  n_features <- as.integer(n_features)
  n_samples <- as.integer(n_samples)

  .with_seed(seed, {
    feature_ids <- sprintf("gene%04d", seq_len(n_features))
    sample_ids <- sprintf("sample%02d", seq_len(n_samples))
    condition <- factor(
      rep(c("control", "treated"), length.out = n_samples),
      levels = c("control", "treated")
    )
    subject_id <- factor(rep(
      sprintf("subject%02d", seq_len(n_samples / 2L)),
      each = 2L
    ))
    batch <- factor(rep(c("batch_a", "batch_a", "batch_b", "batch_b"), length.out = n_samples))
    timepoint <- rep(c(0, 1, 2), length.out = n_samples)
    outcome <- as.integer(condition == "treated")
    risk_score <- as.numeric(scale(0.7 * outcome + stats::rnorm(n_samples)))
    event <- rep(c(1L, 0L, 1L, 1L), length.out = n_samples)
    survival_time <- round(
      pmax(
        5,
        stats::rexp(n_samples, rate = 1 / 500) * exp(-0.3 * risk_score)
      ),
      1L
    )
    sample_data <- data.frame(
      condition = condition,
      batch = batch,
      subject_id = subject_id,
      timepoint = timepoint,
      replicate = subject_id,
      survival_time = survival_time,
      event = event,
      outcome = outcome,
      risk_score = risk_score,
      age = round(stats::runif(n_samples, 40, 80)),
      row.names = sample_ids,
      check.names = FALSE
    )

    baseline <- exp(stats::rnorm(n_features, log(80), 0.6))
    effect <- numeric(n_features)
    effect[seq_len(min(40L, n_features))] <- log(2)
    batch_effect <- numeric(n_features)
    batch_rows <- if (n_features >= 41L) {
      seq.int(41L, min(70L, n_features))
    } else {
      integer()
    }
    batch_effect[batch_rows] <- log(1.3)
    time_effect <- numeric(n_features)
    time_rows <- if (n_features >= 71L) {
      seq.int(71L, min(90L, n_features))
    } else {
      integer()
    }
    time_effect[time_rows] <- log(4)
    interaction_effect <- numeric(n_features)
    interaction_rows <- if (n_features >= 91L) {
      seq.int(91L, min(110L, n_features))
    } else {
      integer()
    }
    interaction_effect[interaction_rows] <- log(4)
    mean_matrix <- exp(
      outer(log(baseline), rep(1, n_samples)) +
        outer(effect, outcome) +
        outer(batch_effect, as.integer(batch == "batch_b")) +
        outer(time_effect, timepoint) +
        outer(interaction_effect, outcome * timepoint)
    )
    counts <- matrix(
      stats::rnbinom(length(mean_matrix), mu = mean_matrix, size = 30),
      nrow = n_features,
      dimnames = list(feature_ids, sample_ids)
    )
    storage.mode(counts) <- "integer"
    library_sizes <- colSums(counts)
    log_expression <- log2(sweep(counts, 2L, library_sizes, "/") * 1e6 + 0.5)
    gene_length <- stats::runif(n_features, 500, 3000)
    per_kb <- counts / (gene_length / 1000)
    tpm <- sweep(per_kb, 2L, colSums(per_kb), "/") * 1e6
    row_data <- data.frame(
      gene_length = gene_length,
      gene_symbol = feature_ids,
      entrez_id = as.character(100000L + seq_len(n_features)),
      gene_id = sprintf("parent%04d", ceiling(seq_len(n_features) / 2)),
      row.names = feature_ids,
      check.names = FALSE
    )
    leaf <- mae_create_experiment(
      assays = list(
        counts = counts,
        log_expression = log_expression,
        tpm = tpm
      ),
      col_data = data.frame(row.names = sample_ids),
      row_data = row_data
    )
    mae_create(stats::setNames(list(leaf), experiment), sample_data)
  })
}
