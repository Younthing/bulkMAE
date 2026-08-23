make_toy_mae <- function(n_features = 3L, n_samples = 4L) {
  counts <- matrix(
    seq_len(n_features * n_samples),
    nrow = n_features,
    dimnames = list(
      paste0("gene", seq_len(n_features)),
      paste0("sample", seq_len(n_samples))
    )
  )
  samples <- data.frame(
    condition = rep(c("a", "b"), length.out = n_samples),
    batch = rep(c("x", "y"), length.out = n_samples),
    row.names = colnames(counts)
  )
  se <- mae_create_experiment(
    assays = list(counts = counts, log_expression = log2(counts + 1)),
    col_data = samples
  )
  mae_create(list(rna = se), samples)
}

make_model_mae <- function(n_features = 100L, n_samples = 8L, seed = 1L) {
  if (n_samples %% 2L != 0L) {
    stop("`n_samples` must be even for the two-group test fixture.")
  }
  set.seed(seed)
  counts <- matrix(
    stats::rnbinom(n_features * n_samples, mu = 100, size = 5),
    nrow = n_features,
    dimnames = list(
      paste0("gene", seq_len(n_features)),
      paste0("sample", seq_len(n_samples))
    )
  )
  samples <- data.frame(
    condition = rep(c("control", "treated"), each = n_samples / 2L),
    batch = rep(c("a", "b"), length.out = n_samples),
    row.names = colnames(counts)
  )
  experiment <- mae_create_experiment(
    assays = list(counts = counts, log_expression = log2(counts + 1)),
    col_data = samples
  )
  mae_create(list(rna = experiment), samples)
}
