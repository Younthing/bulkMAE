#' Filter lowly expressed features with edgeR
#'
#' @inheritParams mae_pull_assay
#' @param group Optional metadata column or grouping vector.
#' @param design Optional design matrix. Supply either `group` or `design`.
#' @param min_count Minimum count passed to [edgeR::filterByExpr()].
#' @param min_total_count Minimum total count passed to
#'   [edgeR::filterByExpr()].
#' @param ... Additional arguments passed to [edgeR::filterByExpr()].
#'
#' @return A named logical vector, one value per feature.
#' @export
filter_expr <- function(
    x,
    experiment,
    assay = "counts",
    group = NULL,
    design = NULL,
    min_count = 10,
    min_total_count = 15,
    ...
) {
  .require_backend("edgeR", "to filter lowly expressed features")
  y <- .make_dge_list(x, experiment, assay)
  data <- mae_samples(x, experiment)

  if (!is.null(group) && !is.null(design)) {
    stop("Supply only one of `group` and `design`.", call. = FALSE)
  }
  if (!is.null(group)) {
    group <- .column_or_vector(group, data, "group")
  }
  if (!is.null(design)) {
    design <- .validate_design_matrix(design, data, "design")
  }

  keep <- edgeR::filterByExpr(
    y,
    group = group,
    design = design,
    min.count = min_count,
    min.total.count = min_total_count,
    ...
  )
  stats::setNames(as.logical(keep), rownames(y))
}

#' Calculate sample-level library metrics
#'
#' @inheritParams mae_pull_assay
#' @param detected_above A feature is detected when its value exceeds this
#'   threshold.
#'
#' @return A data frame with one row per sample.
#' @export
qc_library <- function(
    x,
    experiment,
    assay = "counts",
    detected_above = 0
) {
  matrix <- .pull_matrix(x, experiment, assay)
  if (any(!is.finite(matrix))) {
    stop("The selected assay contains non-finite values.", call. = FALSE)
  }
  if (
    !is.numeric(detected_above) || length(detected_above) != 1L ||
      is.na(detected_above) || !is.finite(detected_above)
  ) {
    stop("`detected_above` must be one finite number.", call. = FALSE)
  }

  data.frame(
    sample = colnames(matrix),
    library_size = colSums(matrix),
    detected_features = colSums(matrix > detected_above),
    zero_fraction = colMeans(matrix == 0),
    row.names = NULL,
    check.names = FALSE
  )
}

#' Calculate sample-to-sample correlations
#'
#' @inheritParams mae_pull_assay
#' @param method Correlation method passed to [stats::cor()].
#' @param use Missing-value policy passed to [stats::cor()].
#' @param features Optional feature subset.
#'
#' @return A sample-by-sample correlation matrix.
#' @export
qc_correlation <- function(
    x,
    experiment,
    assay,
    method = "pearson",
    use = "pairwise.complete.obs",
    features = NULL
) {
  matrix <- .select_features(.pull_matrix(x, experiment, assay), features)
  if (nrow(matrix) < 2L) {
    stop("At least two features are required for sample correlations.", call. = FALSE)
  }
  stats::cor(matrix, method = method, use = use)
}

#' Run sample principal component analysis
#'
#' @inheritParams qc_correlation
#' @param top_n Optional number of most variable features to retain.
#' @param center,scale. Passed to [stats::prcomp()].
#'
#' @return A native `prcomp` object.
#' @export
reduce_pca <- function(
    x,
    experiment,
    assay,
    features = NULL,
    top_n = NULL,
    center = TRUE,
    scale. = FALSE
) {
  matrix <- .select_features(.pull_matrix(x, experiment, assay), features)
  if (any(!is.finite(matrix))) {
    stop("PCA requires a finite assay matrix.", call. = FALSE)
  }
  if (ncol(matrix) < 2L) {
    stop("PCA requires at least two samples.", call. = FALSE)
  }
  matrix <- .top_variable_features(matrix, top_n)
  feature_sd <- apply(matrix, 1L, stats::sd, na.rm = TRUE)
  variable <- is.finite(feature_sd) & feature_sd > 0
  matrix <- matrix[variable, , drop = FALSE]
  if (!nrow(matrix)) {
    stop("No variable features remain for PCA.", call. = FALSE)
  }

  stats::prcomp(t(matrix), center = center, scale. = scale.)
}

#' Run sample multidimensional scaling
#'
#' Uses [limma::plotMDS()] with `plot = FALSE` and returns the native MDS
#' object without opening a graphics device.
#'
#' @inheritParams qc_correlation
#' @param top Number of leading features used by `plotMDS`.
#' @param gene_selection Feature-selection rule accepted by `plotMDS`.
#' @param ... Additional arguments passed to [limma::plotMDS()].
#'
#' @return A native limma MDS object.
#' @export
reduce_mds <- function(
    x,
    experiment,
    assay,
    features = NULL,
    top = 500,
    gene_selection = "pairwise",
    ...
) {
  .require_backend("limma", "to calculate multidimensional scaling")
  matrix <- .select_features(.pull_matrix(x, experiment, assay), features)
  if (any(!is.finite(matrix))) {
    stop("MDS requires a finite assay matrix.", call. = FALSE)
  }
  if (nrow(matrix) < 2L || ncol(matrix) < 2L) {
    stop("MDS requires at least two features and two samples.", call. = FALSE)
  }
  limma::plotMDS(
    matrix,
    top = top,
    gene.selection = gene_selection,
    plot = FALSE,
    ...
  )
}

#' Run UMAP on samples
#'
#' @inheritParams reduce_pca
#' @param neighbors Number of nearest neighbors. When `NULL`, uses the smaller
#'   of 15 and the number of samples minus one.
#' @param components Embedding dimensions.
#' @param metric Distance metric.
#' @param min_distance Minimum embedding distance.
#' @param seed Optional random seed.
#' @param ... Additional arguments passed to [uwot::umap()].
#'
#' @return The native sample-by-component matrix returned by uwot.
#' @export
reduce_umap <- function(
    x,
    experiment,
    assay,
    features = NULL,
    top_n = NULL,
    neighbors = NULL,
    components = 2,
    metric = "euclidean",
    min_distance = 0.01,
    seed = 1,
    ...
) {
  .require_backend("uwot", "to calculate a UMAP embedding")
  matrix <- .select_features(.pull_matrix(x, experiment, assay), features)
  if (any(!is.finite(matrix))) {
    stop("UMAP requires a finite assay matrix.", call. = FALSE)
  }
  matrix <- .top_variable_features(matrix, top_n)
  sample_count <- ncol(matrix)

  if (sample_count < 3L) {
    stop("UMAP requires at least three samples.", call. = FALSE)
  }
  matrix <- .variable_embedding_features(matrix, "UMAP")
  if (is.null(neighbors)) {
    neighbors <- min(15L, sample_count - 1L)
  }
  if (
    !is.numeric(neighbors) || length(neighbors) != 1L || is.na(neighbors) ||
      neighbors != trunc(neighbors) || neighbors < 2L ||
      neighbors >= sample_count
  ) {
    stop(
      "`neighbors` must be an integer from 2 to the number of samples minus one.",
      call. = FALSE
    )
  }
  if (
    !is.numeric(components) || length(components) != 1L || is.na(components) ||
      components != trunc(components) || components < 1L ||
      components >= sample_count
  ) {
    stop(
      "`components` must be a positive integer smaller than the sample count.",
      call. = FALSE
    )
  }
  .assert_scalar_character(metric, "metric")
  if (
    !is.numeric(min_distance) || length(min_distance) != 1L ||
      is.na(min_distance) || !is.finite(min_distance) || min_distance < 0
  ) {
    stop("`min_distance` must be one finite non-negative number.", call. = FALSE)
  }

  embedding <- .with_seed(
    seed,
    uwot::umap(
      t(matrix),
      n_neighbors = as.integer(neighbors),
      n_components = as.integer(components),
      metric = metric,
      min_dist = min_distance,
      ...
    )
  )
  if (is.matrix(embedding)) {
    rownames(embedding) <- colnames(matrix)
  } else if (is.list(embedding) && is.matrix(embedding$embedding)) {
    rownames(embedding$embedding) <- colnames(matrix)
  }
  embedding
}

#' Run t-SNE on samples
#'
#' @inheritParams reduce_pca
#' @param dimensions Embedding dimensions.
#' @param perplexity Perplexity passed to [Rtsne::Rtsne()]. When `NULL`, uses
#'   the largest safe integer no greater than 30.
#' @param seed Optional random seed.
#' @param ... Additional arguments passed to [Rtsne::Rtsne()].
#'
#' @return A native `Rtsne` result object.
#' @export
reduce_tsne <- function(
    x,
    experiment,
    assay,
    features = NULL,
    top_n = NULL,
    dimensions = 2,
    perplexity = NULL,
    seed = 1,
    ...
) {
  .require_backend("Rtsne", "to calculate a t-SNE embedding")
  matrix <- .select_features(.pull_matrix(x, experiment, assay), features)
  if (any(!is.finite(matrix))) {
    stop("t-SNE requires a finite assay matrix.", call. = FALSE)
  }
  matrix <- .top_variable_features(matrix, top_n)
  sample_count <- ncol(matrix)

  if (sample_count < 5L) {
    stop("t-SNE requires at least five samples.", call. = FALSE)
  }
  matrix <- .variable_embedding_features(matrix, "t-SNE")
  if (is.null(perplexity)) {
    perplexity <- min(30L, floor((sample_count - 2L) / 3L))
  }
  if (
    !is.numeric(perplexity) || length(perplexity) != 1L ||
      is.na(perplexity) || !is.finite(perplexity) || perplexity <= 0 ||
      3 * perplexity >= sample_count - 1L
  ) {
    stop(
      "`perplexity` must be positive and satisfy 3 * perplexity < samples - 1.",
      call. = FALSE
    )
  }
  if (
    !is.numeric(dimensions) || length(dimensions) != 1L ||
      is.na(dimensions) || dimensions != trunc(dimensions) ||
      dimensions < 1L || dimensions >= sample_count
  ) {
    stop(
      "`dimensions` must be a positive integer smaller than the sample count.",
      call. = FALSE
    )
  }

  result <- .with_seed(
    seed,
    Rtsne::Rtsne(
      t(matrix),
      dims = as.integer(dimensions),
      perplexity = perplexity,
      ...
    )
  )
  rownames(result$Y) <- colnames(matrix)
  result
}

#' Flag samples with unusually large robust PCA distances
#'
#' This is a transparent screening heuristic, not a test that automatically
#' excludes samples. Distances are calculated from median/MAD-standardized PC
#' scores and compared with a chi-squared quantile.
#'
#' @param pca A `prcomp` object, usually from [reduce_pca()].
#' @param components Principal components to include.
#' @param probability Chi-squared cutoff probability.
#'
#' @return A data frame containing robust distance, cutoff, and flag.
#' @export
qc_outliers <- function(
    pca,
    components = seq_len(min(5L, ncol(pca$x))),
    probability = 0.99
) {
  if (!inherits(pca, "prcomp")) {
    stop("`pca` must be a prcomp object.", call. = FALSE)
  }
  if (!length(components) || any(!components %in% seq_len(ncol(pca$x)))) {
    stop("`components` contains unavailable principal components.", call. = FALSE)
  }
  if (length(probability) != 1L || probability <= 0 || probability >= 1) {
    stop("`probability` must be between zero and one.", call. = FALSE)
  }

  scores <- pca$x[, components, drop = FALSE]
  centers <- apply(scores, 2L, stats::median, na.rm = TRUE)
  scales <- apply(scores, 2L, stats::mad, na.rm = TRUE)
  scales[!is.finite(scales) | scales == 0] <- 1
  standardized <- sweep(scores, 2L, centers, "-")
  standardized <- sweep(standardized, 2L, scales, "/")
  distance <- rowSums(standardized^2)
  cutoff <- stats::qchisq(probability, df = length(components))

  data.frame(
    sample = rownames(scores),
    robust_distance = distance,
    cutoff = cutoff,
    flagged = distance > cutoff,
    row.names = NULL,
    check.names = FALSE
  )
}

.top_variable_features <- function(matrix, top_n) {
  if (is.null(top_n)) {
    return(matrix)
  }
  if (
    !is.numeric(top_n) || length(top_n) != 1L || is.na(top_n) ||
      !is.finite(top_n) || top_n != trunc(top_n) || top_n < 1
  ) {
    stop("`top_n` must be a positive integer.", call. = FALSE)
  }
  variance <- apply(matrix, 1L, stats::var, na.rm = TRUE)
  keep <- order(variance, decreasing = TRUE)[
    seq_len(as.integer(min(top_n, nrow(matrix))))
  ]
  matrix[keep, , drop = FALSE]
}

.variable_embedding_features <- function(matrix, method) {
  variance <- apply(matrix, 1L, stats::var)
  keep <- is.finite(variance) & variance > 0
  matrix <- matrix[keep, , drop = FALSE]
  if (nrow(matrix) < 2L) {
    stop(
      method, " requires at least two features with non-zero variance.",
      call. = FALSE
    )
  }
  matrix
}
