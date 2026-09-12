#' Run consensus clustering
#'
#' @inheritParams mae_pull_assay
#' @param max_k Maximum number of clusters.
#' @param repetitions Number of resampling repetitions.
#' @param item_fraction,feature_fraction Sampling fractions.
#' @param algorithm Clustering algorithm.
#' @param distance Distance metric.
#' @param seed Random seed.
#' @param title Output title or directory used by the backend.
#' @param plot Plot format accepted by ConsensusClusterPlus, or `NULL`.
#' @param ... Additional arguments passed to
#'   [ConsensusClusterPlus::ConsensusClusterPlus()].
#'
#' @return The native ConsensusClusterPlus result list.
#' @export
cluster_consensus <- function(
    x,
    experiment,
    assay,
    max_k = 6,
    repetitions = 1000,
    item_fraction = 0.8,
    feature_fraction = 1,
    algorithm = "hc",
    distance = "pearson",
    seed = 1,
    title = tempdir(),
    plot = NULL,
    ...
) {
  .require_backend("ConsensusClusterPlus", "to run consensus clustering")
  .assert_scalar_character(algorithm, "algorithm")
  .assert_scalar_character(distance, "distance")
  .assert_scalar_character(title, "title")
  if (
    !is.null(plot) &&
      (!is.character(plot) || length(plot) != 1L || is.na(plot) || !nzchar(plot))
  ) {
    stop("`plot` must be NULL or one non-empty string.", call. = FALSE)
  }
  matrix <- .pull_matrix(x, experiment, assay)
  if (nrow(matrix) < 2L || ncol(matrix) < 3L) {
    stop("Consensus clustering needs at least two features and three samples.", call. = FALSE)
  }
  if (any(!is.finite(matrix))) {
    stop("Consensus clustering requires finite assay values.", call. = FALSE)
  }
  feature_variance <- apply(matrix, 1L, stats::var)
  variable <- is.finite(feature_variance) & feature_variance > 0
  removed_features <- rownames(matrix)[!variable]
  if (length(removed_features)) {
    warning(
      "Consensus clustering removed ", length(removed_features),
      " non-variable feature(s).",
      call. = FALSE
    )
    matrix <- matrix[variable, , drop = FALSE]
  }
  if (nrow(matrix) < 2L) {
    stop("Consensus clustering needs at least two variable features.", call. = FALSE)
  }
  .network_assert_whole_number(max_k, "max_k", lower = 2L)
  .network_assert_whole_number(repetitions, "repetitions", lower = 1L)
  if (max_k > ncol(matrix)) {
    stop("`max_k` cannot exceed the number of samples.", call. = FALSE)
  }
  .network_assert_fraction(item_fraction, "item_fraction")
  .network_assert_fraction(feature_fraction, "feature_fraction")
  if (ceiling(ncol(matrix) * item_fraction) < max_k) {
    stop(
      "`item_fraction` samples too few items to fit `max_k` clusters.",
      call. = FALSE
    )
  }
  if (ceiling(nrow(matrix) * feature_fraction) < 1L) {
    stop("`feature_fraction` must retain at least one feature.", call. = FALSE)
  }
  .network_assert_seed(seed)

  arguments <- c(list(
    d = matrix,
    maxK = as.integer(max_k),
    reps = as.integer(repetitions),
    pItem = item_fraction,
    pFeature = feature_fraction,
    clusterAlg = algorithm,
    distance = distance,
    seed = seed,
    title = title,
    plot = plot
  ), list(...))
  result <- .with_seed(
    seed,
    .consensus_cluster_call(arguments, suppress_plot = is.null(plot))
  )
  attr(result, "bulkMAERemovedFeatures") <- removed_features
  result
}

#' Calculate consensus-clustering stability diagnostics
#'
#' @param results Result list returned by [cluster_consensus()].
#' @param title Output title or directory used by ConsensusClusterPlus.
#' @param plot Plot format accepted by [ConsensusClusterPlus::calcICL()], or
#'   `NULL` to suppress plots.
#' @param ... Additional arguments passed to [ConsensusClusterPlus::calcICL()].
#'
#' @return The native list containing cluster and item consensus tables.
#' @export
cluster_consensus_diagnostics <- function(
    results,
    title = tempdir(),
    plot = NULL,
    ...
) {
  .require_backend("ConsensusClusterPlus", "to calculate consensus diagnostics")
  if (!is.list(results) || length(results) < 2L) {
    stop("`results` must be a ConsensusClusterPlus result list.", call. = FALSE)
  }
  .assert_scalar_character(title, "title")
  if (
    !is.null(plot) &&
      (!is.character(plot) || length(plot) != 1L || is.na(plot) || !nzchar(plot))
  ) {
    stop("`plot` must be NULL or one non-empty string.", call. = FALSE)
  }
  if (is.null(plot)) {
    grDevices::pdf(file = NULL)
    on.exit(grDevices::dev.off(), add = TRUE)
  }
  ConsensusClusterPlus::calcICL(results, title = title, plot = plot, ...)
}

#' Extract sample classes from a consensus-clustering solution
#'
#' @param results Result list returned by [cluster_consensus()].
#' @param k Number of clusters to extract.
#'
#' @return A named integer vector of sample classes.
#' @export
cluster_consensus_classes <- function(results, k) {
  if (
    !is.list(results) || !is.numeric(k) || length(k) != 1L || is.na(k) ||
      !is.finite(k) ||
      k != as.integer(k) || k < 2L || k > length(results)
  ) {
    stop("`k` must identify one available consensus solution.", call. = FALSE)
  }
  solution <- results[[as.integer(k)]]
  if (!is.list(solution)) {
    stop("The selected consensus solution is malformed.", call. = FALSE)
  }
  classes <- solution$consensusClass
  if (is.null(classes)) {
    stop("The selected result has no `consensusClass` element.", call. = FALSE)
  }
  classes
}

#' Proportion of ambiguous clustering for each k
#'
#' PAC is the fraction of pairwise consensus values that fall in the open
#' interval (`lower`, `upper`). Lower values indicate a more decisive
#' consensus matrix. This is an extractor, not a new clustering method.
#'
#' @param results Result list returned by [cluster_consensus()].
#' @param lower,upper Open interval bounds used to count ambiguous pairs.
#'
#' @return A named numeric vector of PAC values, one per available `k`.
#' @export
cluster_consensus_pac <- function(results, lower = 0.1, upper = 0.9) {
  .network_assert_range(lower, "lower", lower = 0, upper = 1)
  .network_assert_range(upper, "upper", lower = 0, upper = 1)
  if (lower >= upper) {
    stop("`lower` must be less than `upper`.", call. = FALSE)
  }
  pairs <- .cluster_consensus_pair_values(results)
  vapply(
    pairs,
    function(values) mean(values > lower & values < upper),
    numeric(1)
  )
}

#' Consensus CDF area and relative delta area
#'
#' Integrates the empirical CDF of pairwise consensus values at each `k`.
#' `delta` is the CDF area at `k = 2` and the relative increase
#' `(area_k - area_{k-1}) / area_{k-1}` for larger `k`. This is an extractor,
#' not an automatic k-selection rule.
#'
#' @param results Result list returned by [cluster_consensus()].
#'
#' @return A data frame with `k`, `area`, and `delta` columns.
#' @export
cluster_consensus_delta_area <- function(results) {
  pairs <- .cluster_consensus_pair_values(results)
  k <- as.integer(names(pairs))
  area <- vapply(pairs, .cluster_consensus_cdf_area, numeric(1))
  delta <- area
  if (length(area) > 1L) {
    previous <- area[-length(area)]
    relative <- ifelse(
      previous > 0,
      (area[-1L] - previous) / previous,
      NA_real_
    )
    delta[-1L] <- relative
  }
  data.frame(
    k = k,
    area = unname(area),
    delta = unname(delta),
    row.names = NULL
  )
}

.cluster_consensus_pair_values <- function(results) {
  if (!is.list(results) || length(results) < 2L) {
    stop("`results` must be a ConsensusClusterPlus result list.", call. = FALSE)
  }
  out <- list()
  for (index in seq_along(results)) {
    if (index < 2L) {
      next
    }
    solution <- results[[index]]
    if (!is.list(solution) || is.null(solution$consensusMatrix)) {
      next
    }
    matrix <- as.matrix(solution$consensusMatrix)
    values <- matrix[upper.tri(matrix, diag = FALSE)]
    values <- values[is.finite(values)]
    if (!length(values)) {
      next
    }
    if (any(values < 0 | values > 1)) {
      stop("Consensus values must lie in [0, 1].", call. = FALSE)
    }
    out[[as.character(index)]] <- values
  }
  if (!length(out)) {
    stop("No consensus matrices were available for k >= 2.", call. = FALSE)
  }
  out
}

.cluster_consensus_cdf_area <- function(values) {
  values <- sort(as.numeric(values[is.finite(values)]))
  x <- c(0, values, 1)
  y <- c(0, seq_along(values) / length(values), 1)
  keep <- !duplicated(x, fromLast = TRUE)
  x <- x[keep]
  y <- y[keep]
  sum(diff(x) * y[-length(y)])
}

#' Test differential co-expression between two groups
#'
#' Uses the maintained Bioconductor `diffcoexp` implementation to identify
#' differential links and genes from two condition-specific expression
#' matrices. DCA here means differential co-expression analysis.
#'
#' @inheritParams mae_pull_assay
#' @param group Metadata column or grouping vector.
#' @param contrast Two group labels, in the order passed to `exprs.1` and
#'   `exprs.2`.
#' @param features Optional feature subset.
#' @param correlation Correlation method.
#' @param p_adjust Multiple-testing correction method.
#' @param correlation_threshold,correlation_fdr Thresholds for
#'   condition-specific correlations.
#' @param difference_threshold,difference_fdr Thresholds for changes in
#'   correlation.
#' @param gene_fdr FDR threshold for differential co-expression genes.
#' @param max_pairs Maximum number of feature pairs to analyze. Set to `Inf`
#'   only after considering the quadratic memory and runtime cost.
#'
#' @details Use a normalized, variance-stabilized or log-expression assay.
#' Remove known unwanted variation before this analysis when appropriate.
#' Correlations are estimated separately, so small groups are unstable even
#' when they meet the minimum of four samples; groups below ten emit a warning.
#'
#' @return The native result returned by [diffcoexp::diffcoexp()].
#' @export
coexpr_differential <- function(
    x,
    experiment,
    group,
    contrast,
    assay,
    features = NULL,
    correlation = "pearson",
    p_adjust = "BH",
    correlation_threshold = 0.5,
    correlation_fdr = 0.1,
    difference_threshold = 0.5,
    difference_fdr = 0.1,
    gene_fdr = 0.1,
    max_pairs = 5e6
) {
  .require_backend("diffcoexp", "to test differential co-expression")
  data <- mae_samples(x, experiment)
  group <- .column_or_vector(group, data, "group")
  if (anyNA(group)) {
    stop("`group` cannot contain missing values.", call. = FALSE)
  }
  if (
    (!is.atomic(contrast) && !is.factor(contrast)) ||
      length(contrast) != 2L || anyNA(contrast) ||
      length(unique(as.character(contrast))) != 2L ||
      any(!contrast %in% unique(group))
  ) {
    stop("`contrast` must contain two distinct observed group labels.", call. = FALSE)
  }
  correlation <- match.arg(correlation, c("pearson", "kendall", "spearman"))
  p_adjust <- match.arg(p_adjust, stats::p.adjust.methods)
  .network_assert_range(
    correlation_threshold,
    "correlation_threshold",
    lower = 0,
    upper = 1
  )
  .network_assert_range(correlation_fdr, "correlation_fdr", lower = 0, upper = 1)
  .network_assert_range(
    difference_threshold,
    "difference_threshold",
    lower = 0,
    upper = 2
  )
  .network_assert_range(difference_fdr, "difference_fdr", lower = 0, upper = 1)
  .network_assert_range(gene_fdr, "gene_fdr", lower = 0, upper = 1)
  if (
    !is.numeric(max_pairs) || length(max_pairs) != 1L || is.na(max_pairs) ||
      max_pairs <= 0
  ) {
    stop("`max_pairs` must be positive or `Inf`.", call. = FALSE)
  }

  matrix <- .select_features(.pull_matrix(x, experiment, assay), features)
  if (nrow(matrix) < 2L) {
    stop("Differential co-expression requires at least two features.", call. = FALSE)
  }
  pair_count <- nrow(matrix) * (nrow(matrix) - 1) / 2
  if (pair_count > max_pairs) {
    stop(
      "This analysis would test ", format(pair_count, scientific = FALSE),
      " feature pairs, exceeding `max_pairs = ",
      format(max_pairs, scientific = FALSE),
      "`. Subset `features` or explicitly increase `max_pairs`.",
      call. = FALSE
    )
  }

  selected <- group %in% contrast
  matrix <- matrix[, selected, drop = FALSE]
  group <- group[selected]
  sizes <- table(factor(group, levels = contrast))
  if (any(sizes < 4L)) {
    stop("Differential co-expression requires at least four samples per group.", call. = FALSE)
  }
  if (any(sizes < 10L)) {
    warning(
      "Differential co-expression correlations are unstable with fewer than ",
      "ten samples per group.",
      call. = FALSE
    )
  }
  if (any(!is.finite(matrix))) {
    stop("Differential co-expression requires finite assay values.", call. = FALSE)
  }
  first <- matrix[, group == contrast[[1L]], drop = FALSE]
  second <- matrix[, group == contrast[[2L]], drop = FALSE]
  zero_variance <- apply(first, 1L, stats::var) == 0 |
    apply(second, 1L, stats::var) == 0
  if (any(zero_variance)) {
    removed_features <- rownames(matrix)[zero_variance]
    warning(
      "Differential co-expression removed ", length(removed_features),
      " feature(s) with zero within-group variance.",
      call. = FALSE
    )
    first <- first[!zero_variance, , drop = FALSE]
    second <- second[!zero_variance, , drop = FALSE]
  } else {
    removed_features <- character()
  }
  if (nrow(first) < 2L) {
    stop("Differential co-expression needs at least two variable features.", call. = FALSE)
  }

  result <- diffcoexp::diffcoexp(
    exprs.1 = first,
    exprs.2 = second,
    r.method = correlation,
    q.method = p_adjust,
    rth = correlation_threshold,
    qth = correlation_fdr,
    r.diffth = difference_threshold,
    q.diffth = difference_fdr,
    q.dcgth = gene_fdr
  )
  attr(result, "bulkMAERemovedFeatures") <- removed_features
  result
}

#' Run non-negative matrix factorization
#'
#' @inheritParams mae_pull_assay
#' @param rank Factorization rank or vector of candidate ranks.
#' @param method NMF algorithm.
#' @param runs Number of runs.
#' @param seed Random seed or NMF seed method.
#' @param ... Additional arguments passed to [NMF::nmf()].
#'
#' @return A native NMF fit object.
#' @export
cluster_nmf <- function(
    x,
    experiment,
    assay,
    rank,
    method = "brunet",
    runs = 30,
    seed = 1,
    ...
) {
  .require_backend("NMF", "to run non-negative matrix factorization")
  matrix <- .pull_matrix(x, experiment, assay)
  if (any(!is.finite(matrix)) || any(matrix < 0)) {
    stop("NMF requires finite, non-negative values.", call. = FALSE)
  }
  if (
    !is.numeric(rank) || !length(rank) || anyNA(rank) ||
      any(!is.finite(rank)) || any(rank != trunc(rank)) || any(rank < 1L) ||
      any(rank > min(dim(matrix)))
  ) {
    stop(
      "`rank` must contain positive integers no larger than the smaller ",
      "assay dimension.",
      call. = FALSE
    )
  }
  .network_assert_whole_number(runs, "runs", lower = 1L)
  .assert_scalar_character(method, "method")
  arguments <- c(list(
    x = matrix,
    rank = rank,
    method = method,
    nrun = runs,
    seed = seed
  ), list(...))
  if (is.numeric(seed)) {
    .network_assert_seed(seed)
    return(.with_attached_namespaces(
      "NMF",
      .with_seed(seed, do.call(NMF::nmf, arguments))
    ))
  }
  .with_attached_namespaces("NMF", do.call(NMF::nmf, arguments))
}

#' Detect co-expression modules with WGCNA
#'
#' @inheritParams mae_pull_assay
#' @param power Soft-thresholding power chosen for the dataset.
#' @param top_n Optional number of most variable genes.
#' @param network_type WGCNA network type. Use the same value for module
#'   preservation.
#' @param tom_type Optional topological-overlap type. By default it follows
#'   `network_type` (`unsigned` stays unsigned; signed variants use signed
#'   TOM).
#' @param min_module_size Minimum module size.
#' @param merge_cut_height Module merging threshold.
#' @param numeric_labels Return numeric module labels.
#' @param seed Random seed used by WGCNA without changing the caller's random
#'   number generator state.
#' @param verbose Backend verbosity.
#' @param ... Additional arguments passed to [WGCNA::blockwiseModules()].
#'
#' @return The native WGCNA module result list, augmented with
#'   `bulkMAEQuality`, `bulkMAERemovedSamples`, and `bulkMAERemovedFeatures`
#'   elements documenting automatic quality filtering.
#' @export
coexpr_wgcna <- function(
    x,
    experiment,
    assay,
    power,
    top_n = NULL,
    network_type = "signed",
    tom_type = NULL,
    min_module_size = 30,
    merge_cut_height = 0.25,
    numeric_labels = TRUE,
    seed = 1,
    verbose = 2,
    ...
) {
  .require_backend("WGCNA", "to detect co-expression modules")
  network_type <- match.arg(
    network_type,
    c("unsigned", "signed", "signed hybrid")
  )
  if (is.null(tom_type)) {
    tom_type <- if (identical(network_type, "unsigned")) "unsigned" else "signed"
  }
  .assert_scalar_character(tom_type, "tom_type")
  if (
    !is.numeric(power) || length(power) != 1L || is.na(power) ||
      !is.finite(power) || power <= 0
  ) {
    stop("`power` must be one finite positive number.", call. = FALSE)
  }
  .network_assert_whole_number(
    min_module_size,
    "min_module_size",
    lower = 2L
  )
  .network_assert_range(
    merge_cut_height,
    "merge_cut_height",
    lower = 0,
    upper = 1
  )
  if (
    !is.logical(numeric_labels) || length(numeric_labels) != 1L ||
      is.na(numeric_labels)
  ) {
    stop("`numeric_labels` must be TRUE or FALSE.", call. = FALSE)
  }
  .network_assert_seed(seed)
  matrix <- .top_variable_features(
    .pull_matrix(x, experiment, assay),
    top_n
  )
  data_expression <- t(matrix)
  quality <- WGCNA::goodSamplesGenes(data_expression, verbose = 0)
  removed_samples <- rownames(data_expression)[!quality$goodSamples]
  removed_features <- colnames(data_expression)[!quality$goodGenes]
  if (!quality$allOK) {
    warning(
      "WGCNA removed ", length(removed_samples), " sample(s) and ",
      length(removed_features), " feature(s) failing goodSamplesGenes().",
      call. = FALSE
    )
    data_expression <- data_expression[
      quality$goodSamples,
      quality$goodGenes,
      drop = FALSE
    ]
  }
  if (nrow(data_expression) < 4L || ncol(data_expression) < 3L) {
    stop(
      "Too few samples or features remain after WGCNA quality filtering.",
      call. = FALSE
    )
  }
  if (min_module_size > ncol(data_expression)) {
    stop(
      "`min_module_size` cannot exceed the number of retained features.",
      call. = FALSE
    )
  }
  arguments <- c(list(
    datExpr = data_expression,
    power = power,
    networkType = network_type,
    TOMType = tom_type,
    minModuleSize = min_module_size,
    mergeCutHeight = merge_cut_height,
    numericLabels = numeric_labels,
    randomSeed = seed,
    verbose = verbose
  ), list(...))
  result <- .with_seed(
    seed,
    .wgcna_call("blockwiseModules", arguments)
  )
  result$bulkMAEQuality <- quality
  result$bulkMAERemovedSamples <- removed_samples
  result$bulkMAERemovedFeatures <- removed_features
  result$bulkMAENetworkType <- network_type
  result
}

#' Evaluate candidate WGCNA soft-thresholding powers
#'
#' Runs WGCNA's scale-free topology diagnostics directly from an MAE assay and
#' applies the same automatic sample/feature quality filtering used by
#' [coexpr_wgcna()]. The returned fit indices support, but do not replace, the
#' scientific choice of a soft-thresholding power.
#'
#' @inheritParams coexpr_wgcna
#' @param powers Positive candidate powers.
#' @param r_squared Target scale-free topology fit passed as `RsquaredCut`.
#' @param ... Additional arguments passed to [WGCNA::pickSoftThreshold()].
#'
#' @return The native `pickSoftThreshold()` list, augmented with
#'   `bulkMAEQuality`, `bulkMAERemovedSamples`, and
#'   `bulkMAERemovedFeatures` elements.
#' @export
coexpr_pick_power <- function(
    x,
    experiment,
    assay,
    powers = c(seq_len(10L), seq(12L, 20L, by = 2L)),
    top_n = NULL,
    network_type = "signed",
    r_squared = 0.85,
    verbose = 0,
    ...
) {
  .require_backend("WGCNA", "to evaluate soft-thresholding powers")
  if (
    !is.numeric(powers) || !length(powers) || anyNA(powers) ||
      any(!is.finite(powers)) || any(powers <= 0) || anyDuplicated(powers)
  ) {
    stop("`powers` must contain unique, finite positive values.", call. = FALSE)
  }
  network_type <- match.arg(
    network_type,
    c("unsigned", "signed", "signed hybrid")
  )
  .network_assert_range(r_squared, "r_squared", lower = 0, upper = 1)
  matrix <- .top_variable_features(
    .pull_matrix(x, experiment, assay),
    top_n
  )
  data_expression <- t(matrix)
  quality <- WGCNA::goodSamplesGenes(data_expression, verbose = 0)
  removed_samples <- rownames(data_expression)[!quality$goodSamples]
  removed_features <- colnames(data_expression)[!quality$goodGenes]
  if (!quality$allOK) {
    warning(
      "WGCNA power selection removed ", length(removed_samples),
      " sample(s) and ", length(removed_features), " feature(s).",
      call. = FALSE
    )
    data_expression <- data_expression[
      quality$goodSamples,
      quality$goodGenes,
      drop = FALSE
    ]
  }
  if (nrow(data_expression) < 4L || ncol(data_expression) < 3L) {
    stop(
      "Too few samples or features remain for WGCNA power selection.",
      call. = FALSE
    )
  }
  arguments <- c(list(
    data = data_expression,
    powerVector = powers,
    RsquaredCut = r_squared,
    networkType = network_type,
    verbose = verbose
  ), list(...))
  result <- .wgcna_call("pickSoftThreshold", arguments)
  result$bulkMAEQuality <- quality
  result$bulkMAERemovedSamples <- removed_samples
  result$bulkMAERemovedFeatures <- removed_features
  result$bulkMAENetworkType <- network_type
  result
}

#' Test preservation of reference WGCNA modules in another cohort
#'
#' The two cohorts remain explicit MAE inputs. Features are intersected and
#' ordered identically before constructing the `multiData` and `multiColor`
#' lists required by [WGCNA::modulePreservation()].
#'
#' @param reference,test Reference and test `MultiAssayExperiment` objects.
#' @param reference_experiment,test_experiment Experiment names.
#' @param reference_assay,test_assay Assays on comparable transformed scales.
#' @param module_colors Named vector assigning every reference feature to a
#'   WGCNA module color or label.
#' @param permutations Number of permutations.
#' @param network_type WGCNA network type.
#' @param seed Random seed passed to WGCNA.
#' @param verbose Backend verbosity.
#' @param save_permuted_statistics Save permutation statistics to a file. The
#'   default avoids hidden filesystem output.
#' @param ... Additional arguments passed to [WGCNA::modulePreservation()].
#'
#' @return The native WGCNA module-preservation result list, augmented with a
#'   `bulkMAEInputQuality` element documenting feature intersections and
#'   automatic sample/feature filtering.
#' @export
coexpr_preservation <- function(
    reference,
    reference_experiment,
    test,
    test_experiment,
    module_colors,
    reference_assay,
    test_assay,
    permutations = 200,
    network_type = "signed",
    seed = 1,
    verbose = 2,
    save_permuted_statistics = FALSE,
    ...
) {
  .require_backend("WGCNA", "to test module preservation")
  .network_assert_whole_number(permutations, "permutations", lower = 1L)
  .network_assert_seed(seed)
  if (
    !is.logical(save_permuted_statistics) ||
      length(save_permuted_statistics) != 1L ||
      is.na(save_permuted_statistics)
  ) {
    stop("`save_permuted_statistics` must be TRUE or FALSE.", call. = FALSE)
  }
  reference_matrix <- .pull_matrix(
    reference,
    reference_experiment,
    reference_assay
  )
  test_matrix <- .pull_matrix(test, test_experiment, test_assay)
  if (
    (!is.atomic(module_colors) && !is.factor(module_colors)) ||
      is.null(names(module_colors)) || anyDuplicated(names(module_colors))
  ) {
    stop("`module_colors` must be named by unique reference features.", call. = FALSE)
  }
  if (
    anyNA(names(module_colors)) || any(!nzchar(names(module_colors))) ||
      anyNA(module_colors) || any(!nzchar(as.character(module_colors)))
  ) {
    stop("`module_colors` names and values cannot be missing or empty.", call. = FALSE)
  }
  network_type <- match.arg(
    network_type,
    c("unsigned", "signed", "signed hybrid")
  )

  common_before_quality <- Reduce(
    intersect,
    list(
      rownames(reference_matrix),
      rownames(test_matrix),
      names(module_colors)
    )
  )
  if (length(common_before_quality) < 3L) {
    stop(
      "Reference, test, and `module_colors` need at least three common features.",
      call. = FALSE
    )
  }

  reference_expression <- t(
    reference_matrix[common_before_quality, , drop = FALSE]
  )
  test_expression <- t(test_matrix[common_before_quality, , drop = FALSE])
  reference_quality <- WGCNA::goodSamplesGenes(reference_expression, verbose = 0)
  test_quality <- WGCNA::goodSamplesGenes(test_expression, verbose = 0)
  good_features <- reference_quality$goodGenes & test_quality$goodGenes
  common <- common_before_quality[good_features]
  removed_reference_samples <- rownames(reference_expression)[
    !reference_quality$goodSamples
  ]
  removed_test_samples <- rownames(test_expression)[!test_quality$goodSamples]
  removed_quality_features <- common_before_quality[!good_features]
  if (
    length(removed_reference_samples) || length(removed_test_samples) ||
      length(removed_quality_features)
  ) {
    warning(
      "Module preservation removed ", length(removed_reference_samples),
      " reference sample(s), ", length(removed_test_samples),
      " test sample(s), and ", length(removed_quality_features),
      " feature(s) failing goodSamplesGenes().",
      call. = FALSE
    )
  }
  if (length(common) < 3L) {
    stop("Too few common features remain after WGCNA quality filtering.", call. = FALSE)
  }

  reference_expression <- reference_expression[
    reference_quality$goodSamples,
    good_features,
    drop = FALSE
  ]
  test_expression <- test_expression[
    test_quality$goodSamples,
    good_features,
    drop = FALSE
  ]
  if (nrow(reference_expression) < 4L || nrow(test_expression) < 4L) {
    stop("Each cohort needs at least four samples after quality filtering.", call. = FALSE)
  }

  multi_expression <- list(
    reference = list(data = reference_expression),
    test = list(data = test_expression)
  )
  multi_color <- list(reference = unname(module_colors[common]))
  arguments <- c(list(
    multiData = multi_expression,
    multiColor = multi_color,
    referenceNetworks = 1,
    nPermutations = permutations,
    networkType = network_type,
    randomSeed = seed,
    verbose = verbose,
    savePermutedStatistics = save_permuted_statistics
  ), list(...))
  result <- .with_seed(
    seed,
    .wgcna_call("modulePreservation", arguments)
  )
  result$bulkMAEInputQuality <- list(
    reference = reference_quality,
    test = test_quality,
    removedReferenceSamples = removed_reference_samples,
    removedTestSamples = removed_test_samples,
    excludedReferenceFeatures = setdiff(
      rownames(reference_matrix),
      common_before_quality
    ),
    excludedTestFeatures = setdiff(rownames(test_matrix), common_before_quality),
    excludedModuleFeatures = setdiff(names(module_colors), common_before_quality),
    removedQualityFeatures = removed_quality_features,
    analyzedFeatures = common,
    networkType = network_type
  )
  result
}

#' Infer a gene regulatory network with GENIE3
#'
#' @inheritParams mae_pull_assay
#' @param regulators Optional regulator names or indices.
#' @param targets Optional target names or indices.
#' @param trees Number of trees.
#' @param cores Number of worker cores.
#' @param seed Optional random seed.
#' @param ... Additional arguments passed to [GENIE3::GENIE3()].
#'
#' @return A native GENIE3 weight matrix.
#' @export
network_genie3 <- function(
    x,
    experiment,
    assay,
    regulators = NULL,
    targets = NULL,
    trees = 1000,
    cores = 1,
    seed = NULL,
    ...
) {
  .require_backend("GENIE3", "to infer a regulatory network")
  matrix <- .pull_matrix(x, experiment, assay)
  if (any(!is.finite(matrix))) {
    stop("GENIE3 requires a finite assay matrix.", call. = FALSE)
  }
  .network_assert_whole_number(trees, "trees", lower = 1L)
  .network_assert_whole_number(cores, "cores", lower = 1L)
  if (!is.null(seed)) .network_assert_seed(seed)
  .with_seed(seed, GENIE3::GENIE3(
    exprMatr = matrix,
    regulators = regulators,
    targets = targets,
    nTrees = trees,
    nCores = cores,
    ...
  ))
}

#' Convert GENIE3 weights to a ranked edge list
#'
#' @param weights Weight matrix returned by [network_genie3()].
#' @param threshold Optional minimum edge weight.
#' @param top Optional maximum number of edges.
#'
#' @return The native data frame returned by [GENIE3::getLinkList()].
#' @export
network_genie3_links <- function(weights, threshold = 0, top = NULL) {
  .require_backend("GENIE3", "to extract regulatory links")
  if (
    !is.numeric(threshold) || length(threshold) != 1L ||
      !is.finite(threshold) || threshold < 0
  ) {
    stop("`threshold` must be one finite, non-negative number.", call. = FALSE)
  }
  if (!is.null(top)) {
    .network_assert_whole_number(top, "top", lower = 1L)
  }
  arguments <- list(weightMatrix = weights, threshold = threshold)
  if (!is.null(top)) arguments$reportMax <- top
  do.call(GENIE3::getLinkList, arguments)
}

#' Retrieve STRING interactions for selected assay features
#'
#' @inheritParams mae_pull_assay
#' @param genes Optional feature identifiers; defaults to assay row names.
#' @param species NCBI taxonomy identifier.
#' @param version STRING database version.
#' @param score_threshold Minimum STRING score.
#' @param input_directory Optional STRING download/cache directory.
#' @param remove_unmapped Remove identifiers not mapped by STRING.
#'
#' @return A data frame of STRING interactions.
#' @export
network_string <- function(
    x,
    experiment,
    assay,
    genes = NULL,
    species = 9606,
    version = "12.0",
    score_threshold = 400,
    input_directory = "",
    remove_unmapped = TRUE
) {
  .require_backend("STRINGdb", "to retrieve protein interactions")
  matrix <- .pull_matrix(x, experiment, assay)
  if (is.null(genes)) genes <- rownames(matrix)
  genes <- .clean_gene_ids(genes, "genes")
  .network_assert_whole_number(species, "species", lower = 1L)
  .assert_scalar_character(version, "version")
  .network_assert_range(
    score_threshold,
    "score_threshold",
    lower = 0,
    upper = 1000
  )
  if (
    !is.character(input_directory) || length(input_directory) != 1L ||
      is.na(input_directory)
  ) {
    stop("`input_directory` must be one string.", call. = FALSE)
  }
  if (
    !is.logical(remove_unmapped) || length(remove_unmapped) != 1L ||
      is.na(remove_unmapped)
  ) {
    stop("`remove_unmapped` must be TRUE or FALSE.", call. = FALSE)
  }
  database <- STRINGdb::STRINGdb$new(
    version = version,
    species = species,
    score_threshold = score_threshold,
    input_directory = input_directory
  )
  mapping <- database$map(
    data.frame(feature = unique(as.character(genes))),
    "feature",
    removeUnmappedRows = remove_unmapped,
    takeFirst = TRUE
  )
  string_ids <- unique(mapping$STRING_id)
  string_ids <- string_ids[!is.na(string_ids) & nzchar(string_ids)]
  if (!length(string_ids)) {
    stop("No supplied genes could be mapped to STRING identifiers.", call. = FALSE)
  }
  database$get_interactions(string_ids)
}

.network_assert_whole_number <- function(value, argument, lower = 0L) {
  if (
    !is.numeric(value) || length(value) != 1L || is.na(value) ||
      !is.finite(value) ||
      value != trunc(value) || value < lower ||
      value > .Machine$integer.max
  ) {
    stop(
      "`", argument, "` must be one integer greater than or equal to ",
      lower, ".",
      call. = FALSE
    )
  }
  invisible(value)
}

.network_assert_fraction <- function(value, argument) {
  .network_assert_range(value, argument, lower = 0, upper = 1, open_lower = TRUE)
}

.network_assert_range <- function(
    value,
    argument,
    lower,
    upper,
    open_lower = FALSE
) {
  if (
    !is.numeric(value) || length(value) != 1L || is.na(value) ||
      !is.finite(value)
  ) {
    interval <- if (open_lower) "(" else "["
    stop(
      "`", argument, "` must be in ", interval, lower, ", ", upper, "].",
      call. = FALSE
    )
  }
  invalid_lower <- if (open_lower) value <= lower else value < lower
  if (invalid_lower || value > upper) {
    interval <- if (open_lower) "(" else "["
    stop(
      "`", argument, "` must be in ", interval, lower, ", ", upper, "].",
      call. = FALSE
    )
  }
  invisible(value)
}

.network_assert_seed <- function(seed) {
  if (
    !is.numeric(seed) || length(seed) != 1L || is.na(seed) ||
      !is.finite(seed) || seed < 0 ||
      seed != trunc(seed) || seed > .Machine$integer.max
  ) {
    stop("`seed` must be one non-negative integer.", call. = FALSE)
  }
  invisible(seed)
}

.consensus_cluster_call <- function(arguments, suppress_plot) {
  backend <- getExportedValue(
    "ConsensusClusterPlus",
    "ConsensusClusterPlus"
  )
  compatibility_environment <- new.env(parent = environment(backend))
  compatibility_environment$clusterTrackingPlot <- function(values) {
    if (is.null(dim(values))) values <- base::matrix(values, nrow = 1L)
    tracking_plot <- get(
      "clusterTrackingPlot",
      envir = asNamespace("ConsensusClusterPlus"),
      inherits = FALSE
    )
    tracking_plot(values)
  }
  environment(backend) <- compatibility_environment

  if (suppress_plot) {
    grDevices::pdf(file = NULL)
    on.exit(grDevices::dev.off(), add = TRUE)
  }
  do.call(backend, arguments)
}

.wgcna_call <- function(function_name, arguments) {
  # WGCNA 1.74 resolves a character `corFnc` from its caller in one path.
  # Keep the intended weighted implementation visible without attaching WGCNA.
  cor <- getExportedValue("WGCNA", "cor")
  do.call(getExportedValue("WGCNA", function_name), arguments)
}
