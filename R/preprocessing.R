#' Calculate TMM normalization factors
#'
#' @inheritParams mae_pull_assay
#' @param method Normalization method passed to [edgeR::normLibSizes()].
#' @param ... Additional arguments passed to [edgeR::normLibSizes()].
#'
#' @return A native [edgeR::DGEList] containing normalization factors.
#' @export
normalize_tmm <- function(
    x,
    experiment,
    assay = "counts",
    method = "TMM",
    ...
) {
  .require_backend("edgeR", "to calculate TMM normalization factors")
  y <- .make_dge_list(x, experiment, assay)
  edgeR::normLibSizes(y, method = method, ...)
}

#' Estimate DESeq2 size factors
#'
#' @inheritParams mae_pull_assay
#' @param type Size-factor estimator passed to `DESeq2::estimateSizeFactors()`.
#' @param ... Additional arguments passed to `DESeq2::estimateSizeFactors()`.
#'
#' @return A native `DESeqDataSet` with estimated size factors.
#' @export
normalize_deseq <- function(
    x,
    experiment,
    assay = "counts",
    type = "ratio",
    ...
) {
  .require_backend("DESeq2", "to estimate size factors")
  dds <- .make_deseq_dataset(x, experiment, assay, design = ~1)
  DESeq2::estimateSizeFactors(dds, type = type, ...)
}

#' Convert counts to transcripts per million
#'
#' Gene lengths may be supplied directly, read from one `rowData` column, or
#' read from an aligned assay such as the sample-specific effective-length
#' assay produced by [import_tximport()]. No external preprocessing step is
#' required.
#'
#' @inheritParams mae_pull_assay
#' @param lengths A positive numeric vector, a feature-named numeric vector, a
#'   feature-by-sample matrix, or the name of a `rowData` column or assay.
#'
#' @return A finite feature-by-sample TPM matrix whose columns sum to one
#'   million.
#' @export
normalize_tpm <- function(
    x,
    experiment,
    lengths,
    assay = "counts"
) {
  counts <- .as_count_matrix(x, experiment, assay)
  se <- .pull_se(x, experiment)
  if (is.character(lengths) && length(lengths) == 1L) {
    row_data <- as.data.frame(
      SummarizedExperiment::rowData(se),
      optional = TRUE
    )
    if (lengths %in% names(row_data)) {
      lengths <- row_data[[lengths]]
      names(lengths) <- rownames(se)
    } else if (lengths %in% SummarizedExperiment::assayNames(se)) {
      lengths <- .pull_matrix(x, experiment, lengths)
    } else {
      stop(
        "Character `lengths` must name a `rowData` column or assay.",
        call. = FALSE
      )
    }
  }

  if (is.atomic(lengths) && is.null(dim(lengths))) {
    if (!is.numeric(lengths)) {
      stop("`lengths` must contain numeric values.", call. = FALSE)
    }
    if (!is.null(names(lengths))) {
      .assert_ids(names(lengths), "Names of `lengths`")
      if (!setequal(names(lengths), rownames(counts))) {
        stop("Named `lengths` must contain every feature exactly once.", call. = FALSE)
      }
      lengths <- lengths[rownames(counts)]
    } else if (length(lengths) != nrow(counts)) {
      stop("`lengths` must contain one value per feature.", call. = FALSE)
    }
    lengths <- matrix(
      lengths,
      nrow = nrow(counts),
      ncol = ncol(counts),
      dimnames = dimnames(counts)
    )
  } else {
    lengths <- as.matrix(lengths)
    if (!is.numeric(lengths)) {
      stop("`lengths` must contain numeric values.", call. = FALSE)
    }
    if (is.null(rownames(lengths)) || is.null(colnames(lengths))) {
      stop("Matrix `lengths` must have feature and sample names.", call. = FALSE)
    }
    .assert_ids(rownames(lengths), "`lengths` feature names")
    .assert_ids(colnames(lengths), "`lengths` sample names")
    if (
      !setequal(rownames(lengths), rownames(counts)) ||
        !setequal(colnames(lengths), colnames(counts))
    ) {
      stop("Matrix `lengths` must align exactly with the count assay.", call. = FALSE)
    }
    lengths <- lengths[rownames(counts), colnames(counts), drop = FALSE]
  }
  if (any(!is.finite(lengths)) || any(lengths <= 0)) {
    stop("`lengths` must contain finite, strictly positive values.", call. = FALSE)
  }

  reads_per_kilobase <- counts / (lengths / 1000)
  totals <- colSums(reads_per_kilobase)
  if (any(!is.finite(totals)) || any(totals <= 0)) {
    stop("TPM normalization requires positive per-sample scaled totals.", call. = FALSE)
  }
  sweep(reads_per_kilobase, 2L, totals / 1e6, "/")
}

#' Apply the DESeq2 variance-stabilizing transformation
#'
#' @inheritParams mae_pull_assay
#' @param blind Whether the transformation should ignore the design.
#' @param design Design used when `blind = FALSE`.
#' @param fit_type Dispersion fit type passed to [DESeq2::vst()].
#' @param ... Additional arguments passed to [DESeq2::vst()].
#'
#' @return A native `DESeqTransform` object.
#' @export
transform_vst <- function(
    x,
    experiment,
    assay = "counts",
    blind = TRUE,
    design = ~1,
    fit_type = "parametric",
    ...
) {
  .require_backend("DESeq2", "to apply the VST")
  dds <- .make_deseq_dataset(x, experiment, assay, design = design)
  DESeq2::vst(dds, blind = blind, fitType = fit_type, ...)
}

#' Apply the DESeq2 regularized-log transformation
#'
#' @inheritParams transform_vst
#'
#' @return A native `DESeqTransform` object.
#' @export
transform_rlog <- function(
    x,
    experiment,
    assay = "counts",
    blind = TRUE,
    design = ~1,
    fit_type = "parametric",
    ...
) {
  .require_backend("DESeq2", "to apply the rlog transformation")
  dds <- .make_deseq_dataset(x, experiment, assay, design = design)
  DESeq2::rlog(dds, blind = blind, fitType = fit_type, ...)
}

#' Apply limma-voom transformation and precision weighting
#'
#' @inheritParams mae_pull_assay
#' @param formula Fixed-effect model formula.
#' @param normalize_method edgeR library normalization method.
#' @param plot Whether `voom` should draw its diagnostic plot.
#' @param ... Additional arguments passed to [limma::voom()].
#'
#' @return A native limma `EList`.
#' @export
transform_voom <- function(
    x,
    experiment,
    formula,
    assay = "counts",
    normalize_method = "TMM",
    plot = FALSE,
    ...
) {
  .require_backend("edgeR", "to construct a DGEList")
  .require_backend("limma", "to apply voom")
  data <- mae_samples(x, experiment)
  design <- .model_matrix(formula, data)
  y <- .make_dge_list(x, experiment, assay)
  y <- edgeR::normLibSizes(y, method = normalize_method)
  limma::voom(y, design = design, plot = plot, ...)
}

#' Adjust a continuous assay with ComBat
#'
#' @inheritParams mae_pull_assay
#' @param batch Metadata column name or batch vector.
#' @param preserve Optional formula for biological covariates to preserve.
#' @param parametric Use parametric empirical Bayes priors.
#' @param mean_only Adjust batch means but not scales.
#' @param reference_batch Optional reference batch.
#' @param ... Additional arguments passed to [sva::ComBat()].
#'
#' @return A corrected feature-by-sample matrix.
#' @export
adjust_combat <- function(
    x,
    experiment,
    batch,
    assay,
    preserve = NULL,
    parametric = TRUE,
    mean_only = FALSE,
    reference_batch = NULL,
    ...
) {
  .require_backend("sva", "to apply ComBat")
  matrix <- .pull_matrix(x, experiment, assay)
  if (any(!is.finite(matrix))) {
    stop("ComBat requires a finite continuous assay matrix.", call. = FALSE)
  }
  data <- mae_samples(x, experiment)
  batch <- .column_or_vector(batch, data, "batch")
  .assert_multiple_batches(batch, "batch")
  model <- if (is.null(preserve)) NULL else .model_matrix(preserve, data)
  if (!is.null(reference_batch)) {
    if (length(reference_batch) != 1L || is.na(reference_batch)) {
      stop("`reference_batch` must identify one batch.", call. = FALSE)
    }
    if (!as.character(reference_batch) %in% as.character(batch)) {
      stop("`reference_batch` is not present in `batch`.", call. = FALSE)
    }
  }

  sva::ComBat(
    dat = matrix,
    batch = batch,
    mod = model,
    par.prior = parametric,
    mean.only = mean_only,
    ref.batch = reference_batch,
    ...
  )
}

#' Remove batch effects from continuous expression data with limma
#'
#' This convenience method is intended for visualization or exploratory
#' analysis. Statistical models should usually include batch in their design.
#'
#' @inheritParams mae_pull_assay
#' @param batch Metadata column name or batch vector.
#' @param preserve Formula describing biological effects to preserve.
#' @param batch2 Optional second batch column or vector.
#' @param covariates Optional numeric matrix, or metadata column names, for
#'   continuous nuisance covariates.
#' @param ... Additional arguments passed to [limma::removeBatchEffect()].
#'
#' @return A corrected feature-by-sample matrix.
#' @export
adjust_batch <- function(
    x,
    experiment,
    batch,
    assay,
    preserve = ~1,
    batch2 = NULL,
    covariates = NULL,
    ...
) {
  .require_backend("limma", "to remove batch effects")
  matrix <- .pull_matrix(x, experiment, assay)
  if (any(!is.finite(matrix))) {
    stop(
      "Batch-effect removal requires a finite continuous assay matrix.",
      call. = FALSE
    )
  }
  data <- mae_samples(x, experiment)
  batch <- .column_or_vector(batch, data, "batch")
  .assert_multiple_batches(batch, "batch")
  if (!is.null(batch2)) {
    batch2 <- .column_or_vector(batch2, data, "batch2")
    .assert_multiple_batches(batch2, "batch2")
  }
  if (is.character(covariates)) {
    if (!length(covariates) || anyNA(covariates) || anyDuplicated(covariates)) {
      stop("`covariates` must name unique metadata columns.", call. = FALSE)
    }
    missing <- setdiff(covariates, names(data))
    if (length(missing)) {
      stop(
        "Unknown covariate columns: ", paste(missing, collapse = ", "),
        call. = FALSE
      )
    }
    covariates <- data[, covariates, drop = FALSE]
  }
  if (!is.null(covariates)) {
    covariates <- .sample_matrix(
      covariates,
      data,
      "covariates",
      numeric = TRUE
    )
    variable <- apply(covariates, 2L, stats::sd) > 0
    if (any(!variable)) {
      stop(
        "`covariates` contains constant columns: ",
        paste(colnames(covariates)[!variable], collapse = ", "),
        call. = FALSE
      )
    }
  }
  limma::removeBatchEffect(
    matrix,
    batch = batch,
    batch2 = batch2,
    covariates = covariates,
    design = .model_matrix(preserve, data),
    ...
  )
}

#' Adjust integer counts with ComBat-seq
#'
#' @inheritParams mae_pull_assay
#' @param batch Metadata column name or batch vector.
#' @param group Optional biological group column or vector to preserve.
#' @param covariates Optional formula for additional covariates to preserve.
#' @param ... Additional arguments passed to [sva::ComBat_seq()].
#'
#' @return A corrected count matrix.
#' @export
adjust_combatseq <- function(
    x,
    experiment,
    batch,
    assay = "counts",
    group = NULL,
    covariates = NULL,
    ...
) {
  .require_backend("sva", "to apply ComBat-seq")
  counts <- .as_count_matrix(x, experiment, assay, integer = TRUE)
  data <- mae_samples(x, experiment)
  batch <- .column_or_vector(batch, data, "batch")
  .assert_multiple_batches(batch, "batch")
  if (!is.null(group)) {
    group <- .column_or_vector(group, data, "group")
  }
  covariate_model <- if (is.null(covariates)) {
    NULL
  } else {
    .model_matrix(covariates, data)
  }

  sva::ComBat_seq(
    counts = counts,
    batch = batch,
    group = group,
    covar_mod = covariate_model,
    ...
  )
}

#' Estimate surrogate variables with sva
#'
#' @inheritParams mae_pull_assay
#' @param full Full biological model formula.
#' @param null Null model formula.
#' @param n_surrogates Optional number of surrogate variables. When omitted,
#'   [sva::num.sv()] is used.
#' @param method Method passed to [sva::sva()].
#' @param ... Additional arguments passed to [sva::sva()].
#'
#' @return A native `sva` result list.
#' @export
adjust_sva <- function(
    x,
    experiment,
    assay,
    full,
    null = ~1,
    n_surrogates = NULL,
    method = "irw",
    ...
) {
  .require_backend("sva", "to estimate surrogate variables")
  matrix <- .pull_matrix(x, experiment, assay)
  if (any(!is.finite(matrix))) {
    stop("SVA requires a finite assay matrix.", call. = FALSE)
  }
  data <- mae_samples(x, experiment)
  full_model <- .model_matrix(full, data)
  null_model <- .model_matrix(null, data)
  if (qr(cbind(full_model, null_model))$rank > qr(full_model)$rank) {
    stop("The `null` model must be nested within the `full` model.", call. = FALSE)
  }
  residual_df <- nrow(full_model) - qr(full_model)$rank
  if (is.null(n_surrogates)) {
    n_surrogates <- sva::num.sv(matrix, full_model, method = "leek")
  }
  if (
    !is.numeric(n_surrogates) || length(n_surrogates) != 1L ||
      is.na(n_surrogates) || !is.finite(n_surrogates) ||
      n_surrogates != trunc(n_surrogates) || n_surrogates < 0L ||
      n_surrogates > residual_df
  ) {
    stop(
      "`n_surrogates` must be an integer from zero through the full-model ",
      "residual degrees of freedom (", residual_df, ").",
      call. = FALSE
    )
  }

  sva::sva(
    dat = matrix,
    mod = full_model,
    mod0 = null_model,
    n.sv = n_surrogates,
    method = method,
    ...
  )
}

#' Estimate unwanted factors with RUVg
#'
#' @inheritParams mae_pull_assay
#' @param controls Control-feature names or indices passed as `cIdx`.
#' @param k Number of unwanted factors.
#' @param ... Additional arguments passed to [RUVSeq::RUVg()].
#'
#' @return The native object returned by [RUVSeq::RUVg()].
#' @export
adjust_ruv <- function(
    x,
    experiment,
    controls,
    k = 1,
    assay = "counts",
    ...
) {
  .require_backend("RUVSeq", "to estimate unwanted factors")
  counts <- .as_count_matrix(x, experiment, assay)
  if (is.character(controls)) {
    if (
      !length(controls) || anyNA(controls) || any(!nzchar(controls)) ||
        anyDuplicated(controls)
    ) {
      stop("Character `controls` must contain unique feature names.", call. = FALSE)
    }
    missing <- setdiff(controls, rownames(counts))
    if (length(missing)) {
      stop(
        "Unknown control features: ",
        paste(utils::head(missing, 10L), collapse = ", "),
        call. = FALSE
      )
    }
  } else if (is.logical(controls)) {
    if (length(controls) != nrow(counts) || anyNA(controls) || !any(controls)) {
      stop("Logical `controls` must select at least one assay feature.", call. = FALSE)
    }
  } else if (
    !is.numeric(controls) || !length(controls) || anyNA(controls) ||
      any(!is.finite(controls)) || any(controls != trunc(controls)) ||
      any(controls < 1L) || any(controls > nrow(counts)) ||
      anyDuplicated(controls)
  ) {
    stop("Numeric `controls` must be unique valid feature indices.", call. = FALSE)
  }
  if (
    !is.numeric(k) || length(k) != 1L || is.na(k) || !is.finite(k) ||
      k != trunc(k) || k < 1L || k >= ncol(counts)
  ) {
    stop("`k` must be a positive integer smaller than the sample count.", call. = FALSE)
  }
  RUVSeq::RUVg(counts, cIdx = controls, k = k, ...)
}

#' Partition expression variance among model terms
#'
#' @inheritParams mae_pull_assay
#' @param formula Variance-partition formula. Random effects use `(1|group)`.
#' @param ... Additional arguments passed to
#'   [variancePartition::fitExtractVarPartModel()].
#'
#' @return The native variance-partition result.
#' @export
de_variance <- function(x, experiment, formula, assay, ...) {
  .require_backend("variancePartition", "to partition expression variance")
  matrix <- .pull_matrix(x, experiment, assay)
  data <- mae_samples(x, experiment)
  if (any(!is.finite(matrix))) {
    stop("Variance partitioning requires a finite assay matrix.", call. = FALSE)
  }
  if (!inherits(formula, "formula")) {
    stop("`formula` must be a formula.", call. = FALSE)
  }
  variables <- all.vars(formula)
  missing <- setdiff(variables, names(data))
  if (length(missing)) {
    stop(
      "Unknown variables in `formula`: ", paste(missing, collapse = ", "),
      call. = FALSE
    )
  }
  if (length(variables) && anyNA(data[, variables, drop = FALSE])) {
    stop("Variables in `formula` must not contain missing values.", call. = FALSE)
  }
  variancePartition::fitExtractVarPartModel(
    exprObj = matrix,
    formula = formula,
    data = data,
    ...
  )
}

.make_deseq_dataset <- function(x, experiment, assay, design) {
  se <- .pull_se(x, experiment)
  txi <- .tximport_list(se, assay)
  data_frame <- as.data.frame(
    SummarizedExperiment::colData(se),
    optional = TRUE
  )
  if (inherits(design, "formula")) {
    .model_matrix(design, data_frame)
  } else {
    design <- .validate_design_matrix(design, data_frame, "design")
  }
  data <- S4Vectors::DataFrame(data_frame)

  if (!is.null(txi)) {
    return(DESeq2::DESeqDataSetFromTximport(
      txi = txi,
      colData = data,
      design = design
    ))
  }

  counts <- .as_count_matrix(x, experiment, assay, integer = TRUE)
  DESeq2::DESeqDataSetFromMatrix(
    countData = counts,
    colData = data,
    design = design
  )
}

.make_dge_list <- function(x, experiment, assay = "counts") {
  se <- .pull_se(x, experiment)
  txi <- .tximport_list(se, assay)
  if (!is.null(txi)) {
    if ("DGEListFromTximport" %in% getNamespaceExports("edgeR")) {
      constructor <- getExportedValue("edgeR", "DGEListFromTximport")
      return(constructor(txi))
    }
    return(.dge_list_from_tximport(txi))
  }
  edgeR::DGEList(counts = .as_count_matrix(x, experiment, assay))
}

.dge_list_from_tximport <- function(txi) {
  y <- edgeR::DGEList(counts = txi$counts)
  mode <- txi$countsFromAbundance %||% "no"
  if (!identical(mode, "no")) {
    return(y)
  }

  length_factors <- txi$length /
    exp(rowMeans(log(txi$length)))
  normalized_counts <- txi$counts / length_factors
  effective_libraries <- edgeR::normLibSizes(normalized_counts) *
    colSums(normalized_counts)
  offsets <- log(sweep(length_factors, 2L, effective_libraries, "*"))
  edgeR::scaleOffset(y, offsets)
}

.assert_multiple_batches <- function(batch, argument) {
  if (length(unique(batch)) < 2L) {
    stop("`", argument, "` must contain at least two batches.", call. = FALSE)
  }
  invisible(batch)
}
