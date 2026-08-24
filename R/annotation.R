#' Remove Ensembl version suffixes
#'
#' @param ids Character vector of identifiers.
#'
#' @return A character vector with a terminal dot and digits removed.
#'
#' @examples
#' annotate_ensembl(c("ENSG00000141510.18", "TP53"))
#' @export
annotate_ensembl <- function(ids) {
  if (!is.atomic(ids) || is.list(ids)) {
    stop("`ids` must be an atomic vector of identifiers.", call. = FALSE)
  }
  ids <- as.character(ids)
  if (!length(ids) || anyNA(ids) || any(!nzchar(ids))) {
    stop("`ids` must contain non-missing identifiers.", call. = FALSE)
  }
  sub("\\.[0-9]+$", "", ids)
}

#' Map feature identifiers with an AnnotationDbi database
#'
#' This is a thin wrapper around [AnnotationDbi::mapIds()]. It deliberately
#' returns a mapping table and does not alter row names in the MAE.
#'
#' @param ids Character vector of source identifiers.
#' @param database An installed `AnnotationDb` object, such as `org.Hs.eg.db`.
#' @param from Source key type.
#' @param to Destination column.
#' @param multi_values Policy passed to `multiVals`.
#'
#' @return A long data frame with atomic `source_id` and `target_id` columns.
#'   When `multi_values` returns a list or `CharacterList`, one-to-many mappings
#'   are expanded to repeated source rows so the result can be passed directly
#'   to [annotate_rekey()].
#' @export
annotate_ids <- function(
    ids,
    database,
    from,
    to,
    multi_values = "first"
) {
  .require_backend("AnnotationDbi", "to map feature identifiers")
  ids <- .clean_gene_ids(ids, "ids")
  .assert_scalar_character(from, "from")
  .assert_scalar_character(to, "to")
  mapped <- AnnotationDbi::mapIds(
    x = database,
    keys = ids,
    column = to,
    keytype = from,
    multiVals = multi_values
  )

  source <- names(mapped)
  target <- unname(mapped)
  if (is.list(target) || methods::is(target, "List")) {
    target <- lapply(as.list(target), function(value) {
      value <- as.character(value)
      if (length(value)) value else NA_character_
    })
    source <- rep.int(source, lengths(target))
    target <- unlist(target, use.names = FALSE)
  }
  result <- data.frame(
    source_id = source,
    target_id = target,
    stringsAsFactors = FALSE,
    row.names = NULL
  )
  result[, c("source_id", "target_id"), drop = FALSE]
}

#' Query Ensembl through biomaRt
#'
#' @param values Values supplied to the selected filter.
#' @param attributes Ensembl attributes to return.
#' @param filter Ensembl filter name.
#' @param dataset Ensembl dataset, for example `hsapiens_gene_ensembl`.
#' @param biomart Ensembl BioMart name.
#' @param mirror Optional Ensembl mirror accepted by [biomaRt::useEnsembl()].
#' @param version Optional archived Ensembl version.
#'
#' @return A data frame returned by [biomaRt::getBM()].
#' @export
annotate_biomart <- function(
    values,
    attributes,
    filter,
    dataset,
    biomart = "genes",
    mirror = NULL,
    version = NULL
) {
  .require_backend("biomaRt", "to query Ensembl")
  if (!is.atomic(values) || is.list(values) || !length(values) || anyNA(values)) {
    stop("`values` must contain non-missing filter values.", call. = FALSE)
  }
  if (
    !is.character(attributes) || !length(attributes) || anyNA(attributes) ||
      any(!nzchar(attributes)) || anyDuplicated(attributes)
  ) {
    stop("`attributes` must contain unique, non-empty names.", call. = FALSE)
  }
  .assert_scalar_character(filter, "filter")
  .assert_scalar_character(dataset, "dataset")
  .assert_scalar_character(biomart, "biomart")

  mart_args <- list(biomart = biomart, dataset = dataset)
  if (!is.null(mirror)) {
    mart_args$mirror <- mirror
  }
  if (!is.null(version)) {
    mart_args$version <- version
  }
  mart <- do.call(biomaRt::useEnsembl, mart_args)

  biomaRt::getBM(
    attributes = attributes,
    filters = filter,
    values = values,
    mart = mart
  )
}

#' Retrieve representative transcript lengths from Ensembl
#'
#' This online convenience function queries BioMart for transcript lengths and
#' summarizes multiple annotated transcripts per input gene. These annotation
#' lengths are not sample-specific effective lengths; quantified RNA-seq data
#' should preferentially use the `length` assay from [import_tximport()].
#'
#' @param ids Unique feature identifiers.
#' @param dataset Ensembl dataset, for example `hsapiens_gene_ensembl`.
#' @param id_type BioMart gene identifier attribute and filter.
#' @param summary How to summarize multiple positive transcript lengths per
#'   gene.
#' @param strip_version Remove terminal numeric Ensembl version suffixes before
#'   querying.
#' @param biomart,mirror,version Passed to [annotate_biomart()].
#'
#' @return A feature-aligned data frame with `source_id`, `query_id`,
#'   `gene_length`, and `n_transcripts`. Missing genes retain `NA` lengths.
#' @export
annotate_gene_lengths <- function(
    ids,
    dataset = "hsapiens_gene_ensembl",
    id_type = "ensembl_gene_id",
    summary = c("median", "mean", "max"),
    strip_version = TRUE,
    biomart = "genes",
    mirror = NULL,
    version = NULL
) {
  if (!is.character(ids)) ids <- as.character(ids)
  .assert_ids(ids, "`ids`")
  .assert_scalar_character(dataset, "dataset")
  .assert_scalar_character(id_type, "id_type")
  summary <- match.arg(summary)
  if (
    !is.logical(strip_version) || length(strip_version) != 1L ||
      is.na(strip_version)
  ) {
    stop("`strip_version` must be TRUE or FALSE.", call. = FALSE)
  }
  query_ids <- if (strip_version) annotate_ensembl(ids) else ids
  query <- annotate_biomart(
    values = unique(query_ids),
    attributes = c(id_type, "transcript_length"),
    filter = id_type,
    dataset = dataset,
    biomart = biomart,
    mirror = mirror,
    version = version
  )
  if (!all(c(id_type, "transcript_length") %in% names(query))) {
    stop("BioMart did not return the requested identifier and length fields.", call. = FALSE)
  }
  gene <- as.character(query[[id_type]])
  length_value <- suppressWarnings(as.numeric(query$transcript_length))
  keep <- !is.na(gene) & nzchar(gene) & is.finite(length_value) & length_value > 0
  gene <- gene[keep]
  length_value <- length_value[keep]
  if (!length(gene)) {
    stop("BioMart returned no positive transcript lengths for the supplied IDs.", call. = FALSE)
  }
  values <- split(length_value, gene)
  summarize <- switch(
    summary,
    median = stats::median,
    mean = base::mean,
    max = base::max
  )
  summarized_lengths <- vapply(values, summarize, numeric(1))
  transcript_counts <- base::lengths(values)
  result <- data.frame(
    source_id = ids,
    query_id = query_ids,
    gene_length = unname(summarized_lengths[query_ids]),
    n_transcripts = unname(transcript_counts[query_ids]),
    row.names = ids,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  missing <- is.na(result$gene_length)
  if (any(missing)) {
    warning(
      sum(missing), " feature(s) had no positive Ensembl transcript length.",
      call. = FALSE
    )
  }
  attr(result, "summary") <- summary
  attr(result, "dataset") <- dataset
  result
}
