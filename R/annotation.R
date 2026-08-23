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
#' @return A data frame with `source_id` and `target_id` columns.
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

  target <- unname(mapped)
  if (is.list(target) || methods::is(target, "List")) {
    target <- I(as.list(target))
  }
  result <- data.frame(
    source_id = names(mapped),
    stringsAsFactors = FALSE,
    row.names = NULL
  )
  result$target_id <- target
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
