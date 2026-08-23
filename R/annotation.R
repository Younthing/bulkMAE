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
  sub("\\.[0-9]+$", "", as.character(ids))
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
  mapped <- AnnotationDbi::mapIds(
    x = database,
    keys = unique(as.character(ids)),
    column = to,
    keytype = from,
    multiVals = multi_values
  )

  data.frame(
    source_id = names(mapped),
    target_id = unname(mapped),
    stringsAsFactors = FALSE,
    row.names = NULL
  )
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
