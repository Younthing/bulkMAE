#' Resolve a local organism annotation database
#'
#' @param organism Supported organism name.
#' @param package Optional installed `OrgDb` package name. When omitted,
#'   `org.Hs.eg.db` or `org.Mm.eg.db` is selected from `organism`.
#'
#' @return The installed `OrgDb` object exported by `package`.
#' @export
annotation_orgdb <- function(
    organism = c("human", "mouse"),
    package = NULL
) {
  organism <- match.arg(organism)
  if (is.null(package)) {
    package <- switch(
      organism,
      human = "org.Hs.eg.db",
      mouse = "org.Mm.eg.db"
    )
  } else {
    .resource_assert_scalar_character(package, "package")
  }

  if (!requireNamespace(package, quietly = TRUE)) {
    stop(
      "Annotation package `", package, "` is not installed. Install it with ",
      "pak::pkg_install('bioc::", package, "').",
      call. = FALSE
    )
  }
  database <- tryCatch(
    getExportedValue(package, package),
    error = function(condition) {
      stop(
        "Installed package `", package,
        "` does not export an object named `", package, "`.",
        call. = FALSE
      )
    }
  )
  if (!methods::is(database, "OrgDb")) {
    stop(
      "Exported object `", package, "::", package,
      "` is not an OrgDb database.",
      call. = FALSE
    )
  }
  database
}

#' Apply an identifier mapping to common gene-level values
#'
#' A source identifier that maps to multiple targets is expanded to every
#' target. When multiple input rows or elements then share a target, the
#' explicitly selected `duplicates` rule is applied. Exact duplicate mapping
#' rows are removed before expansion so they cannot multiply a value.
#'
#' @param values A numeric gene-by-sample matrix, a named numeric vector, or a
#'   character vector of identifiers.
#' @param mapping A data frame containing source and target identifier columns.
#' @param source,target Column names in `mapping`.
#' @param duplicates Aggregation rule for numeric values sharing a target.
#' @param drop_unmapped Drop source identifiers without a non-missing target.
#'   When `FALSE`, their original source identifiers are retained.
#'
#' @return A target-keyed numeric matrix, named numeric vector, or character
#'   vector, matching the input value type.
#' @export
annotate_rekey <- function(
    values,
    mapping,
    source = "source_id",
    target = "target_id",
    duplicates = c("sum", "mean", "max_abs", "first"),
    drop_unmapped = TRUE
) {
  duplicates <- match.arg(duplicates)
  .resource_assert_scalar_character(source, "source")
  .resource_assert_scalar_character(target, "target")
  if (identical(source, target)) {
    stop("`source` and `target` must name different columns.", call. = FALSE)
  }
  if (
    !is.logical(drop_unmapped) || length(drop_unmapped) != 1L ||
      is.na(drop_unmapped)
  ) {
    stop("`drop_unmapped` must be TRUE or FALSE.", call. = FALSE)
  }
  mapping <- .rekey_mapping(mapping, source, target)

  if (is.matrix(values)) {
    if (!is.numeric(values)) {
      stop("Matrix `values` must be numeric.", call. = FALSE)
    }
    if (is.null(rownames(values))) {
      stop("Matrix `values` must have source identifiers as row names.", call. = FALSE)
    }
    .resource_assert_ids(rownames(values), "Matrix row names", unique = FALSE)
    if (is.null(colnames(values))) {
      stop("Matrix `values` must have sample column names.", call. = FALSE)
    }
    .resource_assert_ids(colnames(values), "Matrix column names")
    if (any(!is.finite(values))) {
      stop("Numeric `values` must contain only finite values.", call. = FALSE)
    }
    expanded <- .rekey_expand(rownames(values), mapping, drop_unmapped)
    if (!length(expanded$index)) {
      return(values[FALSE, , drop = FALSE])
    }
    expanded_values <- values[expanded$index, , drop = FALSE]
    return(.rekey_aggregate_matrix(
      expanded_values,
      expanded$target,
      duplicates
    ))
  }

  if (is.numeric(values) && is.null(dim(values))) {
    if (is.null(names(values))) {
      stop("Numeric `values` must be named by source identifiers.", call. = FALSE)
    }
    .resource_assert_ids(names(values), "Names of numeric `values`", unique = FALSE)
    if (!length(values) || any(!is.finite(values))) {
      stop("Numeric `values` must contain finite values.", call. = FALSE)
    }
    expanded <- .rekey_expand(names(values), mapping, drop_unmapped)
    if (!length(expanded$index)) {
      return(stats::setNames(numeric(), character()))
    }
    matrix_value <- matrix(
      values[expanded$index],
      ncol = 1L,
      dimnames = list(NULL, "value")
    )
    result <- .rekey_aggregate_matrix(
      matrix_value,
      expanded$target,
      duplicates
    )
    return(stats::setNames(
      as.numeric(result[, 1L]),
      rownames(result)
    ))
  }

  if (is.character(values) && is.null(dim(values))) {
    .resource_assert_ids(values, "Character `values`", unique = FALSE)
    expanded <- .rekey_expand(values, mapping, drop_unmapped)
    return(unique(expanded$target))
  }

  stop(
    "`values` must be a numeric matrix, named numeric vector, or character ",
    "identifier vector.",
    call. = FALSE
  )
}

#' Standardize gene sets as a named list
#'
#' @param x A named list of gene vectors, or a long data frame. For a data
#'   frame, the first two columns are used when `term` and `gene` are omitted.
#' @param term,gene Column name or one-based column index for terms and genes.
#' @param min_size,max_size Inclusive unique-gene size limits.
#'
#' @return A uniquely named list of unique character gene identifiers.
#' @export
gene_sets_prepare <- function(
    x,
    term = NULL,
    gene = NULL,
    min_size = 1,
    max_size = Inf
) {
  .gene_sets_assert_sizes(min_size, max_size)

  if (is.data.frame(x)) {
    if (ncol(x) < 2L) {
      stop("A long gene-set table must contain at least two columns.", call. = FALSE)
    }
    term_column <- .gene_sets_column(x, term, 1L, "term")
    gene_column <- .gene_sets_column(x, gene, 2L, "gene")
    if (identical(term_column, gene_column)) {
      stop("`term` and `gene` must select different columns.", call. = FALSE)
    }
    if (is.list(x[[term_column]]) || is.list(x[[gene_column]])) {
      stop("Gene-set term and gene columns must be atomic vectors.", call. = FALSE)
    }
    terms <- as.character(x[[term_column]])
    genes <- as.character(x[[gene_column]])
    .resource_assert_ids(terms, "Gene-set terms", unique = FALSE)
    .resource_assert_ids(genes, "Gene-set genes", unique = FALSE)
    pairs <- unique(data.frame(
      term = terms,
      gene = genes,
      stringsAsFactors = FALSE
    ))
    term_levels <- unique(pairs$term)
    result <- split(
      pairs$gene,
      factor(pairs$term, levels = term_levels)
    )
    result <- lapply(result, unique)
  } else {
    if (
      !is.list(x) || is.data.frame(x) || !length(x) || is.null(names(x)) ||
        anyNA(names(x)) || any(!nzchar(names(x))) || anyDuplicated(names(x))
    ) {
      stop(
        "`x` must be a non-empty, uniquely named list or a long data frame.",
        call. = FALSE
      )
    }
    result <- lapply(names(x), function(name) {
      value <- x[[name]]
      if (!is.atomic(value) || is.list(value)) {
        stop("Every gene set must be an atomic vector.", call. = FALSE)
      }
      value <- as.character(value)
      if (!length(value)) return(character())
      .resource_assert_ids(
        value,
        paste0("Genes in set `", name, "`"),
        unique = FALSE
      )
      unique(value)
    })
    names(result) <- names(x)
  }

  sizes <- lengths(result)
  keep <- sizes >= min_size & sizes <= max_size
  result <- result[keep]
  if (!length(result)) {
    stop(
      "No gene set remains after `min_size` and `max_size` filtering.",
      call. = FALSE
    )
  }
  result
}

#' Read gene sets from a GMT file
#'
#' @param file Path to a Gene Matrix Transposed (`.gmt`) file.
#' @param ... Size filters passed to [gene_sets_prepare()].
#'
#' @return A named gene-set list. Set descriptions and the normalized source
#'   path are stored in `description` and `source_file` attributes.
#' @export
gene_sets_read_gmt <- function(file, ...) {
  .resource_assert_scalar_character(file, "file")
  if (!file.exists(file) || dir.exists(file)) {
    stop("`file` must identify an existing GMT file.", call. = FALSE)
  }
  lines <- readLines(file, warn = FALSE, encoding = "UTF-8")
  lines <- lines[nzchar(trimws(lines))]
  if (!length(lines)) {
    stop("The GMT file contains no gene sets.", call. = FALSE)
  }
  fields <- strsplit(lines, "\t", fixed = TRUE)
  malformed <- lengths(fields) < 3L
  if (any(malformed)) {
    stop(
      "Every GMT line must contain a term, description, and at least one gene; ",
      "malformed line(s): ",
      paste(which(malformed), collapse = ", "),
      ".",
      call. = FALSE
    )
  }

  terms <- vapply(fields, `[[`, character(1), 1L)
  descriptions <- vapply(fields, `[[`, character(1), 2L)
  .resource_assert_ids(terms, "GMT terms", unique = FALSE)
  genes <- lapply(fields, function(value) value[-c(1L, 2L)])
  invalid_gene <- vapply(
    genes,
    function(value) !length(value) || anyNA(value) || any(!nzchar(value)),
    logical(1)
  )
  if (any(invalid_gene)) {
    stop(
      "GMT genes must be non-missing identifiers; invalid line(s): ",
      paste(which(invalid_gene), collapse = ", "),
      ".",
      call. = FALSE
    )
  }
  description_by_term <- split(descriptions, terms)
  conflicting <- vapply(
    description_by_term,
    function(value) length(unique(value)) > 1L,
    logical(1)
  )
  if (any(conflicting)) {
    stop(
      "Repeated GMT terms must use one description: ",
      paste(names(conflicting)[conflicting], collapse = ", "),
      ".",
      call. = FALSE
    )
  }

  long <- data.frame(
    term = rep(terms, lengths(genes)),
    gene = unlist(genes, use.names = FALSE),
    stringsAsFactors = FALSE
  )
  result <- gene_sets_prepare(long, term = "term", gene = "gene", ...)
  descriptions <- vapply(description_by_term, `[[`, character(1), 1L)
  attr(result, "description") <- descriptions[names(result)]
  attr(result, "source_file") <- normalizePath(file, mustWork = TRUE)
  result
}

#' Retrieve MSigDB gene sets through msigdbr
#'
#' @param species Target species accepted by [msigdbr::msigdbr()].
#' @param collection MSigDB collection code.
#' @param subcollection Optional MSigDB subcollection code.
#' @param id_type Identifier column returned in each set.
#' @param db_species MSigDB source database species.
#' @param ... Additional arguments passed to [msigdbr::msigdbr()].
#'
#' @return A named gene-set list with MSigDB source, version, species,
#'   collection, subcollection, and identifier metadata stored as attributes.
#' @export
gene_sets_msigdb <- function(
    species = "human",
    collection = "H",
    subcollection = NULL,
    id_type = c("gene_symbol", "ncbi_gene", "ensembl_gene"),
    db_species = "HS",
    ...
) {
  .require_backend("msigdbr", "to retrieve MSigDB gene sets")
  .resource_assert_scalar_character(species, "species")
  .resource_assert_scalar_character(collection, "collection")
  if (!is.null(subcollection)) {
    .resource_assert_scalar_character(subcollection, "subcollection")
  }
  .resource_assert_scalar_character(db_species, "db_species")
  id_type <- match.arg(id_type)

  arguments <- list(
    db_species = db_species,
    species = species,
    collection = collection,
    subcollection = subcollection,
    ...
  )
  resource <- do.call(msigdbr::msigdbr, arguments)
  required <- c("gs_name", id_type, "db_version")
  if (!is.data.frame(resource) || any(!required %in% names(resource))) {
    stop(
      "msigdbr returned a table without required columns: ",
      paste(required, collapse = ", "),
      ".",
      call. = FALSE
    )
  }
  identifiers <- as.character(resource[[id_type]])
  terms <- as.character(resource$gs_name)
  valid <- !is.na(terms) & nzchar(terms) &
    !is.na(identifiers) & nzchar(identifiers)
  if (!any(valid)) {
    stop("The selected MSigDB resource contains no usable identifiers.", call. = FALSE)
  }
  if (any(!valid)) {
    warning(
      sum(!valid), " MSigDB row(s) without a term or `", id_type,
      "` identifier were removed.",
      call. = FALSE
    )
  }
  result <- gene_sets_prepare(data.frame(
    term = terms[valid],
    gene = identifiers[valid],
    stringsAsFactors = FALSE
  ))

  attr(result, "source") <- "msigdbr"
  attr(result, "db_version") <- unique(as.character(resource$db_version))
  attr(result, "species") <- species
  attr(result, "db_species") <- db_species
  attr(result, "collection") <- unique(as.character(resource$gs_collection))
  attr(result, "subcollection") <- unique(
    as.character(resource$gs_subcollection)
  )
  attr(result, "id_type") <- id_type
  if ("db_target_species" %in% names(resource)) {
    attr(result, "db_target_species") <- unique(
      as.character(resource$db_target_species)
    )
  }
  result
}

#' List molecular resources available through decoupleR
#'
#' @return The native character vector returned by
#'   [decoupleR::show_resources()].
#' @export
activity_resources <- function() {
  .require_backend("decoupleR", "to list molecular resources")
  decoupleR::show_resources()
}

#' List supported signatureSearch reference databases
#'
#' The databases are distributed through Bioconductor ExperimentHub and are
#' downloaded into its local cache when first requested. This function reports
#' metadata only; it does not download a database.
#'
#' @return A data frame describing supported database identifiers, human Entrez
#'   gene IDs, stored value types, ExperimentHub records, and access mode.
#' @export
drug_lincs_databases <- function() {
  data.frame(
    database = c("cmap", "cmap_expr", "lincs", "lincs_expr", "lincs2"),
    id_type = rep("ENTREZID", 5L),
    value_type = c(
      "log2_fold_change",
      "expression_intensity",
      "moderated_z_score",
      "expression_intensity",
      "moderated_z_score"
    ),
    experiment_hub_id = c("EH3223", "EH3224", "EH3226", "EH3227", "EH7297"),
    access = rep("online_cached", 5L),
    stringsAsFactors = FALSE
  )
}

.resource_assert_scalar_character <- function(value, argument) {
  if (
    !is.character(value) || length(value) != 1L || is.na(value) ||
      !nzchar(value)
  ) {
    stop("`", argument, "` must be one non-empty string.", call. = FALSE)
  }
  invisible(value)
}

.resource_assert_ids <- function(value, label, unique = TRUE) {
  if (
    !length(value) || anyNA(value) || any(!nzchar(as.character(value))) ||
      (unique && anyDuplicated(value))
  ) {
    qualifier <- if (unique) "unique, " else ""
    stop(label, " must be ", qualifier, "non-missing identifiers.", call. = FALSE)
  }
  invisible(value)
}

.rekey_mapping <- function(mapping, source, target) {
  if (!is.data.frame(mapping)) {
    stop("`mapping` must be a data frame.", call. = FALSE)
  }
  missing_columns <- setdiff(c(source, target), names(mapping))
  if (length(missing_columns)) {
    stop(
      "`mapping` is missing column(s): ",
      paste(missing_columns, collapse = ", "),
      ".",
      call. = FALSE
    )
  }
  if (is.list(mapping[[source]]) || is.list(mapping[[target]])) {
    stop("Mapping source and target columns must be atomic vectors.", call. = FALSE)
  }
  source_values <- as.character(mapping[[source]])
  target_values <- as.character(mapping[[target]])
  if (!length(source_values)) {
    return(data.frame(
      source = character(),
      target = character(),
      stringsAsFactors = FALSE
    ))
  }
  .resource_assert_ids(source_values, "Mapping source identifiers", unique = FALSE)
  target_values[is.na(target_values) | !nzchar(target_values)] <- NA_character_
  unique(data.frame(
    source = source_values,
    target = target_values,
    stringsAsFactors = FALSE
  ))
}

.rekey_expand <- function(ids, mapping, drop_unmapped) {
  valid <- !is.na(mapping$target)
  targets <- split(
    mapping$target[valid],
    factor(mapping$source[valid], levels = unique(mapping$source[valid]))
  )
  resolved <- lapply(ids, function(id) {
    mapped <- targets[[as.character(id)]]
    if (is.null(mapped) || !length(mapped)) {
      if (drop_unmapped) return(character())
      mapped <- as.character(id)
    }
    unname(mapped)
  })
  keep <- lengths(resolved) > 0L
  list(
    index = rep.int(which(keep), lengths(resolved)[keep]),
    target = if (any(keep)) {
      unlist(resolved[keep], use.names = FALSE)
    } else {
      character()
    }
  )
}

.rekey_aggregate_matrix <- function(values, targets, duplicates) {
  target_levels <- unique(targets)
  groups <- split(seq_along(targets), factor(targets, levels = target_levels))
  aggregate_group <- switch(
    duplicates,
    sum = function(index) colSums(values[index, , drop = FALSE]),
    mean = function(index) colMeans(values[index, , drop = FALSE]),
    first = function(index) values[index[[1L]], , drop = TRUE],
    max_abs = function(index) {
      current <- values[index, , drop = FALSE]
      apply(current, 2L, function(column) column[[which.max(abs(column))]])
    }
  )
  result <- do.call(rbind, lapply(groups, aggregate_group))
  if (is.null(dim(result))) {
    result <- matrix(result, ncol = ncol(values))
  }
  dimnames(result) <- list(target_levels, colnames(values))
  result
}

.gene_sets_assert_sizes <- function(min_size, max_size) {
  if (
    !is.numeric(min_size) || length(min_size) != 1L || is.na(min_size) ||
      !is.finite(min_size) || min_size != trunc(min_size) || min_size < 1L
  ) {
    stop("`min_size` must be one positive integer.", call. = FALSE)
  }
  if (
    !is.numeric(max_size) || length(max_size) != 1L || is.na(max_size) ||
      max_size < min_size ||
      (is.finite(max_size) && max_size != trunc(max_size))
  ) {
    stop(
      "`max_size` must be an integer at least `min_size`, or `Inf`.",
      call. = FALSE
    )
  }
  invisible(TRUE)
}

.gene_sets_column <- function(x, selector, default, argument) {
  if (is.null(selector)) return(names(x)[[default]])
  if (
    is.character(selector) && length(selector) == 1L && !is.na(selector) &&
      nzchar(selector) && selector %in% names(x)
  ) {
    return(selector)
  }
  if (
    is.numeric(selector) && length(selector) == 1L && !is.na(selector) &&
      is.finite(selector) && selector == trunc(selector) && selector >= 1L &&
      selector <= ncol(x)
  ) {
    return(names(x)[[as.integer(selector)]])
  }
  stop(
    "`", argument, "` must name one gene-set table column or its index.",
    call. = FALSE
  )
}
