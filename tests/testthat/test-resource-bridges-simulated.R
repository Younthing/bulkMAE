test_that("annotation_orgdb resolves installed human and mouse databases", {
  suppressWarnings(skip_if_not_installed("AnnotationDbi"))
  suppressWarnings(skip_if_not_installed("org.Hs.eg.db"))
  suppressWarnings(skip_if_not_installed("org.Mm.eg.db"))

  human <- suppressWarnings(annotation_orgdb("human"))
  mouse <- suppressWarnings(annotation_orgdb("mouse"))
  explicit <- suppressWarnings(annotation_orgdb(package = "org.Hs.eg.db"))

  expect_s4_class(human, "OrgDb")
  expect_s4_class(mouse, "OrgDb")
  expect_identical(explicit, human)
  expect_error(annotation_orgdb("rat"), "arg")
})

test_that("annotation_orgdb reports a reproducible pak installation command", {
  expect_error(
    annotation_orgdb(package = "org.bulkMAE.missing.eg.db"),
    "pak::pkg_install\\('bioc::org.bulkMAE.missing.eg.db'\\)"
  )
  expect_error(annotation_orgdb(package = ""), "package")
})

test_that("annotate_rekey expands mappings and aggregates matrix rows", {
  values <- matrix(
    c(2, -1, -5, 3, 7, 9, 4, 2),
    nrow = 4L,
    byrow = TRUE,
    dimnames = list(c("a", "b", "c", "d"), c("sample1", "sample2"))
  )
  mapping <- data.frame(
    source_id = c("a", "a", "a", "b", "c", "d", "d"),
    target_id = c("X", "Y", "X", "X", NA, "Z", ""),
    stringsAsFactors = FALSE
  )

  summed <- annotate_rekey(values, mapping, duplicates = "sum")
  averaged <- annotate_rekey(values, mapping, duplicates = "mean")
  strongest <- annotate_rekey(values, mapping, duplicates = "max_abs")
  first <- annotate_rekey(values, mapping, duplicates = "first")
  retained <- annotate_rekey(
    values,
    mapping,
    duplicates = "sum",
    drop_unmapped = FALSE
  )

  expect_identical(rownames(summed), c("X", "Y", "Z"))
  expect_identical(colnames(summed), colnames(values))
  expect_equal(unname(summed["X", ]), c(-3, 2))
  expect_equal(unname(summed["Y", ]), c(2, -1))
  expect_equal(unname(summed["Z", ]), c(4, 2))
  expect_equal(unname(averaged["X", ]), c(-1.5, 1))
  expect_equal(unname(strongest["X", ]), c(-5, 3))
  expect_equal(unname(first["X", ]), c(2, -1))
  expect_identical(rownames(retained), c("X", "Y", "c", "Z"))
  expect_equal(unname(retained["c", ]), c(7, 9))
})

test_that("annotate_rekey handles named statistics and identifier vectors", {
  mapping <- data.frame(
    old = c("a", "a", "b", "c", "d"),
    new = c("X", "Y", "X", NA, "Z"),
    stringsAsFactors = FALSE
  )
  statistics <- c(a = 2, b = -5, c = 7, d = 4)

  strongest <- annotate_rekey(
    statistics,
    mapping,
    source = "old",
    target = "new",
    duplicates = "max_abs",
    drop_unmapped = FALSE
  )
  identifiers <- annotate_rekey(
    c("a", "b", "c", "a"),
    mapping,
    source = "old",
    target = "new",
    drop_unmapped = FALSE
  )
  mapped_only <- annotate_rekey(
    c("a", "b", "c"),
    mapping,
    source = "old",
    target = "new"
  )

  expect_equal(strongest, c(X = -5, Y = 2, c = 7, Z = 4))
  expect_identical(identifiers, c("X", "Y", "c"))
  expect_identical(mapped_only, c("X", "Y"))
})

test_that("annotate_rekey rejects ambiguous or malformed inputs", {
  mapping <- data.frame(source_id = "a", target_id = "A")
  unnamed <- c(1, 2)
  non_finite <- c(a = Inf)
  no_names <- matrix(1:4, nrow = 2L)

  expect_error(annotate_rekey(unnamed, mapping), "must be named")
  expect_error(annotate_rekey(non_finite, mapping), "finite")
  expect_error(annotate_rekey(no_names, mapping), "row names")
  expect_error(annotate_rekey(list(a = 1), mapping), "must be a numeric")
  expect_error(annotate_rekey(c("a"), list()), "mapping.*data frame")
  expect_error(
    annotate_rekey(c("a"), data.frame(source_id = "a")),
    "missing column"
  )
  expect_error(
    annotate_rekey(c("a"), mapping, source = "source_id", target = "source_id"),
    "different columns"
  )
  expect_error(annotate_rekey(c("a"), mapping, duplicates = "guess"), "arg")
  expect_error(annotate_rekey(c("a"), mapping, drop_unmapped = NA), "TRUE or FALSE")
})

test_that("annotate_rekey returns correctly shaped empty results", {
  mapping <- data.frame(source_id = "other", target_id = "OTHER")
  matrix_value <- matrix(
    1:4,
    nrow = 2L,
    dimnames = list(c("a", "b"), c("s1", "s2"))
  )

  empty_matrix <- annotate_rekey(matrix_value, mapping)
  empty_vector <- annotate_rekey(c(a = 1, b = 2), mapping)
  empty_ids <- annotate_rekey(c("a", "b"), mapping)

  expect_identical(dim(empty_matrix), c(0L, 2L))
  expect_identical(colnames(empty_matrix), c("s1", "s2"))
  expect_identical(empty_vector, setNames(numeric(), character()))
  expect_identical(empty_ids, character())

  empty_mapping <- data.frame(
    source_id = character(),
    target_id = character()
  )
  expect_identical(
    annotate_rekey(c(a = 1, b = 2), empty_mapping, drop_unmapped = FALSE),
    c(a = 1, b = 2)
  )
  expect_identical(
    annotate_rekey(c(a = 1, b = 2), empty_mapping),
    stats::setNames(numeric(), character())
  )
})

test_that("gene_sets_prepare standardizes lists and arbitrary long tables", {
  list_sets <- gene_sets_prepare(
    list(first = c("g1", "g1", "g2"), second = c("g3", "g4", "g5")),
    min_size = 2L,
    max_size = 2L
  )
  long <- data.frame(
    pathway = c("p2", "p1", "p2", "p1", "p2"),
    weight = seq_len(5L),
    feature = c("g3", "g1", "g2", "g1", "g3"),
    stringsAsFactors = FALSE
  )
  table_sets <- gene_sets_prepare(
    long,
    term = "pathway",
    gene = "feature",
    min_size = 1L,
    max_size = 3L
  )
  indexed <- gene_sets_prepare(long[c("pathway", "feature")], term = 1L, gene = 2L)

  expect_identical(list_sets, list(first = c("g1", "g2")))
  expect_identical(names(table_sets), c("p2", "p1"))
  expect_identical(table_sets$p2, c("g3", "g2"))
  expect_identical(table_sets$p1, "g1")
  expect_identical(indexed, table_sets)
  expect_identical(
    gene_sets_prepare(list(empty = character(), valid = c("g1", "g2"))),
    list(valid = c("g1", "g2"))
  )
})

test_that("gene_sets_prepare validates schemas and size policies", {
  duplicated_names <- structure(list("g1", "g2"), names = c("set", "set"))

  expect_error(gene_sets_prepare(list()), "non-empty")
  expect_error(gene_sets_prepare(duplicated_names), "uniquely named")
  expect_error(
    gene_sets_prepare(data.frame(term = "p", gene = NA_character_)),
    "non-missing"
  )
  expect_error(
    gene_sets_prepare(data.frame(term = "p", gene = "g"), term = "missing"),
    "must name one"
  )
  expect_error(
    gene_sets_prepare(data.frame(term = "p", gene = "g"), term = 1L, gene = 1L),
    "different columns"
  )
  expect_error(gene_sets_prepare(list(set = "g1"), min_size = 0L), "min_size")
  expect_error(
    gene_sets_prepare(list(set = "g1"), min_size = 2L),
    "No gene set remains"
  )
})

test_that("gene_sets_read_gmt parses descriptions and repeated terms", {
  file <- tempfile(fileext = ".gmt")
  on.exit(unlink(file), add = TRUE)
  writeLines(c(
    "pathway_a\tdescription a\tg1\tg2",
    "pathway_b\tdescription b\tg2\tg3\tg4",
    "pathway_a\tdescription a\tg2\tg5"
  ), file)

  sets <- gene_sets_read_gmt(file, min_size = 3L)

  expect_identical(names(sets), c("pathway_a", "pathway_b"))
  expect_identical(sets$pathway_a, c("g1", "g2", "g5"))
  expect_identical(sets$pathway_b, c("g2", "g3", "g4"))
  expect_identical(
    attr(sets, "description"),
    c(pathway_a = "description a", pathway_b = "description b")
  )
  expect_identical(attr(sets, "source_file"), normalizePath(file))
})

test_that("gene_sets_read_gmt rejects missing and malformed files", {
  missing <- tempfile(fileext = ".gmt")
  malformed <- tempfile(fileext = ".gmt")
  conflict <- tempfile(fileext = ".gmt")
  on.exit(unlink(c(malformed, conflict)), add = TRUE)
  writeLines("set\tdescription", malformed)
  writeLines(c("set\tone\tg1", "set\ttwo\tg2"), conflict)

  expect_error(gene_sets_read_gmt(missing), "existing GMT")
  expect_error(gene_sets_read_gmt(malformed), "malformed line")
  expect_error(gene_sets_read_gmt(conflict), "one description")
})

test_that("gene_sets_msigdb returns versioned named MSigDB sets", {
  skip_if_not_installed("msigdbr")
  run_online <- identical(
    tolower(Sys.getenv("BULKMAE_RUN_ONLINE_TESTS", "false")),
    "true"
  )
  skip_if(
    !run_online,
    "MSigDB may download/cache data; set BULKMAE_RUN_ONLINE_TESTS=true."
  )

  sets <- suppressWarnings(gene_sets_msigdb(
    species = "human",
    collection = "H",
    id_type = "gene_symbol",
    db_species = "HS"
  ))

  expect_type(sets, "list")
  expect_gt(length(sets), 0L)
  expect_true(all(nzchar(names(sets))))
  expect_true(all(vapply(sets, is.character, logical(1))))
  expect_identical(attr(sets, "source"), "msigdbr")
  expect_identical(attr(sets, "species"), "human")
  expect_identical(attr(sets, "db_species"), "HS")
  expect_identical(attr(sets, "id_type"), "gene_symbol")
  expect_true(length(attr(sets, "db_version")) >= 1L)
  expect_true(all(nzchar(attr(sets, "db_version"))))
  expect_true("H" %in% attr(sets, "collection"))
  expect_error(gene_sets_msigdb(id_type = "unsupported"), "arg")
})

test_that("activity_resources exposes decoupleR resource names", {
  skip_if_not_installed("decoupleR", minimum_version = "2.0.0")
  run_online <- identical(
    tolower(Sys.getenv("BULKMAE_RUN_ONLINE_TESTS", "false")),
    "true"
  )
  skip_if(
    !run_online,
    "The OmniPath registry is online/cache-backed; set BULKMAE_RUN_ONLINE_TESTS=true."
  )

  resources <- suppressWarnings(suppressMessages(activity_resources()))

  expect_type(resources, "character")
  expect_gt(length(resources), 0L)
  expect_true(all(nzchar(resources)))
  expect_true("PROGENy" %in% resources)
})

test_that("drug_lincs_databases reports static cached-resource metadata", {
  databases <- drug_lincs_databases()

  expect_s3_class(databases, "data.frame")
  expect_identical(
    databases$database,
    c("cmap", "cmap_expr", "lincs", "lincs_expr", "lincs2")
  )
  expect_true(all(databases$id_type == "ENTREZID"))
  expect_true(all(databases$access == "online_cached"))
  expect_true(all(grepl("^EH[0-9]+$", databases$experiment_hub_id)))
  expect_identical(anyDuplicated(databases$database), 0L)
})
