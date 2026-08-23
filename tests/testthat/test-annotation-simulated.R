test_that("annotate_ensembl removes terminal versions without altering IDs", {
  identifiers <- c(
    "ENSG00000141510.18",
    "ENST00000269305.9",
    "TP53",
    "gene.2.extra"
  )

  expect_identical(
    annotate_ensembl(identifiers),
    c("ENSG00000141510", "ENST00000269305", "TP53", "gene.2.extra")
  )
})

test_that("annotate_ids maps stable human symbols with a local OrgDb", {
  suppressWarnings(skip_if_not_installed("AnnotationDbi"))
  suppressWarnings(skip_if_not_installed("org.Hs.eg.db"))
  database <- getExportedValue("org.Hs.eg.db", "org.Hs.eg.db")

  mapped <- annotate_ids(
    c("TP53", "EGFR", "TP53"),
    database = database,
    from = "SYMBOL",
    to = "ENTREZID"
  )

  expect_identical(mapped$source_id, c("TP53", "EGFR"))
  expect_identical(as.character(mapped$target_id), c("7157", "1956"))
  expect_identical(names(mapped), c("source_id", "target_id"))
  expect_error(
    annotate_ids("TP53", database, from = "", to = "ENTREZID"),
    "from.*non-empty string"
  )
})

test_that("annotate_biomart rejects invalid local inputs before any query", {
  suppressWarnings(skip_if_not_installed("biomaRt"))

  expect_error(
    annotate_biomart(
      values = character(),
      attributes = "ensembl_gene_id",
      filter = "hgnc_symbol",
      dataset = "hsapiens_gene_ensembl"
    ),
    "values.*non-missing"
  )
  expect_error(
    annotate_biomart(
      values = "TP53",
      attributes = c("ensembl_gene_id", "ensembl_gene_id"),
      filter = "hgnc_symbol",
      dataset = "hsapiens_gene_ensembl"
    ),
    "attributes.*unique"
  )
  expect_error(
    annotate_biomart(
      values = "TP53",
      attributes = "ensembl_gene_id",
      filter = "",
      dataset = "hsapiens_gene_ensembl"
    ),
    "filter.*non-empty string"
  )
  expect_error(
    annotate_biomart(
      values = "TP53",
      attributes = "ensembl_gene_id",
      filter = "hgnc_symbol",
      dataset = NA_character_
    ),
    "dataset.*non-empty string"
  )
})

test_that("annotate_biomart live Ensembl access is explicitly opt-in", {
  suppressWarnings(skip_if_not_installed("biomaRt"))
  skip_if(
    !identical(
      tolower(Sys.getenv("BULKMAE_RUN_ONLINE_TESTS", "false")),
      "true"
    ),
    "Set BULKMAE_RUN_ONLINE_TESTS=true to run live Ensembl tests."
  )

  mapped <- annotate_biomart(
    values = "TP53",
    attributes = c("hgnc_symbol", "ensembl_gene_id"),
    filter = "hgnc_symbol",
    dataset = "hsapiens_gene_ensembl"
  )

  expect_true(is.data.frame(mapped))
  expect_true("TP53" %in% mapped$hgnc_symbol)
  expect_true(any(nzchar(mapped$ensembl_gene_id)))
})
