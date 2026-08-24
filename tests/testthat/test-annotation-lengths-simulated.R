test_that("gene lengths are summarized and aligned without a live BioMart call", {
  testthat::local_mocked_bindings(
    annotate_biomart = function(values, attributes, filter, dataset, ...) {
      expect_identical(values, c("ENSG1", "ENSG2", "ENSG3"))
      expect_identical(attributes, c("ensembl_gene_id", "transcript_length"))
      expect_identical(filter, "ensembl_gene_id")
      expect_identical(dataset, "hsapiens_gene_ensembl")
      data.frame(
        ensembl_gene_id = c("ENSG1", "ENSG1", "ENSG2", "ENSG2"),
        transcript_length = c(1000, 2000, 600, 800),
        stringsAsFactors = FALSE
      )
    },
    .package = "bulkMAE"
  )

  expect_warning(
    lengths <- annotate_gene_lengths(
      c("ENSG1.2", "ENSG2.4", "ENSG3.1"),
      summary = "median"
    ),
    "1 feature"
  )

  expect_identical(rownames(lengths), c("ENSG1.2", "ENSG2.4", "ENSG3.1"))
  expect_equal(lengths$gene_length, c(1500, 700, NA_real_))
  expect_equal(lengths$n_transcripts, c(2L, 2L, NA_integer_))
  expect_identical(attr(lengths, "summary"), "median")
})

test_that("gene length retrieval validates identifiers and policies locally", {
  expect_error(annotate_gene_lengths(c("ENSG1", "ENSG1")), "unique")
  expect_error(annotate_gene_lengths("ENSG1", summary = "guess"), "arg")
  expect_error(annotate_gene_lengths("ENSG1", strip_version = NA), "TRUE or FALSE")
})

test_that("TPM automatically consumes a tximport-style length assay", {
  counts <- matrix(
    c(10L, 20L, 30L, 40L, 50L, 60L),
    nrow = 3L,
    dimnames = list(paste0("gene", 1:3), c("sample1", "sample2"))
  )
  effective_length <- matrix(
    c(1000, 2000, 500, 900, 1800, 600),
    nrow = 3L,
    dimnames = dimnames(counts)
  )
  samples <- data.frame(group = c("a", "b"), row.names = colnames(counts))
  leaf <- mae_create_experiment(
    list(counts = counts, length = effective_length),
    samples
  )
  mae <- mae_create(list(rna = leaf), samples)

  automatic <- normalize_tpm(mae, "rna")
  explicit <- normalize_tpm(mae, "rna", lengths = "length")

  expect_equal(automatic, explicit)
  expect_equal(colSums(automatic), c(sample1 = 1e6, sample2 = 1e6))
  expect_error(
    normalize_tpm(mae_simulate(40L, 8L), "rna"),
    "lengths.*required"
  )
})
