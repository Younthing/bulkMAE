# Construct a single-cell deconvolution reference

Builds the canonical `SingleCellExperiment` reference accepted directly
by
[`deconv_music()`](https://younthing.github.io/bulkMAE/reference/deconv_music.md)
and
[`deconv_bayesprism()`](https://younthing.github.io/bulkMAE/reference/deconv_bayesprism.md).
This helper only standardizes a user-supplied single-cell count matrix
and its metadata; it does not provide or simulate a biological
reference.

## Usage

``` r
deconv_reference(
  counts,
  cell_data,
  cell_type,
  sample,
  cell_state = NULL,
  assay_name = "counts"
)
```

## Arguments

- counts:

  Raw integer count matrix with genes in rows and cells in columns.
  Dense matrices and sparse `Matrix` objects are accepted without
  densifying the latter. Complete, unique gene and cell names are
  required.

- cell_data:

  Cell-level metadata with row names containing every count matrix cell
  name exactly once. Rows are aligned by cell name.

- cell_type, sample:

  A column name in `cell_data`, or a named vector with one non-missing
  label per cell. Values are stored in the canonical `cell_type` and
  `sample_id` columns.

- cell_state:

  Optional cell-state column name or named vector. When supplied, values
  are stored in the canonical `cell_state` column.

- assay_name:

  Name used for the raw-count assay.

## Value

A `SingleCellExperiment` with aligned cell metadata.
