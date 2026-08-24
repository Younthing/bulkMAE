# Deconvolve bulk counts with BayesPrism

Both bulk and single-cell reference inputs must be raw integer counts.
The MAE assay is unambiguously gene-by-sample and is converted to
BayesPrism's sample-by-gene layout. `reference_orientation` explicitly
declares the supplied reference layout; no orientation is inferred from
dimensions.

## Usage

``` r
deconv_bayesprism(
  x,
  experiment,
  reference,
  cell_type_labels = NULL,
  assay = "counts",
  reference_orientation = c("genes_by_cells", "cells_by_genes"),
  cell_state_labels = NULL,
  key = NULL,
  input_type = "count.matrix",
  cores = 1,
  update_gibbs = TRUE,
  gibbs_control = list(),
  opt_control = list(),
  ...,
  reference_assay = NULL
)
```

## Arguments

- x:

  A `MultiAssayExperiment` object.

- experiment:

  Optional experiment name.

- reference:

  Raw reference count matrix, or a `SingleCellExperiment` created by
  [`deconv_reference()`](https://younthing.github.io/bulkMAE/reference/deconv_reference.md).

- cell_type_labels:

  Cell-type label for each reference cell. For a `SingleCellExperiment`,
  `NULL` uses its `cell_type` column and a single string names another
  `colData` column.

- assay:

  Assay name or one-based assay index.

- reference_orientation:

  Layout of `reference`. The default preserves the package's documented
  gene-by-cell input, while `cells_by_genes` accepts BayesPrism's native
  layout directly.

- cell_state_labels:

  Optional cell-state label for each reference cell. For a
  `SingleCellExperiment`, `NULL` uses `cell_state` when that column is
  present and a single string names another `colData` column.

- key:

  Optional malignant/tumour key accepted by BayesPrism.

- input_type:

  BayesPrism reference input type.

- cores:

  Number of worker cores passed to `run.prism()`.

- update_gibbs:

  Run the final Gibbs update.

- gibbs_control, opt_control:

  Named control lists passed to `BayesPrism::run.prism()`.

- ...:

  Additional arguments passed to `BayesPrism::new.prism()`.

- reference_assay:

  Raw-count assay used when `reference` is a `SingleCellExperiment`.
  `NULL` selects `counts`, or the sole assay when no `counts` assay is
  present.

## Value

The native object returned by `BayesPrism::run.prism()`.
