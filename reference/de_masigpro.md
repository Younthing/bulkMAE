# Fit a maSigPro time-course or dose-response model

Follows the official maSigPro sequence: `make.design.matrix()`,
`p.vector()`, `T.fit()`, and `get.siggenes()`.

## Usage

``` r
de_masigpro(
  x,
  experiment,
  edesign,
  assay,
  degree = 2,
  q = 0.05,
  p_adjust = "BH",
  min_observations = NULL,
  counts = FALSE,
  step_method = "backward",
  alpha = 0.05,
  r_squared = 0.6,
  variables = "groups",
  ...
)
```

## Arguments

- x:

  A `MultiAssayExperiment` object.

- experiment:

  Optional experiment name.

- edesign:

  Experimental design data frame required by maSigPro. Rows must be
  row-named by assay sample ID and include `Time`, `Replicates`, and at
  least one binary group-indicator column. Rows are reordered to assay
  sample order before fitting.

- assay:

  Assay name or one-based assay index.

- degree:

  Polynomial degree for the regression model.

- q, p_adjust:

  Significance threshold and adjustment method for
  `maSigPro::p.vector()`.

- min_observations:

  Optional minimum observations per gene.

- counts:

  Whether the assay contains counts.

- step_method:

  Variable-selection method used by `maSigPro::T.fit()`.

- alpha:

  Significance threshold for variable selection.

- r_squared:

  Minimum model R-squared.

- variables:

  maSigPro variable grouping passed as `vars`.

- ...:

  Additional arguments passed to `maSigPro::p.vector()`.

## Value

The native object returned by `maSigPro::get.siggenes()`.
