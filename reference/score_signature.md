# Score a gene-expression signature

Score a gene-expression signature

## Usage

``` r
score_signature(
  x,
  experiment,
  weights,
  assay,
  center = FALSE,
  scale = FALSE,
  na_rm = FALSE
)
```

## Arguments

- x:

  A `MultiAssayExperiment` object.

- experiment:

  Optional experiment name.

- weights:

  Named numeric feature weights, or an unweighted character vector of
  feature names.

- assay:

  Assay name or one-based assay index.

- center, scale:

  Either a logical value, or a named numeric vector of training-set
  feature centers/scales. `TRUE` estimates values from the current
  cohort and is therefore unsuitable for locked external validation;
  supply frozen training values instead.

- na_rm:

  Remove missing values when summing scores.

## Value

A named numeric vector with one score per sample.
