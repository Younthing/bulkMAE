# Construct limma-compatible contrasts from MAE sample metadata

This helper keeps design-column construction and contrast parsing inside
`bulkMAE`. The returned matrix can be passed directly to
[`de_voom()`](https://younthing.github.io/bulkMAE/reference/de_voom.md)
or
[`de_limma()`](https://younthing.github.io/bulkMAE/reference/de_limma.md).

## Usage

``` r
de_contrast(x, experiment, formula, contrasts)
```

## Arguments

- x:

  A `MultiAssayExperiment` object.

- experiment:

  Optional experiment name.

- formula:

  Model formula evaluated against aligned sample metadata.

- contrasts:

  Character contrast expressions accepted by `limma::makeContrasts()`.

## Value

A numeric contrast matrix whose rows match
[`de_design()`](https://younthing.github.io/bulkMAE/reference/de_design.md)
columns.
