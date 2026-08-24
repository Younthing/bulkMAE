# Simulate a self-contained bulk-transcriptomics MAE

This function is intended for examples, teaching, backend smoke tests,
and reproducible bug reports. Its synthetic annotations and outcomes
have no biological meaning and must not be used as reference data for
real studies.

## Usage

``` r
mae_simulate(n_features = 300L, n_samples = 12L, seed = 1L, experiment = "rna")
```

## Arguments

- n_features:

  Number of simulated expression features.

- n_samples:

  Even number of simulated samples.

- seed:

  Local random seed; the caller's random-number state is restored.

- experiment:

  Experiment name in the returned MAE.

## Value

A
[MultiAssayExperiment::MultiAssayExperiment](https://github.com/waldronlab/MultiAssayExperiment/reference/MultiAssayExperiment.html)
containing integer `counts`, `log_expression`, and `tpm` assays; feature
lengths and synthetic identifiers in row metadata; and common design,
prediction, and survival variables in sample metadata.

## Examples

``` r
example_mae <- mae_simulate(n_features = 100, n_samples = 8, seed = 1)
mae_assays(example_mae, "rna")
#> [1] "counts"         "log_expression" "tpm"           
```
