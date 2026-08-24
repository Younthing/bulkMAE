# Calculate time-dependent ROC curves

Calculate time-dependent ROC curves

## Usage

``` r
surv_roc(x, experiment, time, event, marker, times, cause = 1, ...)
```

## Arguments

- x:

  A `MultiAssayExperiment` object.

- experiment:

  Optional experiment name.

- time, event:

  Metadata columns containing follow-up time and event.

- marker:

  Metadata column or numeric marker vector.

- times:

  Evaluation times.

- cause:

  Event code treated as the outcome of interest.

- ...:

  Additional arguments passed to `timeROC::timeROC()`.

## Value

A native `timeROC` object.
