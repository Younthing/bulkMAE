# Test differential co-expression between two groups

Uses the maintained Bioconductor `diffcoexp` implementation to identify
differential links and genes from two condition-specific expression
matrices. DCA here means differential co-expression analysis.

## Usage

``` r
coexpr_differential(
  x,
  experiment,
  group,
  contrast,
  assay,
  features = NULL,
  correlation = "pearson",
  p_adjust = "BH",
  correlation_threshold = 0.5,
  correlation_fdr = 0.1,
  difference_threshold = 0.5,
  difference_fdr = 0.1,
  gene_fdr = 0.1,
  max_pairs = 5e+06
)
```

## Arguments

- x:

  A `MultiAssayExperiment` object.

- experiment:

  Optional experiment name.

- group:

  Metadata column or grouping vector.

- contrast:

  Two group labels, in the order passed to `exprs.1` and `exprs.2`.

- assay:

  Assay name or one-based assay index.

- features:

  Optional feature subset.

- correlation:

  Correlation method.

- p_adjust:

  Multiple-testing correction method.

- correlation_threshold, correlation_fdr:

  Thresholds for condition-specific correlations.

- difference_threshold, difference_fdr:

  Thresholds for changes in correlation.

- gene_fdr:

  FDR threshold for differential co-expression genes.

- max_pairs:

  Maximum number of feature pairs to analyze. Set to `Inf` only after
  considering the quadratic memory and runtime cost.

## Value

The native result returned by `diffcoexp::diffcoexp()`.

## Details

Use a normalized, variance-stabilized or log-expression assay. Remove
known unwanted variation before this analysis when appropriate.
Correlations are estimated separately, so small groups are unstable even
when they meet the minimum of four samples; groups below ten emit a
warning.
