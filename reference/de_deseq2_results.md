# Extract DESeq2 results

Extract DESeq2 results

## Usage

``` r
de_deseq2_results(
  fit,
  contrast = NULL,
  name = NULL,
  alpha = 0.05,
  independent_filtering = TRUE,
  ...
)
```

## Arguments

- fit:

  A fitted `DESeqDataSet`, usually from
  [`de_deseq2()`](https://younthing.github.io/bulkMAE/reference/de_deseq2.md).

- contrast:

  Optional DESeq2 contrast.

- name:

  Optional coefficient name. Do not supply together with `contrast`.

- alpha:

  Independent-filtering false discovery rate threshold.

- independent_filtering:

  Use DESeq2 independent filtering.

- ...:

  Additional arguments passed to `DESeq2::results()`.

## Value

A native `DESeqResults` object.
