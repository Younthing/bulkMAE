# Run selection-bias-aware enrichment with goseq

Run selection-bias-aware enrichment with goseq

## Usage

``` r
enrich_goseq(
  selected,
  genome,
  id,
  bias = NULL,
  gene_to_category = NULL,
  categories = c("GO:CC", "GO:BP", "GO:MF"),
  method = "Wallenius",
  repetitions = 2000,
  include_uncategorized = FALSE,
  plot_fit = FALSE
)
```

## Arguments

- selected:

  Named binary or logical vector over all tested genes.

- genome:

  Genome identifier accepted by `goseq::nullp()`.

- id:

  Gene identifier type accepted by `goseq::nullp()`.

- bias:

  Optional named or position-aligned numeric bias, usually gene length
  or expression abundance.

- gene_to_category:

  Optional gene-to-category mapping.

- categories:

  GO/KEGG categories passed as `test.cats`.

- method:

  Enrichment method. `Wallenius` is the backend recommendation.

- repetitions:

  Sampling repetitions, used only by the sampling method.

- include_uncategorized:

  Include genes without category annotations.

- plot_fit:

  Plot the probability-weighting function for review.

## Value

The native enrichment data frame returned by `goseq::goseq()`.
