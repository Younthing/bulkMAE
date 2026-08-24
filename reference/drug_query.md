# Prepare an up/down query for LINCS or CMap

Prepare an up/down query for LINCS or CMap

## Usage

``` r
drug_query(statistics, n = 150)
```

## Arguments

- statistics:

  Named numeric differential-expression statistics.

- n:

  Maximum number of genes selected separately from the positive and
  negative statistics.

## Value

A list with `upset` and `downset` gene identifiers.

## Examples

``` r
statistic <- c(gene1 = 3, gene2 = 2, gene3 = -1, gene4 = -4)
drug_query(statistic, n = 2)
#> $upset
#> [1] "gene1" "gene2"
#> 
#> $downset
#> [1] "gene4" "gene3"
#> 
```
