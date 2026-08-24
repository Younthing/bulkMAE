# bulkMAE: stateless bulk-transcriptomics adapters

`bulkMAE` uses `MultiAssayExperiment` as its public data contract,
resolves one `SummarizedExperiment` leaf and its aligned sample
metadata, then returns the native result produced by an established
analysis backend. It does not mutate the input object or keep an
analysis-status registry.

## Details

Start with
[`mae_create()`](https://younthing.github.io/bulkMAE/reference/mae_create.md)
or
[`import_tximport()`](https://younthing.github.io/bulkMAE/reference/import_tximport.md),
inspect alignment with
[`mae_validate()`](https://younthing.github.io/bulkMAE/reference/mae_validate.md)
and
[`mae_samples()`](https://younthing.github.io/bulkMAE/reference/mae_samples.md),
and then call an analysis family such as `qc_*`, `de_*`, `enrich_*`,
`coexpr_*`, or `surv_*`. Statistical engines are optional dependencies
and are checked only when their adapter is called.

## See also

[`vignette("getting-started", package = "bulkMAE")`](https://younthing.github.io/bulkMAE/articles/getting-started.md)

## Author

**Maintainer**: Younthing <fanxingfu3344@gmail.com>
