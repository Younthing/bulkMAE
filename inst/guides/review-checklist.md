# Review checklist

Use this list before accepting a new wrapper or publishing an analysis.

## Package architecture

- [ ] Public assay-level function starts from an MAE and names `experiment`.
- [ ] Public name follows `<family>_<method>()` or
      `<family>_<operation>()`; do not introduce a competing `run_*`,
      `infer_*`, or `prepare_*` alias.
- [ ] Assay access goes through `getWithColData()` and preserves sample order.
- [ ] Function does not modify the input MAE or store run status/history.
- [ ] Return value is the backend's native object where practical.
- [ ] Optional backend is checked only when its wrapper is called.
- [ ] New dependency is documented in `DESCRIPTION` or, for a non-standard
      repository, in the dependency policy.
- [ ] Roxygen comments state input scale, orientation, and return class.
- [ ] Implementation is mapped to an official source in
      `inst/guides/official-sources.md`.

## Data and design

- [ ] Raw counts are non-negative; integer-only backends receive integers.
- [ ] Feature and sample identifiers are unique where required.
- [ ] Reference levels and coefficient names have been inspected.
- [ ] Paired/repeated designs include subject effects explicitly.
- [ ] Batch correction is not substituted for design adjustment without a
      documented reason.
- [ ] SVA/RUV factors are included explicitly in the downstream design.
- [ ] ORA universe is the set of testable genes after filtering.
- [ ] GO/KEGG/Reactome calls name `genes` or `ranks` and set `method`
      explicitly, rather than relying on positional matching.
- [ ] Gene ID type and annotation/database release are recorded.
- [ ] decoupleR method family matches the question: gene-set statistic versus
      signed/weighted regulator-target activity inference.
- [ ] Consensus-clustering `k` is supported by stability and external evidence,
      not selected from one index alone.
- [ ] Differential co-expression uses an outcome-independent feature filter
      and reports all correlation and FDR thresholds.
- [ ] Module preservation keeps reference/test roles and reference module
      labels explicit.
- [ ] WGCNA warnings and returned audit fields identify every sample or feature
      removed by `goodSamplesGenes()`.

## Validation and translation

- [ ] Feature selection occurs inside resampling folds.
- [ ] Hyperparameter tuning and performance evaluation use distinct data.
- [ ] External validation keeps the original model coefficients and cutoff.
- [ ] Survival models check proportional-hazards and calibration assumptions.
- [ ] Deconvolution reference, platform, tissue, and cell labels are recorded.
- [ ] LINCS results are interpreted with cell line, dose, time, and toxicity.
- [ ] Random seeds and backend versions are recorded in the analysis project.

## Release checks

For every release candidate, create a clean package library and install the
declared backends using
[`dependency-installation-zh.md`](dependency-installation-zh.md). Preserve the
installation plan, package versions, session information, and check logs as
release artifacts. This checklist defines required verification; it does not
claim that a particular release candidate has passed until those artifacts are
recorded.

- [ ] Use a mutually compatible R/Bioconductor pair that satisfies every
      minimum version in `DESCRIPTION`; do not reuse a mixed-release library.

- [ ] Replace the provisional maintainer identity in `DESCRIPTION`.
- [ ] Add the project URL and bug tracker after creating the repository.

```r
# Regenerate one standard Rd topic per public function from the roxygen source.
devtools::document("bulkMAE")
devtools::test("bulkMAE")
devtools::check("bulkMAE", args = "--as-cran")
BiocCheck::BiocCheck("bulkMAE")
```

Also exercise one small official-example dataset for every optional backend.
In particular, check API compatibility for edgeR `normLibSizes()`, the GSVA
parameter constructors, decoupleR resource column names, dream contrasts, and
the three non-standard deconvolution packages. Also verify the installed
versions of ConsensusClusterPlus `calcICL()`, diffcoexp, maSigPro, ReactomePA,
singscore, timeROC, uwot, Rtsne, and WGCNA `modulePreservation()`.

Network-, credential-, license-, or externally hosted data-dependent checks
must be reported separately from the default offline suite. Record whether each
such check ran, skipped, or failed; do not treat an unrecorded skip as a pass.
