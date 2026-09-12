# NA

## Agent skills

### Issue tracker

Issues and specs are tracked in GitHub Issues. See
`docs/agents/issue-tracker.md`.

### Domain docs

This is a single-context repository. See `docs/agents/domain.md`.

### Sequence an analysis mainline

Do not start from a lone plot or a lone wrapper. Work one independent
mainline in this order.

1.  Name the mainline. One research question. One vignette. A sequel
    reuses the upstream design, filter, and contrast.
2.  The mainline names the analyses. Draft the vignette next. For each
    chosen analysis, list the figures that analysis commonly carries in
    the literature. That list is the required figure set, not a subset
    the prose can skip. The first draft may fail to render. Treat it as
    the list of missing wrappers and plots. Writing rules live in
    `docs/agents/vignette-writing.md`.
3.  If a chosen analysis has no public function, add a thin wrapper
    only. Do not re-implement the backend. Do not change the backend
    return type.
4.  Draw the full literature figure set for each chosen analysis. Each
    constructor returns one unprinted ggplot from already-computed
    results. Bring every figure in that set to publication quality. Plot
    contracts live in `docs/agents/plotting-functions.md`.
5.  Finish the vignette so every number and figure comes from executed
    code.

Do not pick a plot from another analysis family. Do not drop a standard
figure because the vignette prose can live without it.

### Engineering guides

- For the analysis-mainline order above, stay in this file.
- For plotting-interface or publication-output changes, read
  `docs/agents/plotting-functions.md`.
- For executable vignette, rendering, or article-CI changes, read
  `docs/agents/vignette-writing.md`.

## Cursor Cloud specific instructions

The Cloud Agent environment (`.cursor/environment.json` +
`.cursor/install.sh`) installs R, pandoc, the package hard dependencies,
and the core offline backends (`airway`, `edgeR`, `DESeq2`, `limma`,
`clusterProfiler`, `fgsea`, `org.Hs.eg.db`). The heavy GitHub-only
backends from `full-backend-check` (`MuSiC`, `immunedeconv`,
`BayesPrism`, `WGCNA`, …) are not installed; the default suite skips
them cleanly.

Run the test suite sequentially. The package sets
`Config/testthat/parallel: true`, but `testthat`’s parallel worker
startup fails on this toolchain, so disable it (as `full-backend-check`
does):

``` sh
TESTTHAT_PARALLEL=false BULKMAE_RUN_ONLINE_TESTS=false \
  Rscript -e 'testthat::test_local()'
```
