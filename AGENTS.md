## Agent skills

### Issue tracker

Issues and specs are tracked in GitHub Issues. See `docs/agents/issue-tracker.md`.

### Domain docs

This is a single-context repository. See `docs/agents/domain.md`.

### Sequence an analysis mainline

Do not start from a lone plot or a lone wrapper. Work one independent
mainline in this order.

1. Name the mainline. One research question. One vignette. A sequel reuses
   the upstream design, filter, and contrast.
2. Draft the vignette first. The draft lists the analyses and the figures
   the story needs. The first draft may fail to render. Treat it as the
   list of missing wrappers and plots. Writing rules live in
   `docs/agents/vignette-writing.md`.
3. If a needed analysis has no public function, add a thin wrapper only.
   Do not re-implement the backend. Do not change the backend return type.
4. Draw only the figures the draft named. Each constructor returns one
   unprinted ggplot from already-computed results. Plot contracts live in
   `docs/agents/plotting-functions.md`.
5. Finish the vignette so every number and figure comes from executed code.

Do not add a plot because papers often use that plot type. Add it because
the vignette question needs it.

### Engineering guides

- For the analysis-mainline order above, stay in this file.
- For plotting-interface or publication-output changes, read
  `docs/agents/plotting-functions.md`.
- For executable vignette, rendering, or article-CI changes, read
  `docs/agents/vignette-writing.md`.

## Cursor Cloud specific instructions

The Cloud Agent environment (`.cursor/environment.json` + `.cursor/install.sh`)
installs R, pandoc, the package hard dependencies, and the core offline
backends (`airway`, `edgeR`, `DESeq2`, `limma`, `clusterProfiler`, `fgsea`,
`org.Hs.eg.db`). The heavy GitHub-only backends from `full-backend-check`
(`MuSiC`, `immunedeconv`, `BayesPrism`, `WGCNA`, ...) are not installed; the
default suite skips them cleanly.

Run the test suite sequentially. The package sets
`Config/testthat/parallel: true`, but `testthat`'s parallel worker startup
fails on this toolchain, so disable it (as `full-backend-check` does):

```sh
TESTTHAT_PARALLEL=false BULKMAE_RUN_ONLINE_TESTS=false \
  Rscript -e 'testthat::test_local()'
```
