## Agent skills

### Issue tracker

Issues and specs are tracked in GitHub Issues. See `docs/agents/issue-tracker.md`.

### Domain docs

This is a single-context repository. See `docs/agents/domain.md`.

### Engineering guides

- For plotting-interface or publication-output changes, read
  `docs/agents/plotting-functions.md` first.
- For executable vignette, rendering, or article-CI changes, read
  `docs/agents/vignette-writing.md` first.

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
