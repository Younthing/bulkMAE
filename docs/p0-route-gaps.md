# P0 route gaps

Inventory of existing `bulkMAE` adapters against the three P0 routes in the
direction doc: make those routes runnable and give them Nature/Cell-grade
default plots. Do not stack new methods.

Scan date: 2026-09-12. Sources: `NAMESPACE`, `R/`, `vignettes/`,
`tests/testthat/`, `README.md`, `inst/guides/`.

## Existing inventory

### Co-expression (`coexpr_*`)

| Function | File | Role |
|---|---|---|
| `transform_vst()`, `transform_rlog()`, `transform_voom()` | `R/preprocessing.R` | Continuous assay for WGCNA |
| `coexpr_pick_power()` | `R/network.R` | `WGCNA::pickSoftThreshold()` + `goodSamplesGenes()` |
| `coexpr_wgcna()` | `R/network.R` | `WGCNA::blockwiseModules()`; stores `colors`, `MEs`, QC |
| `coexpr_modules()` | `R/result-inputs.R` | Named feature → module labels for preservation |
| `coexpr_preservation()` | `R/network.R` | Two-MAE `WGCNA::modulePreservation()` |
| `coexpr_differential()` | `R/network.R` | Two-group DCA via `diffcoexp` (not this P0 path) |

WGCNA is optional (`Suggests`). Tests skip when it is missing. This Cloud
image does not install it.

### Molecular subtyping (`cluster_*`, feature selection, scores)

| Function | File | Role |
|---|---|---|
| `.top_variable_features()` | `R/qc.R` | Internal HVG subset used by PCA/UMAP/WGCNA |
| `mae_subset_features()` | `R/access.R` | Row-subset one experiment by explicit IDs |
| `mae_add_experiment()` | `R/workflow-inputs.R` | Add a GSVA/eigengene matrix as a new leaf |
| `score_gsva()`, `score_ssgsea()`, `score_singscore()` | `R/enrichment.R` | Sample-level gene-set scores |
| `cluster_consensus()` | `R/network.R` | `ConsensusClusterPlus`; drops constant features |
| `cluster_consensus_diagnostics()` | `R/network.R` | Native `calcICL()` |
| `cluster_consensus_classes()` | `R/network.R` | Named `consensusClass` for one `k` |
| `cluster_nmf()`, `cluster_nmf_classes()` | `R/network.R`, `R/result-inputs.R` | Alternate factorization path |
| `mae_add_sample_data()` | `R/workflow-inputs.R` | Attach subtype labels to a leaf |
| `de_*`, `de_table()`, `de_selected()` | `R/differential.R` | Subtype characterization |
| `reduce_umap()`, `reduce_tsne()`, `plot_embedding()` | `R/qc.R`, `R/plotting.R` | Exploratory only |

ConsensusClusterPlus, NMF, and GSVA are optional. No public HVG name helper.

### Immune / TME deconvolution (`deconv_*`)

| Function | File | Role |
|---|---|---|
| `deconv_reference()` | `R/deconvolution.R` | User-supplied scRNA → canonical SCE |
| `deconv()` | `R/deconvolution.R` | `immunedeconv::deconvolute()` (EPIC, quanTIseq, xCell, …) |
| `deconv_music()` | `R/deconvolution.R` | `MuSiC::music_prop()` |
| `deconv_bayesprism()` | `R/deconvolution.R` | `BayesPrism::new.prism()` + `run.prism()` |
| `deconv_cibersortx_input()` | `R/deconvolution.R` | Local mixture table only; no upload |

No extractor that turns those native objects into a cell × sample fraction
matrix. MuSiC / immunedeconv / BayesPrism are GitHub-or-heavy Suggests and
are not in the default CI image.

### Plot helpers already present

`theme_bulkmae()`, `plot_save()`, QC/embedding/DE/assay/GSEA/ORA constructors
in `R/plotting.R`, `R/plotting-gsea.R`, `R/plotting-ora.R`. Contract: one
unprinted `ggplot`, no re-analysis, IDs joined by name, recommended cm size
on `bulkmae_dimensions`.

**None** of the P0 required plot types exist. Closest reuse:
`plot_assay_heatmap(..., column_split = subtype)` for a subtype feature
heatmap; `plot_qc_correlation()` is a sample–sample matrix but is correlation,
not consensus.

### Vignettes

`getting-started`, `airway-qc-de`, `enrichment-analysis`, `naming-migration`.
No co-expression, subtyping, or deconvolution article.

## Gap table

| Route | Existing functions | Bridge missing | Default plot fn missing | Suggested smallest fix |
|---|---|---|---|---|
| Co-expression modules | `transform_*` → `coexpr_pick_power` → `coexpr_wgcna` → `coexpr_modules` → optional `coexpr_preservation` | Module–trait correlation from `fit$MEs`; feature–module membership (kME); hub table | Soft-threshold / module-size colours; module–trait heatmap; MM vs GS scatter (prefer over a hub network) | Add `coexpr_module_trait()`, `coexpr_membership()`, `coexpr_hubs()` as extractors on the native WGCNA list. Add `plot_coexpr_power()`, `plot_coexpr_modules()`, `plot_coexpr_trait()`, `plot_coexpr_membership()`. |
| Molecular subtyping | Internal HVG; `mae_subset_features`; `score_gsva`; `cluster_consensus` / `cluster_nmf` + class extractors; `mae_add_sample_data`; `de_*`; `plot_assay_heatmap`; `plot_embedding` | Public HVG name vector so HVG → subset → cluster can be chained without copying `.top_variable_features()` | Consensus heatmap; consensus CDF over `k`; subtype sample counts (optionally split by a clinical column). Feature heatmap: reuse `plot_assay_heatmap` | Add `mae_variable_features()`. Plots: `plot_cluster_consensus()`, `plot_cluster_cdf()`, `plot_cluster_sizes()`. Do not add a new heatmap API. UMAP/t-SNE stay exploratory. |
| Immune / TME deconvolution | `deconv_reference`; `deconv` / `deconv_music` / `deconv_bayesprism`; `deconv_cibersortx_input` (prepare only) | Cell-type × sample fraction matrix from native immunedeconv / MuSiC / matrix results | Stacked fractions; group boxplots; subtype-split heatmap | Add `deconv_fractions()`. Plots: `plot_deconv_stacked()`, `plot_deconv_box()`, `plot_deconv_heatmap()`. Pin the reference in the vignette; do not default to CIBERSORTx upload. |

## What this PR will implement

Highest-leverage missing pieces so **one path per P0** can be chained in a
vignette stub. Prefer extractors and plots over new analysis backends.

1. `mae_variable_features()` — public HVG names (reuses `.top_variable_features()`).
2. `coexpr_module_trait()`, `coexpr_membership()`, `coexpr_hubs()` — WGCNA
   list / MAE extractors; Pearson `cor` / `cor.test`, not a new S3 class.
3. `deconv_fractions()` — cell × sample matrix from immunedeconv tables,
   MuSiC `$Est.prop.weighted`, or an already-aligned matrix.
4. Default `plot_*()` helpers above, all using `theme_bulkmae()` and
   `bulkmae_dimensions`.
5. Three executable vignette stubs on `mae_simulate()` objects. Backend
   chunks run only when the Suggests package is installed; plot helpers
   always run on constructed tables.

## What should NOT be added (reuse-first)

- MEGENA, CEMiTool, or other co-expression shells.
- PPI / molecular-dynamics wrappers (`network_string()` already exists for
  STRING and is out of this P0 path).
- New enrichment database wrappers (GO/KEGG/Reactome/MSigDB already exist).
- P1 activity (`activity_*`), survival (`surv_*`), or LINCS (`drug_*`)
  unless a one-line existing call is needed as a characterization bridge.
- A CIBERSORTx upload / remote-execution default. Keep
  `deconv_cibersortx_input()` as a local table helper only.
- New result classes, `autoplot()` methods, or `theme_set()`.
- A `view =` mega-plot that branches stacked / box / heatmap.
- A second subtype feature heatmap; use `plot_assay_heatmap()`.
- Hub-network / igraph plots; MM vs GS is the chosen hub view.
- Automatic k selection, automatic power selection as a scientific claim,
  or simulated single-cell references presented as biology.
- Survival external-validation claims.
- UMAP/t-SNE as a required subtype figure.

## Next commits after this PR's first slice

See the PR body. Expected follow-ups, not this slice:

1. Optional `top_n` on `cluster_consensus()` so HVG can stay inside one call.
2. Optional BayesPrism `get.fraction()` adapter inside `deconv_fractions()`
   once that backend is on the full-backend image.
3. Real-data (not `mae_simulate()`) tutorials once WGCNA / CCP / a pinned
   reference are in the vignette CI extra-packages.
4. Module-preservation plot helper if a second cohort is in scope.
5. Publication PNG review under ignored `docs/plot-previews/`.
