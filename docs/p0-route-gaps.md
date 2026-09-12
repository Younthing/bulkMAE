# P0 route gaps

Inventory of existing `bulkMAE` adapters against the three P0 routes in the
direction doc: make those routes runnable and give them Nature/Cell-grade
default plots. Do not stack new methods.

Scan date: 2026-09-12. Updated in this PR after the airway implementation.

## Existing inventory

### Co-expression (`coexpr_*`)

| Function | File | Role |
|---|---|---|
| `transform_vst()`, `transform_rlog()`, `transform_voom()` | `R/preprocessing.R` | Continuous assay for WGCNA |
| `coexpr_pick_power()` | `R/network.R` | `WGCNA::pickSoftThreshold()` + `goodSamplesGenes()` |
| `coexpr_wgcna()` | `R/network.R` | `WGCNA::blockwiseModules()`; stores `colors`, `MEs`, QC |
| `coexpr_modules()` | `R/result-inputs.R` | Named feature → module labels for preservation |
| `coexpr_module_trait()` | `R/result-inputs.R` | Eigengene–trait Pearson `cor` / `cor.test` |
| `coexpr_membership()`, `coexpr_hubs()` | `R/result-inputs.R` | kME matrix and intramodular hubs |
| `coexpr_preservation()` | `R/network.R` | Two-MAE `WGCNA::modulePreservation()` |
| `coexpr_differential()` | `R/network.R` | Two-group DCA via `diffcoexp` (not this P0 path) |

### Molecular subtyping (`cluster_*`, feature selection, scores)

| Function | File | Role |
|---|---|---|
| `mae_variable_features()` | `R/access.R` | Public HVG names (row variance) |
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

### Immune / TME deconvolution (`deconv_*`)

| Function | File | Role |
|---|---|---|
| `deconv_reference()` | `R/deconvolution.R` | User-supplied scRNA → canonical SCE |
| `mae_simulate()` | `R/simulate.R` | Existing synthetic bulk MAE when a real atlas is absent |
| `deconv()` | `R/deconvolution.R` | `immunedeconv::deconvolute()` |
| `deconv_music()` | `R/deconvolution.R` | `MuSiC::music_prop()` |
| `deconv_bayesprism()` | `R/deconvolution.R` | `BayesPrism::new.prism()` + `run.prism()` |
| `deconv_fractions()` | `R/deconvolution.R` | Cell-type × sample matrix from native results |
| `deconv_cibersortx_input()` | `R/deconvolution.R` | Local mixture table only; no upload |

### Default P0 plots (this PR)

All return one unprinted `ggplot` with `theme_bulkmae()` and
`bulkmae_dimensions`. They do not re-run analysis.

| Route | Helpers |
|---|---|
| Co-expression | `plot_coexpr_power()`, `plot_coexpr_modules()`, `plot_coexpr_trait()`, `plot_coexpr_membership()`, `plot_coexpr_dendrogram()`, `plot_coexpr_tom()` |
| Subtyping | `plot_cluster_consensus()`, `plot_cluster_cdf()`, `plot_cluster_delta()`, `plot_cluster_pac()`, `plot_cluster_sizes()`; feature heatmap reuses `plot_assay_heatmap()` |
| Deconvolution | `plot_deconv_stacked()`, `plot_deconv_box()`, `plot_deconv_heatmap()` |

## Gap table (after this PR)

| Route | Existing functions | Bridge missing | Default plot fn missing | Suggested smallest fix |
|---|---|---|---|---|
| Co-expression modules | `transform_vst` → `mae_variable_features` → `coexpr_pick_power` → `coexpr_wgcna` → `coexpr_module_trait` / `coexpr_membership` / `coexpr_hubs` → optional `coexpr_preservation` | None for one airway path. Preservation still needs a second cohort. | Dendrogram/TOM/GS-MM hooks are now `plot_coexpr_dendrogram()`, `plot_coexpr_tom()`, and existing `plot_coexpr_membership()`. | Done on `airway`. n = 8 cannot support a scale-free power claim; the 0.8 line is a visible miss. |
| Molecular subtyping | `mae_variable_features` → `mae_subset_features` → `cluster_consensus` → `cluster_consensus_classes` → `mae_add_sample_data` → `plot_assay_heatmap` | Optional `top_n` inside `cluster_consensus()` is still a convenience, not a blocker. GSVA path uses existing `score_gsva` + `mae_add_experiment`. | Delta area and PAC are `cluster_consensus_delta_area()` / `cluster_consensus_pac()` plus matching plots. | Done on `airway`. k is diagnostic, not a subtype discovery. |
| Immune / TME deconvolution | `mae_simulate` → `deconv_reference` → `deconv_fractions` → plots. Real backends (`deconv` / MuSiC / BayesPrism) remain optional. | BayesPrism `get.fraction()` adapter. No real atlas is bundled. | None of the required types. | `airway` has no scRNA reference, so the vignette uses existing `mae_simulate()` and aligns `condition` / `batch`. No new simulator. No CIBERSORTx upload. |

## What this PR implemented

1. Extractors: `mae_variable_features()`, `coexpr_module_trait()`,
   `coexpr_membership()`, `coexpr_hubs()`, `deconv_fractions()`.
2. Default `plot_*()` helpers listed above, including dendrogram/TOM,
   delta-area/PAC, annotation bars, and deconvolution jitter.
3. Co-expression and subtyping vignettes on Bioconductor `airway` (same
   `filter_expr()` / `transform_vst(~ cell + dex)` construction as the
   QC/DE tutorial). Deconvolution uses existing `mae_simulate()` because
   `airway` has no single-cell reference. WGCNA and ConsensusClusterPlus
   chunks run when those Suggests are installed; CI extra-packages now
   include them plus `SingleCellExperiment`.
4. No new deconvolution simulator. Sample grouping stays aligned to
   `mae_samples()`.

## What should NOT be added (reuse-first)

- MEGENA, CEMiTool, or other co-expression shells.
- PPI / molecular-dynamics wrappers.
- New enrichment database wrappers.
- P1 activity (`activity_*`), survival (`surv_*`), or LINCS (`drug_*`)
  unless a one-line existing call is needed as a characterization bridge.
- A CIBERSORTx upload / remote-execution default.
- New result classes, `autoplot()` methods, or `theme_set()`.
- A `view =` mega-plot that branches stacked / box / heatmap.
- A second subtype feature heatmap; use `plot_assay_heatmap()`.
- Hub-network / igraph plots; MM vs GS is the chosen hub view.
- Automatic k selection, automatic power selection as a scientific claim,
  or a second in-package simulator besides `mae_simulate()`.
- Survival external-validation claims.
- UMAP/t-SNE as a required subtype figure.

## Remaining follow-ups (not blockers)

1. Optional `top_n` on `cluster_consensus()`.
2. BayesPrism `get.fraction()` adapter inside `deconv_fractions()`.
3. Module-preservation plot if a second real cohort is in scope.
4. Figure-audit review of the default plots (maintainer; not merged here).
5. A licensed, tissue-matched single-cell reference if a later article
   must show real TME fractions. Do not treat `mae_simulate()` as that
   reference.
