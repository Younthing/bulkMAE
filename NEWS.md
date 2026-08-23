# bulkMAE 0.4.0

- Replaces mixed `run_*`, `infer_*`, `prepare_*`, and backend-specific naming
  with discoverable `<family>_<method>()` / `<family>_<operation>()` APIs.
- Organizes public functions under `mae_*`, `qc_*`, `reduce_*`, `de_*`,
  `enrich_*`, `score_*`, `activity_*`, `cluster_*`, `coexpr_*`, `network_*`,
  `deconv_*`, `surv_*`, `ml_*`, `meta_*`, and `drug_*` families.
- Consolidates GO, KEGG, and Reactome ORA/GSEA combinations behind explicit
  `method` arguments instead of multiplying function names.
- Removes old exported aliases intentionally so autocomplete exposes one
  coherent API; an installed migration map records every breaking rename.

# bulkMAE 0.3.0

- Enforces name-aware sample, design, feature, and annotation alignment across
  stateless adapters, with early failures for missing covariates and aliased
  fixed-effect designs.
- Completes the edgeR quasi-likelihood workflow with dispersion estimation and
  repairs logical-vector handling in goseq.
- Strengthens maSigPro, differential-splicing, DCA, consensus clustering,
  WGCNA, module-preservation, embedding, activity, clinical, and deconvolution
  input contracts.
- Keeps WGCNA's automatic bad-sample/bad-gene filtering while reporting and
  retaining the quality decision in its native result list.
- Localizes wrapper-controlled random seeds so calls do not mutate the user's
  global random-number stream.
- Adds executable local examples, conditional backend tests, installed audit
  guides, and package-level help needed for a standard source build.

# bulkMAE 0.2.0

- Adds a general decoupleR adapter covering enrichment-style and
  network-aware statistics, resource discovery, and consensus scoring.
- Adds KEGG/Reactome ORA and GSEA, CAMERA/FRY/mroast, and singscore.
- Adds consensus-clustering diagnostics, differential co-expression, and
  WGCNA module preservation.
- Adds continuous-expression limma, differential splicing, and maSigPro
  time-course modelling.
- Adds UMAP/t-SNE, limma batch removal, Kaplan-Meier, and time-dependent ROC.
- Expands the tutorial and official-source audit map for every new adapter.

# bulkMAE 0.1.0

- Establishes the stateless MAE-to-native-backend contract.
- Adds import, annotation, QC, preprocessing, differential expression,
  enrichment, network, deconvolution, and clinical adapters.
- Adds a complete walkthrough, official-source audit map, and review checklist.
