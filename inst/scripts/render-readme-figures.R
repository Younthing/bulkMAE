# Render the four README gallery plots from public bulkMAE constructors.
# Run from the repository root:
#   R_LIBS_USER="$HOME/R/library" Rscript inst/scripts/render-readme-figures.R

.libPaths(c(Sys.getenv("R_LIBS_USER", unset = .libPaths()[[1L]]), .libPaths()))
if (!requireNamespace("bulkMAE", quietly = TRUE)) {
  stop("Install bulkMAE before rendering README figures.", call. = FALSE)
}

library(bulkMAE)
library(ggplot2)

set.seed(20260910)
out_dir <- file.path("man", "figures")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

readme_save <- function(filename, plot) {
  path <- file.path(out_dir, filename)
  plot_save(path, plot, scale = 2, dpi = 200, bg = "white")
  path
}

n_features <- 480L
feature_ids <- sprintf("GENE%04d", seq_len(n_features))
marker_ids <- c("IL6", "CXCL8", "TNF", "IFNG", "STAT1", "MX1")
feature_ids[seq_along(marker_ids)] <- marker_ids

effect <- stats::rnorm(n_features, mean = 0, sd = 0.7)
up <- seq_len(28L)
down <- 29L:56L
effect[up] <- stats::runif(length(up), 1.15, 3.1)
effect[down] <- stats::runif(length(down), -3.1, -1.15)
z_score <- effect / 0.45 + stats::rnorm(n_features, mean = 0, sd = 0.85)
p_value <- 2 * stats::pnorm(-abs(z_score))
p_value <- pmin(pmax(p_value, 1e-16), 1)
de_table <- data.frame(
  feature_id = feature_ids,
  effect = effect,
  p_value = p_value,
  adjusted_p_value = stats::p.adjust(p_value, method = "BH"),
  row.names = feature_ids,
  check.names = FALSE
)

volcano <- plot_de_volcano(
  de_table,
  fdr = 0.05,
  min_abs_effect = 1,
  label_features = marker_ids
)
readme_save("readme-volcano.png", volcano)

mae <- mae_simulate(n_features = 160L, n_samples = 12L, seed = 20260910L)
log_expr <- mae_pull_assay(mae, "rna", "log_expression")
feature_var <- apply(log_expr, 1L, stats::sd)
heat_features <- names(sort(feature_var, decreasing = TRUE))[seq_len(32L)]
heatmap <- plot_assay_heatmap(
  mae,
  "rna",
  "log_expression",
  features = heat_features,
  column_split = "condition"
)
readme_save("readme-heatmap.png", heatmap)

ora_sets <- list(
  Inflammatory_response = c(
    "IL6", "CXCL8", "TNF", "IL1B", "NFKB1", "RELA", "STAT3", "JUN", "PTGS2"
  ),
  Interferon_alpha = c(
    "STAT1", "MX1", "ISG15", "IRF7", "OAS1", "IFIT1", "STAT2", "IRF1"
  ),
  TNF_signaling = c(
    "TNF", "NFKB1", "JUN", "MAPK8", "BIRC3", "CXCL8", "TNFAIP3", "IL1B", "RELA"
  ),
  IL6_JAK_STAT3 = c(
    "IL6", "STAT3", "SOCS3", "IL6R", "JAK2", "MYC", "JUN"
  ),
  Apoptosis = c(
    "BIRC3", "TNF", "BCL2A1", "CASP8", "FAS", "BID", "JUN"
  ),
  Hypoxia = c(
    "VEGFA", "ENO1", "LDHA", "HK2", "SLC2A1", "JUN", "NFKB1"
  )
)
ora_universe <- unique(unlist(ora_sets, use.names = FALSE))
ora_effects <- stats::setNames(
  stats::rnorm(length(ora_universe), mean = 0.2, sd = 1.1),
  ora_universe
)
ora_effects[c("IL6", "CXCL8", "TNF", "STAT1", "MX1")] <- c(2.4, 2.1, 1.8, 1.6, 1.4)
ora_effects[c("BCL2A1", "SOCS3")] <- c(-1.7, -1.3)
ora_n <- length(ora_universe)
ora_result <- data.frame(
  ID = names(ora_sets),
  Description = gsub("_", " ", names(ora_sets), fixed = TRUE),
  GeneRatio = vapply(ora_sets, function(members) {
    paste0(length(members), "/", ora_n)
  }, character(1)),
  BgRatio = vapply(ora_sets, function(members) {
    paste0(length(members) + 8L, "/200")
  }, character(1)),
  pvalue = c(1e-6, 4e-5, 2e-4, 8e-4, 0.003, 0.01),
  p.adjust = c(4e-6, 1e-4, 5e-4, 0.002, 0.008, 0.02),
  Count = lengths(ora_sets),
  geneID = vapply(ora_sets, function(members) {
    paste(members, collapse = "/")
  }, character(1)),
  stringsAsFactors = FALSE
)
network <- plot_ora_network(
  ora_result,
  terms = names(ora_sets),
  feature_values = ora_effects,
  label_features = c("IL6", "CXCL8", "TNF", "STAT1", "MX1")
)
readme_save("readme-ora-network.png", network)

n_ranks <- 360L
rank_ids <- sprintf("G%03d", seq_len(n_ranks))
ranks <- stats::setNames(seq(3.2, -3.2, length.out = n_ranks), rank_ids)
ridge_terms <- c(
  "HALLMARK_INTERFERON_GAMMA",
  "HALLMARK_TNFA_SIGNALING",
  "HALLMARK_INFLAMMATORY",
  "HALLMARK_IL6_JAK_STAT3",
  "HALLMARK_OXPHOS",
  "HALLMARK_MYC_TARGETS"
)
ridge_labels <- c(
  HALLMARK_INTERFERON_GAMMA = "Interferon gamma",
  HALLMARK_TNFA_SIGNALING = "TNFA signaling",
  HALLMARK_INFLAMMATORY = "Inflammatory response",
  HALLMARK_IL6_JAK_STAT3 = "IL6 JAK STAT3",
  HALLMARK_OXPHOS = "Oxidative phosphorylation",
  HALLMARK_MYC_TARGETS = "MYC targets"
)
rank_windows <- list(
  seq.int(1L, 90L),
  seq.int(12L, 110L),
  seq.int(24L, 130L),
  seq.int(36L, 150L),
  seq.int(220L, 330L),
  seq.int(250L, 360L)
)
gene_sets <- lapply(rank_windows, function(window) {
  sample(rank_ids[window], size = 36L)
})
names(gene_sets) <- ridge_terms
leading_edge <- lapply(gene_sets, function(members) {
  members[seq_len(20L)]
})
gsea_result <- data.frame(
  pathway = ridge_terms,
  ES = c(0.62, 0.54, 0.48, 0.41, -0.52, -0.47),
  NES = c(2.1, 1.8, 1.6, 1.4, -1.9, -1.6),
  pval = c(1e-5, 4e-4, 0.002, 0.008, 3e-4, 0.006),
  padj = c(4e-5, 0.001, 0.006, 0.02, 8e-4, 0.015),
  qvalue = c(6e-5, 0.002, 0.009, 0.03, 0.001, 0.02),
  size = lengths(gene_sets),
  stringsAsFactors = FALSE
)
gsea_result$leadingEdge <- I(unname(leading_edge))
ridge <- plot_gsea_ridge(
  gsea_result,
  terms = ridge_terms,
  ranks = ranks,
  gene_sets = gene_sets,
  membership = "gene_set",
  term_labels = ridge_labels,
  show_statistics = TRUE
)
readme_save("readme-gsea-ridge.png", ridge)

written <- file.path(
  out_dir,
  c(
    "readme-volcano.png",
    "readme-heatmap.png",
    "readme-ora-network.png",
    "readme-gsea-ridge.png"
  )
)
if (!all(file.exists(written))) {
  stop("Missing rendered README figures.", call. = FALSE)
}
sizes <- file.info(written)$size
if (any(sizes < 10000)) {
  stop("A README figure is unexpectedly small.", call. = FALSE)
}
cat(sprintf("%s\t%d\n", written, sizes), sep = "")
