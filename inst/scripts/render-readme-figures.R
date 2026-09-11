# Compose the README gallery: four near-square panels in one patchwork figure.
# Run from the repository root:
#   R_LIBS_USER="$HOME/R/library" Rscript inst/scripts/render-readme-figures.R

.libPaths(c(Sys.getenv("R_LIBS_USER", unset = .libPaths()[[1L]]), .libPaths()))
needed <- c("bulkMAE", "ggplot2", "patchwork")
missing <- needed[!vapply(needed, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) {
  stop("Install packages first: ", paste(missing, collapse = ", "), call. = FALSE)
}

library(bulkMAE)
library(ggplot2)
library(patchwork)

set.seed(20260911)
out_dir <- file.path("man", "figures")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Adjacent left/right plot.margins are 40 pt + 40 pt = 80 pt of column gap.
# patchwork treats unit(..., "pt") in widths as relative shares, so do not
# insert an 80 pt plot_spacer() column.
square_panel <- function(plot, side = c("left", "right")) {
  side <- match.arg(side)
  panel_margin <- if (identical(side, "left")) {
    margin(6, 40, 6, 6)
  } else {
    margin(6, 6, 6, 40)
  }
  plot +
    theme_bulkmae(base_size = 8) +
    theme(
      plot.margin = panel_margin,
      legend.text = element_text(size = 7),
      legend.title = element_text(size = 8)
    )
}

n_features <- 420L
feature_ids <- sprintf("GENE%04d", seq_len(n_features))
marker_ids <- c("IL6", "CXCL8", "TNF", "STAT1")
feature_ids[seq_along(marker_ids)] <- marker_ids
effect <- stats::rnorm(n_features, mean = 0, sd = 0.65)
up <- seq_len(24L)
down <- 25L:48L
effect[up] <- stats::runif(length(up), 1.2, 2.8)
effect[down] <- stats::runif(length(down), -2.8, -1.2)
z_score <- effect / 0.48 + stats::rnorm(n_features, mean = 0, sd = 0.8)
p_value <- pmin(pmax(2 * stats::pnorm(-abs(z_score)), 1e-16), 1)
de_table <- data.frame(
  feature_id = feature_ids,
  effect = effect,
  p_value = p_value,
  adjusted_p_value = stats::p.adjust(p_value, method = "BH"),
  row.names = feature_ids,
  check.names = FALSE
)
volcano <- square_panel(plot_de_volcano(
  de_table,
  fdr = 0.05,
  min_abs_effect = 1,
  label_features = marker_ids
), side = "left") +
  theme(
    legend.position = "bottom",
    legend.direction = "horizontal"
  )

mae <- mae_simulate(n_features = 60L, n_samples = 8L, seed = 20260911L)
log_expr <- mae_pull_assay(mae, "rna", "log_expression")
heat_features <- names(sort(apply(log_expr, 1L, stats::sd), decreasing = TRUE))[seq_len(8L)]
heatmap <- square_panel(plot_assay_heatmap(
  mae,
  "rna",
  "log_expression",
  features = heat_features,
  cluster_rows = TRUE,
  cluster_columns = TRUE
), side = "right") +
  theme(
    axis.text.x = element_text(size = 7, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 7)
  )

ora_sets <- list(
  Inflammatory = c("IL6", "CXCL8", "TNF", "IL1B", "NFKB1", "RELA", "JUN", "PTGS2"),
  TNF_signaling = c("TNF", "NFKB1", "JUN", "CXCL8", "TNFAIP3", "IL1B", "RELA", "MAPK8"),
  IL6_STAT3 = c("IL6", "STAT3", "SOCS3", "JUN", "IL6R", "JAK2", "MYC"),
  Apoptosis = c("TNF", "BIRC3", "CASP8", "FAS", "BID", "JUN", "BCL2A1"),
  Hypoxia = c("VEGFA", "JUN", "NFKB1", "ENO1", "LDHA", "HK2", "SLC2A1")
)
ora_universe <- unique(unlist(ora_sets, use.names = FALSE))
ora_effects <- stats::setNames(
  stats::rnorm(length(ora_universe), mean = 0.15, sd = 1),
  ora_universe
)
ora_effects[c("IL6", "CXCL8", "TNF", "JUN")] <- c(2.2, 1.9, 1.7, 1.3)
ora_n <- length(ora_universe)
ora_result <- data.frame(
  ID = names(ora_sets),
  Description = gsub("_", " ", names(ora_sets), fixed = TRUE),
  GeneRatio = vapply(ora_sets, function(members) {
    paste0(length(members), "/", ora_n)
  }, character(1)),
  BgRatio = vapply(ora_sets, function(members) {
    paste0(length(members) + 6L, "/180")
  }, character(1)),
  pvalue = c(2e-6, 8e-5, 4e-4, 0.002, 0.008),
  p.adjust = c(6e-6, 2e-4, 9e-4, 0.005, 0.016),
  Count = lengths(ora_sets),
  geneID = vapply(ora_sets, paste, character(1), collapse = "/"),
  stringsAsFactors = FALSE
)
network <- plot_ora_network(
  ora_result,
  terms = names(ora_sets),
  feature_values = ora_effects
) +
  theme_void(base_size = 8) +
  theme(
    legend.position = "none",
    plot.background = element_rect(fill = "white", colour = NA),
    panel.background = element_rect(fill = "white", colour = NA),
    panel.border = element_blank(),
    panel.grid = element_blank(),
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.line = element_blank(),
    plot.margin = margin(6, 40, 6, 6)
  )

n_ranks <- 240L
rank_ids <- sprintf("G%03d", seq_len(n_ranks))
ranks <- stats::setNames(seq(2.8, -2.8, length.out = n_ranks), rank_ids)
ridge_terms <- c(
  "HALLMARK_INTERFERON_GAMMA",
  "HALLMARK_INFLAMMATORY",
  "HALLMARK_OXPHOS",
  "HALLMARK_MYC_TARGETS"
)
ridge_labels <- c(
  HALLMARK_INTERFERON_GAMMA = "IFN gamma",
  HALLMARK_INFLAMMATORY = "Inflammatory",
  HALLMARK_OXPHOS = "OXPHOS",
  HALLMARK_MYC_TARGETS = "MYC targets"
)
rank_windows <- list(
  seq.int(1L, 80L),
  seq.int(16L, 100L),
  seq.int(140L, 220L),
  seq.int(160L, 240L)
)
gene_sets <- lapply(rank_windows, function(window) {
  sample(rank_ids[window], size = 40L)
})
names(gene_sets) <- ridge_terms
gsea_result <- data.frame(
  pathway = ridge_terms,
  ES = c(0.58, 0.46, -0.5, -0.44),
  NES = c(1.9, 1.5, -1.8, -1.5),
  pval = c(2e-5, 0.003, 4e-4, 0.007),
  padj = c(6e-5, 0.008, 0.001, 0.016),
  qvalue = c(8e-5, 0.01, 0.002, 0.02),
  size = lengths(gene_sets),
  stringsAsFactors = FALSE
)
gsea_result$leadingEdge <- I(unname(lapply(gene_sets, function(members) {
  members[seq_len(16L)]
})))
ridge <- square_panel(plot_gsea_ridge(
  gsea_result,
  terms = ridge_terms,
  ranks = ranks,
  gene_sets = gene_sets,
  membership = "gene_set",
  term_labels = ridge_labels,
  show_statistics = FALSE
), side = "right")

gallery <- (volcano + heatmap) / (network + ridge) +
  plot_layout(widths = c(1, 1), heights = c(1, 1)) +
  plot_annotation(tag_levels = "A") &
  theme(
    plot.tag = element_text(size = 11, face = "bold"),
    plot.background = element_rect(fill = "white", colour = NA)
  )

gallery_path <- file.path(out_dir, "readme-gallery.png")
ggplot2::ggsave(
  gallery_path,
  gallery,
  width = 18,
  height = 18,
  units = "cm",
  dpi = 220,
  bg = "white"
)

stale <- file.path(
  out_dir,
  c(
    "readme-volcano.png",
    "readme-heatmap.png",
    "readme-ora-network.png",
    "readme-gsea-ridge.png"
  )
)
unlink(stale)

info <- file.info(gallery_path)
if (is.na(info$size) || info$size < 40000) {
  stop("README gallery figure is missing or unexpectedly small.", call. = FALSE)
}
cat(sprintf("%s\t%d\n", gallery_path, info$size))
