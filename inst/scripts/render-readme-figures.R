# Compose the README gallery: four near-square panels in one patchwork figure.
# Volcano and heatmap come from the airway dexamethasone experiment.
# Run from the repository root:
#   R_LIBS_USER="$HOME/R/library" Rscript inst/scripts/render-readme-figures.R

.libPaths(c(Sys.getenv("R_LIBS_USER", unset = .libPaths()[[1L]]), .libPaths()))
needed <- c("bulkMAE", "ggplot2", "patchwork", "airway", "DESeq2", "edgeR")
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

# Volcano and heatmap use the airway dexamethasone experiment. The ORA and
# GSEA panels stay schematic: a landing-page gallery should not rerun GO.
data("airway", package = "airway")
counts <- SummarizedExperiment::assay(airway, "counts")
sample_data <- as.data.frame(SummarizedExperiment::colData(airway), optional = TRUE)
feature_data <- as.data.frame(SummarizedExperiment::rowData(airway), optional = TRUE)[
  , c("gene_id", "gene_name", "symbol", "gene_biotype"),
  drop = FALSE
]
sample_data$dex <- factor(sample_data$dex, levels = c("untrt", "trt"))
sample_data$cell <- factor(sample_data$cell)
feature_data$display_label <- as.character(feature_data$symbol)
missing_label <- is.na(feature_data$display_label) | !nzchar(feature_data$display_label)
feature_data$display_label[missing_label] <- feature_data$gene_id[missing_label]
feature_data$display_label <- make.unique(feature_data$display_label)
mae <- mae_from_matrix(counts, sample_data, feature_data, "airway", "counts")
keep <- filter_expr(mae, "airway", group = "dex")
mae_filtered <- mae_subset_features(mae, "airway", names(keep)[keep])
vst <- transform_vst(mae_filtered, "airway", blind = FALSE, design = ~ cell + dex)
mae_vst <- mae_add_assay(mae_filtered, "airway", vst, name = "vst")
fit <- de_deseq2(
  mae_filtered, "airway",
  design = ~ cell + dex, fitType = "parametric", quiet = TRUE
)
de_result <- de_deseq2_results(fit, contrast = c("dex", "trt", "untrt"), alpha = 0.05)
crispld2 <- "ENSG00000103196"
nominated <- c(
  CRISPLD2 = crispld2,
  FKBP5 = "ENSG00000096060",
  TSC22D3 = "ENSG00000157514",
  DUSP1 = "ENSG00000120129",
  PER1 = "ENSG00000179094",
  IL6 = "ENSG00000136244"
)
feature_labels <- stats::setNames(
  mae_feature_data(mae_filtered, "airway")$display_label,
  rownames(mae_feature_data(mae_filtered, "airway"))
)
feature_labels[crispld2] <- "CRISPLD2"

volcano <- square_panel(plot_de_volcano(
  de_result,
  fdr = 0.05,
  min_abs_effect = 1,
  label_features = crispld2,
  feature_labels = feature_labels
), side = "left") +
  theme(
    legend.position = "bottom",
    legend.direction = "horizontal"
  )

heatmap <- square_panel(plot_assay_heatmap(
  mae_vst,
  "airway",
  "vst",
  features = unname(nominated),
  scale = "row",
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  column_split = "dex",
  feature_label = "display_label"
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
