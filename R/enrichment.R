#' Run over-representation analysis against arbitrary gene sets
#'
#' @param genes Gene identifiers of interest.
#' @param gene_sets Named list of gene vectors, or a two-column term-to-gene
#'   data frame.
#' @param universe Optional background gene universe.
#' @param ... Additional arguments passed to [clusterProfiler::enricher()].
#'
#' @return A native `enrichResult` object.
#' @export
enrich_ora <- function(genes, gene_sets, universe = NULL, ...) {
  .require_backend("clusterProfiler", "to run over-representation analysis")
  inputs <- .prepare_ora_inputs(genes, universe, "enrich_ora")
  term_to_gene <- .term_to_gene(gene_sets)
  .warn_gene_set_overlap(inputs$genes, term_to_gene$gene, "gene_sets")
  clusterProfiler::enricher(
    gene = inputs$genes,
    universe = inputs$universe,
    TERM2GENE = term_to_gene,
    ...
  )
}

#' Run Gene Ontology enrichment
#'
#' @param genes Gene identifiers of interest. Required only for `method =
#'   "ora"`.
#' @param ranks Named numeric gene-level statistics. Required only for
#'   `method = "gsea"`.
#' @param method Enrichment method: over-representation analysis (`ora`) or
#'   preranked GSEA (`gsea`).
#' @param org_db An installed organism annotation database.
#' @param key_type Identifier type used in `genes` or `ranks`.
#' @param ontology GO ontology: `BP`, `MF`, `CC`, or `ALL`.
#' @param universe Optional background gene universe for ORA.
#' @param ... Additional arguments passed to [clusterProfiler::enrichGO()] or
#'   [clusterProfiler::gseGO()].
#'
#' @return A native `enrichResult` or `gseaResult` object.
#' @export
enrich_go <- function(
    genes = NULL,
    ranks = NULL,
    method = c("ora", "gsea"),
    org_db,
    key_type = "ENTREZID",
    ontology = "BP",
    universe = NULL,
    ...
) {
  .require_backend("clusterProfiler", "to run GO enrichment")
  method <- match.arg(method)
  if (method == "ora") {
    if (is.null(genes)) {
      stop("`genes` is required when `method = \"ora\"`.", call. = FALSE)
    }
    inputs <- .prepare_ora_inputs(genes, universe, "enrich_go")
    return(clusterProfiler::enrichGO(
      gene = inputs$genes,
      OrgDb = org_db,
      keyType = key_type,
      ont = ontology,
      universe = inputs$universe,
      ...
    ))
  }
  if (is.null(ranks)) {
    stop("`ranks` is required when `method = \"gsea\"`.", call. = FALSE)
  }
  .assert_named_numeric(ranks, "ranks")
  clusterProfiler::gseGO(
    geneList = sort(ranks, decreasing = TRUE),
    OrgDb = org_db,
    keyType = key_type,
    ont = ontology,
    ...
  )
}

#' Run selection-bias-aware enrichment with goseq
#'
#' @param selected Named binary or logical vector over all tested genes.
#' @param genome Genome identifier accepted by [goseq::nullp()].
#' @param id Gene identifier type accepted by [goseq::nullp()].
#' @param bias Optional named or position-aligned numeric bias, usually gene
#'   length or expression abundance.
#' @param gene_to_category Optional gene-to-category mapping.
#' @param categories GO/KEGG categories passed as `test.cats`.
#' @param method Enrichment method. `Wallenius` is the backend recommendation.
#' @param repetitions Sampling repetitions, used only by the sampling method.
#' @param include_uncategorized Include genes without category annotations.
#' @param plot_fit Plot the probability-weighting function for review.
#'
#' @return The native enrichment data frame returned by [goseq::goseq()].
#' @export
enrich_goseq <- function(
    selected,
    genome,
    id,
    bias = NULL,
    gene_to_category = NULL,
    categories = c("GO:CC", "GO:BP", "GO:MF"),
    method = "Wallenius",
    repetitions = 2000,
    include_uncategorized = FALSE,
    plot_fit = FALSE
) {
  .require_backend("goseq", "to run selection-bias-aware enrichment")
  if (is.logical(selected)) {
    selected_names <- names(selected)
    selected <- as.integer(selected)
    names(selected) <- selected_names
  }
  if (
    !is.numeric(selected) || is.null(names(selected)) || anyNA(names(selected)) ||
      any(!nzchar(names(selected))) || any(!selected %in% c(0, 1)) ||
      anyDuplicated(names(selected))
  ) {
    stop("`selected` must be a named binary or logical vector.", call. = FALSE)
  }
  if (length(unique(selected)) != 2L) {
    stop("`selected` must contain both selected and unselected genes.", call. = FALSE)
  }
  if (!is.null(bias)) {
    if (!is.null(names(bias))) {
      if (
        anyNA(names(bias)) || any(!nzchar(names(bias))) ||
          anyDuplicated(names(bias))
      ) {
        stop("Named `bias` must have unique, non-empty gene names.", call. = FALSE)
      }
      missing <- setdiff(names(selected), names(bias))
      if (length(missing)) {
        stop("Named `bias` is missing tested genes.", call. = FALSE)
      }
      bias <- bias[names(selected)]
    }
    if (
      !is.numeric(bias) || length(bias) != length(selected) ||
        any(!is.finite(bias)) || any(bias < 0)
    ) {
      stop(
        "`bias` must have one finite, non-negative value per tested gene.",
        call. = FALSE
      )
    }
    if (length(unique(bias)) < 2L) {
      stop("`bias` must vary across tested genes.", call. = FALSE)
    }
  }
  if (
    !is.character(method) || length(method) != 1L || is.na(method) ||
      !nzchar(method)
  ) {
    stop("`method` must name one goseq enrichment method.", call. = FALSE)
  }
  if (
    !is.numeric(repetitions) || length(repetitions) != 1L ||
      !is.finite(repetitions) || repetitions < 1L ||
      repetitions != as.integer(repetitions)
  ) {
    stop("`repetitions` must be a positive integer.", call. = FALSE)
  }
  if (
    !is.logical(include_uncategorized) || length(include_uncategorized) != 1L ||
      is.na(include_uncategorized) || !is.logical(plot_fit) ||
      length(plot_fit) != 1L || is.na(plot_fit)
  ) {
    stop("`include_uncategorized` and `plot_fit` must be logical scalars.", call. = FALSE)
  }
  probability_weights <- goseq::nullp(
    selected,
    genome = genome,
    id = id,
    bias.data = bias,
    plot.fit = plot_fit
  )
  goseq::goseq(
    probability_weights,
    genome = genome,
    id = id,
    gene2cat = gene_to_category,
    test.cats = categories,
    method = method,
    repcnt = repetitions,
    use_genes_without_cat = include_uncategorized
  )
}

#' Run KEGG enrichment
#'
#' @param genes Gene identifiers of interest. Required only for `method =
#'   "ora"`.
#' @param ranks Named numeric gene-level statistics. Required only for
#'   `method = "gsea"`.
#' @param method Enrichment method: over-representation analysis (`ora`) or
#'   preranked GSEA (`gsea`).
#' @param organism KEGG organism code, for example `hsa` or `mmu`.
#' @param key_type Identifier type accepted by the clusterProfiler backend.
#' @param universe Optional background gene universe for ORA.
#' @param ... Additional arguments passed to [clusterProfiler::enrichKEGG()]
#'   or [clusterProfiler::gseKEGG()].
#'
#' @return A native `enrichResult` or `gseaResult` object.
#' @export
enrich_kegg <- function(
    genes = NULL,
    ranks = NULL,
    method = c("ora", "gsea"),
    organism = "hsa",
    key_type = "kegg",
    universe = NULL,
    ...
) {
  .require_backend("clusterProfiler", "to run KEGG enrichment")
  method <- match.arg(method)
  if (method == "ora") {
    if (is.null(genes)) {
      stop("`genes` is required when `method = \"ora\"`.", call. = FALSE)
    }
    inputs <- .prepare_ora_inputs(genes, universe, "enrich_kegg")
    return(clusterProfiler::enrichKEGG(
      gene = inputs$genes,
      organism = organism,
      keyType = key_type,
      universe = inputs$universe,
      ...
    ))
  }
  if (is.null(ranks)) {
    stop("`ranks` is required when `method = \"gsea\"`.", call. = FALSE)
  }
  .assert_named_numeric(ranks, "ranks")
  clusterProfiler::gseKEGG(
    geneList = sort(ranks, decreasing = TRUE),
    organism = organism,
    keyType = key_type,
    ...
  )
}

#' Run Reactome enrichment
#'
#' @param genes Entrez gene identifiers of interest. Required only for
#'   `method = "ora"`.
#' @param ranks Named numeric Entrez-level statistics. Required only for
#'   `method = "gsea"`.
#' @param method Enrichment method: over-representation analysis (`ora`) or
#'   preranked GSEA (`gsea`).
#' @param organism Organism accepted by [ReactomePA::enrichPathway()].
#' @param universe Optional background gene universe for ORA.
#' @param ... Additional arguments passed to [ReactomePA::enrichPathway()] or
#'   [ReactomePA::gsePathway()].
#'
#' @return A native `enrichResult` or `gseaResult` object.
#' @export
enrich_reactome <- function(
    genes = NULL,
    ranks = NULL,
    method = c("ora", "gsea"),
    organism = "human",
    universe = NULL,
    ...
) {
  .require_backend("ReactomePA", "to run Reactome enrichment")
  method <- match.arg(method)
  if (method == "ora") {
    if (is.null(genes)) {
      stop("`genes` is required when `method = \"ora\"`.", call. = FALSE)
    }
    inputs <- .prepare_ora_inputs(genes, universe, "enrich_reactome")
    return(ReactomePA::enrichPathway(
      gene = inputs$genes,
      organism = organism,
      universe = inputs$universe,
      ...
    ))
  }
  if (is.null(ranks)) {
    stop("`ranks` is required when `method = \"gsea\"`.", call. = FALSE)
  }
  .assert_named_numeric(ranks, "ranks")
  ReactomePA::gsePathway(
    geneList = sort(ranks, decreasing = TRUE),
    organism = organism,
    ...
  )
}

#' Run preranked gene-set enrichment with fgsea
#'
#' @param ranks Named numeric gene-level statistics.
#' @param pathways Named list of gene sets.
#' @param min_size,max_size Gene-set size limits.
#' @param ... Additional arguments passed to [fgsea::fgseaMultilevel()].
#'
#' @return A native fgsea result table.
#' @export
enrich_fgsea <- function(
    ranks,
    pathways,
    min_size = 15,
    max_size = 500,
    ...
) {
  .require_backend("fgsea", "to run preranked enrichment")
  .assert_named_numeric(ranks, "ranks")
  ranks <- sort(ranks, decreasing = TRUE)
  fgsea::fgseaMultilevel(
    pathways = pathways,
    stats = ranks,
    minSize = min_size,
    maxSize = max_size,
    ...
  )
}

#' Run preranked enrichment against arbitrary gene sets with clusterProfiler
#'
#' @param ranks Named numeric gene-level statistics.
#' @param gene_sets Named list of gene sets, or a term-to-gene data frame.
#' @param ... Additional arguments passed to [clusterProfiler::GSEA()].
#'
#' @return A native `gseaResult` object.
#' @export
enrich_gsea <- function(ranks, gene_sets, ...) {
  .require_backend("clusterProfiler", "to run gene-set enrichment")
  .assert_named_numeric(ranks, "ranks")
  clusterProfiler::GSEA(
    geneList = sort(ranks, decreasing = TRUE),
    TERM2GENE = .term_to_gene(gene_sets),
    ...
  )
}

#' Run CAMERA competitive gene-set testing
#'
#' The selected assay must contain an approximately homoscedastic log-scale
#' expression matrix. For voom analyses, supply the corresponding observation
#' weights through `weights`; omitting them changes the intended mean-variance
#' model.
#'
#' @inheritParams mae_pull_assay
#' @param gene_sets Named list of gene sets.
#' @param formula Model formula used to construct the design matrix.
#' @param contrast Contrast accepted by [limma::camera()].
#' @param weights Optional observation-level precision-weight matrix, or the
#'   name of an assay containing one. It must have the same dimensions and
#'   dimnames as `assay`. When omitted, no voom precision weights are used.
#' @param ... Additional arguments passed to [limma::camera()].
#'
#' @return The native CAMERA result data frame.
#' @export
enrich_camera <- function(
    x,
    experiment,
    gene_sets,
    formula,
    contrast,
    assay,
    weights = NULL,
    ...
) {
  .require_backend("limma", "to run CAMERA")
  matrix <- .pull_matrix(x, experiment, assay)
  if (any(!is.finite(matrix))) {
    stop("CAMERA requires a finite expression matrix.", call. = FALSE)
  }
  design <- .model_matrix(formula, mae_samples(x, experiment))
  index <- limma::ids2indices(gene_sets, rownames(matrix))
  index <- .assert_nonempty_indices(index, "gene_sets")
  expression <- .limma_expression_input(x, experiment, matrix, weights)
  limma::camera(expression, index, design = design, contrast = contrast, ...)
}

#' Run FRY rotation gene-set testing
#'
#' @inheritParams enrich_camera
#'
#' @return The native FRY result data frame.
#' @export
enrich_fry <- function(
    x,
    experiment,
    gene_sets,
    formula,
    contrast,
    assay,
    weights = NULL,
    ...
) {
  .require_backend("limma", "to run FRY")
  matrix <- .pull_matrix(x, experiment, assay)
  if (any(!is.finite(matrix))) {
    stop("FRY requires a finite expression matrix.", call. = FALSE)
  }
  design <- .model_matrix(formula, mae_samples(x, experiment))
  index <- limma::ids2indices(gene_sets, rownames(matrix))
  index <- .assert_nonempty_indices(index, "gene_sets")
  expression <- .limma_expression_input(x, experiment, matrix, weights)
  limma::fry(expression, index, design = design, contrast = contrast, ...)
}

#' Run mroast rotation gene-set testing
#'
#' @inheritParams enrich_camera
#' @param rotations Number of rotations passed as `nrot`.
#' @param seed Optional random seed.
#' @param ... Additional arguments passed to [limma::mroast()].
#'
#' @return The native mroast result data frame.
#' @export
enrich_roast <- function(
    x,
    experiment,
    gene_sets,
    formula,
    contrast,
    assay,
    weights = NULL,
    rotations = 1999,
    seed = NULL,
    ...
) {
  .require_backend("limma", "to run mroast")
  if (
    !is.numeric(rotations) || length(rotations) != 1L ||
      !is.finite(rotations) || rotations < 1L ||
      rotations != as.integer(rotations)
  ) {
    stop("`rotations` must be a positive integer.", call. = FALSE)
  }
  matrix <- .pull_matrix(x, experiment, assay)
  if (any(!is.finite(matrix))) {
    stop("mroast requires a finite expression matrix.", call. = FALSE)
  }
  design <- .model_matrix(formula, mae_samples(x, experiment))
  index <- limma::ids2indices(gene_sets, rownames(matrix))
  index <- .assert_nonempty_indices(index, "gene_sets")
  expression <- .limma_expression_input(x, experiment, matrix, weights)
  .with_seed(
    seed,
    limma::mroast(
      expression,
      index,
      design = design,
      contrast = contrast,
      nrot = rotations,
      ...
    )
  )
}

#' Calculate GSVA scores
#'
#' Uses the current parameter-object API: [GSVA::gsvaParam()] followed by
#' [GSVA::gsva()].
#'
#' @inheritParams mae_pull_assay
#' @param gene_sets Named list of gene sets.
#' @param kcdf Kernel used by GSVA.
#' @param min_size,max_size Gene-set size limits.
#' @param verbose Show backend progress.
#' @param ... Additional arguments passed to [GSVA::gsvaParam()].
#'
#' @return The native matrix-like object returned by [GSVA::gsva()].
#' @export
score_gsva <- function(
    x,
    experiment,
    gene_sets,
    assay,
    kcdf = "auto",
    min_size = 1,
    max_size = Inf,
    verbose = FALSE,
    ...
) {
  .require_backend("GSVA", "to calculate GSVA scores")
  if (!is.logical(verbose) || length(verbose) != 1L || is.na(verbose)) {
    stop("`verbose` must be TRUE or FALSE.", call. = FALSE)
  }
  matrix <- .pull_matrix(x, experiment, assay)
  if (any(!is.finite(matrix))) {
    stop("GSVA requires a finite assay matrix.", call. = FALSE)
  }
  gene_sets <- .prepare_sample_gene_sets(
    gene_sets,
    rownames(matrix),
    min_size,
    max_size
  )
  parameter <- GSVA::gsvaParam(
    matrix,
    gene_sets,
    kcdf = kcdf,
    minSize = min_size,
    maxSize = max_size,
    ...
  )
  GSVA::gsva(parameter, verbose = verbose)
}

#' Calculate single-sample GSEA scores
#'
#' Uses the current [GSVA::ssgseaParam()] parameter-object API.
#'
#' @inheritParams score_gsva
#' @param normalize Normalize ssGSEA scores by their range.
#' @param alpha Tail-weight exponent used by ssGSEA.
#'
#' @return The native matrix-like object returned by [GSVA::gsva()].
#' @export
score_ssgsea <- function(
    x,
    experiment,
    gene_sets,
    assay,
    min_size = 1,
    max_size = Inf,
    normalize = TRUE,
    alpha = 0.25,
    verbose = FALSE,
    ...
) {
  .require_backend("GSVA", "to calculate ssGSEA scores")
  if (
    !is.logical(normalize) || length(normalize) != 1L || is.na(normalize) ||
      !is.logical(verbose) || length(verbose) != 1L || is.na(verbose)
  ) {
    stop("`normalize` and `verbose` must be TRUE or FALSE.", call. = FALSE)
  }
  if (
    !is.numeric(alpha) || length(alpha) != 1L || is.na(alpha) ||
      !is.finite(alpha) || alpha < 0
  ) {
    stop("`alpha` must be one finite, non-negative number.", call. = FALSE)
  }
  matrix <- .pull_matrix(x, experiment, assay)
  if (any(!is.finite(matrix))) {
    stop("ssGSEA requires a finite assay matrix.", call. = FALSE)
  }
  gene_sets <- .prepare_sample_gene_sets(
    gene_sets,
    rownames(matrix),
    min_size,
    max_size
  )
  parameter <- GSVA::ssgseaParam(
    matrix,
    gene_sets,
    minSize = min_size,
    maxSize = max_size,
    normalize = normalize,
    alpha = alpha,
    ...
  )
  GSVA::gsva(parameter, verbose = verbose)
}

#' Score samples with rank-based singscore
#'
#' @inheritParams mae_pull_assay
#' @param up_set Features expected to be up-regulated.
#' @param down_set Optional features expected to be down-regulated.
#' @param ties_method Ranking method passed to [singscore::rankGenes()].
#' @param known_direction Whether signature direction is known.
#' @param ... Additional arguments passed to [singscore::simpleScore()].
#'
#' @return The native data frame returned by [singscore::simpleScore()].
#' @export
score_singscore <- function(
    x,
    experiment,
    up_set,
    down_set = NULL,
    assay,
    ties_method = "min",
    known_direction = TRUE,
    ...
) {
  .require_backend("singscore", "to calculate rank-based signature scores")
  matrix <- .pull_matrix(x, experiment, assay)
  if (any(!is.finite(matrix))) {
    stop("singscore input must contain finite expression values.", call. = FALSE)
  }
  up_set <- .prepare_singscore_set(up_set, rownames(matrix), "up_set")
  if (!is.null(down_set)) {
    down_set <- .prepare_singscore_set(
      down_set,
      rownames(matrix),
      "down_set"
    )
    overlap <- intersect(up_set, down_set)
    if (length(overlap)) {
      stop(
        "`up_set` and `down_set` must not overlap; shared features include: ",
        paste(utils::head(overlap, 10L), collapse = ", "),
        call. = FALSE
      )
    }
  }
  ranked <- singscore::rankGenes(matrix, tiesMethod = ties_method)
  singscore::simpleScore(
    rankData = ranked,
    upSet = up_set,
    downSet = down_set,
    knownDirection = known_direction,
    ...
  )
}

#' List methods available through decoupleR
#'
#' @return The native method table returned by [decoupleR::show_methods()].
#' @export
activity_methods <- function() {
  .require_backend("decoupleR", "to list activity-inference methods")
  decoupleR::show_methods()
}

#' Retrieve a network or gene-set resource through decoupleR
#'
#' @param name Resource name accepted by [decoupleR::get_resource()].
#' @param organism Organism name accepted by [decoupleR::get_resource()].
#' @param ... Additional arguments passed to [decoupleR::get_resource()].
#'
#' @return The native resource table returned by decoupleR.
#' @export
activity_resource <- function(name, organism = "human", ...) {
  .require_backend("decoupleR", "to retrieve a molecular resource")
  decoupleR::get_resource(name = name, organism = organism, ...)
}

#' Run one or more decoupleR methods
#'
#' This general adapter deliberately exposes both enrichment-style methods
#' (`aucell`, `fgsea`, `gsva`, and `ora`) and network-aware activity methods
#' (`mlm`, `ulm`, `viper`, `wmean`, and `wsum`). With `statistics = NULL`,
#' decoupleR runs its documented default methods and can calculate a consensus.
#'
#' @inheritParams mae_pull_assay
#' @param network Long-format network or gene-set table.
#' @param statistics Method names accepted by [decoupleR::decouple()], or
#'   `NULL` for its defaults.
#' @param method_args List of method-specific argument lists.
#' @param consensus Calculate a consensus score.
#' @param consensus_statistics Optional decoupleR result-statistic names used
#'   for consensus scoring, for example `norm_mlm` and `norm_ulm`. `NULL`
#'   preserves the backend default.
#' @param source,target Column names identifying regulators/sets and targets.
#' @param mor Optional column containing signed interaction weights. It is
#'   mapped to decoupleR's conventional `mor` column without discarding the
#'   original column.
#' @param min_size Minimum number of targets per source.
#' @param include_time Include execution time in the result.
#' @param ... Additional arguments passed to [decoupleR::decouple()].
#'
#' @return The native long-format table returned by [decoupleR::decouple()].
#' @export
activity_decouple <- function(
    x,
    experiment,
    network,
    assay,
    statistics = NULL,
    method_args = list(NULL),
    consensus = TRUE,
    consensus_statistics = NULL,
    source = "source",
    target = "target",
    mor = NULL,
    min_size = 5,
    include_time = FALSE,
    ...
) {
  .require_backend("decoupleR", "to run enrichment or activity inference")
  matrix <- .pull_matrix(x, experiment, assay)
  if (any(!is.finite(matrix))) {
    stop("decoupleR requires a finite assay matrix.", call. = FALSE)
  }
  statistics <- .expand_decouple_statistics(statistics)
  network <- .prepare_decouple_network(
    network,
    rownames(matrix),
    source,
    target,
    mor,
    min_size
  )
  .assert_decouple_mor(network, statistics)
  method_args <- .prepare_decouple_arguments(statistics, method_args)
  .validate_decouple_consensus(
    statistics,
    consensus,
    consensus_statistics
  )
  arguments <- c(
    list(
      mat = matrix,
      network = network,
      .source = source,
      .target = target,
      statistics = statistics,
      args = method_args,
      consensus_score = consensus,
      include_time = include_time,
      minsize = min_size
    ),
    list(...)
  )
  if (!is.null(consensus_statistics)) {
    arguments$consensus_stats <- consensus_statistics
  }
  do.call(decoupleR::decouple, arguments)
}

#' Infer pathway activity with PROGENy and decoupleR
#'
#' @inheritParams mae_pull_assay
#' @param organism Organism accepted by [decoupleR::get_progeny()].
#' @param top Number of responsive genes retained per pathway.
#' @param min_size Minimum regulon size passed to [decoupleR::run_mlm()].
#' @param ... Additional arguments passed to [decoupleR::run_mlm()].
#'
#' @return The native decoupleR result table.
#' @export
activity_progeny <- function(
    x,
    experiment,
    assay,
    organism = "human",
    top = 500,
    min_size = 5,
    ...
) {
  .require_backend("decoupleR", "to infer PROGENy pathway activity")
  matrix <- .pull_matrix(x, experiment, assay)
  if (any(!is.finite(matrix))) {
    stop("PROGENy activity inference requires a finite assay matrix.", call. = FALSE)
  }
  network <- decoupleR::get_progeny(organism = organism, top = top)
  network <- .prepare_decouple_network(
    network,
    rownames(matrix),
    "source",
    "target",
    "weight",
    min_size
  )
  decoupleR::run_mlm(
    mat = matrix,
    network = network,
    .source = "source",
    .target = "target",
    .mor = "mor",
    minsize = min_size,
    ...
  )
}

#' Infer transcription-factor activity with decoupleR
#'
#' @inheritParams mae_pull_assay
#' @param organism Organism accepted by the selected network resource.
#' @param resource TF network: current CollecTRI or DoRothEA.
#' @param min_size Minimum regulon size passed to [decoupleR::run_ulm()].
#' @param ... Additional arguments passed to [decoupleR::run_ulm()].
#'
#' @return The native decoupleR result table.
#' @export
activity_tf <- function(
    x,
    experiment,
    assay,
    organism = "human",
    resource = c("collectri", "dorothea"),
    min_size = 5,
    ...
) {
  .require_backend("decoupleR", "to infer transcription-factor activity")
  resource <- match.arg(resource)
  matrix <- .pull_matrix(x, experiment, assay)
  if (any(!is.finite(matrix))) {
    stop("TF activity inference requires a finite assay matrix.", call. = FALSE)
  }
  network <- if (resource == "collectri") {
    decoupleR::get_collectri(organism = organism, split_complexes = FALSE)
  } else {
    decoupleR::get_dorothea(organism = organism)
  }
  network <- .prepare_decouple_network(
    network,
    rownames(matrix),
    "source",
    "target",
    "mor",
    min_size
  )
  decoupleR::run_ulm(
    mat = matrix,
    network = network,
    .source = "source",
    .target = "target",
    .mor = "mor",
    minsize = min_size,
    ...
  )
}

.term_to_gene <- function(gene_sets) {
  if (is.data.frame(gene_sets)) {
    if (ncol(gene_sets) < 2L) {
      stop("A term-to-gene table needs at least two columns.", call. = FALSE)
    }
    result <- gene_sets[, 1:2, drop = FALSE]
    names(result) <- c("term", "gene")
  } else {
    if (
      !is.list(gene_sets) || is.null(names(gene_sets)) ||
        any(!nzchar(names(gene_sets))) || anyDuplicated(names(gene_sets))
    ) {
      stop(
        "`gene_sets` must be a uniquely named list or term-to-gene table.",
        call. = FALSE
      )
    }
    gene_sets <- lapply(gene_sets, as.character)
    result <- data.frame(
      term = rep(names(gene_sets), lengths(gene_sets)),
      gene = unlist(gene_sets, use.names = FALSE),
      stringsAsFactors = FALSE
    )
  }

  result$term <- as.character(result$term)
  result$gene <- as.character(result$gene)
  invalid <- is.na(result$term) | !nzchar(result$term) |
    is.na(result$gene) | !nzchar(result$gene)
  if (any(invalid)) {
    stop("Gene-set terms and identifiers must be non-missing strings.", call. = FALSE)
  }
  unique(result)
}

.assert_named_numeric <- function(x, argument) {
  if (
    !is.numeric(x) || !length(x) || is.null(names(x)) || anyNA(names(x)) ||
      any(!nzchar(names(x)))
  ) {
    stop("`", argument, "` must be a named numeric vector.", call. = FALSE)
  }
  if (any(!is.finite(x)) || anyDuplicated(names(x))) {
    stop("`", argument, "` must have finite values and unique names.", call. = FALSE)
  }
  invisible(x)
}

.prepare_ora_inputs <- function(genes, universe, caller) {
  genes <- .clean_gene_ids(genes, "genes")
  if (is.null(universe)) {
    warning(
      "`", caller, "()` received no `universe`; the backend annotation or ",
      "gene-set collection will define the background. For RNA-seq ORA, pass ",
      "all genes that were actually eligible for testing.",
      call. = FALSE
    )
    return(list(genes = genes, universe = NULL))
  }

  universe <- .clean_gene_ids(universe, "universe")
  outside <- setdiff(genes, universe)
  if (length(outside)) {
    warning(
      length(outside), " selected gene(s) are absent from `universe`. Check ",
      "that gene identifiers and ID types match; the enrichment backend may ",
      "discard them.",
      call. = FALSE
    )
  }
  list(genes = genes, universe = universe)
}

.clean_gene_ids <- function(ids, argument) {
  if (!is.atomic(ids) || is.list(ids)) {
    stop("`", argument, "` must be an atomic vector of gene identifiers.", call. = FALSE)
  }
  ids <- as.character(ids)
  if (!length(ids) || anyNA(ids) || any(!nzchar(ids))) {
    stop("`", argument, "` must contain non-missing gene identifiers.", call. = FALSE)
  }
  unique(ids)
}

.warn_gene_set_overlap <- function(genes, set_genes, argument) {
  overlap <- intersect(genes, unique(as.character(set_genes)))
  if (!length(overlap)) {
    stop(
      "No identifiers in `genes` overlap `", argument,
      "`; check the identifier namespace.",
      call. = FALSE
    )
  }
  invisible(overlap)
}

.assert_nonempty_indices <- function(index, argument) {
  if (!is.list(index) || !length(index) || !any(lengths(index) > 0L)) {
    stop(
      "No features in `", argument,
      "` overlap the selected assay identifiers.",
      call. = FALSE
    )
  }
  empty <- sum(lengths(index) == 0L)
  if (empty) {
    warning(empty, " gene set(s) have no assay overlap and will be ignored.", call. = FALSE)
  }
  index[lengths(index) > 0L]
}

.limma_expression_input <- function(x, experiment, matrix, weights) {
  if (is.null(weights)) {
    return(matrix)
  }
  if (is.character(weights) && length(weights) == 1L) {
    weights <- .pull_matrix(x, experiment, weights)
  } else {
    weights <- as.matrix(weights)
    storage.mode(weights) <- "double"
  }
  if (is.null(rownames(weights)) || is.null(colnames(weights))) {
    stop("`weights` must have feature and sample names.", call. = FALSE)
  }
  if (
    !setequal(rownames(weights), rownames(matrix)) ||
      !setequal(colnames(weights), colnames(matrix))
  ) {
    stop("`weights` must describe exactly the selected assay entries.", call. = FALSE)
  }
  weights <- weights[rownames(matrix), colnames(matrix), drop = FALSE]
  if (any(!is.finite(weights)) || any(weights < 0)) {
    stop("`weights` must be finite and non-negative.", call. = FALSE)
  }
  expression <- methods::new("EList")
  expression$E <- matrix
  expression$weights <- weights
  expression
}

.prepare_sample_gene_sets <- function(gene_sets, features, min_size, max_size) {
  if (
    length(min_size) != 1L || !is.numeric(min_size) || is.na(min_size) ||
      min_size < 1 || min_size != as.integer(min_size)
  ) {
    stop("`min_size` must be a positive integer.", call. = FALSE)
  }
  if (
    length(max_size) != 1L || !is.numeric(max_size) || is.na(max_size) ||
      max_size < min_size ||
      (is.finite(max_size) && max_size != as.integer(max_size))
  ) {
    stop("`max_size` must be an integer at least `min_size`, or `Inf`.", call. = FALSE)
  }
  if (
    !is.list(gene_sets) || is.null(names(gene_sets)) || !length(gene_sets) ||
      any(!nzchar(names(gene_sets))) || anyDuplicated(names(gene_sets))
  ) {
    stop("`gene_sets` must be a non-empty, uniquely named list.", call. = FALSE)
  }

  overlapped <- lapply(gene_sets, function(set) {
    set <- as.character(set)
    set <- unique(set[!is.na(set) & nzchar(set)])
    intersect(set, features)
  })
  keep <- lengths(overlapped) >= min_size & lengths(overlapped) <= max_size
  if (!any(keep)) {
    stop(
      "No gene set remains after assay overlap and `min_size`/`max_size` filtering.",
      call. = FALSE
    )
  }
  if (any(!keep)) {
    warning(
      sum(!keep), " gene set(s) were removed after assay-overlap size filtering.",
      call. = FALSE
    )
  }
  overlapped[keep]
}

.prepare_singscore_set <- function(set, features, argument) {
  if (!is.atomic(set) || is.list(set)) {
    stop("`", argument, "` must be a vector of feature identifiers.", call. = FALSE)
  }
  set <- unique(as.character(set))
  if (!length(set) || anyNA(set) || any(!nzchar(set))) {
    stop("`", argument, "` must contain non-missing identifiers.", call. = FALSE)
  }
  overlap <- intersect(set, features)
  if (!length(overlap)) {
    stop("`", argument, "` has no overlap with the selected assay.", call. = FALSE)
  }
  if (length(overlap) < length(set)) {
    warning(
      length(set) - length(overlap), " feature(s) from `", argument,
      "` are absent from the assay and were removed.",
      call. = FALSE
    )
  }
  overlap
}

.prepare_decouple_network <- function(
    network,
    features,
    source,
    target,
    mor,
    min_size
) {
  if (!is.data.frame(network)) {
    stop("`network` must be a data frame.", call. = FALSE)
  }
  for (value in list(source = source, target = target)) {
    if (!is.character(value) || length(value) != 1L || is.na(value) || !nzchar(value)) {
      stop("`source` and `target` must each name one network column.", call. = FALSE)
    }
  }
  if (identical(source, target) || any(!c(source, target) %in% names(network))) {
    stop("`source` and `target` must name distinct network columns.", call. = FALSE)
  }
  if (
    length(min_size) != 1L || !is.numeric(min_size) || is.na(min_size) ||
      min_size < 1 || min_size != as.integer(min_size)
  ) {
    stop("`min_size` must be a positive integer.", call. = FALSE)
  }

  network[[source]] <- as.character(network[[source]])
  network[[target]] <- as.character(network[[target]])
  invalid <- is.na(network[[source]]) | !nzchar(network[[source]]) |
    is.na(network[[target]]) | !nzchar(network[[target]])
  if (any(invalid)) {
    stop("Network source and target identifiers must be non-missing strings.", call. = FALSE)
  }
  if (!is.null(mor)) {
    if (
      !is.character(mor) || length(mor) != 1L || is.na(mor) ||
        !mor %in% names(network)
    ) {
      stop("`mor` must name a network weight column.", call. = FALSE)
    }
    if (!identical(mor, "mor") && "mor" %in% c(source, target)) {
      stop("A mapped `mor` column must not overwrite `source` or `target`.", call. = FALSE)
    }
    values <- network[[mor]]
    if (!is.numeric(values) || any(!is.finite(values))) {
      stop("The `mor` network column must contain finite numeric weights.", call. = FALSE)
    }
    if (
      !identical(mor, "mor") && "mor" %in% names(network) &&
        !identical(network[["mor"]], values)
    ) {
      stop("The network already has a different `mor` column.", call. = FALSE)
    }
    network[["mor"]] <- values
  }
  if ("mor" %in% names(network)) {
    if (!is.numeric(network[["mor"]]) || any(!is.finite(network[["mor"]]))) {
      stop("The normalized `mor` column must contain finite numeric weights.", call. = FALSE)
    }
  }

  keep <- network[[target]] %in% features
  if (!any(keep)) {
    stop(
      "No network targets overlap the selected assay; check identifier types.",
      call. = FALSE
    )
  }
  if (any(!keep)) {
    warning(sum(!keep), " network edge(s) with absent targets were removed.", call. = FALSE)
  }
  network <- network[keep, , drop = FALSE]
  sizes <- vapply(
    split(network[[target]], network[[source]]),
    function(value) length(unique(value)),
    integer(1)
  )
  if (!length(sizes) || max(sizes) < min_size) {
    stop("No network source reaches `min_size` after assay overlap.", call. = FALSE)
  }
  network
}

.assert_decouple_mor <- function(network, statistics) {
  methods <- if (is.null(statistics)) {
    c("mlm", "ulm", "wsum")
  } else {
    sub("^run_", "", tolower(statistics))
  }
  weighted <- c("mlm", "ulm", "viper", "wmean", "wsum")
  if (any(methods %in% weighted) && !"mor" %in% names(network)) {
    stop(
      "The requested weighted decoupleR method(s) require a normalized `mor` ",
      "column. Supply `mor` as the name of the signed network-weight column.",
      call. = FALSE
    )
  }
  invisible(network)
}

.prepare_decouple_arguments <- function(statistics, method_args) {
  if (!is.null(statistics)) {
    if (
      !is.character(statistics) || !length(statistics) || anyNA(statistics) ||
        any(!nzchar(statistics)) || anyDuplicated(statistics)
    ) {
      stop("`statistics` must contain unique method names, or be `NULL`.", call. = FALSE)
    }
  }
  if (!is.list(method_args) || !length(method_args)) {
    stop("`method_args` must be a non-empty list.", call. = FALSE)
  }
  valid_element <- vapply(
    method_args,
    function(value) is.null(value) || is.list(value),
    logical(1)
  )
  if (!all(valid_element)) {
    stop("Every `method_args` element must be a list or `NULL`.", call. = FALSE)
  }
  if (is.null(statistics)) {
    if (length(method_args) != 1L) {
      stop(
        "With `statistics = NULL`, supply only the default `list(NULL)` ",
        "method arguments.",
        call. = FALSE
      )
    }
    return(method_args)
  }

  if (length(method_args) == 1L && is.null(method_args[[1L]])) {
    method_args <- rep(list(NULL), length(statistics))
  }
  if (length(method_args) != length(statistics)) {
    stop("`method_args` must have one element per requested statistic.", call. = FALSE)
  }
  argument_names <- names(method_args)
  if (!is.null(argument_names) && any(nzchar(argument_names))) {
    if (any(!nzchar(argument_names)) || !setequal(argument_names, statistics)) {
      stop("Named `method_args` must be named exactly as `statistics`.", call. = FALSE)
    }
    method_args <- method_args[statistics]
  } else {
    names(method_args) <- statistics
  }
  method_args
}

.expand_decouple_statistics <- function(statistics) {
  if (
    !is.null(statistics) && length(statistics) == 1L &&
      is.character(statistics) && !is.na(statistics) &&
      identical(tolower(statistics), "all")
  ) {
    methods <- decoupleR::show_methods()
    function_column <- if ("Function" %in% names(methods)) {
      methods$Function
    } else {
      methods[[1L]]
    }
    statistics <- sub("^run_", "", tolower(as.character(function_column)))
    statistics <- setdiff(statistics, "consensus")
  }
  statistics
}

.validate_decouple_consensus <- function(
    statistics,
    consensus,
    consensus_statistics
) {
  if (!is.logical(consensus) || length(consensus) != 1L || is.na(consensus)) {
    stop("`consensus` must be `TRUE` or `FALSE`.", call. = FALSE)
  }
  if (!is.null(consensus_statistics)) {
    if (
      !is.character(consensus_statistics) || !length(consensus_statistics) ||
        anyNA(consensus_statistics) || any(!nzchar(consensus_statistics)) ||
        anyDuplicated(consensus_statistics)
    ) {
      stop(
        "`consensus_statistics` must contain unique statistic names.",
        call. = FALSE
      )
    }
    if (!consensus) {
      stop("`consensus_statistics` is only meaningful when `consensus = TRUE`.", call. = FALSE)
    }
    requested <- sub("^run_", "", tolower(statistics))
    consensus_methods <- sub(
      "^norm_",
      "",
      sub("^run_", "", tolower(consensus_statistics))
    )
    if (!is.null(statistics) && any(!consensus_methods %in% requested)) {
      stop(
        "Methods named by `consensus_statistics` must be requested in ",
        "`statistics`.",
        call. = FALSE
      )
    }
  }
  if (!consensus || is.null(statistics)) {
    return(invisible(TRUE))
  }

  selected <- if (is.null(consensus_statistics)) statistics else consensus_statistics
  selected <- sub("^norm_", "", sub("^run_", "", tolower(selected)))
  enrichment <- c("aucell", "fgsea", "gsva", "ora")
  regulatory <- c("mdt", "mlm", "udt", "ulm", "viper", "wmean", "wsum")
  if (any(selected %in% enrichment) && any(selected %in% regulatory)) {
    stop(
      "Do not form one consensus across enrichment-style and network-aware ",
      "decoupleR statistics. Set `consensus = FALSE`, or explicitly choose ",
      "same-family `consensus_statistics`.",
      call. = FALSE
    )
  }
  invisible(TRUE)
}
