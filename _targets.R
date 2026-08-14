# =============================================================================
# _targets.R
# T-cell Activation Dynamics — Bulk RNA-Seq Time-Course Pipeline
# Run the entire analysis: targets::tar_make()
# Visualise the DAG:        targets::tar_visnetwork()
# =============================================================================

library(targets)
library(tarchetypes)

tar_option_set(
  packages = c(
    # Bioconductor
    "DESeq2", "GEOquery", "SummarizedExperiment", "Biobase",
    "BiocParallel", "org.Hs.eg.db", "AnnotationDbi",
    "clusterProfiler", "glmGamPoi", "apeglm", "ashr",
    # CRAN
    "dplyr", "tidyr", "purrr", "tibble", "stringr", "readr",
    "ggplot2", "ggrepel", "patchwork", "pheatmap", "RColorBrewer",
    "scales", "Mfuzz", "NMF", "mclust", "splines",
    "edgeR", "matrixStats"
  ),
  memory           = "transient",   # free RAM between targets
  garbage_collection = TRUE,
  seed             = 2024
)

# Source all function scripts
source("scripts/01_load_data.R")
source("scripts/02_deseq2_lrt.R")
source("scripts/03_clustering.R")
source("scripts/04_pathways.R")
source("scripts/05_visualization.R")

# =============================================================================
list(

  # ── 1. Configuration ────────────────────────────────────────────────────────
  tar_target(gse_id,       "GSE197067"),
  tar_target(data_dir,     "data/"),
  tar_target(figures_dir,  "figures/"),
  tar_target(results_dir,  "results/"),

  # ── 2. Data Ingestion ────────────────────────────────────────────────────────
  tar_target(
    raw_data,
    download_geo_data(gse_id, data_dir),
    format = "rds"
  ),
  tar_target(
    counts_raw,
    extract_count_matrix(raw_data),
    format = "rds"
  ),
  tar_target(
    sample_meta,
    build_sample_metadata(raw_data),
    format = "rds"
  ),
  tar_target(
    counts_filtered,
    filter_low_counts(counts_raw, sample_meta, min_cpm = 1, min_samples = 3),
    format = "rds"
  ),

  # ── 3. DESeq2 — Base object for VST ─────────────────────────────────────────
  tar_target(
    dds_base,
    build_deseq2_base(counts_filtered, sample_meta),
    format = "rds"
  ),
  tar_target(
    vst_matrix,
    compute_vst(dds_base),
    format = "rds"
  ),

  # ── 4. DESeq2 — Primary LRT: condition × time interaction ────────────────────
  tar_target(
    dds_lrt_interaction,
    run_lrt_interaction(counts_filtered, sample_meta),
    format = "rds"
  ),
  tar_target(
    res_interaction,
    extract_lrt_results(dds_lrt_interaction),
    format = "rds"
  ),

  # ── 5. DESeq2 — Secondary LRT: activated samples only ───────────────────────
  tar_target(
    dds_lrt_activated,
    run_lrt_activated_only(counts_filtered, sample_meta),
    format = "rds"
  ),
  tar_target(
    res_activated,
    extract_lrt_results(dds_lrt_activated),
    format = "rds"
  ),

  # ── 6. DEG classification ───────────────────────────────────────────────────
  tar_target(
    deg_categories,
    classify_deg_categories(res_interaction, res_activated),
    format = "rds"
  ),

  # ── 7. Tertiary: pairwise contrasts at each time point ───────────────────────
  tar_target(
    pairwise_contrasts,
    run_pairwise_contrasts(counts_filtered, sample_meta),
    format = "rds"
  ),

  # ── 8. Temporal clustering ──────────────────────────────────────────────────
  tar_target(
    top_deg_genes,
    get_top_interaction_genes(res_interaction, fdr_cutoff = 0.05, lfc_cutoff = 1),
    format = "rds"
  ),
  tar_target(
    cluster_mfuzz,
    run_mfuzz_clustering(vst_matrix, top_deg_genes, sample_meta, n_clusters = 5),
    format = "rds"
  ),
  tar_target(
    cluster_nmf,
    run_nmf_clustering(vst_matrix, top_deg_genes, sample_meta, n_clusters = 5),
    format = "rds"
  ),
  tar_target(
    cluster_kmeans_spline,
    run_spline_kmeans(vst_matrix, top_deg_genes, sample_meta, n_clusters = 5),
    format = "rds"
  ),
  tar_target(
    cluster_comparison,
    compare_clustering_solutions(cluster_mfuzz, cluster_nmf, cluster_kmeans_spline),
    format = "rds"
  ),
  tar_target(
    final_clusters,
    build_consensus_clusters(cluster_mfuzz, cluster_comparison),
    format = "rds"
  ),

  # ── 9. Pathway enrichment ───────────────────────────────────────────────────
  tar_target(
    pathway_go,
    run_go_enrichment_per_cluster(final_clusters, ont = "BP"),
    format = "rds"
  ),
  tar_target(
    pathway_kegg,
    run_kegg_enrichment_per_cluster(final_clusters),
    format = "rds"
  ),
  tar_target(
    temporal_pathway_matrix,
    build_temporal_pathway_matrix(pathway_go, pairwise_contrasts, sample_meta),
    format = "rds"
  ),

  # ── 10. Figures ─────────────────────────────────────────────────────────────
  tar_target(
    fig_pca,
    plot_pca(vst_matrix, sample_meta, figures_dir),
    format = "file"
  ),
  tar_target(
    fig_sample_dist,
    plot_sample_distance_heatmap(vst_matrix, sample_meta, figures_dir),
    format = "file"
  ),
  tar_target(
    fig_volcano,
    plot_volcano_panels(pairwise_contrasts, figures_dir),
    format = "file"
  ),
  tar_target(
    fig_ma,
    plot_ma_interaction(res_interaction, figures_dir),
    format = "file"
  ),
  tar_target(
    fig_trajectories,
    plot_gene_trajectories(vst_matrix, final_clusters, sample_meta, figures_dir),
    format = "file"
  ),
  tar_target(
    fig_cluster_heatmap,
    plot_cluster_heatmap(vst_matrix, final_clusters, sample_meta, figures_dir),
    format = "file"
  ),
  tar_target(
    fig_pathway_dotplot,
    plot_pathway_dotplots(pathway_go, figures_dir),
    format = "file"
  ),
  tar_target(
    fig_temporal_pathway,
    plot_temporal_pathway_heatmap(temporal_pathway_matrix, figures_dir),
    format = "file"
  ),
  tar_target(
    fig_ari,
    plot_cluster_concordance(cluster_comparison, figures_dir),
    format = "file"
  ),

  # ── 11. Result tables ───────────────────────────────────────────────────────
  tar_target(
    table_deg,
    save_deg_table(res_interaction, deg_categories, results_dir),
    format = "file"
  ),
  tar_target(
    table_clusters,
    save_cluster_table(final_clusters, results_dir),
    format = "file"
  ),
  tar_target(
    table_pathways,
    save_pathway_table(pathway_go, pathway_kegg, results_dir),
    format = "file"
  )
)
