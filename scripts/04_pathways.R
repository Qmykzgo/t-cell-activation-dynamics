# =============================================================================
# scripts/04_pathways.R
# PURPOSE : Pathway interpretation of temporal gene clusters.
#
#   1. GO Biological Process enrichment — per cluster (clusterProfiler::enrichGO)
#   2. KEGG pathway enrichment         — per cluster (clusterProfiler::enrichKEGG)
#   3. Temporal pathway matrix         — re-run enrichment at each time point
#      using DEGs from pairwise contrasts → produces a pathway × time-point
#      –log10(padj) matrix for Figure 8.
#
# GENE ID STRATEGY:
#   Gene IDs in the count matrix may be ENSEMBL or SYMBOL.
#   gene_to_entrez() tries ENSEMBL first, then SYMBOL.
#   Entrez IDs required for KEGG; enrichGO also accepts SYMBOL/ENSEMBL
#   but Entrez is most reliable.
#
# INPUT  : final_clusters, pairwise_contrasts, sample_meta
# OUTPUT : pathway_go (list), pathway_kegg (list), temporal_pathway_matrix
# =============================================================================

library(clusterProfiler)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(dplyr)
library(purrr)
library(tibble)
library(stringr)

# -----------------------------------------------------------------------------
# gene_to_entrez()
# Converts a character vector of gene IDs (ENSEMBL or SYMBOL) to Entrez IDs.
# Drops unmapped IDs. Caches the mapping for the whole organism DB.
# -----------------------------------------------------------------------------
gene_to_entrez <- function(gene_ids, db = org.Hs.eg.db) {
  gene_ids <- unique(gene_ids[!is.na(gene_ids) & gene_ids != ""])

  # ── Try ENSEMBL first ─────────────────────────────────────────────────────
  is_ensembl <- grepl("^ENSG[0-9]+", gene_ids[1], ignore.case = TRUE)

  if (is_ensembl) {
    # Strip Ensembl version suffixes ("ENSG00000000003.14" -> "ENSG00000000003")
    gene_ids <- sub("\\.[0-9]+$", "", gene_ids)
    map <- AnnotationDbi::mapIds(db, keys = gene_ids,
                                  keytype = "ENSEMBL", column = "ENTREZID",
                                  multiVals = "first")
  } else {
    map <- AnnotationDbi::mapIds(db, keys = gene_ids,
                                  keytype = "SYMBOL", column = "ENTREZID",
                                  multiVals = "first")
  }

  # ── Fallback to SYMBOL if ENSEMBL gave mostly NAs ─────────────────────────
  if (mean(is.na(map)) > 0.5 && is_ensembl) {
    sym_ids <- AnnotationDbi::mapIds(db, keys = gene_ids,
                                      keytype = "ENSEMBL", column = "SYMBOL",
                                      multiVals = "first")
    sym_ids <- sym_ids[!is.na(sym_ids)]
    map2 <- AnnotationDbi::mapIds(db, keys = as.character(sym_ids),
                                   keytype = "SYMBOL", column = "ENTREZID",
                                   multiVals = "first")
    map[names(map2)] <- map2
  }

  map <- map[!is.na(map)]
  as.character(map)
}


# -----------------------------------------------------------------------------
# run_go_enrichment_per_cluster()
# Runs GO:BP enrichment on each cluster using the full tested gene set as
# universe (recommended over genome-wide universe for focused analyses).
# Returns a named list of enrichResult objects.
# -----------------------------------------------------------------------------
run_go_enrichment_per_cluster <- function(final_clusters, ont = "BP",
                                           pval_cutoff = 0.05,
                                           qval_cutoff = 0.10,
                                           min_gs      = 10,
                                           max_gs      = 500) {
  cl_hard  <- final_clusters$cluster_hard
  cl_labels <- final_clusters$trajectory_labels
  cl_ids   <- sort(unique(cl_hard))

  # Universe = all tested (clustered) genes
  all_entrez <- gene_to_entrez(names(cl_hard))
  message("\n>>> GO:", ont, " enrichment — ", length(all_entrez),
          " background genes, ", length(cl_ids), " clusters")

  results <- purrr::map(cl_ids, function(cl) {
    g_ids   <- names(cl_hard)[cl_hard == cl]
    g_ent   <- gene_to_entrez(g_ids)
    label   <- cl_labels[as.character(cl)]

    if (length(g_ent) < 5) {
      message("  Cluster ", cl, " [", label, "]: too few genes (", length(g_ent), ") — skip")
      return(NULL)
    }

    message("  Cluster ", cl, " [", label, "]: n=", length(g_ent), " genes")

    tryCatch(
      clusterProfiler::enrichGO(
        gene          = g_ent,
        universe      = all_entrez,
        OrgDb         = org.Hs.eg.db,
        ont           = ont,
        pAdjustMethod = "BH",
        pvalueCutoff  = pval_cutoff,
        qvalueCutoff  = qval_cutoff,
        minGSSize     = min_gs,
        maxGSSize     = max_gs,
        readable      = TRUE   # convert Entrez back to gene symbols in output
      ),
      error = function(e) {
        warning("GO enrichment failed cluster ", cl, ": ", conditionMessage(e))
        NULL
      }
    )
  })

  names(results) <- paste0("C", cl_ids, "_", cl_labels[as.character(cl_ids)])

  n_enriched <- sum(!sapply(results, is.null))
  message(">>> GO enrichment done: ", n_enriched, "/", length(cl_ids), " clusters returned results")

  results
}


# -----------------------------------------------------------------------------
# run_kegg_enrichment_per_cluster()
# KEGG enrichment using Entrez IDs (required for KEGG organism code "hsa").
# Note: enrichKEGG queries the KEGG REST API — needs internet connectivity.
# -----------------------------------------------------------------------------
run_kegg_enrichment_per_cluster <- function(final_clusters,
                                              pval_cutoff = 0.05,
                                              qval_cutoff = 0.20) {
  cl_hard   <- final_clusters$cluster_hard
  cl_labels <- final_clusters$trajectory_labels
  cl_ids    <- sort(unique(cl_hard))
  all_entrez <- gene_to_entrez(names(cl_hard))

  message("\n>>> KEGG enrichment — ", length(cl_ids), " clusters")

  results <- purrr::map(cl_ids, function(cl) {
    g_ent <- gene_to_entrez(names(cl_hard)[cl_hard == cl])
    label <- cl_labels[as.character(cl)]

    if (length(g_ent) < 5) return(NULL)

    message("  Cluster ", cl, " [", label, "]: n=", length(g_ent))

    tryCatch(
      clusterProfiler::enrichKEGG(
        gene          = g_ent,
        universe      = all_entrez,
        organism      = "hsa",
        pAdjustMethod = "BH",
        pvalueCutoff  = pval_cutoff,
        qvalueCutoff  = qval_cutoff
      ),
      error = function(e) {
        warning("KEGG failed cluster ", cl, ": ", conditionMessage(e))
        NULL
      }
    )
  })

  names(results) <- paste0("C", cl_ids, "_", cl_labels[as.character(cl_ids)])
  results
}


# -----------------------------------------------------------------------------
# build_temporal_pathway_matrix()
# For each non-zero time point, runs a quick GO:BP enrichment on the
# significant DEGs from the pairwise contrast at that time.
# Assembles the top n_top_pathways into a matrix of –log10(padj) values
# across time, suitable for the temporal heatmap (Figure 8).
# -----------------------------------------------------------------------------
build_temporal_pathway_matrix <- function(pathway_go,
                                           pairwise_contrasts,
                                           sample_meta,
                                           n_top_pathways = 25,
                                           fdr_cutoff     = 0.05) {
  message("\n>>> Building temporal pathway enrichment matrix ...")

  time_keys <- names(pairwise_contrasts)   # e.g. "T6h_activated_vs_control"
  time_pts  <- as.numeric(str_extract(time_keys, "[0-9]+"))
  names(time_keys) <- paste0("T", time_pts, "h")

  # ── Run quick enrichment per time point ───────────────────────────────────
  tp_enrich <- purrr::imap(time_keys, function(key, tp_name) {
    tp_res <- pairwise_contrasts[[key]]
    if (is.null(tp_res)) return(NULL)

    sig_genes <- dplyr::filter(tp_res, !is.na(padj), padj < fdr_cutoff) |>
      dplyr::pull(gene_id)
    all_genes <- dplyr::filter(tp_res, !is.na(padj)) |> dplyr::pull(gene_id)

    sig_ent <- gene_to_entrez(sig_genes)
    all_ent <- gene_to_entrez(all_genes)

    if (length(sig_ent) < 10) return(NULL)

    message("  ", tp_name, ": ", length(sig_ent), " significant genes")

    tryCatch(
      clusterProfiler::enrichGO(
        gene          = sig_ent,
        universe      = all_ent,
        OrgDb         = org.Hs.eg.db,
        ont           = "BP",
        pAdjustMethod = "BH",
        pvalueCutoff  = 0.10,    # relaxed for discovery
        minGSSize     = 10,
        readable      = TRUE
      ),
      error = function(e) NULL
    )
  })

  # ── Collect top pathways across all time points ───────────────────────────
  top_paths <- purrr::map_dfr(tp_enrich, function(enr) {
    if (is.null(enr) || nrow(enr@result) == 0) return(NULL)
    dplyr::arrange(enr@result, p.adjust) |>
      dplyr::slice_head(n = 12) |>
      dplyr::select(Description)
  }) |>
    dplyr::distinct(Description) |>
    dplyr::slice_head(n = n_top_pathways)

  if (nrow(top_paths) == 0) {
    warning("No pathway results — temporal matrix will be empty.")
    return(matrix(0, nrow = 0, ncol = 0))
  }

  path_names <- top_paths$Description

  # ── Build –log10(padj) matrix ─────────────────────────────────────────────
  pw_matrix <- matrix(
    0,
    nrow     = length(path_names),
    ncol     = length(tp_enrich),
    dimnames = list(path_names, names(tp_enrich))
  )

  for (tp_name in names(tp_enrich)) {
    enr <- tp_enrich[[tp_name]]
    if (is.null(enr)) next
    df <- enr@result

    for (pth in path_names) {
      hit <- which(df$Description == pth)
      if (length(hit) > 0) {
        pw_matrix[pth, tp_name] <- -log10(pmax(df$p.adjust[hit[1]], 1e-10))
      }
    }
  }

  # Trim pathways with no signal anywhere
  pw_matrix <- pw_matrix[rowSums(pw_matrix) > 0, , drop = FALSE]

  message(">>> Temporal pathway matrix: ", nrow(pw_matrix), " pathways × ",
          ncol(pw_matrix), " time points")
  pw_matrix
}
