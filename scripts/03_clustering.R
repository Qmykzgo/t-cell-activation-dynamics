# =============================================================================
# scripts/03_clustering.R
# PURPOSE : Temporal pattern analysis of top interaction DEGs.
#
#   Three independent clustering methods applied to the same gene set:
#   1. Mfuzz      — fuzzy c-means on z-scored mean time-course expression
#   2. NMF        — non-negative matrix factorisation (brunet algorithm)
#   3. Spline k-means — k-means on cubic spline coefficient vectors
#
#   Concordance measured with Adjusted Rand Index (ARI).
#   Mfuzz used as primary assignment (soft membership interpretable biologically).
#   Clusters labelled by trajectory shape against expected T-cell biology:
#     "Early Innate Burst"      — peak ≤ 6h, rapid decline
#     "Sustained Inflammatory"  — monotonic rise, plateau ≥ 24h
#     "Delayed Proliferative"   — peak 48–72h (cell cycle genes)
#     "Transient Metabolic"     — peak 12–24h, return to baseline
#     "Repressed Quiescence"    — sustained downregulation (FOXO targets)
#
# INPUT  : vst_matrix (genes×samples), top_deg_genes (character vector),
#          sample_meta (data.frame)
# OUTPUT : cluster_mfuzz, cluster_nmf, cluster_kmeans_spline, cluster_comparison,
#          final_clusters (all list objects saved as RDS by targets)
# =============================================================================

library(Mfuzz)
library(NMF)
library(mclust)
library(splines)
library(dplyr)
library(tidyr)
library(purrr)
library(tibble)
library(Biobase)

# Human-readable trajectory labels (index matches cluster integer)
TRAJECTORY_LABELS <- c(
  "Early Innate Burst",
  "Sustained Inflammatory",
  "Delayed Proliferative",
  "Transient Metabolic",
  "Repressed Quiescence"
)


# -----------------------------------------------------------------------------
# get_top_interaction_genes()
# Selects high-confidence interaction DEGs for clustering.
# LFC threshold guards against statistically significant but tiny effects.
# Cap at max_genes to keep Mfuzz and NMF memory-efficient.
# -----------------------------------------------------------------------------
get_top_interaction_genes <- function(res_interaction,
                                       fdr_cutoff = 0.05,
                                       lfc_cutoff = 1.0,
                                       max_genes  = 2000) {
  top <- res_interaction |>
    dplyr::filter(
      !is.na(padj),
      padj < fdr_cutoff,
      !is.na(log2FoldChange),
      abs(log2FoldChange) >= lfc_cutoff
    ) |>
    dplyr::arrange(padj) |>
    dplyr::slice_head(n = max_genes) |>
    dplyr::pull(gene_id)

  message(">>> ", length(top), " genes selected for clustering ",
          "(FDR < ", fdr_cutoff, " | |LFC| ≥ ", lfc_cutoff, " | top ", max_genes, ")")
  top
}


# -----------------------------------------------------------------------------
# .make_mean_time_matrix()  [internal]
# Computes per-time-point mean VST expression for a given condition.
# Rows = genes, cols = time points (named "T0h", "T6h", …).
# -----------------------------------------------------------------------------
.make_mean_time_matrix <- function(vst_matrix, gene_ids, sample_meta,
                                    condition_use = "activated") {

  act_samples   <- rownames(sample_meta)[sample_meta$condition == condition_use]
  common_genes  <- intersect(gene_ids, rownames(vst_matrix))
  common_samples <- intersect(act_samples, colnames(vst_matrix))

  expr   <- vst_matrix[common_genes, common_samples, drop = FALSE]
  meta_s <- sample_meta[common_samples, , drop = FALSE]
  tps    <- sort(unique(meta_s$time_hours))

  mean_mat <- vapply(tps, function(tp) {
    tp_cols <- rownames(meta_s)[meta_s$time_hours == tp]
    rowMeans(expr[, tp_cols, drop = FALSE])
  }, numeric(length(common_genes)))

  colnames(mean_mat) <- paste0("T", tps, "h")
  rownames(mean_mat) <- common_genes
  mean_mat
}


# =============================================================================
# METHOD 1: Mfuzz soft fuzzy c-means
# =============================================================================
run_mfuzz_clustering <- function(vst_matrix, top_deg_genes, sample_meta,
                                  n_clusters = 5, m_param = 1.25) {
  message("\n>>> [Clustering 1/3] Mfuzz (c=", n_clusters, ", m=", m_param, ") ...")

  mean_mat <- .make_mean_time_matrix(vst_matrix, top_deg_genes, sample_meta)

  # Mfuzz requires ExpressionSet
  eset <- Biobase::ExpressionSet(assayData = mean_mat)
  eset <- Mfuzz::filter.NA(eset,  thres = 0.25)
  eset <- Mfuzz::fill.NA(eset,    mode  = "mean")
  eset <- Mfuzz::filter.std(eset, min.std = 0.001)
  eset_std <- Mfuzz::standardise(eset)

  set.seed(2024)
  cl <- Mfuzz::mfuzz(eset_std, c = n_clusters, m = m_param)

  # Hard assignment = cluster with maximum membership score
  cluster_hard <- apply(cl$membership, 1, which.max)
  names(cluster_hard) <- rownames(Biobase::exprs(eset_std))

  message(">>> Mfuzz cluster sizes: ", paste(table(cluster_hard), collapse = " | "))

  list(
    method         = "mfuzz",
    n_clusters     = n_clusters,
    cluster_hard   = cluster_hard,          # named int: gene → cluster
    membership     = cl$membership,          # soft membership matrix
    cluster_centers = cl$centers,
    cluster_object = cl,
    mean_matrix    = Biobase::exprs(eset_std),   # standardised (for plotting)
    raw_mean_mat   = mean_mat                # VST scale (for labelling)
  )
}


# =============================================================================
# METHOD 2: NMF (non-negative matrix factorisation)
# =============================================================================
run_nmf_clustering <- function(vst_matrix, top_deg_genes, sample_meta,
                                 n_clusters = 5, n_runs = 50) {
  message("\n>>> [Clustering 2/3] NMF (rank=", n_clusters, ", runs=", n_runs, ") ...")

  mean_mat <- .make_mean_time_matrix(vst_matrix, top_deg_genes, sample_meta)

  # NMF requires non-negative values: shift so minimum = ε
  min_val      <- min(mean_mat, na.rm = TRUE)
  mean_mat_nn  <- mean_mat - min_val + 1e-6

  set.seed(2024)
  nmf_res <- NMF::nmf(
    mean_mat_nn,
    rank     = n_clusters,
    method   = "brunet",   # Lee & Seung 2001, robust for gene expression
    nrun     = n_runs,
    seed     = "nndsvd",   # deterministic initialisation, aids convergence
    .options = "v"
  )

  W <- NMF::basis(nmf_res)           # genes × rank
  H <- NMF::coef(nmf_res)            # rank × time_points

  # Gene assigned to the NMF component in which it has the highest loading
  cluster_hard <- apply(W, 1, which.max)
  names(cluster_hard) <- rownames(mean_mat)

  message(">>> NMF cluster sizes: ", paste(table(cluster_hard), collapse = " | "))

  list(
    method        = "nmf",
    n_clusters    = n_clusters,
    cluster_hard  = cluster_hard,
    W_matrix      = W,
    H_matrix      = H,
    nmf_object    = nmf_res,
    mean_matrix   = mean_mat
  )
}


# =============================================================================
# METHOD 3: k-means on cubic spline coefficients
# =============================================================================
# Rationale: fitting a cubic spline to each gene's time course reduces its
# temporal profile to 4 numbers (intercept + 3 basis coefficients).
# Clustering in coefficient space groups genes by trajectory SHAPE, not
# just by expression level — a more biologically meaningful distance.
# =============================================================================
run_spline_kmeans <- function(vst_matrix, top_deg_genes, sample_meta,
                               n_clusters = 5, nstart = 200) {
  message("\n>>> [Clustering 3/3] Spline k-means (k=", n_clusters,
          ", nstart=", nstart, ") ...")

  mean_mat   <- .make_mean_time_matrix(vst_matrix, top_deg_genes, sample_meta)
  time_pts   <- as.numeric(gsub("[^0-9]", "", colnames(mean_mat)))

  # Fit cubic spline per gene → extract 4 coefficients
  fit_spline_coefs <- function(expr_vec) {
    df_g <- data.frame(t = time_pts, y = as.numeric(expr_vec))
    fit  <- lm(y ~ splines::ns(t, df = 3), data = df_g)
    as.numeric(coef(fit))     # length 4: intercept + spl1 + spl2 + spl3
  }

  message("    Fitting splines to ", nrow(mean_mat), " genes ...")
  coef_mat <- t(apply(mean_mat, 1, fit_spline_coefs))
  colnames(coef_mat) <- c("intercept", "spl1", "spl2", "spl3")

  # Scale each coefficient dimension before clustering
  coef_scaled <- scale(coef_mat)

  set.seed(2024)
  km <- kmeans(coef_scaled, centers = n_clusters,
               nstart = nstart, iter.max = 1000)

  cluster_hard <- km$cluster
  names(cluster_hard) <- rownames(mean_mat)

  message(">>> Spline k-means cluster sizes: ",
          paste(table(cluster_hard), collapse = " | "))

  list(
    method        = "kmeans_spline",
    n_clusters    = n_clusters,
    cluster_hard  = cluster_hard,
    kmeans_object = km,
    coef_matrix   = coef_mat,
    mean_matrix   = mean_mat
  )
}


# =============================================================================
# compare_clustering_solutions()
# Pairwise ARI between all three methods.
# ARI = 1 → perfect agreement; ARI = 0 → chance level.
# =============================================================================
compare_clustering_solutions <- function(cluster_mfuzz, cluster_nmf,
                                          cluster_kmeans_spline) {
  message("\n>>> Comparing clustering solutions (Adjusted Rand Index) ...")

  common <- Reduce(
    intersect,
    list(names(cluster_mfuzz$cluster_hard),
         names(cluster_nmf$cluster_hard),
         names(cluster_kmeans_spline$cluster_hard))
  )
  message("    Common genes for comparison: ", length(common))

  cl_mf <- cluster_mfuzz$cluster_hard[common]
  cl_nm <- cluster_nmf$cluster_hard[common]
  cl_km <- cluster_kmeans_spline$cluster_hard[common]

  ari_mf_nm <- mclust::adjustedRandIndex(cl_mf, cl_nm)
  ari_mf_km <- mclust::adjustedRandIndex(cl_mf, cl_km)
  ari_nm_km <- mclust::adjustedRandIndex(cl_nm, cl_km)

  ari_mat <- matrix(
    c(1,        ari_mf_nm, ari_mf_km,
      ari_mf_nm, 1,        ari_nm_km,
      ari_mf_km, ari_nm_km, 1),
    nrow     = 3,
    dimnames = list(c("Mfuzz","NMF","Kmeans_Spline"),
                    c("Mfuzz","NMF","Kmeans_Spline"))
  )

  message(">>> ARI matrix:")
  print(round(ari_mat, 3))

  list(
    ari_matrix    = ari_mat,
    common_genes  = common,
    labels_mfuzz  = cl_mf,
    labels_nmf    = cl_nm,
    labels_kmeans = cl_km,
    crosstab_mf_nm = table(Mfuzz = cl_mf, NMF      = cl_nm),
    crosstab_mf_km = table(Mfuzz = cl_mf, KMeans   = cl_km)
  )
}


# =============================================================================
# build_consensus_clusters()
# Adopts Mfuzz as the canonical solution (highest biological interpretability
# via soft memberships); assigns trajectory labels automatically by
# analysing the shape of each cluster's mean time-course.
# =============================================================================
build_consensus_clusters <- function(cluster_mfuzz, cluster_comparison) {
  message("\n>>> Building consensus cluster labels ...")

  cl_hard    <- cluster_mfuzz$cluster_hard
  raw_mat    <- cluster_mfuzz$raw_mean_mat
  time_pts   <- as.numeric(gsub("[^0-9]", "", colnames(raw_mat)))
  cl_ids     <- sort(unique(cl_hard))

  # Compute per-cluster mean trajectory (in original VST scale)
  cl_means <- vapply(cl_ids, function(cl) {
    g <- names(cl_hard)[cl_hard == cl]
    g <- intersect(g, rownames(raw_mat))
    colMeans(raw_mat[g, , drop = FALSE])
  }, numeric(ncol(raw_mat)))
  rownames(cl_means) <- colnames(raw_mat)
  colnames(cl_means) <- as.character(cl_ids)

  # Heuristic trajectory labelling ──────────────────────────────────────────
  label_trajectory <- function(traj, tps) {
    t0_val    <- traj[which(tps == 0)]
    if (length(t0_val) == 0) t0_val <- traj[1]
    peak_tp   <- tps[which.max(traj)]
    trough_tp <- tps[which.min(traj)]
    early_mean <- mean(traj[tps <= 12])
    late_mean  <- mean(traj[tps >= 48])

    overall_up <- mean(traj) > t0_val

    if (!overall_up) {
      "Repressed Quiescence"
    } else if (peak_tp <= 6) {
      "Early Innate Burst"
    } else if (peak_tp >= 48 && late_mean > early_mean) {
      "Delayed Proliferative"
    } else if (peak_tp %in% c(12, 24) && late_mean < early_mean * 1.1) {
      "Transient Metabolic"
    } else {
      "Sustained Inflammatory"
    }
  }

  traj_labels <- setNames(
    vapply(cl_ids, function(cl) {
      label_trajectory(cl_means[, as.character(cl)], time_pts)
    }, character(1)),
    as.character(cl_ids)
  )

  message(">>> Trajectory labels:")
  for (cl in cl_ids) {
    n_g <- sum(cl_hard == cl)
    message("    Cluster ", cl, " (n=", n_g, "): ", traj_labels[as.character(cl)])
  }

  list(
    cluster_hard      = cl_hard,
    trajectory_labels = traj_labels,
    cluster_means     = cl_means,
    mean_matrix       = raw_mat,
    membership        = cluster_mfuzz$membership
  )
}
