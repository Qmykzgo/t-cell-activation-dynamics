# =============================================================================
# scripts/02_deseq2_lrt.R
# PURPOSE : Likelihood Ratio Tests for T-cell activation dynamics.
#
#   Model 1 — Interaction LRT (PRIMARY)
#     Full:    ~ donor + condition + ns(time,df=3) + condition:ns(time,df=3)
#     Reduced: ~ donor + condition + ns(time,df=3)
#     Tests: does the SHAPE of the time trajectory differ between conditions?
#
#   Model 2 — Activated-only LRT (SECONDARY)
#     Full:    ~ donor + ns(time,df=3)
#     Reduced: ~ donor
#     Tests: which genes change over time within activated T-cells?
#
#   Model 3 — Pairwise contrasts (TERTIARY)
#     Identifies early (6h) vs sustained (24h, 72h) responders with
#     lfcShrink(type="ashr") for accurate effect size estimates.
#
# STATISTICAL JUSTIFICATION:
#   - Spline basis (df=3) captures non-linear temporal dynamics with fewer
#     parameters than a saturated factor model (Love et al., 2014).
#   - LRT is preferred over Wald tests for multi-parameter hypotheses
#     (Schulz et al., 2018, PMC5826275).
#   - glmGamPoi fitting speeds up NB dispersion estimation ~10× on laptops.
#
# INPUT  : counts_filtered (list $counts/$meta from 01_load_data.R)
# OUTPUT : dds objects (RDS), results tibbles (RDS)
# =============================================================================

library(DESeq2)
library(splines)
library(glmGamPoi)
library(dplyr)
library(tibble)
library(purrr)

# -----------------------------------------------------------------------------
# build_deseq2_base()
# Creates a DESeqDataSet with a simple formula design.
# Used ONLY for VST normalisation — NOT for the LRT hypothesis tests.
# -----------------------------------------------------------------------------
build_deseq2_base <- function(counts_filtered, sample_meta) {
  counts <- counts_filtered$counts
  meta   <- counts_filtered$meta

  meta$donor      <- droplevels(meta$donor)
  meta$condition  <- droplevels(meta$condition)
  meta$time_factor <- factor(meta$time_hours,
                              levels = sort(unique(meta$time_hours)))

  dds <- DESeqDataSetFromMatrix(
    countData = counts,
    colData   = meta,
    design    = ~ donor + condition + time_factor
  )
  dds <- estimateSizeFactors(dds)

  message(">>> Base DESeq2 object: ",
          nrow(dds), " genes × ", ncol(dds), " samples")
  dds
}


# -----------------------------------------------------------------------------
# compute_vst()
# Variance-stabilising transformation (blind=FALSE uses the base design).
# Returns a genes × samples numeric matrix for PCA, heatmaps, clustering.
# -----------------------------------------------------------------------------
compute_vst <- function(dds_base) {
  message(">>> Computing VST ...")
  vst_obj <- DESeq2::vst(dds_base, blind = FALSE)
  assay(vst_obj)
}


# -----------------------------------------------------------------------------
# run_lrt_interaction()
# PRIMARY ANALYSIS: tests condition × time interaction via LRT.
#
# Technical note: DESeq2 accepts a pre-built model.matrix as the `design`
# argument. We pass both full and reduced matrices to DESeq(test="LRT").
# The LRT statistic has df = ncol(full) - ncol(reduced) = 3 (one spline
# column per interaction term, activated × spl1/spl2/spl3).
# -----------------------------------------------------------------------------
run_lrt_interaction <- function(counts_filtered, sample_meta) {
  counts <- counts_filtered$counts
  meta   <- counts_filtered$meta

  meta$donor     <- droplevels(meta$donor)
  meta$condition <- droplevels(meta$condition)

  mm_full <- model.matrix(
    ~ donor + condition + ns(time_hours, df = 3) +
      condition:ns(time_hours, df = 3),
    data = meta
  )

  mm_reduced <- model.matrix(
    ~ donor + condition + ns(time_hours, df = 3),
    data = meta
  )

  message(">>> Interaction LRT: full model (", ncol(mm_full), " params) vs ",
          "reduced (", ncol(mm_reduced), " params) — ",
          ncol(mm_full) - ncol(mm_reduced), " interaction df")

  dds <- DESeqDataSetFromMatrix(
    countData = counts,
    colData   = meta,
    design    = mm_full
  )
  dds <- estimateSizeFactors(dds)

  message(">>> Running DESeq2 LRT — interaction model ...")
  dds <- DESeq(
    dds,
    test     = "LRT",
    reduced  = mm_reduced,
    fitType  = "glmGamPoi",
    parallel = FALSE
  )

  message(">>> Interaction LRT complete.")
  dds
}


# -----------------------------------------------------------------------------
# run_lrt_activated_only()
# SECONDARY ANALYSIS: activated samples only.
# Identifies genes that change over time within the activation program,
# regardless of control trajectory.
# -----------------------------------------------------------------------------
run_lrt_activated_only <- function(counts_filtered, sample_meta) {
  counts <- counts_filtered$counts
  meta   <- counts_filtered$meta

  act_idx    <- which(meta$condition == "activated")
  counts_act <- counts[, act_idx, drop = FALSE]
  meta_act   <- meta[act_idx, , drop = FALSE]
  meta_act$donor <- droplevels(meta_act$donor)

  message(">>> Activated-only analysis: ", ncol(counts_act), " samples across ",
          length(unique(meta_act$time_hours)), " time points")

  mm_full_act    <- model.matrix(~ donor + ns(time_hours, df = 3), data = meta_act)
  mm_reduced_act <- model.matrix(~ donor,                           data = meta_act)

  dds_act <- DESeqDataSetFromMatrix(
    countData = counts_act,
    colData   = meta_act,
    design    = mm_full_act
  )
  dds_act <- estimateSizeFactors(dds_act)

  message(">>> Running DESeq2 LRT — activated-only model ...")
  dds_act <- DESeq(
    dds_act,
    test     = "LRT",
    reduced  = mm_reduced_act,
    fitType  = "glmGamPoi",
    parallel = FALSE
  )

  message(">>> Activated-only LRT complete.")
  dds_act
}


# -----------------------------------------------------------------------------
# extract_lrt_results()
# Wraps DESeq2::results() and returns a sorted tibble.
# Note: log2FoldChange from model.matrix LRT reflects the last coefficient;
# it is used here only for filtering/sorting, NOT as a definitive effect size.
# Accurate per-timepoint LFCs come from run_pairwise_contrasts().
# -----------------------------------------------------------------------------
extract_lrt_results <- function(dds_lrt, alpha = 0.05) {
  res <- DESeq2::results(dds_lrt, alpha = alpha)

  res_tbl <- res |>
    as.data.frame() |>
    tibble::rownames_to_column("gene_id") |>
    tibble::as_tibble() |>
    dplyr::arrange(padj, dplyr::desc(abs(stat)))

  n_sig <- sum(!is.na(res_tbl$padj) & res_tbl$padj < alpha)
  message(">>> Significant genes (FDR < ", alpha, "): ", n_sig)

  res_tbl
}


# -----------------------------------------------------------------------------
# classify_deg_categories()
# Cross-references the two LRT gene lists to classify dynamics:
#
#   STIMULATION_SPECIFIC — interaction LRT sig, activated-only NOT sig
#     Interpretation: trajectory shape differs between conditions;
#     the change is not simply driven by time but requires stimulation.
#
#   GENERAL_ACTIVATION   — significant in BOTH LRTs
#     Interpretation: strong temporal dynamics exist within activation;
#     interaction LRT confirms they are not merely mirrored by control.
#
#   TEMPORAL_ONLY        — activated-only sig, interaction NOT sig
#     Interpretation: temporal dynamics in activated cells mirror control;
#     may reflect circadian/culture effects rather than stimulation.
# -----------------------------------------------------------------------------
classify_deg_categories <- function(res_interaction, res_activated,
                                     fdr_cutoff = 0.05) {
  sig_int <- dplyr::filter(res_interaction, padj < fdr_cutoff) |>
    dplyr::pull(gene_id)
  sig_act <- dplyr::filter(res_activated,    padj < fdr_cutoff) |>
    dplyr::pull(gene_id)

  stim_specific  <- setdiff(sig_int, sig_act)
  general_activ  <- intersect(sig_int, sig_act)
  temporal_only  <- setdiff(sig_act, sig_int)

  message(">>> DEG classification (FDR < ", fdr_cutoff, "):")
  message("    STIMULATION_SPECIFIC: ", length(stim_specific))
  message("    GENERAL_ACTIVATION:   ", length(general_activ))
  message("    TEMPORAL_ONLY:        ", length(temporal_only))

  list(
    stimulation_specific = stim_specific,
    general_activation   = general_activ,
    temporal_only        = temporal_only,
    all_interaction_sig  = sig_int
  )
}


# -----------------------------------------------------------------------------
# run_pairwise_contrasts()
# TERTIARY ANALYSIS: activated vs control at each time point.
# Uses a group-level factor design (~ donor + group) where
# group = paste(condition, time_hours, sep="_").
#
# LFC shrinkage: type="ashr" (Stephens 2017) is used because it:
#   1. Works with any contrast (not tied to a named coefficient)
#   2. Is appropriate when the prior is unknown
#   3. Produces s-values (local false sign rate) alongside s-values
#
# Returns: named list of tibbles, one per time point.
# -----------------------------------------------------------------------------
run_pairwise_contrasts <- function(counts_filtered, sample_meta) {
  counts <- counts_filtered$counts
  meta   <- counts_filtered$meta

  meta$donor     <- droplevels(meta$donor)
  meta$condition <- droplevels(meta$condition)
  meta$group     <- factor(paste0(meta$condition, "_", meta$time_hours, "h"))

  dds_pw <- DESeqDataSetFromMatrix(
    countData = counts,
    colData   = meta,
    design    = ~ donor + group
  )
  dds_pw <- estimateSizeFactors(dds_pw)

  message(">>> Running DESeq2 for pairwise contrasts ...")
  dds_pw <- DESeq(dds_pw, fitType = "glmGamPoi", parallel = FALSE)

  # All non-zero time points
  time_points <- sort(unique(meta$time_hours))

  contrasts_list <- list()

  for (tp in time_points) {
    grp_act  <- paste0("activated_", tp, "h")
    grp_ctrl <- paste0("control_",   tp, "h")

    # Check both levels exist (0h control may be the implicit reference)
    if (!grp_act %in% levels(meta$group) || !grp_ctrl %in% levels(meta$group)) {
      message("    Skipping ", tp, "h — group levels not found")
      next
    }

    key <- paste0("T", tp, "h_activated_vs_control")
    message("    Contrast: ", grp_act, " vs ", grp_ctrl)

    tryCatch({
      res_raw <- DESeq2::results(
        dds_pw,
        contrast = c("group", grp_act, grp_ctrl),
        alpha    = 0.05
      )

      res_shrunk <- DESeq2::lfcShrink(
        dds_pw,
        contrast = c("group", grp_act, grp_ctrl),
        res      = res_raw,
        type     = "ashr",
        quiet    = TRUE
      )

      contrasts_list[[key]] <- res_shrunk |>
        as.data.frame() |>
        tibble::rownames_to_column("gene_id") |>
        tibble::as_tibble() |>
        dplyr::mutate(time_hours = tp, contrast = key) |>
        dplyr::arrange(padj)

    }, error = function(e) {
      warning("Contrast failed at ", tp, "h: ", conditionMessage(e))
    })
  }

  n_contrasts <- length(contrasts_list)
  n_sig_total <- sum(
    purrr::map_int(contrasts_list, ~ sum(!is.na(.x$padj) & .x$padj < 0.05))
  )
  message(">>> Pairwise contrasts: ", n_contrasts, " completed, ",
          n_sig_total, " total significant genes across time points")

  contrasts_list
}
