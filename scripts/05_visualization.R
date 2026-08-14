# =============================================================================
# scripts/05_visualization.R
# PURPOSE : Generate all publication-quality figures (300 DPI PNG + PDF).
#
#   Fig 01 — PCA coloured by condition × time with confidence ellipses
#   Fig 02 — Sample distance heatmap (time | condition | donor annotations)
#   Fig 03 — Volcano plots at 6h, 24h, 72h (multi-panel)
#   Fig 04 — MA plot for interaction LRT
#   Fig 05 — Gene trajectory plots per cluster (mean ± SEM + faint gene lines)
#   Fig 06 — Cluster heatmap (genes ordered by cluster, samples by time)
#   Fig 07 — Pathway enrichment dot plots faceted by cluster
#   Fig 08 — Temporal pathway enrichment heatmap
#   Fig 09 — ARI concordance matrix (clustering method comparison)
#
# STYLE CONVENTIONS:
#   - theme_tcell() extends theme_classic() — clean, publication-ready
#   - All colour palettes are colourblind-safe (ColorBrewer / manual)
#   - Fonts ≥ 10 pt; axis labels clearly described
#   - Every figure saved as both PNG (embedded in README/Quarto) and PDF/SVG
#     (for journal submission at arbitrary resolution)
#
# INPUT  : targets objects (vst_matrix, sample_meta, results, clusters, pathways)
# OUTPUT : PNG + PDF files written to figures/
# =============================================================================

library(ggplot2)
library(ggrepel)
library(patchwork)
library(pheatmap)
library(RColorBrewer)
library(dplyr)
library(tidyr)
library(purrr)
library(tibble)
library(stringr)
library(scales)
library(matrixStats)
library(grid)

# ── Shared design constants ────────────────────────────────────────────────────
CONDITION_COLORS <- c(activated = "#E64B35FF", control = "#4DBBD5FF")

TIME_COLORS <- setNames(
  RColorBrewer::brewer.pal(6, "YlOrRd"),
  c("0", "6", "12", "24", "48", "72")
)

CLUSTER_COLORS <- c(
  "Early Innate Burst"      = "#F77F00",
  "Sustained Inflammatory"  = "#D62828",
  "Delayed Proliferative"   = "#7B2D8B",
  "Transient Metabolic"     = "#2DC653",
  "Repressed Quiescence"    = "#1A6AB8"
)

DONOR_COLORS <- setNames(
  RColorBrewer::brewer.pal(8, "Set2"),
  paste0("D", 1:8)
)

TIME_SHAPES <- setNames(
  c(16L, 17L, 15L, 18L, 3L, 4L),
  c("0", "6", "12", "24", "48", "72")
)

# ── Publication theme ──────────────────────────────────────────────────────────
theme_tcell <- function(base_size = 12) {
  theme_classic(base_size = base_size) +
    theme(
      plot.title       = element_text(face = "bold", size = base_size + 2,
                                       margin = margin(b = 4)),
      plot.subtitle    = element_text(colour = "grey40", size = base_size - 1,
                                       margin = margin(b = 6)),
      plot.caption     = element_text(colour = "grey55", size = base_size - 3,
                                       hjust = 0),
      axis.title       = element_text(size = base_size),
      axis.text        = element_text(size = base_size - 1),
      legend.title     = element_text(face = "bold", size = base_size - 1),
      legend.text      = element_text(size = base_size - 2),
      strip.background = element_rect(fill = "grey92", colour = NA),
      strip.text       = element_text(face = "bold", size = base_size - 1),
      panel.grid.major = element_line(colour = "grey93", linewidth = 0.4),
      panel.grid.minor = element_blank()
    )
}

# ── I/O helpers ───────────────────────────────────────────────────────────────
.save_gg <- function(p, stem, dir, w = 10, h = 7, dpi = 300) {
  dir.create(dir, showWarnings = FALSE, recursive = TRUE)
  png_f <- file.path(dir, paste0(stem, ".png"))
  pdf_f <- file.path(dir, paste0(stem, ".pdf"))
  ggplot2::ggsave(png_f, p, width = w, height = h, dpi = dpi, bg = "white")
  ggplot2::ggsave(pdf_f, p, width = w, height = h,
                  device = grDevices::cairo_pdf, bg = "white")
  message(">>> Saved: ", basename(png_f))
  png_f
}

.save_ph <- function(ph_obj, stem, dir, w = 10, h = 8, dpi = 300) {
  dir.create(dir, showWarnings = FALSE, recursive = TRUE)
  png_f <- file.path(dir, paste0(stem, ".png"))
  pdf_f <- file.path(dir, paste0(stem, ".pdf"))

  grDevices::png(png_f, width = w, height = h, units = "in", res = dpi)
  grid.newpage(); grid.draw(ph_obj$gtable)
  grDevices::dev.off()

  grDevices::pdf(pdf_f, width = w, height = h)
  grid.newpage(); grid.draw(ph_obj$gtable)
  grDevices::dev.off()

  message(">>> Saved: ", basename(png_f))
  png_f
}

.align_samples <- function(vst_matrix, sample_meta) {
  shared <- intersect(colnames(vst_matrix), rownames(sample_meta))
  list(vst  = vst_matrix[, shared, drop = FALSE],
       meta = sample_meta[shared, , drop = FALSE])
}


# =============================================================================
# FIG 01 — PCA
# =============================================================================
plot_pca <- function(vst_matrix, sample_meta, figures_dir) {
  message("\n>>> [Fig 01] PCA ...")
  d <- .align_samples(vst_matrix, sample_meta)

  rv     <- matrixStats::rowVars(d$vst)
  top_n  <- min(500L, length(rv))
  pca    <- stats::prcomp(t(d$vst[order(rv, decreasing = TRUE)[seq_len(top_n)], ]),
                           scale. = TRUE, center = TRUE)
  pct    <- round(100 * summary(pca)$importance[2, 1:4], 1)

  pca_df <- as.data.frame(pca$x[, 1:3]) |>
    tibble::rownames_to_column("sample_id") |>
    dplyr::left_join(d$meta, by = "sample_id") |>
    dplyr::mutate(
      time_chr = as.character(time_hours),
      size_val = rescale(time_hours, to = c(2, 6))
    )

  make_pca_panel <- function(x_pc, y_pc, x_pct, y_pct, show_ellipse = TRUE) {
    p <- ggplot(pca_df,
                aes(x = .data[[x_pc]], y = .data[[y_pc]],
                    colour = condition, shape = time_chr)) +
      geom_point(aes(size = size_val), alpha = 0.85) +
      scale_colour_manual(name = "Condition", values = CONDITION_COLORS) +
      scale_shape_manual(name = "Time (h)", values = TIME_SHAPES) +
      scale_size_identity() +
      labs(x = paste0(x_pc, " (", x_pct, "%)"),
           y = paste0(y_pc, " (", y_pct, "%)")) +
      theme_tcell(11)

    if (show_ellipse)
      p <- p + stat_ellipse(aes(group = condition), type = "norm",
                             level = 0.90, linetype = "dashed",
                             linewidth = 0.6, alpha = 0.5)
    p + geom_text_repel(aes(label = donor), colour = "grey30",
                        size = 2.5, max.overlaps = 20,
                        box.padding = 0.3, segment.size = 0.25)
  }

  p12 <- make_pca_panel("PC1","PC2", pct[1], pct[2]) +
    labs(title = "PCA — PC1 vs PC2",
         subtitle = paste0("Top ", top_n, " variable genes; ellipses = 90% CI per condition"))
  p13 <- make_pca_panel("PC1","PC3", pct[1], pct[3], show_ellipse = FALSE) +
    labs(title = "PCA — PC1 vs PC3")

  p <- (p12 | p13) +
    patchwork::plot_annotation(
      caption = paste0("VST-transformed counts; n = ", ncol(d$vst), " samples"),
      theme   = theme_tcell(11)
    )
  .save_gg(p, "fig01_pca", figures_dir, w = 14, h = 6)
}


# =============================================================================
# FIG 02 — Sample Distance Heatmap
# =============================================================================
plot_sample_distance_heatmap <- function(vst_matrix, sample_meta, figures_dir) {
  message("\n>>> [Fig 02] Sample distance heatmap ...")
  d      <- .align_samples(vst_matrix, sample_meta)
  dists  <- stats::dist(t(d$vst))
  dmat   <- as.matrix(dists)

  ann <- data.frame(
    Condition = d$meta$condition,
    Time_h    = factor(d$meta$time_hours),
    Donor     = d$meta$donor,
    row.names = rownames(d$meta)
  )
  n_donors <- nlevels(d$meta$donor)
  ann_col  <- list(
    Condition = CONDITION_COLORS,
    Time_h    = setNames(RColorBrewer::brewer.pal(6, "YlOrRd"),
                          levels(ann$Time_h)),
    Donor     = setNames(
      RColorBrewer::brewer.pal(max(n_donors, 3), "Set2")[seq_len(n_donors)],
      levels(d$meta$donor))
  )

  ph <- pheatmap::pheatmap(
    dmat,
    annotation_col    = ann,
    annotation_row    = ann,
    annotation_colors = ann_col,
    color             = colorRampPalette(
      rev(RColorBrewer::brewer.pal(9, "Blues")))(100),
    clustering_distance_rows = dists,
    clustering_distance_cols = dists,
    show_colnames     = FALSE,
    fontsize_row      = 7,
    main              = "Sample-to-Sample Euclidean Distance (VST)",
    silent            = TRUE
  )
  .save_ph(ph, "fig02_sample_distance_heatmap", figures_dir, w = 11, h = 9)
}


# =============================================================================
# FIG 03 — Volcano Plots (multi-panel)
# =============================================================================
plot_volcano_panels <- function(pairwise_contrasts, figures_dir,
                                 time_show  = c(6L, 24L, 72L),
                                 lfc_cut    = 1,
                                 fdr_cut    = 0.05,
                                 n_label    = 14) {
  message("\n>>> [Fig 03] Volcano plots at ", paste(time_show, collapse=", "), "h ...")

  one_volcano <- function(df, tp) {
    df <- df |>
      dplyr::filter(!is.na(padj), !is.na(log2FoldChange)) |>
      dplyr::mutate(
        nlp    = pmin(-log10(padj + 1e-300), 55),
        status = dplyr::case_when(
          padj < fdr_cut & log2FoldChange >  lfc_cut ~ "Up",
          padj < fdr_cut & log2FoldChange < -lfc_cut ~ "Down",
          TRUE ~ "NS"
        )
      )

    top_lab <- dplyr::bind_rows(
      dplyr::filter(df, status == "Up")   |> dplyr::arrange(padj) |> dplyr::slice_head(n = ceiling(n_label/2)),
      dplyr::filter(df, status == "Down") |> dplyr::arrange(padj) |> dplyr::slice_head(n = floor(n_label/2))
    )

    n_up   <- sum(df$status == "Up")
    n_down <- sum(df$status == "Down")

    ggplot(df, aes(log2FoldChange, nlp)) +
      geom_point(aes(colour = status), size = 0.7, alpha = 0.45, stroke = 0) +
      geom_point(data = top_lab, aes(colour = status), size = 1.4) +
      geom_text_repel(data = top_lab, aes(label = gene_id), size = 2.3,
                      max.overlaps = 18, box.padding = 0.35,
                      segment.size = 0.25, segment.alpha = 0.5) +
      geom_vline(xintercept = c(-lfc_cut, lfc_cut),
                 linetype = "dashed", colour = "grey55", linewidth = 0.45) +
      geom_hline(yintercept = -log10(fdr_cut),
                 linetype = "dashed", colour = "grey55", linewidth = 0.45) +
      scale_colour_manual(values = c(Up="#E64B35",Down="#4DBBD5",NS="grey72"),
                          guide = "none") +
      annotate("text",
               x = max(df$log2FoldChange, na.rm=TRUE) * 0.88,
               y = max(df$nlp, na.rm=TRUE) * 0.96,
               label = paste0("\u2191 ", n_up),
               colour = "#E64B35", fontface = "bold", size = 3.5) +
      annotate("text",
               x = min(df$log2FoldChange, na.rm=TRUE) * 0.88,
               y = max(df$nlp, na.rm=TRUE) * 0.96,
               label = paste0("\u2193 ", n_down),
               colour = "#4DBBD5", fontface = "bold", size = 3.5) +
      labs(title = paste0("Activated vs Control — ", tp, "h"),
           x = expression(log[2]~"Fold Change  (ashr-shrunk)"),
           y = expression(-log[10]~"(adj. p-value)")) +
      theme_tcell(10)
  }

  panels <- list()
  for (tp in time_show) {
    key <- paste0("T", tp, "h_activated_vs_control")
    if (key %in% names(pairwise_contrasts))
      panels[[as.character(tp)]] <- one_volcano(pairwise_contrasts[[key]], tp)
  }
  if (length(panels) == 0) {
    warning("No contrasts available for volcano plots")
    return(invisible(file.path(figures_dir, "fig03_volcano.png")))
  }

  p <- patchwork::wrap_plots(panels, nrow = 1) +
    patchwork::plot_annotation(
      title    = "Differential Expression: Activated vs Unactivated T-cells",
      subtitle = paste0("FDR < ", fdr_cut, "  |  |LFC| > ", lfc_cut,
                        "  |  LFC shrinkage: ashr"),
      theme    = theme_tcell(12)
    )
  .save_gg(p, "fig03_volcano_panels", figures_dir, w = 15, h = 6)
}


# =============================================================================
# FIG 04 — MA Plot (interaction LRT)
# =============================================================================
plot_ma_interaction <- function(res_interaction, figures_dir, fdr_cut = 0.05) {
  message("\n>>> [Fig 04] MA plot ...")

  df <- res_interaction |>
    dplyr::filter(!is.na(baseMean), baseMean > 0) |>
    dplyr::mutate(
      lmean  = log10(baseMean),
      sig    = !is.na(padj) & padj < fdr_cut,
      col_g  = dplyr::case_when(
        sig & log2FoldChange > 0  ~ "Up",
        sig & log2FoldChange <= 0 ~ "Down",
        TRUE                       ~ "NS"
      )
    ) |>
    dplyr::arrange(!sig)

  n_up   <- sum(df$col_g == "Up")
  n_down <- sum(df$col_g == "Down")

  p <- ggplot(df, aes(lmean, log2FoldChange)) +
    geom_point(aes(colour = col_g), size = 0.65, alpha = 0.35, stroke = 0) +
    geom_hline(yintercept = 0, colour = "grey30", linewidth = 0.5) +
    geom_smooth(colour = "navy", method = "loess", se = FALSE,
                linewidth = 0.9, formula = y ~ x) +
    scale_colour_manual(
      values = c(Up = "#E64B35", Down = "#4DBBD5", NS = "grey72"),
      name   = paste0("FDR < ", fdr_cut)
    ) +
    annotate("text", x = Inf, y = Inf, hjust = 1.1, vjust = 1.5,
             label = paste0("\u2191 ", n_up, "    \u2193 ", n_down),
             size = 3.5, colour = "grey30") +
    labs(
      title    = "MA Plot — Interaction LRT",
      subtitle = "condition \u00d7 time spline interaction (full vs reduced model)",
      x        = expression(log[10]~"Mean Normalised Count"),
      y        = expression(log[2]~"Fold Change")
    ) +
    theme_tcell()
  .save_gg(p, "fig04_ma_interaction", figures_dir, w = 9, h = 6)
}


# =============================================================================
# FIG 05 — Gene Trajectory Plots
# =============================================================================
plot_gene_trajectories <- function(vst_matrix, final_clusters, sample_meta,
                                    figures_dir) {
  message("\n>>> [Fig 05] Gene trajectory plots ...")

  cl_hard <- final_clusters$cluster_hard
  cl_lab  <- final_clusters$trajectory_labels

  # ── Activated samples only, z-scored per gene ─────────────────────────────
  act_s  <- intersect(
    rownames(sample_meta)[sample_meta$condition == "activated"],
    colnames(vst_matrix)
  )
  meta_a <- sample_meta[act_s, ]
  vst_a  <- vst_matrix[intersect(names(cl_hard), rownames(vst_matrix)),
                        act_s, drop = FALSE]
  vst_z  <- t(scale(t(vst_a)))

  # Long format — all genes
  long_all <- vst_z |>
    as.data.frame() |>
    tibble::rownames_to_column("gene_id") |>
    tidyr::pivot_longer(-gene_id, names_to = "sample_id", values_to = "z") |>
    dplyr::left_join(
      meta_a |>
        dplyr::select(sample_id, time_hours, donor),
      by = "sample_id"
    ) |>
    dplyr::mutate(
      cluster_id    = cl_hard[gene_id],
      cluster_label = cl_lab[as.character(cluster_id)]
    ) |>
    dplyr::filter(!is.na(cluster_label))

  # Summary: mean ± SEM per cluster × time
  cl_summary <- long_all |>
    dplyr::group_by(cluster_label, time_hours) |>
    dplyr::summarise(
      mean_z = mean(z, na.rm = TRUE),
      sem_z  = sd(z, na.rm = TRUE) / sqrt(dplyr::n()),
      .groups = "drop"
    )

  # Gene-level mean (for faint lines) — random 35 per cluster
  set.seed(42)
  gene_subset <- long_all |>
    dplyr::distinct(gene_id, cluster_label) |>
    dplyr::group_by(cluster_label) |>
    dplyr::slice_sample(n = 35) |>
    dplyr::pull(gene_id)

  gene_means <- long_all |>
    dplyr::filter(gene_id %in% gene_subset) |>
    dplyr::group_by(gene_id, cluster_label, time_hours) |>
    dplyr::summarise(gm = mean(z, na.rm = TRUE), .groups = "drop")

  # Cluster size annotation
  cl_sizes <- long_all |>
    dplyr::distinct(gene_id, cluster_label) |>
    dplyr::count(cluster_label) |>
    dplyr::mutate(label = paste0(cluster_label, "\n(n=", n, ")"))
  size_map  <- setNames(cl_sizes$label, cl_sizes$cluster_label)
  gene_means$panel <- size_map[gene_means$cluster_label]
  cl_summary$panel <- size_map[cl_summary$cluster_label]

  p <- ggplot() +
    geom_line(data  = gene_means,
              aes(time_hours, gm, group = gene_id),
              colour = "grey78", linewidth = 0.2, alpha = 0.4) +
    geom_ribbon(data  = cl_summary,
                aes(time_hours,
                    ymin = mean_z - sem_z, ymax = mean_z + sem_z,
                    fill = cluster_label),
                alpha = 0.22) +
    geom_line(data  = cl_summary,
              aes(time_hours, mean_z, colour = cluster_label),
              linewidth = 1.2) +
    geom_point(data = cl_summary,
               aes(time_hours, mean_z, colour = cluster_label),
               size = 2.5) +
    scale_colour_manual(values = CLUSTER_COLORS, guide = "none") +
    scale_fill_manual(values   = CLUSTER_COLORS, guide = "none") +
    scale_x_continuous(breaks  = c(0,6,12,24,48,72),
                       labels  = c("0h","6h","12h","24h","48h","72h")) +
    facet_wrap(~ panel, nrow = 2, scales = "free_y") +
    labs(
      title    = "Temporal Gene Expression Clusters — Activated T-cells",
      subtitle = "Mfuzz soft clustering on z-scored VST counts",
      x        = "Time post-activation",
      y        = "Standardised Expression (z-score)",
      caption  = "Bold line = cluster mean \u00b1 SEM  |  Faint lines = random gene sample (n\u226435)"
    ) +
    theme_tcell()
  .save_gg(p, "fig05_gene_trajectories", figures_dir, w = 13, h = 8)
}


# =============================================================================
# FIG 06 — Cluster Heatmap
# =============================================================================
plot_cluster_heatmap <- function(vst_matrix, final_clusters, sample_meta,
                                  figures_dir, max_per_cl = 80) {
  message("\n>>> [Fig 06] Cluster heatmap ...")

  cl_hard <- final_clusters$cluster_hard
  cl_lab  <- final_clusters$trajectory_labels

  # Activated samples ordered by time then donor
  act_s  <- intersect(
    rownames(sample_meta)[sample_meta$condition == "activated"],
    colnames(vst_matrix)
  )
  meta_a <- sample_meta[act_s, ] |> dplyr::arrange(time_hours, donor)
  act_s  <- rownames(meta_a)

  # Select top-membership genes per cluster
  sel_genes <- if (!is.null(final_clusters$membership)) {
    unlist(lapply(sort(unique(cl_hard)), function(cl) {
      g <- names(cl_hard)[cl_hard == cl]
      m <- final_clusters$membership[g, cl]
      head(names(sort(m, decreasing = TRUE)), max_per_cl)
    }))
  } else {
    names(cl_hard)
  }
  sel_genes <- intersect(sel_genes, rownames(vst_matrix))

  # z-score and cap
  expr_z <- t(scale(t(vst_matrix[sel_genes, act_s, drop = FALSE])))
  expr_z <- pmin(pmax(expr_z, -3), 3)

  # Annotations
  gene_labels <- cl_lab[as.character(cl_hard[sel_genes])]
  ann_row <- data.frame(
    Cluster = factor(gene_labels, levels = names(CLUSTER_COLORS)),
    row.names = sel_genes
  )
  ann_col <- data.frame(
    Time_h    = factor(meta_a$time_hours),
    Donor     = meta_a$donor,
    row.names = rownames(meta_a)
  )
  n_donors <- nlevels(meta_a$donor)
  ann_colors <- list(
    Cluster = CLUSTER_COLORS,
    Time_h  = setNames(RColorBrewer::brewer.pal(6,"YlOrRd"), levels(ann_col$Time_h)),
    Donor   = setNames(RColorBrewer::brewer.pal(max(n_donors,3),"Set2")[seq_len(n_donors)],
                        levels(meta_a$donor))
  )

  row_ord <- order(match(ann_row$Cluster, names(CLUSTER_COLORS)))

  ph <- pheatmap::pheatmap(
    expr_z[row_ord, ],
    annotation_col    = ann_col,
    annotation_row    = ann_row,
    annotation_colors = ann_colors,
    color             = colorRampPalette(c("#2166AC","white","#D6604D"))(100),
    cluster_rows      = FALSE,
    cluster_cols      = FALSE,
    show_rownames     = FALSE,
    fontsize_col      = 7,
    main              = "Temporal Cluster Heatmap — Activated T-cells (z-score, capped \u00b13)",
    silent            = TRUE
  )
  .save_ph(ph, "fig06_cluster_heatmap", figures_dir, w = 12, h = 10)
}


# =============================================================================
# FIG 07 — Pathway Dot Plots
# =============================================================================
plot_pathway_dotplots <- function(pathway_go, figures_dir, n_top = 10) {
  message("\n>>> [Fig 07] Pathway dot plots ...")

  df <- purrr::map_dfr(names(pathway_go), function(cl_name) {
    enr <- pathway_go[[cl_name]]
    if (is.null(enr) || nrow(enr@result) == 0) return(NULL)
    enr@result |>
      dplyr::arrange(p.adjust) |>
      dplyr::slice_head(n = n_top) |>
      dplyr::mutate(
        cluster     = cl_name,
        nlq         = -log10(p.adjust),
        gene_ratio  = purrr::map_dbl(GeneRatio, function(gr) {
          p <- strsplit(gr, "/")[[1]]; as.numeric(p[1]) / as.numeric(p[2])
        }),
        Description = stringr::str_trunc(Description, 55)
      )
  })

  if (nrow(df) == 0) {
    warning("No GO results for dot plots")
    return(invisible(file.path(figures_dir, "fig07_pathway_dotplot.png")))
  }

  p <- ggplot(df, aes(cluster, reorder(Description, nlq),
                       size = gene_ratio, colour = nlq)) +
    geom_point() +
    scale_size_continuous(name = "Gene Ratio", range = c(2, 7)) +
    scale_colour_gradient(name = expression(-log[10]~"(adj. p)"),
                          low = "#FED976", high = "#BD0026") +
    labs(
      title    = "GO Biological Process Enrichment by Temporal Cluster",
      subtitle = paste0("Top ", n_top, " terms per cluster  |  BH-adjusted p < 0.05"),
      x        = NULL,
      y        = "GO Biological Process"
    ) +
    theme_tcell() +
    theme(axis.text.x = element_text(angle = 35, hjust = 1, size = 9),
          axis.text.y = element_text(size = 8))
  .save_gg(p, "fig07_pathway_dotplots", figures_dir, w = 14, h = 10)
}


# =============================================================================
# FIG 08 — Temporal Pathway Heatmap
# =============================================================================
plot_temporal_pathway_heatmap <- function(temporal_pathway_matrix, figures_dir) {
  message("\n>>> [Fig 08] Temporal pathway heatmap ...")

  if (is.null(temporal_pathway_matrix) ||
      nrow(temporal_pathway_matrix) == 0) {
    warning("Empty temporal pathway matrix — skipping Fig 08")
    return(invisible(file.path(figures_dir, "fig08_temporal_pathway_heatmap.png")))
  }

  ph <- pheatmap::pheatmap(
    temporal_pathway_matrix,
    color         = colorRampPalette(c("white","#FED976","#D62828"))(100),
    cluster_cols  = FALSE,
    cluster_rows  = TRUE,
    show_rownames = TRUE,
    show_colnames = TRUE,
    fontsize_row  = 8,
    fontsize_col  = 9,
    border_color  = NA,
    main          = expression("Temporal GO:BP Enrichment  –"*log[10]*"(BH-adj. p)"),
    silent        = TRUE
  )
  .save_ph(ph, "fig08_temporal_pathway_heatmap", figures_dir, w = 12, h = 9)
}


# =============================================================================
# FIG 09 — ARI Concordance Matrix
# =============================================================================
plot_cluster_concordance <- function(cluster_comparison, figures_dir) {
  message("\n>>> [Fig 09] ARI concordance plot ...")

  ari_mat <- cluster_comparison$ari_matrix
  df <- as.data.frame(ari_mat) |>
    tibble::rownames_to_column("Method_A") |>
    tidyr::pivot_longer(-Method_A, names_to = "Method_B", values_to = "ARI") |>
    dplyr::mutate(
      Method_A = factor(Method_A, levels = rownames(ari_mat)),
      Method_B = factor(Method_B, levels = colnames(ari_mat))
    )

  p <- ggplot(df, aes(Method_A, Method_B, fill = ARI)) +
    geom_tile(colour = "white", linewidth = 1.5) +
    geom_text(aes(label = round(ARI, 3)), size = 5.5, fontface = "bold") +
    scale_fill_gradient2(low = "#4DBBD5", mid = "white", high = "#E64B35",
                         midpoint = 0.5, limits = c(0, 1), name = "ARI") +
    coord_equal() +
    labs(
      title    = "Clustering Method Concordance",
      subtitle = "Adjusted Rand Index: Mfuzz · NMF · Spline k-means",
      x        = NULL, y = NULL,
      caption  = "ARI = 1.0: perfect agreement  |  ARI = 0: chance-level agreement"
    ) +
    theme_tcell() +
    theme(panel.grid = element_blank(),
          axis.text  = element_text(size = 11))
  .save_gg(p, "fig09_ari_concordance", figures_dir, w = 7, h = 6)
}
