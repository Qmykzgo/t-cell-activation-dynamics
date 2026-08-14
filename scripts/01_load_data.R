# =============================================================================
# scripts/01_load_data.R
# PURPOSE : Download GSE197067 from NCBI GEO; parse count matrix and sample
#           metadata; apply CPM-based low-count filter.
# INPUT   : GEO accession string, local data directory
# OUTPUT  : count matrix (integer genes×samples), sample_meta data.frame,
#           counts_filtered list($counts, $meta)
# DATASET : Rade et al. 2023, Genome Biology (PMID 37974140)
#           Pan T-cell Verification Set — GSE197067
#           44 samples | 6 time points (0–72 h) | 2 conditions | 4 donors
# =============================================================================

library(GEOquery)
library(edgeR)
library(dplyr)
library(stringr)
library(readr)
library(tibble)

# -----------------------------------------------------------------------------
# download_geo_data()
# Downloads the GSE series matrix (metadata) and supplementary count files.
# Returns: list($gse, $supp_files, $supp_dir)
# -----------------------------------------------------------------------------
download_geo_data <- function(gse_id, data_dir) {
  message("\n>>> [01] Downloading ", gse_id, " from NCBI GEO ...")

  geo_raw_dir <- file.path(data_dir, "geo_raw")
  dir.create(geo_raw_dir, showWarnings = FALSE, recursive = TRUE)

  # --- Series matrix (pData / metadata) -------------------------------------
  gse <- tryCatch(
    GEOquery::getGEO(
      GEO      = gse_id,
      destdir  = geo_raw_dir,
      GSEMatrix = TRUE,
      AnnotGPL  = FALSE,
      getGPL    = FALSE
    ),
    error = function(e) stop("GEO download failed: ", conditionMessage(e))
  )

  # --- Supplementary count files --------------------------------------------
  supp_dir <- file.path(geo_raw_dir, gse_id)
  if (!dir.exists(supp_dir) ||
      length(list.files(supp_dir, pattern = "\\.gz$|\\.txt$|\\.csv$")) == 0) {
    message(">>> Downloading supplementary files ...")
    GEOquery::getGEOSuppFiles(
      GEO           = gse_id,
      makeDirectory = TRUE,
      baseDir       = geo_raw_dir,
      fetch_files   = TRUE
    )
  } else {
    message(">>> Supplementary files already cached.")
  }

  supp_files <- list.files(supp_dir, full.names = TRUE)
  message(">>> Supplementary files found: ",
          paste(basename(supp_files), collapse = " | "))

  list(gse = gse, supp_files = supp_files, supp_dir = supp_dir)
}


# -----------------------------------------------------------------------------
# extract_count_matrix()
# Identifies and parses the raw count matrix from GEO supplementary files.
# Handles .txt.gz, .csv.gz, .tsv.gz with tab or comma delimiters.
# Returns: integer matrix, rownames = gene IDs, colnames = GSM IDs
# -----------------------------------------------------------------------------
extract_count_matrix <- function(raw_data) {
  supp_files <- raw_data$supp_files

  # Priority: files explicitly named "count" or "raw"; fallback to any gz/txt
  count_file <- supp_files[
    grepl("count|raw|expression|matrix", basename(supp_files), ignore.case = TRUE) &
    grepl("\\.gz$|\\.txt$|\\.csv$|\\.tsv$", basename(supp_files))
  ]

  if (length(count_file) == 0) {
    count_file <- supp_files[grepl("\\.gz$|\\.txt$", basename(supp_files))]
  }

  if (length(count_file) == 0 || is.na(count_file[1])) {
    stop(
      "Cannot locate count matrix in: ", raw_data$supp_dir, "\n",
      "Available files: ", paste(basename(raw_data$supp_files), collapse = ", ")
    )
  }

  target_file <- count_file[1]
  message(">>> Reading count matrix: ", basename(target_file))

  # Auto-detect delimiter from first two lines
  header_lines <- readLines(target_file, n = 2)
  sep_char     <- if (grepl("\t", header_lines[2])) "\t" else ","

  counts_df <- readr::read_delim(
    target_file,
    delim          = sep_char,
    show_col_types = FALSE,
    progress       = FALSE
  )

  gene_col     <- colnames(counts_df)[1]
  gene_ids     <- counts_df[[gene_col]]
  count_matrix <- as.matrix(counts_df[, -1, drop = FALSE])
  rownames(count_matrix) <- gene_ids

  # Sanity cleaning: remove NA gene IDs, all-zero rows
  keep         <- !is.na(rownames(count_matrix)) &
                  rownames(count_matrix) != "" &
                  rowSums(count_matrix, na.rm = TRUE) > 0
  count_matrix <- count_matrix[keep, , drop = FALSE]
  storage.mode(count_matrix) <- "integer"

  message(">>> Count matrix dimensions: ",
          nrow(count_matrix), " genes × ", ncol(count_matrix), " samples")
  count_matrix
}


# -----------------------------------------------------------------------------
# build_sample_metadata()
# Parses pData from the GEO series matrix to extract:
#   donor, condition (activated / control), time_hours, time_factor
# Uses regex to handle varied GEO characteristic string formats.
# Returns: data.frame, rownames = GEO sample IDs (GSMxxxxxx)
# -----------------------------------------------------------------------------
build_sample_metadata <- function(raw_data) {
  gse   <- raw_data$gse
  if (is.list(gse)) gse <- gse[[1]]          # handle multi-platform GSE
  pheno <- Biobase::pData(gse)

  message(">>> Building sample metadata from ", nrow(pheno), " samples")

  # Helper: search pData columns by pattern, return first match column name
  find_col <- function(df, patterns) {
    for (p in patterns) {
      hit <- grep(p, colnames(df), ignore.case = TRUE, value = TRUE)[1]
      if (!is.na(hit)) return(hit)
    }
    NA_character_
  }

  # Collapse all characteristics columns into one searchable string per sample
  char_cols  <- grep("characteristics", colnames(pheno), value = TRUE)
  char_blob  <- if (length(char_cols) > 0)
    apply(pheno[, char_cols, drop = FALSE], 1, paste, collapse = " | ")
  else
    pheno$title

  # ── Time extraction ────────────────────────────────────────────────────────
  # Priority 1: dedicated column (e.g. "time point:ch1")
  time_col <- find_col(pheno, "^time point")
  if (!is.na(time_col)) {
    time_hours <- suppressWarnings(as.numeric(
      stringr::str_extract(as.character(pheno[[time_col]]), "[0-9]+")))
  } else {
    # Priority 2: matches "0h", "6 h", "72 hours", "timepoint: 24h" etc.
    time_raw   <- str_extract(char_blob, "(?i)(?:time[^:]*:\\s*)?([0-9]+)\\s*h(?:ours?)?")
    time_hours <- as.numeric(str_extract(time_raw, "[0-9]+"))
  }

  # Fallback: parse from sample title
  if (all(is.na(time_hours))) {
    time_hours <- as.numeric(str_extract(pheno$title, "[0-9]+(?=\\s*h)"))
  }

  # ── Condition extraction ───────────────────────────────────────────────────
  # Activated = stimulation; control = unstimulated/not activated.
  # NB: "Not activated" ALSO contains "activ", so negation must be checked first.
  # Priority 1: dedicated column (e.g. "activation:ch1")
  cond_col <- find_col(pheno, "activation")
  if (!is.na(cond_col)) {
    cond_raw <- tolower(as.character(pheno[[cond_col]]))
    condition <- ifelse(
      grepl("not activ|without|unstim|untreat|negativ|control", cond_raw),
      "control",
      ifelse(grepl("activ|stimul|cd3|cd28|treat", cond_raw), "activated", NA_character_)
    )
  } else {
    # Priority 2: source_name / title ("after activation" vs "without any treatments")
    source_text <- paste(as.character(pheno$source_name), as.character(pheno$title))
    condition <- ifelse(
      grepl("without|not activ|unstim|untreat|control", source_text, ignore.case = TRUE),
      "control",
      ifelse(grepl("activ|stimul|cd3|cd28|treat", source_text, ignore.case = TRUE),
             "activated", NA_character_)
    )
  }

  # ── Donor extraction ───────────────────────────────────────────────────────
  # Priority 1: dedicated column (e.g. "biological replicate:ch1")
  rep_col <- find_col(pheno, "biological replicate|donor|subject|individual")
  if (!is.na(rep_col)) {
    donor <- paste0("D", suppressWarnings(as.numeric(
      stringr::str_extract(as.character(pheno[[rep_col]]), "[0-9]+"))))
  } else {
    # Priority 2: regex over characteristics blob
    donor_raw <- str_extract(char_blob, "(?i)(donor|subject|individual|patient|source)\\s*[0-9]+")
    donor <- if (all(is.na(donor_raw))) {
      # Priority 3: sample title suffix "_D1", "_D2", ...
      donor_from_title <- str_extract(pheno$title, "(?i)_?D[0-9]+$")
      if (all(!is.na(donor_from_title))) {
        toupper(str_replace(donor_from_title, "^_", ""))
      } else {
        paste0("D", as.integer(factor(pheno$source_name)))
      }
    } else {
      str_replace_all(donor_raw, "\\s+", "_")
    }
  }

  # ── Assemble metadata data.frame ──────────────────────────────────────────
  meta <- data.frame(
    sample_id    = rownames(pheno),
    geo_title    = as.character(pheno$title),
    time_hours   = time_hours,
    time_factor  = factor(time_hours, levels = sort(unique(time_hours))),
    condition    = factor(condition, levels = c("control", "activated")),
    donor        = factor(donor),
    stringsAsFactors = FALSE,
    row.names    = rownames(pheno)
  )

  # ── Validation report ─────────────────────────────────────────────────────
  expected_n_times    <- 6    # 0, 6, 12, 24, 48, 72
  expected_n_cond     <- 2
  expected_n_donors   <- 4

  message(">>> Time points detected (h): ",
          paste(sort(unique(meta$time_hours)), collapse = ", "))
  message(">>> Conditions: ",
          paste(levels(meta$condition), collapse = ", "))
  message(">>> Donors: ",
          paste(levels(meta$donor), collapse = ", "))

  if (any(is.na(meta$time_hours))) {
    warning(sum(is.na(meta$time_hours)),
            " sample(s) have NA time — inspect metadata manually.\n",
            "  Affected samples: ",
            paste(meta$sample_id[is.na(meta$time_hours)], collapse = ", "))
  }

  meta
}


# -----------------------------------------------------------------------------
# filter_low_counts()
# Aligns count matrix columns to sample metadata rows (handles GEO ID mismatches).
# Applies CPM > min_cpm in >= min_samples filter (standard bulk RNA-Seq practice).
# Returns: list($counts = filtered integer matrix, $meta = aligned metadata)
# -----------------------------------------------------------------------------
filter_low_counts <- function(count_matrix, sample_meta,
                               min_cpm = 1, min_samples = 3) {

  # --- Column-to-row alignment ----------------------------------------------
  shared <- intersect(colnames(count_matrix), rownames(sample_meta))

  if (length(shared) == 0) {
    # GEO sometimes uses "GSMxxxxxx" in counts but different IDs in pData;
    # fall back to positional alignment when dimensions match.
    if (ncol(count_matrix) == nrow(sample_meta)) {
      message(">>> Sample IDs differ — using positional alignment.")
      colnames(count_matrix) <- rownames(sample_meta)
      shared <- rownames(sample_meta)
    } else {
      stop(
        "Cannot align count matrix (", ncol(count_matrix), " cols) ",
        "to metadata (", nrow(sample_meta), " rows).\n",
        "Inspect colnames(count_matrix) vs rownames(sample_meta)."
      )
    }
  }

  count_matrix <- count_matrix[, shared, drop = FALSE]
  sample_meta  <- sample_meta[shared, , drop = FALSE]
  message(">>> Aligned ", length(shared), " samples")

  # --- CPM filter -----------------------------------------------------------
  cpm_mat <- edgeR::cpm(count_matrix)
  keep    <- rowSums(cpm_mat >= min_cpm) >= min_samples

  message(">>> Genes pre-filter:  ", nrow(count_matrix))
  message(">>> Genes post-filter: ", sum(keep),
          "  (CPM >= ", min_cpm, " in >= ", min_samples, " samples)")

  list(
    counts = count_matrix[keep, , drop = FALSE],
    meta   = sample_meta
  )
}


# -----------------------------------------------------------------------------
# save_deg_table()  /  save_cluster_table()  /  save_pathway_table()
# Convenience wrappers called by targets for results export.
# -----------------------------------------------------------------------------
save_deg_table <- function(res_interaction, deg_categories, results_dir) {
  dir.create(file.path(results_dir, "tables"), showWarnings = FALSE, recursive = TRUE)
  out <- file.path(results_dir, "tables", "deg_interaction_lrt.csv")

  res_interaction |>
    dplyr::mutate(
      deg_category = dplyr::case_when(
        gene_id %in% deg_categories$stimulation_specific ~ "stimulation_specific",
        gene_id %in% deg_categories$general_activation   ~ "general_activation",
        gene_id %in% deg_categories$temporal_only        ~ "temporal_only",
        TRUE ~ "not_significant"
      )
    ) |>
    readr::write_csv(out)

  message(">>> Saved: ", out)
  out
}

save_cluster_table <- function(final_clusters, results_dir) {
  dir.create(file.path(results_dir, "tables"), showWarnings = FALSE, recursive = TRUE)
  out <- file.path(results_dir, "tables", "cluster_membership.csv")

  data.frame(
    gene_id       = names(final_clusters$cluster_hard),
    cluster_id    = final_clusters$cluster_hard,
    cluster_label = final_clusters$trajectory_labels[
      as.character(final_clusters$cluster_hard)],
    stringsAsFactors = FALSE
  ) |>
    dplyr::arrange(cluster_id) |>
    readr::write_csv(out)

  message(">>> Saved: ", out)
  out
}

save_pathway_table <- function(pathway_go, pathway_kegg, results_dir) {
  dir.create(file.path(results_dir, "tables"), showWarnings = FALSE, recursive = TRUE)
  out <- file.path(results_dir, "tables", "pathway_enrichment.csv")

  extract_results <- function(enrich_list, db_name) {
    purrr::map_dfr(names(enrich_list), function(cl_name) {
      enr <- enrich_list[[cl_name]]
      if (is.null(enr) || nrow(enr@result) == 0) return(NULL)
      enr@result |>
        dplyr::mutate(cluster = cl_name, database = db_name)
    })
  }

  dplyr::bind_rows(
    extract_results(pathway_go,   "GO_BP"),
    extract_results(pathway_kegg, "KEGG")
  ) |>
    dplyr::filter(p.adjust < 0.05) |>
    dplyr::arrange(cluster, p.adjust) |>
    readr::write_csv(out)

  message(">>> Saved: ", out)
  out
}
