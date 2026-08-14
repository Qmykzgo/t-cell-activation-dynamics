# =============================================================================
# Makefile — T-cell Activation Dynamics Pipeline
# Usage:
#   make all        Run pipeline + render report
#   make pipeline   Run targets pipeline only
#   make report     Render Quarto report (HTML + PDF)
#   make clean      Remove targets store and rendered outputs
#   make dag        Print pipeline DAG to browser
#   make status     Show which targets are outdated
#   make env        Create conda environment
# =============================================================================

.PHONY: all pipeline report clean dag status env install check

# ── Main entry point ──────────────────────────────────────────────────────────
all: pipeline report

# ── 1. Run the full targets pipeline ─────────────────────────────────────────
pipeline:
	@echo "==> Running targets pipeline ..."
	Rscript -e "targets::tar_make()"

# ── 2. Render Quarto report ───────────────────────────────────────────────────
report: analysis_report.qmd
	@echo "==> Rendering Quarto report (HTML + PDF) ..."
	quarto render analysis_report.qmd --to html
	quarto render analysis_report.qmd --to pdf
	@echo "==> Report written to: analysis_report.html / analysis_report.pdf"

# ── 3. Interactive pipeline DAG ───────────────────────────────────────────────
dag:
	Rscript -e "targets::tar_visnetwork()"

# ── 4. Check which targets need re-running ────────────────────────────────────
status:
	Rscript -e "targets::tar_outdated()"

# ── 5. Install R packages via renv ───────────────────────────────────────────
install:
	@echo "==> Restoring R package environment with renv ..."
	Rscript -e "if (!requireNamespace('renv', quietly=TRUE)) install.packages('renv'); renv::restore()"

# ── 6. Create conda environment ───────────────────────────────────────────────
env:
	conda env create -f environment.yml
	@echo "==> Activate with: conda activate tcell-activation"

# ── 7. Check R package versions ───────────────────────────────────────────────
check:
	Rscript -e "
	pkgs <- c('targets','DESeq2','GEOquery','Mfuzz','NMF','mclust',
	          'clusterProfiler','ggplot2','patchwork','glmGamPoi','ashr')
	installed <- rownames(installed.packages())
	missing   <- setdiff(pkgs, installed)
	if (length(missing)==0) cat('All packages OK\n') else
	  cat('Missing packages:', paste(missing, collapse=', '), '\n')
	"

# ── 8. Clean generated outputs ────────────────────────────────────────────────
clean:
	@echo "==> Removing targets store, figures, results ..."
	Rscript -e "targets::tar_destroy(ask=FALSE)" 2>/dev/null || true
	rm -rf figures/*.png figures/*.pdf
	rm -rf results/tables/*.csv
	rm -f analysis_report.html analysis_report.pdf
	@echo "==> Clean complete. Run 'make all' to regenerate."

# ── 9. Deep clean (also removes downloaded data) ──────────────────────────────
distclean: clean
	@echo "==> Also removing downloaded GEO data ..."
	rm -rf data/geo_raw/
	@echo "==> Full clean complete."
