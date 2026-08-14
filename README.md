# T-cell Activation Dynamics: Bulk RNA-Seq Time-Course Analysis

> **Multi-method temporal analysis of the human T-cell transcriptome during
> the first 72 hours of anti-CD3/CD28 activation** — DESeq2 LRT with natural
> splines, three clustering methods (Mfuzz / NMF / spline k-means), and
> clusterProfiler pathway enrichment, orchestrated end-to-end with `targets`.

---

## Biological Question

Which genes show **stimulation-specific** temporal dynamics during T-cell
activation, as opposed to passive temporal drift present in both activated
and unactivated cells? And how do those genes organise into biologically
coherent trajectory clusters?

Understanding the temporal architecture of T-cell activation is directly
relevant to:

- **CAR-T manufacturing** — activation protocols define the epigenetic memory
  of the final product
- **Checkpoint immunotherapy** — exhaustion begins during early activation
- **Vaccine design** — timing of adjuvant delivery relative to antigen shapes
  the quality of the T-cell response

---

## Dataset

| Field | Value |
|---|---|
| GEO Accession | [GSE197067](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE197067) |
| Publication | Rade *et al.* (2023) *Genome Biology*, [PMID 37974140](https://pubmed.ncbi.nlm.nih.gov/37974140/) |
| Samples | 44 (4 donors × 6 time points × 2 conditions) |
| Time Points | 0h, 6h, 12h, 24h, 48h, 72h |
| Conditions | Anti-CD3/CD28 bead activated · Unactivated control |
| Organism | *Homo sapiens* (peripheral blood T-cells) |
| Data type | Bulk RNA-Seq raw counts (from GEO supplementary files) |

---

## Pipeline Architecture

```
┌─────────────────────────────────────────────────────────────────┐
│              T-cell Activation Dynamics Pipeline                 │
│                   targets::tar_make()                            │
└─────────────────────────────────────────────────────────────────┘

GEO: GSE197067
      │  GEOquery::getGEO() + getGEOSuppFiles()
      ▼
┌─────────────────┐
│  01_load_data   │  Count matrix │ Sample metadata │ CPM ≥ 1 filter
└────────┬────────┘
         │  counts_filtered  ·  sample_meta
         ▼
┌──────────────────────────────────────────────────────────────────┐
│                    02_deseq2_lrt                                  │
│                                                                   │
│  ┌─────────────────────────────────────────────────────────────┐ │
│  │  MODEL 1 — Interaction LRT (PRIMARY)                         │ │
│  │  Full:    ~ donor + condition + ns(time,df=3) +             │ │
│  │                               condition:ns(time,df=3)       │ │
│  │  Reduced: ~ donor + condition + ns(time,df=3)               │ │
│  │  Tests: does trajectory SHAPE differ between conditions?    │ │
│  └─────────────────────────────────────────────────────────────┘ │
│                                                                   │
│  ┌─────────────────────────────────────────────────────────────┐ │
│  │  MODEL 2 — Activated-only LRT (SECONDARY)                   │ │
│  │  ~ donor + ns(time,df=3)  vs  ~ donor                      │ │
│  └─────────────────────────────────────────────────────────────┘ │
│                                                                   │
│  ┌─────────────────────────────────────────────────────────────┐ │
│  │  MODEL 3 — Pairwise contrasts (TERTIARY)                    │ │
│  │  Activated vs Control at each time point                    │ │
│  │  LFC shrinkage: lfcShrink(type="ashr")                      │ │
│  └─────────────────────────────────────────────────────────────┘ │
│                                                                   │
│  DEG CATEGORIES:                                                  │
│    [Model 1 only] → STIMULATION-SPECIFIC                         │
│    [Both models]  → GENERAL ACTIVATION                           │
│    [Model 2 only] → TEMPORAL ONLY (likely batch/culture effect)  │
└──────────────────────────────────────────┬───────────────────────┘
                                            │
                              top DEGs (FDR<0.05, |LFC|≥1)
                                            │
                ┌───────────────────────────┼───────────────────────┐
                │                           │                       │
                ▼                           ▼                       ▼
        ┌──────────────┐          ┌──────────────┐        ┌──────────────────┐
        │    Mfuzz     │          │     NMF      │        │  Spline k-means  │
        │ Fuzzy c-means│          │ Brunet algo  │        │ k-means on cubic │
        │ z-score VST  │          │ rank = 5     │        │ spline coefs     │
        └──────┬───────┘          └──────┬───────┘        └────────┬─────────┘
               └──────────────────┬──────┘                         │
                                  │     Adjusted Rand Index        │
                                  └────────────────────────────────┘
                                             │
                                    Consensus clusters
                                    (Mfuzz as primary)
                                             │
                    ┌──────────────────────────────────────────┐
                    │           Trajectory Labels               │
                    │  Early Innate Burst      (peak ≤ 6h)    │
                    │  Sustained Inflammatory  (monotonic ↑)  │
                    │  Delayed Proliferative   (peak 48–72h)  │
                    │  Transient Metabolic     (peak 12–24h)  │
                    │  Repressed Quiescence    (sustained ↓)  │
                    └──────────────────────────────────────────┘
                                             │
                              ┌──────────────┴──────────────┐
                              ▼                              ▼
                    ┌──────────────────┐        ┌──────────────────────┐
                    │  GO:BP per       │        │  Temporal pathway    │
                    │  cluster         │        │  enrichment matrix   │
                    │  clusterProfiler │        │  (per time point)    │
                    └──────────────────┘        └──────────────────────┘
                                             │
                              ┌──────────────┴──────────────┐
                              ▼                              ▼
                     figures/ (PNG+PDF)           results/tables/ (CSV)
```

---

## Key Results

> *(Placeholder — fill in after running the pipeline)*

- **~N stimulation-specific DEGs** identified by the interaction LRT
  (FDR < 0.05), of which ~N% are not detectable by simple pairwise comparison
  at any single time point alone.

- **Five trajectory clusters** with high concordance across methods
  (ARI_Mfuzz_NMF ≈ X.XX) — validating the biological robustness of the
  temporal patterns.

- **Pathway timing:** TCR signalling peaks at 6h → Warburg metabolism at
  12–24h → cell cycle entry at 48–72h, consistent with known T-cell
  activation biology and providing a reference framework for CAR-T product
  characterisation.

### Figure Gallery

| Figure | Description |
|---|---|
| ![PCA](figures/fig01_pca.png) | PCA: condition × time structure |
| ![Trajectories](figures/fig05_gene_trajectories.png) | Cluster mean trajectories |
| ![Heatmap](figures/fig06_cluster_heatmap.png) | Cluster heatmap |
| ![Pathways](figures/fig08_temporal_pathway_heatmap.png) | Temporal pathway activation |

---

## Repository Structure

```
tcell-activation-dynamics/
│
├── _targets.R                  ← Pipeline DAG (one command to run all)
├── analysis_report.qmd         ← Quarto report (rendered → GitHub Pages)
├── renv.lock                   ← Exact R package versions (renv)
├── environment.yml             ← Conda environment (Python/system deps)
├── Makefile                    ← Convenience targets (make report, make clean)
│
├── scripts/
│   ├── 01_load_data.R          ← GEO download · count parsing · CPM filter
│   ├── 02_deseq2_lrt.R         ← DESeq2 LRT · spline models · contrasts
│   ├── 03_clustering.R         ← Mfuzz · NMF · spline k-means · ARI
│   ├── 04_pathways.R           ← clusterProfiler GO/KEGG · temporal matrix
│   └── 05_visualization.R      ← All ggplot2 + pheatmap figures (300 DPI)
│
├── data/                       ← Downloaded from GEO (git-ignored)
├── results/
│   ├── rds/                    ← targets object store (git-ignored)
│   └── tables/                 ← CSV exports (DEGs, clusters, pathways)
├── figures/                    ← PNG + PDF outputs
├── notebooks/                  ← Exploratory Rmd scratchpads
└── docs/
    └── portfolio_summary.md    ← One-page summary for personal website
```

---

## How to Reproduce

### Prerequisites

- R ≥ 4.3.0
- Internet connection (for GEO download and KEGG API)
- RAM ≥ 8 GB (16 GB recommended for NMF)
- Storage ≥ 2 GB

### Step-by-step

```bash
# 1. Clone the repository
git clone https://github.com/YOUR_USERNAME/tcell-activation-dynamics
cd tcell-activation-dynamics

# 2. Restore exact R package environment
Rscript -e "install.packages('renv'); renv::restore()"

# 3. Run the entire pipeline (downloads data, runs all analyses, saves all figures)
Rscript -e "targets::tar_make()"

# 4. Render the Quarto report
quarto render analysis_report.qmd

# OR: use Make
make all
```

### Check pipeline status without running

```r
library(targets)
tar_visnetwork()   # interactive DAG in browser
tar_outdated()     # list targets that need re-running
tar_manifest()     # full target manifest
```

### Run a single target for development

```r
targets::tar_make(names = "fig_trajectories")
```

---

## Methods Summary

| Step | Tool | Version | Notes |
|---|---|---|---|
| GEO download | `GEOquery` | 2.70+ | GSEMatrix + supplementary files |
| Low-count filter | `edgeR::cpm` | 4.0+ | CPM ≥ 1 in ≥ 3 samples |
| Normalisation | `DESeq2::vst` | 1.42+ | `blind=FALSE` |
| Interaction LRT | `DESeq2` + `glmGamPoi` | 1.42+ | model.matrix with `ns(time,df=3)` |
| LFC shrinkage | `ashr` | 2.2+ | Works with any contrast |
| Soft clustering | `Mfuzz` | 2.66+ | c=5, m=1.25 |
| NMF | `NMF` | 0.28+ | brunet, nrun=50 |
| Spline k-means | `stats::kmeans` | base | Coefficients from `ns(time,df=3)` |
| Concordance | `mclust::adjustedRandIndex` | 6.0+ | ARI between 3 methods |
| Pathway enrichment | `clusterProfiler` | 4.10+ | GO:BP + KEGG per cluster |
| Annotation | `org.Hs.eg.db` | 3.18+ | ENSEMBL/SYMBOL → Entrez |
| Pipeline | `targets` | 1.6+ | `tar_make()` orchestration |
| Environment | `renv` | 1.0+ | Lock file reproducibility |
| Report | `Quarto` | 1.4+ | HTML + PDF output |

**Statistical justification:** The spline-basis LRT approach is preferred
over pairwise static comparisons for time-series data (Schulz *et al.* 2018,
[PMC5826275](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5826275/)) because
it (1) uses all time points jointly, increasing power; (2) formally tests
trajectory *shape*, not just magnitude; and (3) reduces multiple testing
burden by testing one model per gene rather than one contrast per time point
per gene.

---

## System Requirements

| Resource | Minimum | Recommended |
|---|---|---|
| RAM | 8 GB | 16 GB |
| Storage | 700 MB | 2 GB |
| CPU cores | 1 | 4+ (NMF parallelises) |
| R version | 4.2.0 | 4.3.0+ |
| OS | Linux / macOS | Linux / macOS |

> **Windows users:** `glmGamPoi` and some Bioconductor packages may require
> Rtools. WSL2 is the recommended path on Windows.

---

## Limitations

1. **n = 4 donors** — insufficient power for donor × condition interactions.
2. **Bulk RNA-Seq** — CD4⁺/CD8⁺ subset heterogeneity is averaged out.
3. **In vitro beads** — supraphysiological stimulus; in vivo dynamics differ.
4. **No proteomics** — translational lag not captured.
5. **Cross-sectional donors** — each donor measured at all time points, but
   PBMC batches may introduce subtle between-donor technical variation.

---

## What I Would Do With More Resources

| Scale-up | Scientific value |
|---|---|
| scRNA-seq (10x Chromium) of the same activation | Resolve CD4⁺ vs CD8⁺ vs Treg trajectories; identify rare early-responding subsets |
| CITE-seq (RNA + surface protein) | Quantify translational lag; validate mRNA→protein synchrony |
| CRISPR perturbation screen (pooled, activated T-cells) | Functionally validate top *Sustained Inflammatory* hub genes |
| CAR-T infusion product scRNA-seq | Test whether GSE197067 activation clusters predict persistence and efficacy in vivo |
| Phospho-proteomics (IMAC-MS) at 5–30 min post-activation | Capture the immediate signalling layer upstream of the transcriptional response |

---

## Citation

If you use this analysis or adapt any code, please cite the original dataset:

> Rade M, Böttcher M, Weinberger T, *et al.* (2023). Pan-T cell
> activation transcriptomics for T cell subset verification in bulk
> RNA-Seq. *Genome Biology*, 24, 256.
> [https://doi.org/10.1186/s13059-023-03104-7](https://doi.org/10.1186/s13059-023-03104-7)

And the statistical framework:

> Schulz MH, Devanny WE, Gitter A, Zhong S, Ernst J, Bar-Joseph Z (2012).
> DREM 2.0: Improved reconstruction of dynamic regulatory networks from
> time-series expression data. *Frontiers in Genetics*.
> [PMC5826275](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5826275/)

---

## License

MIT License — see [LICENSE](LICENSE) for details.

Code is free to use and adapt with attribution.
Data belongs to the original authors (Rade *et al.* 2023, CC BY 4.0).
