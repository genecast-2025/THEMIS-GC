# THEMIS-GC · Source Code Appendix

Reproducible R analysis scripts for all figures in:

> **Multi-modal cell-free DNA profiling enables early detection and prognosis of gastric cancer**

THEMIS-GC is an ensemble classifier integrating four cfDNA modalities — methylation fragment ratio (MFR), fragment size index (FSI), fragment end motif (FEM), and chromosomal arm fragment frequency (CAFF) — for gastric cancer detection and prognosis.

---

## Prerequisites

### Software

| Tool | Minimum version |
|------|-----------------|
| R | 4.1 |
| Rscript (CLI) | 4.1 |

### R packages

```r
# CRAN
install.packages(c(
  "tidyverse",     # ggplot2, dplyr, tidyr, readr, purrr, stringr
  "ggpubr", "ggsci", "rlang", "svglite",
  "RColorBrewer", "gridExtra", "ggsignif", "ggforce",
  "pROC",          # ROC curves and DeLong test
  "rstatix",       # pairwise Wilcoxon with BH correction
  "rmda",          # decision curve analysis
  "survival", "survminer",   # Kaplan-Meier and Cox regression
  "rms",           # cph() for Chi-Square contribution
  "timeROC",       # time-dependent ROC
  "pheatmap",      # heatmaps
  "here",          # relative path management
  "systemfonts"    # optional Arial font registration
))

# Bioconductor
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(c("ComplexHeatmap", "circlize"))
```

> **Optional font:** Place `arial.ttf` in a `fonts/` sub-directory alongside `data/` to render plots with Arial. The scripts fall back to the system default font if the file is absent.

---

## Data setup

Each script reads from `data/`, writes figures to `figures/`, and writes result tables to `results/`. Create these directories once from your working directory:

```bash
mkdir -p data figures results
```

Export the following sheets from `source_data/Table S1-S9.xlsx` as **tab-separated** files:

| File to create | Excel sheet | Scripts that read it |
|----------------|-------------|----------------------|
| `data/Table_S3.tsv` | Table S3 | Figure2, Figure3, FigureS1, FigureS2, FigureS3, FigureS4, FigureS5 |
| `data/Table_S5.tsv` | Table S5 | Figure4, FigureS6, FigureS7, FigureS8 |
| `data/Table_S8.tsv` | Table S8 | Figure5, FigureS9 |

Convert the two model-score workbooks in `source_data/paitents_with_protein_data/` to TSV:

| File to create | Source workbook |
|----------------|----------------|
| `data/model1_cfDNA_features.tsv` | `model1_cfDNA_features.new.data.df.xls` |
| `data/model2_cfDNA_protein_features.tsv` | `model2_combiend_cfDNA_features_with protein biomarkers_AFP+CA125+CA199+CEA.new.data.df.xls` |

> **FigureS1 and FigureS2** require large per-sample feature matrices (MFR, FSI, FEM, CAFF) that are not redistributed here. Placeholder filenames (`data/*.demo.tsv`) are documented at the top of each script and must be replaced with the real matrices before those figures can be generated.

---

## Quick Start

Reproduce **Figure 2** (THEMIS-GC diagnostic performance, Training and Validation cohorts) after placing `data/Table_S3.tsv` in your working directory:

```bash
mkdir -p data figures results
Rscript R_scripts/Figure2.R
```

Expected outputs in `figures/`:

```
Figure_2A.pdf
Figure_2_Training_ROC_multi_set.pdf
Figure_2_Validation_ROC_multi_set.pdf
Figure2D_THEMIS-GC_Stage_Cancer_spf0.95_senci.hist.svg
Training.dca_final_matched.pdf
Validation.dca_final_matched.pdf
```

---

## Repository structure

```
.
├── R_scripts/                          # Compiled, publication-ready scripts
│   ├── Figure2.R                       # THEMIS-GC performance — Training & Validation
│   ├── Figure3.R                       # THEMIS-GC performance — external cohorts
│   ├── Figure4.R                       # cfDNA vs cfDNA+protein models (ROC + landscape)
│   ├── Figure5.R                       # Prognosis 
│   ├── FigureS1.R                      # cfDNA feature distributions & PCA  [demo data]
│   ├── FigureS2.R                      # Chromosome-arm CNA heatmap          [demo data]
│   ├── FigureS3.R                      # Individual feature boxplots — Training/Validation
│   ├── FigureS4.R                      # THEMIS-GC sensitivity by clinical subgroup
│   ├── FigureS5.R                      # Individual feature boxplots — external cohorts
│   ├── FigureS6.R                      # Score/biomarker distributions — protein sub-cohort
│   ├── FigureS7.R                      # ROC: THEMIS-GC vs protein biomarkers
│   ├── FigureS8.R                      # cfDNA vs cfDNA+protein — Testing Dataset 1
│   └── FigureS9.R                      # Prognosis 
│
├── source_data/
│   ├── Table S1-S9.xlsx                # All supplementary tables
│   └── paitents_with_protein_data/
│       ├── model1_cfDNA_features.new.data.df.xls
│       └── model2_combiend_cfDNA_features_with protein biomarkers_AFP+CA125+CA199+CEA.new.data.df.xls


```

### Script dependency graph

Two supplemental scripts inherit shared functions from their paired main-figure script via `source()`:

```
Figure4.R  ──source()──▶  FigureS8.R
Figure5.R  ──source()──▶  FigureS9.R
```

Running `FigureS8.R` or `FigureS9.R` will also regenerate the Figure 4 and Figure 5 outputs respectively.

---

## Input data summary

| Supplementary table | Rows | Cohorts | Used for |
|---------------------|------|---------|----------|
| Table S3 | 1,281 | Training, Validation, Testing 1, MONITOR | Figures 2, 3, S1–S5 |
| Table S5 | 367 | Xijing-Training, Xijing-Validation, Testing 1 | Figures 4, S6, S7, S8 |
| Table S8 | 125 | Xijing (surgery patients) | Figures 5, S9 |

---

## Session info

Developed and tested under R 4.x on macOS. Run `sessionInfo()` after loading packages to record your exact environment for reproducibility reporting.

---

## Citation

> To be completed upon publication.
