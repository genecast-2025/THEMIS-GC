# =============================================================================
# Script:      FigureS2.R
# Description: Genome-wide copy number alteration of chromosome arms.
#              Heatmap of z-score normalised mean coverage per chromosome arm
#              in Training and Validation datasets (Xijing cohort), ordered by
#              disease stage and overall copy-number alteration burden.
#
# --- Required Input Data ------------------------------------------------------
# File 1:  data/Table_S3.tsv                    (real data — Supplementary Table S3)
#          Columns: ID, Diagnosis, Laurens' subtypes, Tumor stage, Cohort,
#                   Gender, Age, MFR, FSI, FEM, CAFF, THEMIS-GC
#
# File 2:  data/CAFF_arm_coverage.demo.tsv       (demo placeholder)
#          Pre-computed mean DELFI coverage per chromosome arm per sample.
#          Rows = chromosome arms (e.g. 1p, 1q, … 22q).
#          Columns = arm (row label) + one column per sample ID.
#          Values = mean coverage_corrected_delfi across bins in that arm.
#          Normalisation (column-wise then row-wise z-score) is performed
#          inside this script; the demo file should contain raw arm means.
#
# --- Output -------------------------------------------------------------------
# Figure:   figures/Figure_S2_CAFF_heatmap.pdf
# =============================================================================

library(pheatmap)
library(RColorBrewer)
library(here)
library(systemfonts)

font_path <- here("fonts", "arial.ttf")
if (file.exists(font_path)) register_font(name = "Arial", plain = font_path)

# =============================================================================
# Load & harmonise sample metadata (Table S3)
# =============================================================================
sample.df <- read.table(here("data", "Table_S3.tsv"), header = TRUE, sep = "\t",
                          check.names = FALSE, stringsAsFactors = FALSE)
sample.df <- subset(sample.df, Cohort %in% c("Xijing-Training", "Xijing-Validation"))

sample.df$Sample <- sample.df$ID
sample.df$Stage  <- ifelse(sample.df$Diagnosis %in% c("HEALTHY", "BENIGN"),
                             sample.df$Diagnosis, sample.df$`Tumor stage`)
sample.df$Group  <- ifelse(sample.df$Diagnosis %in% c("HEALTHY", "BENIGN"),
                             sample.df$Diagnosis, sample.df$`Laurens' subtypes`)
sample.df$Cohort <- gsub("Xijing-", "", sample.df$Cohort)

anno.df <- data.frame(
  Stage = factor(sample.df$Stage,
                  levels = c("HEALTHY", "BENIGN", "I", "II", "III", "IV", "Unknown")),
  Group = factor(sample.df$Group,
                  levels = c("HEALTHY", "BENIGN", "DT", "IT", "MT", "Unknown")),
  row.names = sample.df$Sample,
  stringsAsFactors = FALSE
)

baseline.v <- sample.df$Sample[sample.df$Stage %in% c("HEALTHY", "BENIGN") &
                                   sample.df$Cohort == "Training"]

anno_colors <- list(
  Group = c(HEALTHY = brewer.pal(9, "Blues")[8],
             BENIGN  = brewer.pal(9, "Blues")[7],
             DT      = brewer.pal(6, "Dark2")[3],
             IT      = brewer.pal(6, "Dark2")[4],
             MT      = brewer.pal(6, "Dark2")[5],
             Unknown = brewer.pal(6, "Dark2")[6]),
  Stage = c(HEALTHY = brewer.pal(9, "Blues")[8],
             BENIGN  = brewer.pal(9, "Blues")[7],
             I       = brewer.pal(8, "Reds")[5],
             II      = brewer.pal(8, "Reds")[6],
             III     = brewer.pal(8, "Reds")[7],
             IV      = brewer.pal(8, "Reds")[8],
             Unknown = brewer.pal(8, "Reds")[4])
)

# =============================================================================
# Load pre-computed arm coverage matrix (demo)
# =============================================================================
cov.df    <- read.table(here("data", "CAFF_arm_coverage.demo.tsv"),
                          header = TRUE, check.names = FALSE, sep = "\t")
chr_cov.m <- as.matrix(cov.df[, -1])
rownames(chr_cov.m) <- cov.df[[1]]
chr_cov.m <- chr_cov.m[, sample.df$Sample]

# =============================================================================
# Normalisation
# (1) Column-wise: express each arm as fraction of total arm coverage per sample
# (2) Row-wise:    z-score relative to Training healthy/benign baseline
# =============================================================================
chr_cov.m <- apply(chr_cov.m, 2, function(x) x / sum(x))

chr_cov.m <- t(apply(chr_cov.m, 1, function(x) {
  baseline_vals <- x[baseline.v]
  (x - mean(baseline_vals)) / sd(baseline_vals)
}))

# =============================================================================
# Column ordering: Stage, then descending CNA burden (sum |z|)
# =============================================================================
z.v <- apply(abs(chr_cov.m), 2, sum) * (-1)
temp_anno.df <- data.frame(anno.df,
                             z = z.v[match(rownames(anno.df), names(z.v))])
temp_anno.df <- temp_anno.df[order(temp_anno.df$Stage, temp_anno.df$z),
                                c("Group", "Stage"), drop = FALSE]

# =============================================================================
# Figure S2 - CAFF chromosome-arm heatmap
# =============================================================================
my_colors <- colorRampPalette(c("#5C88DAE5", "#FFFFFFFF", "#CC0C00E5"))(100)

pdf(here("figures", "Figure_S2_CAFF_heatmap.pdf"))
pheatmap(
  chr_cov.m[, match(rownames(temp_anno.df), colnames(chr_cov.m))],
  cluster_rows    = FALSE,
  cluster_cols    = FALSE,
  annotation_col  = temp_anno.df,
  show_rownames   = TRUE,
  show_colnames   = FALSE,
  breaks          = seq(-15, 15, length.out = 101),
  cellwidth       = 0.52,
  cellheight      = 10.4,
  annotation_colors = anno_colors,
  color           = my_colors
)
dev.off()
