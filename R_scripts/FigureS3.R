# =============================================================================
# Script:      FigureS3.R
# Description: Distribution of individual cfDNA features in the Training and
#              Validation datasets from the Xijing cohort.
#              Boxplots comparing score distributions across healthy/benign
#              individuals and GC patients at different disease stages:
#              (A) MFR  (B) FSI  (C) FEM  (D) CAFF
#              The centre line = median; box limits = upper/lower quartiles.
#
# --- Required Input Data ------------------------------------------------------
# File:        data/Table_S3.tsv
# Source:      Manuscript Supplementary Table S3
#              ("Cancer prediction scores by individual modalities and THEMIS-GC")
# Format:      Tab-separated, columns:
#              ID, Diagnosis, Laurens' subtypes, Tumor stage, Cohort, Gender,
#              Age, MFR, FSI, FEM, CAFF, THEMIS-GC
#
# --- Output -------------------------------------------------------------------
# Figures:     figures/Figure_S3A.pdf   (MFR score)
#              figures/Figure_S3B.pdf   (FSI score)
#              figures/Figure_S3C.pdf   (FEM score)
#              figures/Figure_S3D.pdf   (CAFF score)
#
# =============================================================================

library(ggplot2)
library(ggsci)
library(rlang)
library(here)
library(systemfonts)

font_path <- here("fonts", "arial.ttf")
if (file.exists(font_path)) register_font(name = "Arial", plain = font_path)

# =============================================================================
# Load & harmonise Table S3
# =============================================================================
df <- read.table(here("data", "Table_S3.tsv"), header = TRUE, sep = "\t",
                  check.names = FALSE, stringsAsFactors = FALSE)

# Cohort: Table S3 uses "Xijing-Training" / "Xijing-Validation"
df$Cohort <- gsub("Xijing-", "", df$Cohort)
df$Cohort <- gsub("Testing Dataset 1", "testing1", df$Cohort)
df$Cohort <- gsub("Testing Dataset 2", "testing2", df$Cohort)

# Stage (x-axis): Table S3 splits disease status into Diagnosis
# (HEALTHY/BENIGN/Cancer) and `Tumor stage` (I-IV/Unknown, NA for non-cancer).
# Recombine into a single Stage factor for plotting.
df$Stage <- ifelse(df$Diagnosis %in% c("HEALTHY", "BENIGN"),
                    df$Diagnosis, df$`Tumor stage`)

# =============================================================================
# Figures S3A-D - Individual cfDNA feature scores by stage, Training vs Validation
# =============================================================================
# Panel labels and corresponding Table S3 column names
panels <- list(
  A = "MFR",
  B = "FSI",
  C = "FEM",
  D = "CAFF"
)

plot.df <- subset(df, Cohort %in% c("Training", "Validation"))
plot.df$Cohort <- factor(plot.df$Cohort, levels = c("Training", "Validation"))
plot.df$Stage  <- factor(plot.df$Stage,
                          levels = c("HEALTHY", "BENIGN", "I", "II", "III", "IV", "Unknown"))

colr <- pal_startrek("uniform", alpha = 0.9)(2)
names(colr) <- c("Training", "Validation")

for (panel in names(panels)) {
  model <- panels[[panel]]

  p <- ggplot(plot.df, aes(x = Stage, y = !!sym(model), color = Cohort)) +
    geom_boxplot() +
    labs(title = paste0(model, " score")) +
    coord_cartesian(ylim = c(0, 1)) +
    theme_bw() +
    theme(legend.position = "bottom",
          plot.title = element_text(hjust = 0.5),
          panel.background = element_rect(fill = "transparent"),
          panel.border = element_rect(fill = NA, colour = "black"),
          panel.grid.major = element_line(colour = NA),
          panel.grid.minor = element_line(colour = NA)) +
    scale_colour_manual(values = colr)

  ggsave(here("figures", paste0("Figure_S3", panel, ".pdf")), plot = p, width = 7, height = 5)
}
