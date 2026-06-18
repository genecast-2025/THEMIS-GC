# =============================================================================
# Script:      FigureS5.R
# Description: Distribution of individual cfDNA features in the two external
#              datasets (Testing Dataset 1 and Testing Dataset 2 / MONITOR).
#              (A) MFR  (B) FSI  (C) FEM  (D) CAFF
#              Box plots across healthy/benign individuals and GC patients at
#              different disease stages; centre line = median, box limits =
#              upper/lower quartiles.
#
# --- Required Input Data ------------------------------------------------------
# File:        data/Table_S3.tsv
# Source:      Manuscript Supplementary Table S3
# Format:      Tab-separated, columns:
#              ID, Diagnosis, Laurens' subtypes, Tumor stage, Cohort, Gender,
#              Age, MFR, FSI, FEM, CAFF, THEMIS-GC
#
# --- Output -------------------------------------------------------------------
# Figures:     figures/Figure_S5A.pdf   (MFR)
#              figures/Figure_S5B.pdf   (FSI)
#              figures/Figure_S5C.pdf   (FEM)
#              figures/Figure_S5D.pdf   (CAFF)
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

df$Cohort <- gsub("Xijing-", "", df$Cohort)
df$Cohort <- gsub("Testing Dataset 1", "testing1", df$Cohort)
df$Cohort <- gsub("Testing Dataset 2", "testing2", df$Cohort)
df$Cohort <- gsub("testing2", "MONITOR", df$Cohort)

df$Stage <- ifelse(df$Diagnosis %in% c("HEALTHY", "BENIGN"),
                    df$Diagnosis, df$`Tumor stage`)

# =============================================================================
# Figures S5A-D - Individual cfDNA features by stage, Testing Dataset 1 vs MONITOR
# =============================================================================
panels <- list(A = "MFR", B = "FSI", C = "FEM", D = "CAFF")

plot.df <- subset(df, Cohort %in% c("testing1", "MONITOR"))
plot.df$Cohort <- factor(plot.df$Cohort, levels = c("testing1", "MONITOR"))
plot.df$Stage  <- factor(plot.df$Stage,
                          levels = c("HEALTHY", "BENIGN", "I", "II", "III", "IV", "Unknown"))

colr <- pal_startrek("uniform", alpha = 0.9)(7)[c(4, 6)]
names(colr) <- c("testing1", "MONITOR")

for (panel in names(panels)) {
  model <- panels[[panel]]

  p <- ggplot(plot.df, aes(x = Stage, y = !!sym(model), color = Cohort)) +
    geom_boxplot() +
    labs(title = paste0(model, " score")) +
    coord_cartesian(ylim = c(0, 1)) +
    theme_bw() +
    facet_wrap(vars(Cohort), ncol = 2, scales = "free") +
    theme(legend.position = "bottom",
          plot.title = element_text(hjust = 0.5),
          panel.background = element_rect(fill = "transparent"),
          panel.border = element_rect(fill = NA, colour = "black"),
          panel.grid.major = element_line(colour = NA),
          panel.grid.minor = element_line(colour = NA)) +
    scale_colour_manual(values = colr)

  ggsave(here("figures", paste0("Figure_S5", panel, ".pdf")), plot = p, width = 8, height = 5)
}
