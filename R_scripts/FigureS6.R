# =============================================================================
# Script:      FigureS6.R
# Description: Distribution of THEMIS-GC scores and protein biomarker
#              concentrations in patients with paired NGS and protein data.
#              (A) Box plots for Xijing sub-cohort (Training + Validation)
#              (B) Box plots for Testing Dataset 1 sub-cohort
#              Each panel shows THEMIS-GC score and AFP, CEA, CA125, CA72-4,
#              CA19-9 concentrations by disease stage (BENIGN / I / II / III / IV).
#
# --- Required Input Data ------------------------------------------------------
# File:   data/Table_S5.tsv
# Source: Supplementary Table S5
# Format: Tab-separated, columns:
#         ID, Diagnosis, Laurens' subtypes, Tumor stage, Cohort, Gender, Age,
#         AFP, CEA, CA125, CA72-4, CA19-9, THEMIS-GC
#
# --- Output -------------------------------------------------------------------
# Figures:  figures/Figure_S6A_XJ_subcohort_boxplot.pdf
#           figures/Figure_S6B_Testing1_subcohort_boxplot.pdf
# =============================================================================

library(gridExtra)
library(ggplot2)
library(grid)
library(rlang)
library(ggsci)
library(here)
library(systemfonts)

font_path <- here("fonts", "arial.ttf")
if (file.exists(font_path)) register_font(name = "Arial", plain = font_path)

# =============================================================================
# Load & harmonise Table S5
# =============================================================================
df <- read.table(here("data", "Table_S5.tsv"), header = TRUE, sep = "\t",
                  check.names = FALSE, stringsAsFactors = FALSE)

colnames(df)[colnames(df) == "THEMIS-GC"]   <- "Ensemble"
colnames(df)[colnames(df) == "Tumor stage"]  <- "Tumor_stages"
colnames(df)[colnames(df) == "CA72-4"]       <- "CA724"
colnames(df)[colnames(df) == "CA19-9"]       <- "CA199"

# Derive Stage: Tumor stage for cancer patients, Diagnosis for BENIGN
df_cancer     <- subset(df, Diagnosis == "Cancer")
df_cancer$Stage <- df_cancer$Tumor_stages
df_non_cancer <- subset(df, Diagnosis != "Cancer")
df_non_cancer$Stage <- df_non_cancer$Diagnosis
df <- rbind(df_cancer, df_non_cancer)
df <- subset(df, Stage != "Unknown")
df$Stage <- factor(df$Stage, levels = c("HEALTHY", "BENIGN", "I", "II", "III", "IV"))

markers <- c("Ensemble", "AFP", "CEA", "CA125", "CA724", "CA199")
ybreaks <- c(1, 10, 100, 200, 500, 800)

# =============================================================================
# Shared boxplot function
# =============================================================================
cohort_boxplot <- function(data, cohort_values, cohort_label, color, outfile) {
  plot.df <- subset(data, Cohort %in% cohort_values)
  plot.df$Cohort <- cohort_label

  colr        <- color
  names(colr) <- cohort_label

  plot.l <- vector("list", length(markers))
  for (i in seq_along(markers)) {
    model <- markers[i]
    if (model == "Ensemble") {
      plot.l[[i]] <- ggplot(plot.df, aes(x = Stage, y = !!sym(model), color = Cohort)) +
        geom_boxplot() +
        labs(title = model) +
        coord_cartesian(ylim = c(0, 1)) +
        theme_bw() +
        theme(legend.position    = "bottom",
              plot.title         = element_text(hjust = 0.5),
              panel.background   = element_rect(fill = "transparent"),
              panel.border       = element_rect(fill = NA, colour = "black"),
              panel.grid.major   = element_line(colour = NA),
              panel.grid.minor   = element_line(colour = NA)) +
        scale_colour_manual(values = colr)
    } else {
      plot.l[[i]] <- ggplot(plot.df, aes(x = Stage, y = !!sym(model), color = Cohort)) +
        geom_boxplot() +
        labs(title = model) +
        scale_y_continuous(trans = "log2", breaks = ybreaks, labels = ybreaks) +
        theme_bw() +
        theme(legend.position    = "bottom",
              plot.title         = element_text(hjust = 0.5),
              panel.background   = element_rect(fill = "transparent"),
              panel.border       = element_rect(fill = NA, colour = "black"),
              panel.grid.major   = element_line(colour = NA),
              panel.grid.minor   = element_line(colour = NA)) +
        scale_colour_manual(values = colr)
    }
  }

  pdf(outfile, width = 16, height = 10)
  do.call(grid.arrange, c(plot.l, ncol = 3))
  dev.off()
}

# =============================================================================
# Figure S6A - Xijing sub-cohort (Training + Validation)
# =============================================================================
cohort_boxplot(
  data          = df,
  cohort_values = c("Xijing-Training", "Xijing-Validation"),
  cohort_label  = "XJ_subset",
  color         = pal_startrek("uniform", alpha = 0.9)(7)[1],
  outfile       = here("figures", "Figure_S6A_XJ_subcohort_boxplot.pdf")
)

# =============================================================================
# Figure S6B - Testing Dataset 1 sub-cohort
# =============================================================================
cohort_boxplot(
  data          = df,
  cohort_values = "Testing Dataset 1",
  cohort_label  = "Testing Dataset 1",
  color         = pal_startrek("uniform", alpha = 0.9)(7)[4],
  outfile       = here("figures", "Figure_S6B_Testing1_subcohort_boxplot.pdf")
)
