# =============================================================================
# Script:      FigureS7.R
# Description: ROC curve comparison of THEMIS-GC vs protein biomarkers
#              (AFP, CEA, CA125, CA72-4, CA19-9) for GC detection.
#              (A) Xijing sub-cohort (Training + Validation combined)
#              (B) Testing Dataset 1 sub-cohort
#
# --- Required Input Data ------------------------------------------------------
# File:   data/Table_S5.tsv
# Source: Supplementary Table S5
# Format: Tab-separated, columns:
#         ID, Diagnosis, Laurens' subtypes, Tumor stage, Cohort, Gender, Age,
#         AFP, CEA, CA125, CA72-4, CA19-9, THEMIS-GC
#
# --- Output -------------------------------------------------------------------
# Figures:  figures/Figure_S7A_XJ_subcohort_ROC.pdf
#           figures/Figure_S7B_Testing1_ROC.pdf
# =============================================================================

library(pROC)
library(ggplot2)
library(ggpubr)
library(ggsci)
library(dplyr)
library(here)
library(systemfonts)

font_path <- here("fonts", "arial.ttf")
if (file.exists(font_path)) register_font(name = "Arial", plain = font_path)

# =============================================================================
# Load & harmonise Table S5
# =============================================================================
df <- read.table(here("data", "Table_S5.tsv"), header = TRUE, sep = "\t",
                  check.names = FALSE, stringsAsFactors = FALSE)

colnames(df)[colnames(df) == "THEMIS-GC"] <- "Ensemble"
colnames(df)[colnames(df) == "CA72-4"]    <- "CA724"
colnames(df)[colnames(df) == "CA19-9"]    <- "CA199"

df$GC_label <- ifelse(df$Diagnosis == "Cancer", "GC", "Non_GC")

markers <- c("Ensemble", "AFP", "CEA", "CA125", "CA724", "CA199")

# =============================================================================
# ROC function: all markers on one plot, one colour per marker
# =============================================================================
ROC_PLOT <- function(data, title_key, outfile) {
  cancer   <- "GC"
  roc_list <- list()
  auc_list <- list()
  label    <- ""

  for (marker in markers) {
    newdata   <- data[data$GC_label %in% c(cancer, "Non_GC"), ]
    samplenum <- sum(newdata$GC_label == cancer)
    rocobj    <- roc(newdata$GC_label, newdata[[marker]],
                      levels = c("Non_GC", cancer))
    roc_list[[marker]] <- rocobj
    auc_list[[marker]] <- auc(rocobj)[1]
    label <- paste0(label, marker, ": AUC=", round(auc(rocobj)[1], 2),
                     " (", round(ci(rocobj)[1], 2), "-", round(ci(rocobj)[3], 2), ")\n")
  }

  colr <- pal_startrek("uniform", alpha = 0.9)(length(roc_list))
  names(colr) <- names(roc_list)

  title_n <- paste0("GC=", samplenum, " Control=", sum(data$GC_label == "Non_GC"))

  p <- ggroc(roc_list) +
    ggtitle(paste0(title_key, " (", title_n, ")")) +
    theme_bw() +
    theme(plot.title       = element_text(hjust = 0.5),
          panel.background = element_rect(fill = "transparent"),
          panel.border     = element_rect(fill = NA, colour = "black"),
          panel.grid.major = element_line(colour = NA),
          panel.grid.minor = element_line(colour = NA),
          legend.position  = "none") +
    scale_colour_manual(values = colr)

  label_lines <- Filter(nchar, strsplit(label, "\n", fixed = TRUE)[[1]])
  for (i in seq_along(label_lines)) {
    p <- p + annotate(geom = "text", x = 0.35, y = 0.045 * i,
                       label = label_lines[i], color = colr[i])
  }

  pdf(outfile, width = 5, height = 5)
  print(p)
  dev.off()
}

# =============================================================================
# Figure S7A - Xijing sub-cohort (Training + Validation)
# =============================================================================
xj_data <- subset(df, Cohort %in% c("Xijing-Training", "Xijing-Validation"))

ROC_PLOT(xj_data, "Xijing sub-cohort",
          here("figures", "Figure_S7A_XJ_subcohort_ROC.pdf"))

# =============================================================================
# Figure S7B - Testing Dataset 1
# =============================================================================
test1_data <- subset(df, Cohort == "Testing Dataset 1")

ROC_PLOT(test1_data, "Testing Dataset 1",
          here("figures", "Figure_S7B_Testing1_ROC.pdf"))
