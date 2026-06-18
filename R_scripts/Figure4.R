# =============================================================================
# Script:      Figure4.R
# Description: Comparison of a cfDNA-only model (THEMIS-GC) against a combined
#              cfDNA + protein biomarker model in the Xijing cohort.
#              (A) ROC: cfDNA feature vs cfDNA+AFP+CA125+CA199+CEA — Training
#              (B) Clinical landscape — Training dataset
#              (C) ROC — Validation dataset
#              (D) Clinical landscape — Validation dataset
#
# --- Required Input Data ------------------------------------------------------
# File 1:  data/model2_cfDNA_protein_features.tsv
#          Combined model scores (cfDNA + protein biomarkers); Ensemble column
#          = model2 score.
# File 2:  data/model1_cfDNA_features.tsv
#          cfDNA-only model scores; Ensemble column = THEMIS-GC.
# Source:  Both files correspond to participants listed in Supplementary
#          Table S5. Cohorts: Xijing-Training (n=156), Xijing-Validation (n=72).
# Format:  Tab-separated. Shared key columns:
#          GSA_individual_name, Diagnosis, Tumor_stages, Cohort, Stage,
#          AFP, CEA, CA125, CA199 (binary 0/1), Ensemble (model score)
#
# --- Output -------------------------------------------------------------------
# Figures:  figures/Figure_4A_Training_ROC.pdf
#           figures/Figure_4B_Training_landscape.pdf
#           figures/Figure_4C_Validation_ROC.pdf
#           figures/Figure_4D_Validation_landscape.pdf
# Results:  results/Figure_4A_Training.delong.stat.tsv
#           results/Figure_4C_Validation.delong.stat.tsv
# =============================================================================
library(readxl)
library(tidyverse)
library(pROC)
library(ggsci)
library(ComplexHeatmap)
library(circlize)
library(here)
library(systemfonts)

font_path <- here("fonts", "arial.ttf")
if (file.exists(font_path)) register_font(name = "Arial", plain = font_path)
project_dir=here("figures")
if (!dir.exists(project_dir)) {
  dir.create(project_dir)
}
# =============================================================================
# Load & harmonise data
# =============================================================================
dat_model2 <- read.table(here("data", "cfDNA_features+AFP+CA125+CA199+CEA.new.data.df.xls"),
                           header = TRUE, sep = "\t", check.names = FALSE,
                           stringsAsFactors = FALSE)
dat_model2$`cfDNA+AFP+CA125+CA199+CEA` <- dat_model2$Ensemble

dat_model1 <- read.table(here("data", "cfDNA_features.new.data.df.xls"),
                           header = TRUE, sep = "\t", check.names = FALSE,
                           stringsAsFactors = FALSE)
dat_model1$`cfDNA feature` <- dat_model1$Ensemble
dat_model1 <- dat_model1[, c("GSA_individual_name", "cfDNA feature")]

df <- merge(dat_model2, dat_model1, by = "GSA_individual_name")
df <- na.omit(df)

df$Cohort <- gsub("Testing Dataset 1", "testing1", df$Cohort)
df$GC_label <- ifelse(df$Diagnosis == "Cancer", "GC", "Non_GC")

roc_markers <- c("cfDNA feature", "cfDNA+AFP+CA125+CA199+CEA")
new_colors   <- pal_startrek("uniform", alpha = 0.9)(7)

# =============================================================================
# ROC function with DeLong pairwise test
# =============================================================================
ROC_PLOT <- function(data, response, predictor, title_key, fig_prefix) {
  cancer   <- "GC"
  roc_list <- list()
  auc_list <- list()
  label    <- ""

  for (marker in predictor) {
    newdata   <- data[data[[response]] %in% c(cancer, "Non_GC"), ]
    samplenum <- sum(newdata[[response]] == cancer)
    rocobj    <- roc(newdata[[response]], newdata[[marker]], levels = c("Non_GC", cancer))
    roc_list[[marker]] <- rocobj
    auc_list[[marker]] <- auc(rocobj)[1]
    label <- paste0(label, marker, ": AUC=", round(auc(rocobj)[1], 2),
                     " (", round(ci(rocobj)[1], 2), "-", round(ci(rocobj)[3], 2), ")\n")
  }

  result_df <- data.frame(Marker1 = character(), Marker2 = character(),
                            Z_value = numeric(), P_value = numeric(),
                            AUC_diff = numeric(), CI_lower = numeric(), CI_upper = numeric(),
                            stringsAsFactors = FALSE)
  marker_pairs <- combn(names(roc_list), 2)
  for (i in 1:ncol(marker_pairs)) {
    m1 <- marker_pairs[1, i]; m2 <- marker_pairs[2, i]
    tst <- roc.test(roc_list[[m1]], roc_list[[m2]], method = "delong")
    result_df <- bind_rows(result_df, data.frame(
      Marker1 = m1, Marker2 = m2, Z_value = tst$statistic, P_value = tst$p.value,
      AUC_diff = auc_list[[m1]] - auc_list[[m2]],
      CI_lower = tst$conf.int[1], CI_upper = tst$conf.int[2]))
  }

  pval_numeric <- result_df$P_value[1]
  pval_label   <- ifelse(pval_numeric < 0.001, "<0.001", paste0("= ", round(pval_numeric, 3)))

  result_df <- result_df %>% mutate(
    P_value  = ifelse(P_value < 1e-20, "<1e-20", as.character(round(P_value, 3))),
    AUC_diff = round(AUC_diff, 3), CI_lower = round(CI_lower, 3), CI_upper = round(CI_upper, 3))
  write.table(result_df,
              here("figures", paste0(fig_prefix, ".delong.stat.tsv")),
              sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

  colr    <- pal_startrek("uniform", alpha = 0.9)(7)
  names(colr) <- factor(predictor)
  title_n <- paste0("GC=", samplenum, " Control=", sum(data[[response]] == "Non_GC"))

  p <- ggroc(roc_list) +
    ggtitle(paste0(title_key, " (", title_n, ")")) +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5),
          panel.background = element_rect(fill = "transparent"),
          panel.border = element_rect(fill = NA, colour = "black"),
          panel.grid.major = element_line(colour = NA),
          panel.grid.minor = element_line(colour = NA),
          legend.position = "none") +
    scale_colour_manual(values = colr) +
    annotate(geom = "text", x = 0.35, y = 0.4,
             label = paste0("DeLong p-value ", pval_label))

  label_lines <- strsplit(label, "\n", fixed = TRUE)[[1]]
  label_lines <- label_lines[nchar(label_lines) > 0]
  for (i in seq_along(label_lines)) {
    p <- p + annotate(geom = "text", x = 0.35, y = 0.045 * i,
                       label = label_lines[i], color = colr[i])
  }

  pdf(here("figures", paste0(fig_prefix, "_ROC.pdf")), width = 5, height = 5)
  print(p)
  dev.off()
  return(auc_list)
}

# =============================================================================
# Landscape function (ComplexHeatmap annotation)
# =============================================================================
make_landscape <- function(dat_baseline, outfile) {
  dat_baseline$AFP   <- as.character(dat_baseline$AFP)
  dat_baseline$CEA   <- as.character(dat_baseline$CEA)
  dat_baseline$CA199 <- as.character(dat_baseline$CA199)
  dat_baseline$CA125 <- as.character(dat_baseline$CA125)

  ht_list <- HeatmapAnnotation(
    Diagnosis        = dat_baseline$Diagnosis,
    `Tumor stages`   = dat_baseline$Tumor_stages,
    `cfDNA feature`  = dat_baseline$`cfDNA feature`,
    `cfDNA+protein`  = dat_baseline$`cfDNA+AFP+CA125+CA199+CEA`,
    AFP              = dat_baseline$AFP,
    CEA              = dat_baseline$CEA,
    `CA19-9`         = dat_baseline$CA199,
    CA125            = dat_baseline$CA125,
    annotation_legend_param = list(direction = "horizontal", nrow = 2),
    col = list(
      Diagnosis      = c(BENIGN = new_colors[2], Cancer = new_colors[1]),
      `Tumor stages` = c(I = new_colors[3], II = new_colors[7], III = new_colors[4],
                          IV = new_colors[1], Unknown = new_colors[5], BENIGN = new_colors[2]),
      `cfDNA feature` = colorRamp2(c(0, 0.5, 1),
                                    c("#5C88DAE5", "#FFFFFFFF", "#CC0C00E5")),
      `cfDNA+protein` = colorRamp2(c(0, 0.5, 1),
                                    c("#5C88DAE5", "#FFFFFFFF", "#CC0C00E5")),
      AFP     = c("1" = new_colors[1], "0" = new_colors[5], None = new_colors[7]),
      CEA     = c("1" = new_colors[1], "0" = new_colors[5], None = new_colors[7]),
      `CA19-9` = c("1" = new_colors[1], "0" = new_colors[5], None = new_colors[7]),
      CA125   = c("1" = new_colors[1], "0" = new_colors[5], None = new_colors[7])
    )
  ) %v% NULL

  pdf(outfile, onefile = FALSE)
  draw(ht_list, merge_legend = TRUE,
       annotation_legend_side = "bottom", heatmap_legend_side = "bottom")
  dev.off()
}


main <- function() {
# =============================================================================
# Figure 4A - ROC Training dataset
# =============================================================================
ROC_PLOT(df[df$Cohort == "Training", ], "GC_label", roc_markers,
          "Training", "Figure_4A_Training")

# =============================================================================
# Figure 4B - Landscape Training dataset
# =============================================================================
dat_train <- subset(df, Cohort == "Training")
dat_train$Stage <- factor(dat_train$Stage,
                            levels = c("BENIGN", "I", "II", "III", "IV", "Unknown"))
dat_train <- dat_train[order(dat_train$Stage, dat_train$`cfDNA+AFP+CA125+CA199+CEA`), ]
make_landscape(dat_train, here("figures", "Figure_4B_Training_landscape.pdf"))

# =============================================================================
# Figure 4C - ROC Validation dataset
# =============================================================================
ROC_PLOT(df[df$Cohort == "Validation", ], "GC_label", roc_markers,
          "Validation", "Figure_4C_Validation")

# =============================================================================
# Figure 4D - Landscape Validation dataset
# =============================================================================
dat_valid <- subset(df, Cohort == "Validation")
dat_valid$Stage <- factor(dat_valid$Stage,
                            levels = c("BENIGN", "I", "II", "III", "IV", "Unknown"))
dat_valid <- dat_valid[order(dat_valid$Stage, dat_valid$`cfDNA+AFP+CA125+CA199+CEA`), ]
make_landscape(dat_valid, here("figures", "Figure_4D_Validation_landscape.pdf"))
      }

if (sys.nframe() == 0) {
          main()
        }

