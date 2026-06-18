# =============================================================================
# Script:      FigureS4.R
# Description: Detection sensitivity of THEMIS-GC in clinical subgroups.
#              (A-B) Sensitivity at 95% specificity stratified by tumour stage
#                    (I/II/III/IV) x age (<= 65 / >65) — Training (A) and
#                    Validation (B) shown as side-by-side facets.
#              (C)   Sensitivity at 95% specificity by Laurens' subtype
#                    (DT/IT/MT) — Training and Validation together.
#              (D)   Sensitivity at 95% specificity by gender (female/male)
#                    — Training and Validation together.
#              95% CIs calculated by Clopper-Pearson method.
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
# Figures:     figures/Figure_S4AB.pdf   (Stage x Age, Training / Validation)
#              figures/Figure_S4C.svg    (Laurens' subtype, Training + Validation)
#              figures/Figure_S4D.svg    (Gender, Training + Validation)
#
# Notes:
#  - Only THEMIS-GC is evaluated here (renamed from "Ensemble" in the internal
#    analysis pipeline to match Table S3 column names).
#  - Data harmonisation is identical to Figure2.R and FigureS3.R; all three
#    scripts read from the same data/Table_S3.tsv.
# =============================================================================

library(tidyverse)
library(ggsci)
library(svglite)
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

# Stage: recombine Diagnosis + `Tumor stage` into a single factor
df$Stage <- ifelse(df$Diagnosis %in% c("HEALTHY", "BENIGN"),
                    df$Diagnosis, df$`Tumor stage`)

# Age group
df$Age_group <- ifelse(df$Age > 65, ">65", "<=65")

# Compound group for panels A-B: stage x age
df$Stage_age <- paste(df$Stage, df$Age_group, sep = ".")

# Group column: mirrors original pipeline logic —
# HEALTHY/BENIGN controls labelled "HEALTHY"; cancer rows get Laurens' subtype
df$Group <- ifelse(df$Diagnosis %in% c("HEALTHY", "BENIGN"),
                    "HEALTHY", df$`Laurens' subtypes`)

# Rename Laurens' subtypes to Laurens_subtype for use as a grouping variable
df$Laurens_subtype <- df$`Laurens' subtypes`

# Restrict to Training / Validation cohorts throughout
d <- subset(df, Cohort %in% c("Training", "Validation"))

# =============================================================================
# Shared helpers
# =============================================================================
negative_group  <- "HEALTHY"
marker          <- "THEMIS-GC"

# Clopper-Pearson 95% CI for a proportion
ci_calc <- function(total_number, positive_number, spf = 0.95) {
  n  <- total_number
  n1 <- positive_number
  alpha <- 0.05
  f1 <- qf(1 - alpha / 2, 2 * n1, 2 * (n - n1 + 1), lower.tail = FALSE)
  f2 <- qf(alpha / 2, 2 * (n1 + 1), 2 * (n - n1), lower.tail = FALSE)
  pl <- (1 + (n - n1 + 1) / (n1 * f1)) ^ (-1)
  pu <- (1 + (n - n1) / ((n1 + 1) * f2)) ^ (-1)
  pl <- ifelse(is.na(pl), 0, pl)
  pu <- ifelse(is.na(pu), 1, pu)
  paste0(round(pl * 100, 2), "-", round(pu * 100, 2))
}

# Derive per-95%-specificity cutoff from Training healthy controls
discovery_healthy <- d[d$Cohort == "Training" & d$Group == negative_group, ]
spf <- 0.95
sample_number <- round(nrow(discovery_healthy) * spf) + 1
discovery_healthy_sorted <- discovery_healthy[order(discovery_healthy[[marker]]), ]
cutoff <- discovery_healthy_sorted[[marker]][sample_number - 1]

# Generic sensitivity table builder
build_sen_table <- function(data, group_col, marker_col, cutoff, negative_group) {
  d_tmp <- data[, c("Cohort", "Stage", "Diagnosis", group_col, marker_col)]
  colnames(d_tmp)[4] <- "group"
  colnames(d_tmp)[5] <- "marker"
  d_tmp <- na.omit(d_tmp)

  cancer <- d_tmp[d_tmp$Diagnosis != negative_group, ]
  cancer$group <- "Overall"
  d_tmp <- rbind(d_tmp, cancer)

  d_tmp$marker <- as.numeric(d_tmp$marker)
  d_tmp$pred   <- ifelse(d_tmp$marker < cutoff, "negative", "positive")

  sen <- d_tmp %>%
    group_by(Cohort, Diagnosis, group, pred) %>%
    summarise(count = n(), .groups = "drop") %>%
    spread(pred, count)

  if (!"positive" %in% colnames(sen)) sen$positive <- 0
  if (!"negative" %in% colnames(sen)) sen$negative <- 0
  sen$positive <- ifelse(is.na(sen$positive), 0, sen$positive)
  sen$negative <- ifelse(is.na(sen$negative), 0, sen$negative)

  sen <- sen %>%
    mutate(total = negative + positive,
           Sens  = round((positive / total) * 100, 2),
           CI    = ci_calc(total, positive, spf),
           Spef  = spf)
  sen
}

# =============================================================================
# Panels A-B - THEMIS-GC sensitivity by Stage x Age, faceted by Cohort
# (Training = A, Validation = B)
# =============================================================================
ci_plot_by_cohort <- function(sen_table, title, group_name, outfile) {
  d <- sen_table
  d <- d[!grepl("Unknown", d$group), ]
  d <- d[!grepl("Overall", d$group), ]
  d <- d[d$Diagnosis != "HEALTHY", ]
  d <- d[!is.na(d$group), ]

  complete_combinations <- do.call(rbind, lapply(unique(d$Cohort), function(cohort) {
    expand.grid(Cohort = cohort, group = unique(d$group), Diagnosis = unique(d$Diagnosis))
  }))
  d <- merge(complete_combinations, d, by = c("Cohort", "group", "Diagnosis"), all.x = TRUE)
  d <- d %>%
    mutate(across(where(is.numeric), ~ ifelse(is.na(.), 0, .)),
           CI = ifelse(is.na(CI), "0-0", as.character(CI)))

  d$CI_low  <- as.numeric(gsub("-.*", "", d$CI))
  d$CI_high <- as.numeric(gsub(".*-", "", d$CI))
  d$Cohort  <- factor(d$Cohort, levels = c("Training", "Validation"))
  d <- d %>% group_by(group, Diagnosis) %>%
    mutate(group_number = paste0(total, collapse = "|"))
  d$Stage_new <- paste0(d$group, "\n(", d$group_number, ")")
  d$label <- paste0(d$Sens, "%")

  p <- ggplot(d, aes(x = Stage_new, y = Sens, ymin = CI_low, ymax = CI_high,
                      fill = Cohort, colour = Cohort)) +
    geom_errorbar(position = position_dodge(width = 0.7), width = 0.4, color = "black", size = 0.4) +
    geom_point(position = position_dodge(width = 0.7)) +
    scale_fill_manual(values  = c(pal_npg("nrc")(3), "black"), labels = c("Training", "Validation")) +
    scale_color_manual(values = c(pal_npg("nrc")(3), "black"), labels = c("Training", "Validation")) +
    theme_bw() +
    facet_wrap(vars(Cohort), ncol = 2, scales = "free") +
    theme(strip.background  = element_blank(),
          strip.text         = element_text(face = "bold", size = 12),
          plot.title         = element_text(size = 15, face = "bold"),
          panel.border       = element_rect(fill = NA, colour = "black"),
          panel.grid.major   = element_line(colour = NA),
          panel.grid.minor   = element_line(colour = NA)) +
    scale_y_continuous(limits = c(-10, 110),
                        labels = function(x) paste0(x, "%"),
                        breaks = seq(0, 100, by = 10)) +
    ylab("Sensitivity") + xlab(group_name) + ggtitle(title) +
    geom_text(aes(x = Stage_new, y = -5, label = label), size = 12 / .pt,
              position = position_dodge(width = 0.9))

  pdf(outfile, width = 16, height = 9)
  print(p)
  dev.off()
}

sen_stage_age <- build_sen_table(d, "Stage_age", marker, cutoff, negative_group)
write.table(sen_stage_age,
            here("figures", "Figure_S4AB_THEMIS-GC_Stage_age_sensitivity.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
ci_plot_by_cohort(sen_stage_age,
                  title      = "THEMIS-GC sensitivity by Stage x Age",
                  group_name = "Stage x Age group",
                  outfile    = here("figures", "Figure_S4AB.pdf"))

# =============================================================================
# Panels C-D - THEMIS-GC sensitivity by Laurens' subtype and Gender,
# Training and Validation shown together, faceted by Diagnosis
# =============================================================================
ci_plot_by_diagnosis <- function(sen_table, title, group_name, outfile) {
  d <- sen_table
  d <- d[d$Diagnosis != "HEALTHY", ]
  d <- d[d$group != "Unknown", ]
  d <- d[!is.na(d$group), ]

  complete_combinations <- do.call(rbind, lapply(unique(d$Cohort), function(cohort) {
    expand.grid(Cohort = cohort, group = unique(d$group), Diagnosis = unique(d$Diagnosis))
  }))
  d <- merge(complete_combinations, d, by = c("Cohort", "group", "Diagnosis"), all.x = TRUE)
  d <- d %>%
    mutate(across(where(is.numeric), ~ ifelse(is.na(.), 0, .)),
           CI = ifelse(is.na(CI), "0-0", as.character(CI)))

  d$CI_low  <- as.numeric(gsub("-.*", "", d$CI))
  d$CI_high <- as.numeric(gsub(".*-", "", d$CI))
  d$Cohort  <- factor(d$Cohort, levels = c("Training", "Validation"))
  d <- d %>% group_by(group, Diagnosis) %>%
    mutate(group_number = paste0(total, collapse = "|"))
  d$Stage_new <- paste0(d$group, "\n(", d$group_number, ")")
  d$label <- paste0(d$Sens, "%")

  p <- ggplot(d, aes(x = Stage_new, y = Sens, ymin = CI_low, ymax = CI_high,
                      fill = Cohort, colour = Cohort)) +
    geom_errorbar(position = position_dodge(width = 0.7), width = 0.4, color = "black", size = 0.4) +
    geom_point(position = position_dodge(width = 0.7)) +
    scale_fill_manual(values  = c(pal_npg("nrc")(3), "black"), labels = c("Training", "Validation")) +
    scale_color_manual(values = c(pal_npg("nrc")(3), "black"), labels = c("Training", "Validation")) +
    theme_bw() +
    facet_wrap(vars(Diagnosis), ncol = 1, scales = "free_x") +
    theme(strip.background = element_blank(),
          strip.text        = element_text(face = "bold", size = 12),
          plot.title        = element_text(size = 15, face = "bold")) +
    scale_y_continuous(limits = c(-10, 110),
                        labels = function(x) paste0(x, "%"),
                        breaks = seq(0, 100, by = 10)) +
    ylab("Sensitivity") + xlab(group_name) + ggtitle(title) +
    geom_text(aes(x = Stage_new, y = -5, label = label), size = 12 / .pt,
              position = position_dodge(width = 0.9))

  svglite(outfile, width = 10, height = 8, pointsize = 12, fix_text_size = FALSE)
  print(p)
  dev.off()
}

# Panel C - Laurens' subtype
sen_laurens <- build_sen_table(d, "Laurens_subtype", marker, cutoff, negative_group)
write.table(sen_laurens,
            here("figures", "Figure_S4C_THEMIS-GC_Laurens_subtype_sensitivity.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
ci_plot_by_diagnosis(sen_laurens,
                     title      = "THEMIS-GC sensitivity by Laurens' subtype",
                     group_name = "Laurens' subtype",
                     outfile    = here("figures", "Figure_S4C.svg"))

# Panel D - Gender
sen_gender <- build_sen_table(d, "Gender", marker, cutoff, negative_group)
write.table(sen_gender,
            here("figures", "Figure_S4D_THEMIS-GC_Gender_sensitivity.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
ci_plot_by_diagnosis(sen_gender,
                     title      = "THEMIS-GC sensitivity by Gender",
                     group_name = "Gender",
                     outfile    = here("figures", "Figure_S4D.svg"))
