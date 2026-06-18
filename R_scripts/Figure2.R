# =============================================================================
# Script:      Figure2.R
# Description: Development of the ensemble model (THEMIS-GC) for cancer
#              detection in the Xijing cohort.
#              (A) Box plot of THEMIS-GC score by stage, Training/Validation
#              (B-C) ROC curves for THEMIS-GC and individual cfDNA features
#                    (MFR, FSI, FEM, CAFF) - Training (B) and Validation (C)
#              (D) Sensitivity at 95% specificity by tumor stage and age
#                  group, Training and Validation
#              (E-F) Decision curve analysis (DCA) - Training (E) and
#                    Validation (F)
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
# Figures:     figures/Figure_2A.pdf
#              figures/Figure_2_Training_ROC_multi_set.pdf   (Figure 2B)
#              figures/Figure_2_Validation_ROC_multi_set.pdf (Figure 2C)
#              figures/Figure2D_<marker>_<Stage|Age_group>_Cancer_spf0.95_senci.hist.svg
#              figures/Training.dca_final_matched.pdf        (Figure 2E)
#              figures/Validation.dca_final_matched.pdf       (Figure 2F)
#
# =============================================================================
library(readxl)
library(tidyverse)
library(ggsci)
library(ggpubr)
library(rlang)
library(pROC)
library(rmda)
library(svglite)
library(systemfonts)
library(here)

font_path <- here("fonts", "arial.ttf")
if (file.exists(font_path)) register_font(name = "Arial", plain = font_path)
project_dir=here("figures")
if (!dir.exists(project_dir)) {
  dir.create(project_dir)
}
# =============================================================================
# Load & harmonise Table S3
# =============================================================================
df <- read_excel(
  here("data", "Table_S1-S9.xlsx"),
  sheet = "Table S3",
  skip = 1
)

df$`THEMIS-GC` = df$`OMNI-GC`

# Cohort: Table S3 uses "Xijing-Training" / "Xijing-Validation" / "Testing Dataset 1/2"
df$Cohort <- gsub("Xijing-", "", df$Cohort)
df$Cohort <- gsub("Testing Dataset 1", "testing1", df$Cohort)
df$Cohort <- gsub("Testing Dataset 2", "testing2", df$Cohort)

# Stage (HEALTHY/BENIGN/I-IV/Unknown): Table S3 splits this into Diagnosis
# (HEALTHY/BENIGN/Cancer) and `Tumor stage` (I-IV/Unknown, NA for non-cancer)
df$Stage <- ifelse(df$Diagnosis %in% c("HEALTHY", "BENIGN"),
                    df$Diagnosis, df$`Tumor stage`)

# Age group for panel D
df$Age_group <- ifelse(df$Age > 65, ">65", "<=65")

# Binary cancer indicator for DCA (panels E-F)
df$GROUP <- ifelse(df$Diagnosis %in% c("HEALTHY", "BENIGN"), 0, 1)

# GC / Non_GC label for ROC (panels B-C)
df$GC_label <- ifelse(df$Diagnosis == "Cancer", "GC", "Non_GC")


# =============================================================================
# Figure 2A - THEMIS-GC score by stage, Training vs Validation
# =============================================================================
model <- "THEMIS-GC"

fig2a.df <- subset(df, Cohort %in% c("Training", "Validation"))
fig2a.df$Cohort <- factor(fig2a.df$Cohort, levels = c("Training", "Validation"))
fig2a.df$Stage  <- factor(fig2a.df$Stage,
                           levels = c("HEALTHY", "BENIGN", "I", "II", "III", "IV", "Unknown"))

colr_a <- pal_startrek("uniform", alpha = 0.9)(2)
names(colr_a) <- c("Training", "Validation")

p_2a <- ggplot(fig2a.df, aes(x = Stage, y = !!sym(model), color = Cohort)) +
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
  scale_colour_manual(values = colr_a)

ggsave(here("figures", "Figure_2A.pdf"), plot = p_2a, width = 7, height = 5)


# =============================================================================
# Figure 2B-C - ROC curves (THEMIS-GC, MFR, FSI, FEM, CAFF)
# =============================================================================
roc_markers <- c("THEMIS-GC", "MFR", "FSI", "FEM", "CAFF")

AUC_LIST <- function(auclist) {
  result <- data.frame()
  for (index in seq_along(auclist)) {
    cancer_name <- names(auclist)[index]
    cancer_auc <- round(auclist[[index]], 3)
    result <- bind_rows(result, as.data.frame(t(c(cancer_name, cancer_auc))))
  }
  colnames(result) <- c("Cancer_type", "AUC")
  return(result)
}


ROC_PLOT <- function(data,
                     response,
                     predictor,
                     title_key,
                     outfile) {

  cancer <- "GC"

  ## 保留GC和对照
  newdata <- data[data[[response]] %in% c(cancer, "Non_GC"), ]

  ## 样本数
  samplenum  <- sum(newdata[[response]] == cancer, na.rm = TRUE)
  controlnum <- sum(newdata[[response]] == "Non_GC", na.rm = TRUE)

  ## 存储对象
  roc_list <- list()
  auc_list <- list()
  label <- ""

  ## ROC计算
  for (marker in predictor) {

    rocobj <- pROC::roc(
      response  = newdata[[response]],
      predictor = as.numeric(newdata[[marker]]),
      levels = c("Non_GC", cancer),
      quiet = TRUE
    )

    roc_list[[marker]] <- rocobj

    auc_value <- as.numeric(pROC::auc(rocobj))
    auc_list[[marker]] <- auc_value

    ci_auc <- pROC::ci.auc(rocobj)

    label <- paste0(
      label,
      marker,
      ": AUC=",
      sprintf("%.2f", auc_value),
      " (",
      sprintf("%.2f", ci_auc[1]),
      "-",
      sprintf("%.2f", ci_auc[3]),
      ")\n"
    )
  }

  ## 颜色
  colr <- pal_startrek("uniform", alpha = 0.9)(length(predictor))
  names(colr) <- predictor

  ## 标题
  title <- paste0(
    "GC=", samplenum,
    "  Control=", controlnum
  )

  ## ROC图
  p <- ggroc(roc_list, size = 1.2) +
    ggtitle(
      paste0(
        title_key,
        " (",
        title,
        ")"
      )
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(
        hjust = 0.5,
        size = 14,
        face = "bold"
      ),
      panel.background = element_rect(fill = "transparent"),
      panel.border = element_rect(
        fill = NA,
        colour = "black"
      ),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "none"
    ) +
    scale_colour_manual(values = colr) +
    labs(
      x = "1 - Specificity",
      y = "Sensitivity"
    )

  ## 添加AUC文字
  label_lines <- strsplit(
    trimws(label),
    "\n",
    fixed = TRUE
  )[[1]]

  for (i in seq_along(label_lines)) {

    p <- p +
      annotate(
        "text",
        x = 0.75,
        y = 0.05 * i,
        label = label_lines[i],
        color = colr[i],
        hjust = 0,
        size = 4
      )
  }

  ## 输出PDF
  pdf(outfile,
      width = 5,
      height = 5)

  print(p)

  dev.off()

  return(auc_list)
}

train_data <- subset(df, Cohort == "Training")

auclist_training <- ROC_PLOT(
  train_data,
  response = "GC_label",
  predictor = roc_markers,
  title_key = "Training",
  outfile = here(
    "figures",
    "Figure_2_Training_ROC_multi_set.pdf"
  )
)

validation_data <- subset(df, Cohort == "Validation")

auclist_validation <- ROC_PLOT(
  validation_data,
  response = "GC_label",
  predictor = roc_markers,
  title_key = "Validation",
  outfile = here(
    "figures",
    "Figure_2_Validation_ROC_multi_set.pdf"
  )
)


result_training   <- AUC_LIST(auclist_training)
result_validation <- AUC_LIST(auclist_validation)

# =============================================================================
# Figure 2D - Sensitivity at 95% specificity by Stage and Age group
# =============================================================================
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
  pl <- round(pl * 100, 2)
  pu <- round(pu * 100, 2)
  paste0(as.character(pl), "-", as.character(pu))
}

ci_plot <- function(sen_table, title, group_name, outfile) {
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

  d <- d %>% group_by(group, Diagnosis) %>% mutate(group_number = paste0(total, collapse = "|"))
  d$Stage_new <- paste0(d$group, "\n", "(", d$group_number, ")")
  d$label <- paste0(d$Sens, "%")

  p <- ggplot(data = d, aes(x = Stage_new, y = Sens, ymin = CI_low, ymax = CI_high,
                             fill = Cohort, colour = Cohort)) +
    geom_errorbar(position = position_dodge(width = 0.7), width = 0.4, color = "black", size = 0.4) +
    geom_point(position = position_dodge(width = 0.7)) +
    scale_fill_manual(values = c(pal_npg("nrc")(3), "black"), labels = c("Training", "Validation")) +
    scale_color_manual(values = c(pal_npg("nrc")(3), "black"), labels = c("Training", "Validation")) +
    theme_bw() +
    facet_wrap(vars(Diagnosis), ncol = 1, scales = "free_x") +
    theme(strip.background = element_blank(), strip.text = element_text(face = "bold", size = 12),
          plot.title = element_text(size = 15, face = "bold")) +
    scale_y_continuous(limits = c(-10, 110), labels = function(x) paste0(x, "%"), breaks = seq(0, 100, by = 10)) +
    ylab("Sensitivity") +
    xlab(group_name) +
    ggtitle(title) +
    geom_text(aes(x = Stage_new, y = -5, label = label), size = 12 / .pt,
              position = position_dodge(width = 0.9))

  svglite(outfile, width = 10, height = 8, pointsize = 12, fix_text_size = FALSE)
  print(p)
  dev.off()
}

sens_markers <- c("MFR", "FSI", "CAFF", "FEM", "THEMIS-GC")
sens_groups  <- c("Stage", "Age_group")
negative_group <- "HEALTHY"

d <- subset(df, Cohort %in% c("Training", "Validation"))
d[d=="BENIGN"] = negative_group
discovery_healthy <- d[d$Cohort == "Training" & d$Diagnosis == negative_group, ]

for (marker in sens_markers) {
  for (group in sens_groups) {
    spf <- 0.95
    sample_number <- round(nrow(discovery_healthy) * spf) + 1
    discovery_healthy_tmp <- discovery_healthy[order(discovery_healthy[, marker]), ]
    cutoff <- discovery_healthy_tmp[[marker]][sample_number - 1]

    d_tmp <- d[, c("Cohort", "Stage", "Diagnosis", group, marker)]
    colnames(d_tmp)[4] <- "group"
    colnames(d_tmp)[5] <- "marker"
    d_tmp <- na.omit(d_tmp)

    cancer <- d_tmp[d_tmp$Diagnosis != negative_group, ]
    cancer$group <- "Overall"
    d_tmp <- rbind(d_tmp, cancer)

    d_tmp$marker <- as.numeric(d_tmp$marker)
    d_tmp$pred <- ifelse(d_tmp$marker < cutoff, "negative", "positive")

    sen_result <- d_tmp %>%
      group_by(Cohort, Diagnosis, group, pred) %>%
      summarise(count = n(), .groups = "drop") %>%
      spread(pred, count)

    if (!"positive" %in% colnames(sen_result)) sen_result$positive <- 0
    if (!"negative" %in% colnames(sen_result)) sen_result$negative <- 0
    sen_result$positive <- ifelse(is.na(sen_result$positive), 0, sen_result$positive)
    sen_result$negative <- ifelse(is.na(sen_result$negative), 0, sen_result$negative)

    sen_result <- sen_result %>%
      mutate(total = negative + positive,
             Sens = round((positive / total) * 100, 2),
             CI = ci_calc(total, positive, spf),
             Spef = spf)

    out_prefix <- paste0("Figure2D_", marker, "_", group, "_Cancer_spf", spf, "_senci")
    write.table(sen_result, here("figures", paste0(out_prefix, ".tsv")), sep = "\t", quote = FALSE, row.names = FALSE)
    ci_plot(sen_result, out_prefix, gsub("_", " ", group),
            here("figures", paste0(out_prefix, ".hist.svg")))
  }
}


# =============================================================================
# Figure 2E-F - Decision curve analysis (Training / Validation)
# =============================================================================
dca_markers <- c("MFR", "FSI", "CAFF", "FEM", "THEMIS-GC")
dca_cols <- (pal_startrek("uniform", alpha = 0.9)(7))[c(1, 4, 5, 7, 6)]

dca_cohorts <- list("Training" = subset(df, Cohort == "Training"),
                     "Validation" = subset(df, Cohort == "Validation"))

for (n in names(dca_cohorts)) {
  dca_data <- dca_cohorts[[n]][, c("GROUP", dca_markers)]
  dca_data$GROUP <- as.numeric(as.character(dca_data$GROUP))
  dca_data <- na.omit(dca_data)

  dca_results <- list()
  for (marker in dca_markers) {
    set.seed(123456)
    # marker names containing "-" (e.g. "THEMIS-GC") need backticks to form a
    # valid formula
    dca_results[[marker]] <- decision_curve(
      formula = as.formula(paste0("GROUP ~ `", marker, "`")),
      data = dca_data,
      study.design = "cohort",
      thresholds = seq(0, 1, by = 0.01),
      bootstraps = 0,
      confidence.intervals = FALSE
    )
  }

  for (marker in dca_markers) {
    dca_table <- data.frame(
      threshold = dca_results[[marker]]$derived.data$thresholds,
      net_benefit = dca_results[[marker]]$derived.data$NB,
      marker = marker
    )
    write.table(dca_table, here("figures", paste(n, marker, "dca_results.tsv", sep = ".")),
                row.names = FALSE, sep = "\t")
  }

  pdf(here("figures", paste(n, "dca_final_matched.pdf", sep = ".")), width = 8, height = 6)
  plot_decision_curve(
    dca_results,
    curve.names = dca_markers,
    col = dca_cols,
    confidence.intervals = FALSE,
    cost.benefit.axis = FALSE,
    legend.pos = "topright",
    xlab = "Threshold Probability",
    ylab = "Net Benefit",
    main = paste0("DCA - ", n)
  )
  dev.off()
}
