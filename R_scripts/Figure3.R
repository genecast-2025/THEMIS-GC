# =============================================================================
# Script:      Figure3.R
# Description: Performance validation of THEMIS-GC for cancer detection in
#              the two external cohorts.
#              (A-B) THEMIS-GC score distribution — Testing Dataset 1 (A) and
#                    Testing Dataset 2 / MONITOR cohort (B)
#              (C)   ROC curves for THEMIS-GC in Testing Dataset 1 and
#                    Testing Dataset 2 on the same plot
#              (D-E) THEMIS-GC sensitivity at 95% specificity by tumour stage
#                    — Testing Dataset 1 (D) and Testing Dataset 2 (E)
#              (F-G) Decision curve analysis for THEMIS-GC — Testing Dataset 1
#                    (F) and Testing Dataset 2 / MONITOR cohort (G)
#
# --- Required Input Data ------------------------------------------------------
# File:        data/Table_S3.tsv
# Source:      Manuscript Supplementary Table S3
# Format:      Tab-separated, columns:
#              ID, Diagnosis, Laurens' subtypes, Tumor stage, Cohort, Gender,
#              Age, MFR, FSI, FEM, CAFF, THEMIS-GC
# =============================================================================
library(readxl)
library(tidyverse)
library(ggsci)
library(rlang)
library(pROC)
library(rmda)
library(here)
library(systemfonts)

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


df$Cohort <- gsub("Xijing-", "", df$Cohort)
df$Cohort <- gsub("Testing Dataset 1", "testing1", df$Cohort)
df$Cohort <- gsub("Testing Dataset 2", "testing2", df$Cohort)
df$Cohort <- gsub("testing2", "MONITOR", df$Cohort)

df$Stage <- ifelse(df$Diagnosis %in% c("HEALTHY", "BENIGN"),
                    df$Diagnosis, df$`Tumor stage`)

df$GC_label <- ifelse(df$Diagnosis == "Cancer", "GC", "Non_GC")
df$GROUP    <- ifelse(df$Diagnosis %in% c("HEALTHY", "BENIGN"), 0, 1)
df$Group    <- ifelse(df$Diagnosis %in% c("HEALTHY", "BENIGN"), "HEALTHY",
                       df$`Laurens' subtypes`)

colr_ext <- pal_startrek("uniform", alpha = 0.9)(7)[c(4, 6)]
names(colr_ext) <- c("testing1", "MONITOR")

# =============================================================================
# Figures 3A-B - THEMIS-GC score by stage, Testing Dataset 1 vs MONITOR
# (testing1 = A, MONITOR = B shown as side-by-side facets)
# =============================================================================
box.df <- subset(df, Cohort %in% c("testing1", "MONITOR"))
box.df$Cohort <- factor(box.df$Cohort, levels = c("testing1", "MONITOR"))
box.df$Stage  <- factor(box.df$Stage,
                         levels = c("HEALTHY", "BENIGN", "I", "II", "III", "IV", "Unknown"))

p_3ab <- ggplot(box.df, aes(x = Stage, y = `THEMIS-GC`, color = Cohort)) +
  geom_boxplot() +
  labs(title = "THEMIS-GC score") +
  coord_cartesian(ylim = c(0, 1)) +
  theme_bw() +
  facet_wrap(vars(Cohort), ncol = 2, scales = "free") +
  theme(legend.position = "bottom",
        plot.title = element_text(hjust = 0.5),
        panel.background = element_rect(fill = "transparent"),
        panel.border = element_rect(fill = NA, colour = "black"),
        panel.grid.major = element_line(colour = NA),
        panel.grid.minor = element_line(colour = NA)) +
  scale_colour_manual(values = colr_ext)

ggsave(here("figures", "Figure_3AB.pdf"), plot = p_3ab, width = 8, height = 5)

# =============================================================================
# Figure 3C - ROC curves for THEMIS-GC, Testing Dataset 1 and MONITOR
# (one curve per cohort on the same axes)
# =============================================================================
roc_cohorts <- c("testing1", "MONITOR")
marker_roc  <- "THEMIS-GC"

roc_list  <- list()
auc_label <- ""

for (cc in roc_cohorts) {
  newdata    <- subset(df, Cohort == cc & GC_label %in% c("GC", "Non_GC"))
  cohort_num <- nrow(newdata)
  samplenum  <- sum(newdata$GC_label == "GC")
  rocobj     <- roc(newdata$GC_label, newdata[[marker_roc]], levels = c("Non_GC", "GC"))
  roc_list[[cc]] <- rocobj
  auc_label <- paste0(auc_label, cc, " (n=", cohort_num, "): AUC=",
                       round(auc(rocobj)[1], 2),
                       " (", round(ci(rocobj)[1], 2), "-", round(ci(rocobj)[3], 2), ")\n")
}

p_3c <- ggroc(roc_list) +
  ggtitle(marker_roc) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5),
        panel.background = element_rect(fill = "transparent"),
        panel.border = element_rect(fill = NA, colour = "black"),
        panel.grid.major = element_line(colour = NA),
        panel.grid.minor = element_line(colour = NA),
        legend.position = "none") +
  scale_colour_manual(values = colr_ext)

label_lines <- strsplit(auc_label, "\n", fixed = TRUE)[[1]]
for (i in seq_along(label_lines)) {
  p_3c <- p_3c + annotate(geom = "text", x = 0.45, y = 0.045 * i,
                            label = label_lines[i], color = colr_ext[i])
}

pdf(here("figures", "Figure_3C_THEMIS-GC_ROC_multi_set.pdf"), width = 5, height = 5)
print(p_3c)
dev.off()

# =============================================================================
# Figures 3D-E - THEMIS-GC sensitivity at 95% specificity by Stage
# (Testing Dataset 1 = D, MONITOR = E, shown as side-by-side facets)
# =============================================================================
ci_calc <- function(total_number, positive_number, spf = 0.95) {
  n <- total_number; n1 <- positive_number; alpha <- 0.05
  f1 <- qf(1 - alpha / 2, 2 * n1, 2 * (n - n1 + 1), lower.tail = FALSE)
  f2 <- qf(alpha / 2, 2 * (n1 + 1), 2 * (n - n1), lower.tail = FALSE)
  pl <- (1 + (n - n1 + 1) / (n1 * f1)) ^ (-1)
  pu <- (1 + (n - n1) / ((n1 + 1) * f2)) ^ (-1)
  pl <- ifelse(is.na(pl), 0, pl); pu <- ifelse(is.na(pu), 1, pu)
  paste0(round(pl * 100, 2), "-", round(pu * 100, 2))
}

negative_group <- "HEALTHY"
df[df=="BENIGN"] = negative_group

spf <- 0.95
discovery_healthy <- df[df$Cohort == "Training" & df$Group == "HEALTHY", ]
sample_number <- round(nrow(discovery_healthy) * spf) + 1
cutoff <- sort(discovery_healthy[["THEMIS-GC"]])[sample_number - 1]

d_ext <- subset(df, Cohort %in% c("testing1", "MONITOR"))

d_tmp <- d_ext[, c("Cohort", "Stage", "Diagnosis", "Stage", "THEMIS-GC")]
d_tmp <- d_ext[, c("Cohort", "Stage", "Diagnosis", "THEMIS-GC")]
colnames(d_tmp)[4] <- "marker"
d_tmp <- na.omit(d_tmp)

cancer_rows <- d_tmp[d_tmp$Diagnosis != "HEALTHY", ]
cancer_rows$group <- "Overall"
d_tmp$group <- d_tmp$Stage
d_tmp <- rbind(d_tmp, cancer_rows[, c("Cohort", "Stage", "Diagnosis", "marker", "group")])

d_tmp$marker <- as.numeric(d_tmp$marker)
d_tmp$pred   <- ifelse(d_tmp$marker < cutoff, "negative", "positive")

sen_stage <- d_tmp %>%
  group_by(Cohort, Diagnosis, group, pred) %>%
  summarise(count = n(), .groups = "drop") %>%
  spread(pred, count)

if (!"positive" %in% colnames(sen_stage)) sen_stage$positive <- 0
if (!"negative" %in% colnames(sen_stage)) sen_stage$negative <- 0
sen_stage$positive <- ifelse(is.na(sen_stage$positive), 0, sen_stage$positive)
sen_stage$negative <- ifelse(is.na(sen_stage$negative), 0, sen_stage$negative)
sen_stage <- sen_stage %>%
  mutate(total = negative + positive,
         Sens  = round((positive / total) * 100, 1),
         CI    = ci_calc(total, positive, spf),
         Spef  = spf)

write.table(sen_stage, here("figures", "Figure_3DE_THEMIS-GC_Stage_sensitivity.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

plot_de <- sen_stage
plot_de <- plot_de[plot_de$Diagnosis != "HEALTHY", ]
plot_de <- plot_de[!grepl("Unknown", plot_de$group), ]
plot_de <- plot_de[!is.na(plot_de$group), ]

complete_comb <- do.call(rbind, lapply(unique(plot_de$Cohort), function(co) {
  expand.grid(Cohort = co, group = unique(plot_de$group), Diagnosis = unique(plot_de$Diagnosis))
}))
plot_de <- merge(complete_comb, plot_de, by = c("Cohort", "group", "Diagnosis"), all.x = TRUE)
plot_de <- plot_de %>%
  mutate(across(where(is.numeric), ~ ifelse(is.na(.), 0, .)),
         CI = ifelse(is.na(CI), "0-0", as.character(CI)))

plot_de$CI_low  <- as.numeric(gsub("-.*", "", plot_de$CI))
plot_de$CI_high <- as.numeric(gsub(".*-", "", plot_de$CI))
plot_de$Cohort  <- factor(plot_de$Cohort, levels = c("testing1", "MONITOR"))
plot_de <- plot_de %>%
  group_by(group, Diagnosis) %>%
  mutate(group_number = paste0(total, collapse = "|"))
plot_de$Stage_new <- paste0(plot_de$group, "\n(", plot_de$group_number, ")")
plot_de$label <- paste0(plot_de$Sens, "%")

stage_order <- c("Overall", "I", "II", "III", "IV", "Unknown")
plot_de$group <- factor(plot_de$group, levels = stage_order[stage_order %in% plot_de$group])
plot_de <- plot_de %>% arrange(group)
plot_de$Stage_new <- factor(plot_de$Stage_new, levels = unique(plot_de$Stage_new))

p_3de <- ggplot(plot_de, aes(x = Stage_new, y = Sens, ymin = CI_low, ymax = CI_high,
                               fill = Cohort, colour = Cohort)) +
  geom_errorbar(position = position_dodge(width = 0.7), width = 0.4, color = "black", size = 0.4) +
  geom_point(position = position_dodge(width = 0.7)) +
  scale_fill_manual(values  = colr_ext, labels = c("testing1", "MONITOR")) +
  scale_color_manual(values = colr_ext, labels = c("testing1", "MONITOR")) +
  theme_bw() +
  facet_wrap(vars(Cohort), ncol = 2, scales = "free") +
  theme(strip.background = element_blank(),
        strip.text = element_text(face = "bold", size = 12),
        plot.title = element_text(size = 15, face = "bold"),
        panel.background = element_rect(fill = "transparent"),
        panel.border = element_rect(fill = NA, colour = "black"),
        panel.grid.major = element_line(colour = NA),
        panel.grid.minor = element_line(colour = NA)) +
  scale_y_continuous(limits = c(-10, 110), labels = function(x) paste0(x, "%"),
                      breaks = seq(0, 100, by = 25)) +
  ylab("Sensitivity") + xlab("Tumour Stage") +
  ggtitle("THEMIS-GC sensitivity by Stage") +
  geom_text(aes(x = Stage_new, y = -5, label = label), size = 12 / .pt,
            position = position_dodge(width = 0.9))

pdf(here("figures", "Figure_3DE_THEMIS-GC_Stage_sensitivity.pdf"), width = 9, height = 3)
print(p_3de)
dev.off()

# =============================================================================
# Figures 3F-G - DCA for THEMIS-GC, Testing Dataset 1 (F) and MONITOR (G)
# =============================================================================
dca_marker <- "THEMIS-GC"
dca_col    <- pal_startrek("uniform", alpha = 0.9)(7)[4]

dca_cohorts <- list(
  "testing1" = subset(df, Cohort == "testing1"),
  "MONITOR"  = subset(df, Cohort == "MONITOR")
)

for (n in names(dca_cohorts)) {
  dca_data <- dca_cohorts[[n]][, c("GROUP", dca_marker)]
  dca_data$GROUP <- as.numeric(as.character(dca_data$GROUP))
  dca_data <- na.omit(dca_data)

  if (nrow(dca_data) < 10) { warning(paste("Cohort", n, "too small, skipping")); next }

  set.seed(123456)
  dca_res <- decision_curve(
    formula = as.formula(paste0("GROUP ~ `", dca_marker, "`")),
    data = dca_data,
    study.design = "cohort",
    thresholds = seq(0, 1, by = 0.01),
    bootstraps = 0,
    confidence.intervals = FALSE
  )
  dca_results <- list(dca_res)
  names(dca_results) <- dca_marker

  dca_table <- data.frame(
    threshold   = dca_res$derived.data$thresholds,
    net_benefit = dca_res$derived.data$NB,
    marker = dca_marker
  )
  write.table(dca_table, here("figures", paste(n, dca_marker, "dca_results.tsv", sep = ".")),
              row.names = FALSE, sep = "\t")

  pdf(here("figures", paste0("Figure_3", ifelse(n == "testing1", "F", "G"),
                              "_", n, "_dca.pdf")), width = 8, height = 6)
  plot_decision_curve(
    dca_results,
    curve.names = dca_marker,
    col = dca_col,
    confidence.intervals = FALSE,
    cost.benefit.axis = FALSE,
    legend.pos = "topright",
    xlab = "Threshold Probability",
    ylab = "Net Benefit",
    main = paste0("DCA - ", n)
  )
  dev.off()
}
