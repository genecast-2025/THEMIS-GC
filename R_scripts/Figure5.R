# =============================================================================
# Script:      Figure5.R
# Description: Prognostic performance of THEMIS-GC risk group (optimal cutoff)
#              in 125 GC patients with available follow-up data.
#              (A) Kaplan-Meier DFS — high- vs low-risk (optimal cutoff 0.762)
#              (B) Kaplan-Meier OS  — high- vs low-risk (optimal cutoff 0.762)
#              (C) KM DFS within stage I, II, and III subgroups
#              (D) KM OS  within stage I, II, and III subgroups
#              (E) Multivariate Cox forest plot for DFS
#                  (THEMIS-GC risk + Age + Gender + Tumor stage + Laurens subtype)
#              (F) Chi-Square contribution of each factor in the DFS joint model
#
# --- Required Input Data ------------------------------------------------------
# File:   data/Table_S8.tsv
# Source: Supplementary Table S8
# Format: Tab-separated, columns:
#         ID, Gender, Age, Laurens' subtypes, Tumor stage, OS months,
#         OS status, DFS months, DFS status, Cohort, THEMIS-GC,
#         Risk group by THSMIS-GC (optimal cutoff at 0.762),
#         Risk group by THSMIS-GC (median cutoff at 0.695)
#
# --- Output -------------------------------------------------------------------
# Figures:  figures/Figure_5A_DFS_KM.pdf
#           figures/Figure_5B_OS_KM.pdf
#           figures/Figure_5C_DFS_Stage_I_KM.pdf   (and Stage_II, Stage_III)
#           figures/Figure_5D_OS_Stage_I_KM.pdf    (and Stage_II, Stage_III)
#           figures/Figure_5E_DFS_multivariate_forest.pdf
#           figures/Figure_5F_DFS_contribution_pie.pdf
# Results:  results/Figure_5E_DFS_multivariate_cox.tsv
#           results/Figure_5F_DFS_contribution_pie.tsv
# =============================================================================
library(readxl)
library(survival)
library(survminer)
library(rms)
library(ggforce)
library(ggplot2)
library(ggsci)
library(dplyr)
library(here)
library(systemfonts)
library(litedown)
font_path <- here("fonts", "arial.ttf")
if (file.exists(font_path)) register_font(name = "Arial", plain = font_path)
project_dir=here("figures")
if (!dir.exists(project_dir)) {
  dir.create(project_dir)
}

# =============================================================================
# Load & harmonise Table S8
# =============================================================================
#df <- read.table(here("data", "Table_S8.tsv"), header = TRUE, sep = "\t",
#                  check.names = FALSE, stringsAsFactors = FALSE)

df <- read_excel(
  here("data", "Table_S1-S9.xlsx"),
  sheet = "Table S8",
  skip = 1
)

df$`THEMIS-GC` = df$`OMNI-GC`
df$`Risk group by THSMIS-GC (optimal cutoff at 0.762 )` = df$`Risk group by OMNI-GC (cut-off at 0.762 )`

df$PID             <- df$ID
df$OS              <- df$`OS months`
df$re_OS           <- df$`OS status`
df$DFS             <- df$`DFS months`
df$re_DFS          <- df$`DFS status`
df$Ensemble        <- df$`THEMIS-GC`
df$Laurens_subtypes  <- df$`Laurens' subtypes`
df$Laurens_subtypes2 <- ifelse(df$`Laurens' subtypes` == "IT", "IT", "MT_DT")
df$Tumor_stages      <- ifelse(df$`Tumor stage` == "III", "III", "I_II")
df$Age               <- factor(df$Age, levels = c("≤65", ">65"))

df$Stage <- df$`Tumor stage`

df$Ensemble_risk_optimal <- gsub("-", "_",
  df$`Risk group by THSMIS-GC (optimal cutoff at 0.762 )`)
#df$Ensemble_risk_median  <- gsub("-", "_",
#  df$`Risk group by THSMIS-GC (median cutoff at 0.695 )`)

df$Ensemble_risk <- factor(df$Ensemble_risk_optimal, levels = c("low_risk", "high_risk"))

df$Ensemble_risk_vs_Stage = paste(df$Ensemble_risk,df$Stage,sep="_")
clinical_vars <- c("Age", "Gender", "Tumor_stages", "Laurens_subtypes2")
predictor     <- "Ensemble_risk"

# =============================================================================
# Kaplan-Meier function (2-group: high_risk vs low_risk)
# =============================================================================
km_plot <- function(data, time_col, status_col, group_col, ylab_str, outfile) {
  d <- data[, c("PID", group_col, time_col, status_col)]
  d <- na.omit(d)
  #d[[group_col]] <- factor(d[[group_col]], levels = c("low_risk", "high_risk"))
# 自动将分组列转为因子，按出现顺序或字母顺序排序
d[[group_col]] <- as.factor(d[[group_col]])
  #fit  <- survfit(as.formula(paste0("Surv(", time_col, ",", status_col, ") ~ ", group_col)),
  #                 data = d)
  #fit3 <- coxph(as.formula(paste0("Surv(", time_col, ",", status_col, ") ~ ", group_col)),
  #               data = d)
print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
formula_surv <- as.formula(
    paste0(
        "Surv(",
        time_col,
        ",",
        status_col,
        ") ~ ",
        group_col
    )
)
print(formula_surv)
fit  <- survfit(formula_surv, data=d)
svv_diy_labs = c()
    for(i in rownames(as.matrix(fit$strata))){
      i = strsplit(i,'=')[[1]][2]
      svv_diy_labs = c(svv_diy_labs,i)
    }

print("#######")
  # 修复 'symbol' is not subsettable 错误：将 call 中记录的符号替换为实际的公式对象
  fit$call$formula <- formula_surv 

fit3 <- coxph(formula_surv, data=d)
  s    <- summary(fit3)
print(s)
  pval <- round(s$sctest["pvalue"], 4)
  hr   <- round(s$conf.int[, "exp(coef)"], 2)
  ci   <- paste0(round(s$conf.int[, "lower .95"], 2), "-", round(s$conf.int[, "upper .95"], 2))
  pv_label <- paste0("P = ", pval, "\nHR = ", hr, "\n95%CI: ", ci)

  strata_labs <- gsub(paste0(group_col, "="), "", names(fit$strata))
if(length(svv_diy_labs)==9){
     cololr_list_diy = c("#6495ED","#0000FF","#F08080","#FF0000", "#4DBBD5FF",  "#7E6148FF", "black","grey40","#663366")
     linetype = c(2,2,2,2,2,1,1,1,1)

    }
    if(length(svv_diy_labs)==8){
      cololr_list_diy = c("#6495ED","#0000FF","#F08080","#FF0000", "#4DBBD5FF",  "#7E6148FF", "black","grey40")
      linetype = c(2,2,2,2,1,1,1,1)
    }
   if(length(svv_diy_labs)==7){
      cololr_list_diy = c("#6495ED","#0000FF","#F08080","#FF0000", "#4DBBD5FF",  "#7E6148FF", "black")
      linetype = c(2,2,2,1,1,1,1)
    }
    if(length(svv_diy_labs)==6){
      cololr_list_diy = c("#6495ED","#0000FF","#F08080","#FF0000", "#4DBBD5FF",  "#7E6148FF", "#8491B4FF")
      cololr_list_diy = c("#F08080","#FF0000", "darkred","#4DBBD5FF","#0000FF","cyan")
  cololr_list_diy = c("#F08080","darkred","#FF0000", "#4DBBD5FF","cyan","#0000FF")
cololr_list_diy = c("#F08080","darkred","#FF0000", "olivedrab","cyan","#0000FF")
      linetype = c(1,1,1,1,1,1)
    }
  if(length(svv_diy_labs)==5){
      cololr_list_diy = c("#6495ED","#0000FF","#F08080","#FF0000", "#4DBBD5FF",  "#7E6148FF", "#8491B4FF")
      linetype = c(1,1,1,1,1)
    }
if(length(svv_diy_labs)==4){
  cololr_list_diy  = c("#374E5599","#374E55FF","#DF8F4499","#DF8F44FF")
    cololr_list_diy  = c("#6495ED","#0000FF","#F08080","#FF0000")
  linetype=c(1,1,1,1)
}
if(length(svv_diy_labs)==3){
  cololr_list_diy = c("#374E55FF","#DF8F44E5","#B24745FF")
  cololr_list_diy = c("#0000FF","#7E6148FF","#FF0000")
  linetype = c(1,1,1)
}
if(length(svv_diy_labs)==2){
  cololr_list_diy = c("#0000FF","#FF0000")
 # cololr_list_diy = c("#374E55FF","#DF8F44FF")
  linetype = c(1,1)
}

print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
  k2 <- ggsurvplot(
    fit, data = d,
    censor.shape = "|", censor.size = 1.5,
    break.time.by = 6,
    linetype = linetype,
    legend.title = "", legend.labs = strata_labs,
    pval = pv_label, pval.size = 3,
    legend = "none",
    palette = cololr_list_diy,#c("#0000FF", "#FF0000"),
    ylab = ylab_str, xlab = "Time since surgery (months)",
    surv.scale = "percent",
surv.median.line = "hv",
    risk.table = TRUE, risk.table.fontsize = 2, risk.table.col = "strata",
    ggtheme = theme_classic(), risk.table.height = 0.2, lwd = 0.5
  )
  #k2$table <- k2$table +
  #  theme(axis.title.y = element_blank(), axis.title.x = element_blank(),
   #       text = element_text(size = 4), plot.margin = margin(0, 1, 1, 1, "pt"),
    #      axis.text.y = element_text(size = 4), axis.text.x = element_text(size = 4))

  pdf(outfile, onefile = FALSE, width = 6, height = 5)
  print(k2)
  dev.off()
}

# =============================================================================
# Stage-stratified KM: loops over I / II / III subsets
# =============================================================================
km_plot_stages <- function(data, time_col, status_col, group_col, ylab_str, fig_prefix) {
  for (s in c("I", "II", "III")) {
    sub <- subset(data, Stage == s)
    if (nrow(sub) < 5) next
    km_plot(sub, time_col, status_col, group_col, ylab_str,
            here("figures", paste0(fig_prefix, "_Stage_", s, "_KM.pdf")))
  }
}

# =============================================================================
# Multivariate Cox forest plot
# =============================================================================
cox_multi <- function(data, time_col, status_col, clinical_vars, predictor, outprefix) {
  vars     <- c(predictor, clinical_vars)
  data_sub <- data %>% select(all_of(c(time_col, status_col, vars)))
  data_sub[data_sub == ""] <- NA
  data_sub <- na.omit(data_sub)

  formula <- as.formula(paste0("Surv(", time_col, ",", status_col, ") ~ ."))
  fit <- coxph(formula, data = data_sub)

  sumary <- summary(fit)
  res_df <- data.frame(
    term     = rownames(sumary$coefficients),
    p.value  = ifelse(sumary$coefficients[, "Pr(>|z|)"] < 0.001, "<0.001",
                       as.character(round(sumary$coefficients[, "Pr(>|z|)"], 3))),
    HR       = round(sumary$conf.int[, "exp(coef)"], 2),
    CI_lower = round(sumary$conf.int[, "lower .95"], 2),
    CI_upper = round(sumary$conf.int[, "upper .95"], 2),
    check.names = FALSE
  )
  res_df$`HR(95%CI)` <- paste0(res_df$HR, " (", res_df$CI_lower, "-", res_df$CI_upper, ")")
  write.table(res_df, here("figures", paste0(outprefix, "_multivariate_cox.tsv")),
              sep = "\t", row.names = FALSE, quote = FALSE)

  pdf(here("figures", paste0(outprefix, "_forest.pdf")), width = 8, height = 6,onefile=F)
  print(ggforest(fit, data = data_sub))
  dev.off()
}

# =============================================================================
# Chi-Square contribution pie chart (uses rms::cph + ggforce::geom_arc_bar)
# =============================================================================
looc_sta_fun <- function(data, r_len) {
  data$start <- NA; data$end <- NA; pre_list <- list()
  for (i in data$factorname) {
    data[data$factorname == i, ]$start <- sum(data[data$factorname %in% pre_list, ]$Chi_propotion) * 360 / 100
    pre_list <- c(unlist(pre_list), i)
    data[data$factorname == i, ]$end <- sum(data[data$factorname %in% pre_list, ]$Chi_propotion) * 360 / 100
  }
  data$mark  <- sprintf("%s\n(%.2f%%)", data$factorname, round(data$Chi_propotion * 100, 2))
  data$alpha <- data$end - ((data$end - data$start) / 2)
  data$y     <- cos(data$alpha * 100 * pi / 180) * r_len
  data$x     <- sin(data$alpha * 100 * pi / 180) * r_len
  data
}

pie_fun <- function(data, time_col, status_col, clinical_vars, predictor, outprefix) {
  vars     <- c(predictor, clinical_vars)
  data_sub <- data %>% select(all_of(c(time_col, status_col, vars)))
  data_sub[data_sub == ""] <- NA
  data_sub <- na.omit(data_sub)
  #dd <- datadist(data_sub); options(datadist = "dd")
	  dd <- datadist(data_sub)
	  assign("dd", dd, envir = .GlobalEnv)
	  options(datadist = "dd")
  formula <- as.formula(paste0("Surv(", time_col, ",", status_col, ") ~ ",
                                 paste(vars, collapse = " + ")))
  print(formula)
  fit_rms  <- cph(formula, data = data_sub)
  anova_df <- as.data.frame(anova(fit_rms))
  anova_sub <- anova_df[vars[vars %in% rownames(anova_df)], , drop = FALSE]
  anova_sub$Chi_propotion <- anova_sub$`Chi-Square` / sum(anova_sub$`Chi-Square`, na.rm = TRUE)
  anova_sub$factorname    <- rownames(anova_sub)
  anova_sub <- anova_sub[order(anova_sub$Chi_propotion, decreasing = TRUE), ]
  anova_sub$factorname <- factor(anova_sub$factorname, levels = anova_sub$factorname)

  write.table(anova_sub[, c("factorname", "Chi-Square", "Chi_propotion")],
              here("figures", paste0(outprefix, "_contribution_pie.tsv")),
              sep = "\t", row.names = FALSE, quote = FALSE)

  plot_data <- looc_sta_fun(anova_sub, 2.2)

  p <- ggplot() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          aspect.ratio = 1, axis.ticks = element_blank(),
          legend.position = "bottom", legend.key.height = unit(1, "pt"),
          legend.key.width = unit(2, "pt"), axis.text = element_blank(),
          panel.border = element_blank(), panel.background = element_blank()) +
    xlab("") + ylab("") +
    labs(title = "Proportion of each factor") +
    guides(fill = guide_legend(nrow = 3, byrow = TRUE)) +
    scale_fill_manual(values = pal_startrek("uniform", alpha = 0.9)(7)) +
    scale_x_continuous(limits = c(-2.5, 2.5), breaks = NULL) +
    scale_y_continuous(limits = c(-2.5, 2.5), breaks = NULL) +
    geom_arc_bar(data = plot_data, stat = "pie",
                 aes(x0 = 0, y0 = 0, r0 = 0, r = 1.6,
                     amount = Chi_propotion, fill = factorname),
                 color = "white") +
    annotate("text", x = plot_data$x, y = plot_data$y, label = plot_data$mark,
             hjust = 0.5, size = 10 / .pt)

  pdf(here("figures", paste0(outprefix, "_pie.pdf")), width = 6, height = 6)
  print(p)
  dev.off()
}

main <- function() {
# =============================================================================
# Figure 5A - KM DFS (optimal cutoff)
# =============================================================================
km_plot(df, "DFS", "re_DFS", "Ensemble_risk",
         "Disease-free survival",
         here("figures", "Figure_5A_DFS_KM.pdf"))

# =============================================================================
# Figure 5B - KM OS (optimal cutoff)
# =============================================================================
km_plot(df, "OS", "re_OS", "Ensemble_risk",
         "Overall survival",
         here("figures", "Figure_5B_OS_KM.pdf"))

# =============================================================================
# Figure 5C - KM DFS within stage I, II, III subgroups (optimal cutoff)
# =============================================================================
km_plot_stages(df, "DFS", "re_DFS", "Ensemble_risk",
               "Disease-free survival", "Figure_5C")
#Ensemble_risk_vs_Stage
print(table(df$Ensemble_risk_vs_Stage))
km_plot(df, "DFS", "re_DFS", "Ensemble_risk_vs_Stage",
         "Disease-free survival",
         here("figures", "Figure_5C_DFS_KM.pdf"))


# =============================================================================
# Figure 5D - KM OS within stage I, II, III subgroups (optimal cutoff)
# =============================================================================
km_plot_stages(df, "OS", "re_OS", "Ensemble_risk",
               "Overall survival", "Figure_5D")
km_plot(df, "OS", "re_OS", "Ensemble_risk_vs_Stage",
         "Overall survival",
         here("figures", "Figure_5D_OS_KM.pdf"))


# =============================================================================
# Figure 5E - Multivariate Cox forest plot for DFS
# =============================================================================
cox_multi(df, "DFS", "re_DFS", clinical_vars, predictor, "Figure_5E")

# =============================================================================
# Figure 5F - Chi-Square contribution pie for DFS joint model
# =============================================================================
pie_fun(df, "DFS", "re_DFS", clinical_vars, predictor, "Figure_5F")
      }

if (sys.nframe() == 0) {
          main()
        }

