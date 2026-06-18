# =============================================================================
# Script:      FigureS9.R
# Description: Prognostic performance of THEMIS-GC using the median cutoff
#              (0.695), stage subgroup KM, and time-dependent ROC analysis.
#              (A) Kaplan-Meier DFS — high- vs low-risk (median cutoff 0.695)
#              (B) Kaplan-Meier OS  — high- vs low-risk (median cutoff 0.695)
#              (C) KM DFS within stage I, II, and III subgroups (median cutoff)
#              (D) KM OS  within stage I, II, and III subgroups (median cutoff)
#              (E) Multivariate Cox forest plot for OS
#              (F) Chi-Square contribution of each factor in the OS joint model
#              (G) Time-dependent ROC for DFS joint model (12- and 18-month)
#              (H) Time-dependent ROC for OS  joint model (12- and 18-month)
#
# This script sources Figure5.R to reuse km_plot(), km_plot_stages(),
# cox_multi(), looc_sta_fun(), pie_fun(), and the harmonised data frame (df).
# Running FigureS9.R will also regenerate Figure 5 outputs.
#
# --- Output -------------------------------------------------------------------
# Figures:  figures/Figure_S9A_DFS_KM.pdf
#           figures/Figure_S9B_OS_KM.pdf
#           figures/Figure_S9C_DFS_Stage_I_KM.pdf  (and Stage_II, Stage_III)
#           figures/Figure_S9D_OS_Stage_I_KM.pdf   (and Stage_II, Stage_III)
#           figures/Figure_S9E_OS_multivariate_forest.pdf
#           figures/Figure_S9F_OS_contribution_pie.pdf
#           figures/Figure_S9G_DFS_timeROC.pdf
#           figures/Figure_S9H_OS_timeROC.pdf
# Results:  results/Figure_S9E_OS_multivariate_cox.tsv
#           results/Figure_S9F_OS_contribution_pie.tsv
# =============================================================================

library(here)
library(timeROC)
source(here("Figure5.R"))

# Switch risk group to median cutoff for all S9 panels
df$Ensemble_risk <- factor(df$Ensemble_risk_median, levels = c("low_risk", "high_risk"))

# =============================================================================
# Time-dependent ROC function (joint Cox model → predicted risk → timeROC)
# =============================================================================
time_roc_plot <- function(data, time_col, status_col, clinical_vars, predictor,
                           time_points, outfile) {
  vars     <- c(predictor, clinical_vars)
  data_sub <- data[, c(time_col, status_col, vars)]
  data_sub[data_sub == ""] <- NA
  data_sub <- na.omit(data_sub)

  formula <- as.formula(paste0("Surv(", time_col, ",", status_col, ") ~ ",
                                 paste(vars, collapse = " + ")))
  model <- coxph(formula, data = data_sub)

  time_roc_res <- timeROC(
    T       = data_sub[[time_col]],
    delta   = data_sub[[status_col]],
    marker  = predict(model, type = "risk"),
    cause   = 1,
    weighting = "marginal",
    times   = time_points,
    ROC     = TRUE,
    iid     = TRUE
  )

  colr <- pal_startrek("uniform", alpha = 0.9)(7)[c(1, 4, 6)][seq_along(time_points)]

  pdf(outfile, width = 5, height = 5)
  par(mar = c(5, 4, 4, 6))
  plot(time_roc_res, time = time_points[1], col = colr[1], lwd = 2,
       xlab = "1 - Specificity", ylab = "Sensitivity", bty = "n")
  for (i in seq_along(time_points)[-1]) {
    plot(time_roc_res, time = time_points[i], col = colr[i], add = TRUE, lwd = 2)
  }
  legend("bottomright",
         legend = sapply(seq_along(time_points), function(i)
           paste0(time_points[i], " months (AUC: ", round(time_roc_res$AUC[i], 2), ")")),
         col = colr, lwd = 2, bty = "n")
  dev.off()
}

# =============================================================================
# Figure S9A - KM DFS (median cutoff)
# =============================================================================
km_plot(df, "DFS", "re_DFS", "Ensemble_risk",
         "Disease-free survival",
         here("figures", "Figure_S9A_DFS_KM.pdf"))

# =============================================================================
# Figure S9B - KM OS (median cutoff)
# =============================================================================
km_plot(df, "OS", "re_OS", "Ensemble_risk",
         "Overall survival",
         here("figures", "Figure_S9B_OS_KM.pdf"))

# =============================================================================
# Figure S9C - KM DFS within stage I, II, III subgroups (median cutoff)
# =============================================================================
km_plot_stages(df, "DFS", "re_DFS", "Ensemble_risk",
               "Disease-free survival", "Figure_S9C")

# =============================================================================
# Figure S9D - KM OS within stage I, II, III subgroups (median cutoff)
# =============================================================================
km_plot_stages(df, "OS", "re_OS", "Ensemble_risk",
               "Overall survival", "Figure_S9D")

# =============================================================================
# Figure S9E - Multivariate Cox forest plot for OS
# =============================================================================
cox_multi(df, "OS", "re_OS", clinical_vars, predictor, "Figure_S9E")

# =============================================================================
# Figure S9F - Chi-Square contribution pie for OS joint model
# =============================================================================
pie_fun(df, "OS", "re_OS", clinical_vars, predictor, "Figure_S9F")

# =============================================================================
# Figure S9G - Time-dependent ROC for DFS joint model
# =============================================================================
time_roc_plot(df, "DFS", "re_DFS", clinical_vars, predictor,
               time_points = c(12, 18),
               outfile = here("figures", "Figure_S9G_DFS_timeROC.pdf"))

# =============================================================================
# Figure S9H - Time-dependent ROC for OS joint model
# =============================================================================
time_roc_plot(df, "OS", "re_OS", clinical_vars, predictor,
               time_points = c(12, 18),
               outfile = here("figures", "Figure_S9H_OS_timeROC.pdf"))
