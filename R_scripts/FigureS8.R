# =============================================================================
# Script:      FigureS8.R
# Description: Comparison of a cfDNA-only model (THEMIS-GC) against a combined
#              cfDNA + protein biomarker model in the external Testing Dataset 1.
#              (A) ROC: cfDNA feature vs cfDNA+AFP+CA125+CA199+CEA — Testing 1
#              (B) Clinical landscape — Testing Dataset 1
#
# This script sources Figure4.R to reuse ROC_PLOT(), make_landscape(), and the
# merged data frame (df). Running FigureS8.R will also regenerate Figure 4 outputs.
#
# --- Output -------------------------------------------------------------------
# Figures:  figures/Figure_S8A_Testing1_ROC.pdf
#           figures/Figure_S8B_Testing1_landscape.pdf
# Results:  results/Figure_S8A_Testing1.delong.stat.tsv
# =============================================================================

library(here)
source(here("Figure4.R"))

# =============================================================================
# Figure S8A - ROC Testing Dataset 1
# =============================================================================
ROC_PLOT(df[df$Cohort == "testing1", ], "GC_label", roc_markers,
          "Testing Dataset 1", "Figure_S8A_Testing1")

# =============================================================================
# Figure S8B - Landscape Testing Dataset 1
# =============================================================================
dat_test <- subset(df, Cohort == "testing1")
dat_test$Stage <- factor(dat_test$Stage,
                           levels = c("BENIGN", "I", "II", "III", "IV", "Unknown"))
dat_test <- dat_test[order(dat_test$Stage, dat_test$`cfDNA+AFP+CA125+CA199+CEA`), ]
make_landscape(dat_test, here("figures", "Figure_S8B_Testing1_landscape.pdf"))
