# =============================================================================
# Script:      FigureS1.R
# Description: cfDNA feature distributions in Training and Validation datasets.
#              (A) Boxplot of median MFR per sample across HEALTHY, BENIGN, GC
#              (B) PCA of MFR profiles
#              (C) Boxplot of FSI Pearson correlation with Training baseline
#              (D) PCA of FSI profiles
#              (E) PCA of FEM profiles
#              Pairwise comparisons use Wilcoxon rank sum test with
#              Benjamini-Hochberg correction (*: p<=0.05, **: p<=0.01,
#              ***: p<=0.001, ns: non-significant).
#
# --- Required Input Data ------------------------------------------------------
# File 1:  data/Table_S3.tsv              (real data — Supplementary Table S3)
#          Columns: ID, Diagnosis, Laurens' subtypes, Tumor stage, Cohort,
#                   Gender, Age, MFR, FSI, FEM, CAFF, THEMIS-GC
#
# File 2:  data/MFR_matrix.demo.tsv       (demo placeholder)
#          Rows = genomic windows ("chrN:start-end"); columns = Window + sample IDs
#
# File 3:  data/FSI_matrix.demo.tsv       (demo placeholder)
#          Rows = FSI bins; columns = sample_id + sample IDs
#
# File 4:  data/FSI_bins.demo.tsv         (demo placeholder)
#          Rows = FSI bins; columns = seqnames, arm, start
#
# File 5:  data/FEM_matrix.demo.tsv       (demo placeholder)
#          Rows = FEM features; columns = Sample + sample IDs
#
# --- Output -------------------------------------------------------------------
# Figures:  figures/Figure_S1A_MFR_boxplot.svg
#           figures/Figure_S1B_MFR_PCA.svg
#           figures/Figure_S1C_FSI_boxplot.svg
#           figures/Figure_S1D_FSI_PCA.svg
#           figures/Figure_S1E_FEM_PCA.svg
# =============================================================================

library(RColorBrewer)
library(ggplot2)
library(ggsignif)
library(svglite)
library(rstatix)
library(ggpubr)
library(dplyr)
library(here)
library(systemfonts)

font_path <- here("fonts", "arial.ttf")
if (file.exists(font_path)) register_font(name = "Arial", plain = font_path)

# =============================================================================
# Load & harmonise sample metadata (Table S3)
# =============================================================================
sample.df <- read.table(here("data", "Table_S3.tsv"), header = TRUE, sep = "\t",
                          check.names = FALSE, stringsAsFactors = FALSE)
sample.df <- subset(sample.df, Cohort %in% c("Xijing-Training", "Xijing-Validation"))

sample.df$Sample <- sample.df$ID
sample.df$Stage  <- ifelse(sample.df$Diagnosis %in% c("HEALTHY", "BENIGN"),
                             sample.df$Diagnosis, sample.df$`Tumor stage`)
sample.df$Group  <- ifelse(sample.df$Diagnosis %in% c("HEALTHY", "BENIGN"),
                             sample.df$Diagnosis, sample.df$`Laurens' subtypes`)
sample.df$Cohort <- gsub("Xijing-", "", sample.df$Cohort)

anno.df <- data.frame(
  Stage = factor(sample.df$Stage,  levels = c("HEALTHY","BENIGN","I","II","III","IV","Unknown")),
  Group = factor(sample.df$Group,  levels = c("HEALTHY","BENIGN","DT","IT","MT","Unknown")),
  Cohort = sample.df$Cohort,
  row.names = sample.df$Sample,
  stringsAsFactors = FALSE
)

baseline.v <- sample.df$Sample[sample.df$Stage %in% c("HEALTHY","BENIGN") &
                                   sample.df$Cohort == "Training"]

# Shared colour scheme
group_colors <- c(
  HEALTHY = brewer.pal(9, "Blues")[8],
  BENIGN  = brewer.pal(9, "Blues")[7],
  GC      = brewer.pal(8, "Reds")[6]
)

# =============================================================================
# Shared PCA plot function
# =============================================================================
pca_plot <- function(pca_obj, anno, outfile) {
  pc.explain <- round(pca_obj$sdev^2 / sum(pca_obj$sdev^2) * 100, 2)

  pca.df <- data.frame(anno)
  pca.df$Group <- as.character(pca.df$Group)
  pca.df$PC1 <- pca_obj$x[match(rownames(anno), rownames(pca_obj$x)), "PC1"]
  pca.df$PC2 <- pca_obj$x[match(rownames(anno), rownames(pca_obj$x)), "PC2"]
  pca.df$Group[!pca.df$Group %in% c("HEALTHY", "BENIGN")] <- "GC"
  pca.df$Group <- factor(pca.df$Group, levels = c("HEALTHY", "BENIGN", "GC"))

  p <- ggplot(pca.df, aes(PC1, PC2, color = Group, shape = Group)) +
    geom_point() +
    labs(x = paste0("PC1 (", pc.explain[1], "%)"),
         y = paste0("PC2 (", pc.explain[2], "%)")) +
    theme(plot.title   = element_text(hjust = 0.5, face = "bold"),
          panel.background = element_rect(fill = "white", color = "black"),
          panel.border  = element_rect(fill = NA, color = "black", linewidth = 0.5),
          legend.position = "bottom") +
    scale_shape_manual(values = c(20, 23, 17)) +
    scale_color_manual(values = group_colors)

  svglite(outfile, width = 5, height = 5, fix_text_size = FALSE)
  plot(p)
  dev.off()
}

# =============================================================================
# Shared significance boxplot helper
# =============================================================================
sig_boxplot <- function(plot.df, y_col, ylab_str, outfile) {
  plot.df$Stage <- ifelse(plot.df$Stage %in% c("HEALTHY", "BENIGN"),
                           plot.df$Stage, "GC")
  plot.df$Stage <- factor(plot.df$Stage, levels = c("HEALTHY", "BENIGN", "GC"))

  kruskal_result <- kruskal.test(as.formula(paste(y_col, "~ Stage")), data = plot.df)
  pairwise_test  <- plot.df %>%
    pairwise_wilcox_test(as.formula(paste(y_col, "~ Stage")),
                          p.adjust.method = "BH") %>%
    add_xy_position(x = "Stage")

  p <- ggplot(plot.df, aes(x = Stage, y = .data[[y_col]], colour = Stage)) +
    geom_boxplot(outlier.shape = NA, width = 0.5) +
    geom_point(position = position_jitterdodge(jitter.width = 1),
               alpha = 0.8, size = 0.7) +
    stat_pvalue_manual(pairwise_test, label = "p.adj.signif",
                        tip.length = 0.005, step.increase = 0.09) +
    labs(x = "Stage", y = ylab_str,
         subtitle = paste0("Kruskal-Wallis p = ",
                            signif(kruskal_result$p.value, 3))) +
    scale_color_manual(values = rep("black", 3)) +
    theme_bw() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          plot.title = element_text(hjust = 0.5, face = "bold"),
          legend.position = "none")

  svglite(outfile, width = 5, height = 5, fix_text_size = FALSE)
  plot(p)
  dev.off()
}

# =============================================================================
# Panels A & B — MFR boxplot and PCA
# =============================================================================
mfr.df <- read.table(here("data", "MFR_matrix.demo.tsv"),
                       header = TRUE, check.names = FALSE)
mfr.df  <- mfr.df[, c("Window", rownames(anno.df))]

chr.v   <- sapply(strsplit(as.character(mfr.df$Window), ":"), `[`, 1)
start.v <- as.numeric(sapply(strsplit(
                         sapply(strsplit(as.character(mfr.df$Window), ":"), `[`, 2),
                         "-"), `[`, 1))
chr.v   <- factor(chr.v, levels = paste0("chr", 1:22))
idx     <- order(chr.v, start.v)
mfr.m   <- as.matrix(mfr.df[idx, -1])

# Panel A — median MFR per sample
mfr.v      <- apply(mfr.m, 2, median)
mfr_plot.df <- data.frame(anno.df,
                            mfr = mfr.v[match(rownames(anno.df), names(mfr.v))],
                            stringsAsFactors = FALSE)
mfr_plot.df$Stage <- as.character(mfr_plot.df$Stage)
sig_boxplot(mfr_plot.df, "mfr", "Median MFR",
             here("figures", "Figure_S1A_MFR_boxplot.svg"))

# Panel B — MFR PCA
mfr.pca <- prcomp(t(mfr.m), center = TRUE, scale. = FALSE)
pca_plot(mfr.pca, anno.df, here("figures", "Figure_S1B_MFR_PCA.svg"))

# =============================================================================
# Panels C & D — FSI correlation boxplot and PCA
# =============================================================================
fsi_cv_cut <- 0.5   # fraction of features to retain (lowest CV)

fsi.df <- read.table(here("data", "FSI_matrix.demo.tsv"),
                       header = TRUE, check.names = FALSE)
fsi.df  <- fsi.df[, c("sample_id", rownames(anno.df))]

fsi.m2  <- as.matrix(fsi.df[, -1])
rownames(fsi.m2) <- fsi.df$sample_id

cv_fun <- function(mat) apply(mat, 1, function(x) sd(x) / mean(x))
cv.v   <- cv_fun(fsi.m2)
n_keep <- round(nrow(fsi.m2) * fsi_cv_cut)
keep   <- names(sort(cv.v))[1:n_keep]

bin.df  <- read.table(here("data", "FSI_bins.demo.tsv"),
                        header = TRUE, check.names = FALSE)
bin.df  <- bin.df[as.numeric(sub("bin", "", keep)), ]
idx     <- order(as.numeric(sub("chr", "", bin.df$seqnames)),
                  bin.df$arm, bin.df$start)
fsi.m   <- fsi.m2[keep[idx], ]

# Panel C — Pearson correlation with Training baseline mean
baseline.fsi   <- fsi.m[, colnames(fsi.m) %in% baseline.v]
baseline.mean  <- rowMeans(baseline.fsi)
cor.v          <- apply(fsi.m, 2, function(x) cor(baseline.mean, x, method = "pearson"))

temp_anno.df   <- anno.df
fsi_plot.df    <- data.frame(temp_anno.df,
                               Correlation = cor.v[match(rownames(temp_anno.df), names(cor.v))],
                               stringsAsFactors = FALSE)
fsi_plot.df$Stage <- as.character(fsi_plot.df$Stage)
sig_boxplot(fsi_plot.df, "Correlation", "FSI Pearson correlation",
             here("figures", "Figure_S1C_FSI_boxplot.svg"))

# Panel D — FSI PCA
fsi.pca <- prcomp(t(fsi.m), center = FALSE, scale. = FALSE)
pca_plot(fsi.pca, anno.df, here("figures", "Figure_S1D_FSI_PCA.svg"))

# =============================================================================
# Panel E — FEM PCA
# =============================================================================
FEM.df <- read.table(here("data", "FEM_matrix.demo.tsv"),
                       header = TRUE, check.names = FALSE)
FEM.df  <- FEM.df[, c("Sample", rownames(anno.df))]
FEM.m   <- as.matrix(FEM.df[, -1])

FEM.pca <- prcomp(t(FEM.m), center = TRUE, scale. = TRUE)
pca_plot(FEM.pca, anno.df, here("figures", "Figure_S1E_FEM_PCA.svg"))
