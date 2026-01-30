####------------------------------------------------------- Figure S5 -------------------------------------------------------
rm(list=ls())
library(ComplexHeatmap)
library(circlize)
library(grid)

# Prevalence
cor <- read.csv("E:/2025.8-2026.7/iScience_revision/14_pheno_analysis/analysis_prevalence/cor_all_final.csv", header = TRUE, row.names = 1)
pval <- read.csv("E:/2025.8-2026.7/iScience_revision/14_pheno_analysis/analysis_prevalence/qval_all_final.csv", header = TRUE, row.names = 1)
type <- read.csv("E:/2025.8-2026.7/iScience_revision/14_pheno_analysis/pheno_final_cate_46.csv", header = TRUE)
unique(type$Type)

rownames(cor) <- gsub("_", " ", rownames(cor))
rownames(pval) <- gsub("_", " ", rownames(pval))
species_order <- gsub("_", " ", c(
  "Pseudoalteromonas_sp._3J6",
  "Prescottella_equi",
  "Bdellovibrio_bacteriovorus",
  "Parabacteroides_distasonis",
  "Enterobacter_hormaechei",
  "Xylella_fastidiosa",
  "Wolbachia_pipientis",
  "Yersinia_enterocolitica",
  "Elizabethkingia_anophelis",
  "Treponema_pallidum",
  "Klebsiella_michiganensis",
  "Ochrobactrum_quorumnocens"
))

pval <- pval[, match(colnames(cor), colnames(pval))]

colnames(cor) <- type$Pheno[match(colnames(cor), type$pheno)]
colnames(pval) <- type$Pheno[match(colnames(pval), type$pheno)]

if (!all(colnames(cor) == colnames(pval))) {
  stop("Column names of cor and pval are not consistent after renaming.")
}

cor <- cor[, order(match(colnames(cor), type$Pheno))]
pval <- pval[, order(match(colnames(pval), type$Pheno))]

annotation_col <- data.frame(Type = type$Type[match(colnames(cor), type$Pheno)])
rownames(annotation_col) <- colnames(cor)

cor <- cor[species_order, , drop = FALSE]
pval <- pval[species_order, , drop = FALSE]

significance_matrix <- ifelse(is.na(pval), "", 
                              ifelse(pval < 0.001, "***",
                                     ifelse(pval < 0.01, "**",
                                            ifelse(pval < 0.05, "*", ""))))
rownames(significance_matrix) <- rownames(pval)
colnames(significance_matrix) <- colnames(pval)

phenotype_order <- c("Host", "Pregnancy", "Diet", "Disease/Medication",
                     "Adaptive", "Cognitive", "Behavior", "Social/Communication")
annotation_col$Type <- factor(annotation_col$Type, levels = phenotype_order)

if (!all(colnames(cor) %in% rownames(annotation_col))) {
  stop("Error: Some columns in 'cor' do not have matching rows in 'annotation_col'.")
}

cor1 <- cor

quantile(as.matrix(cor1), 
         probs = c(0, 0.01, 0.05, 0.25, 0.5, 0.75, 0.95, 0.99, 1),
         na.rm = TRUE)

# 0%          1%          5%         25%         50%         75%         95%         99%        100% 
# -14.1251956 -12.3240436  -6.6892446  -0.6087086  -0.1151275   0.1179545   0.9063407   3.6434932   7.5385913
col_fun <- colorRamp2(c(-0.6, 0, 0.6), c("#0f86a9", "white", "#FC8452"))

ht1 <- Heatmap(
  cor1,
  name = "Correlation",
  col = col_fun,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_names = TRUE,
  show_column_names = FALSE,
  column_split = annotation_col$Type,
  border = TRUE,
  rect_gp = gpar(col = "white", lwd = 1.5),
  border_gp = gpar(col = "#0f86a9", lty = 2, lwd = 1.2),
  cell_fun = function(j, i, x, y, width, height, fill) {
    grid.text(significance_matrix[i, j], x, y, gp = gpar(fontsize = 10))
  },
  row_names_gp = gpar(fontsize = 13, fontface = "italic"),
  column_names_rot = 45,
  column_names_gp = gpar(fontsize = 13),
  column_title_gp = gpar(fontsize = 14, fontface = "bold"),
  heatmap_legend_param = list(
    title = "Correlation",
    title_gp = gpar(fontsize = 14, fontface = "bold"),
    labels_gp = gpar(fontsize = 12),
    at = c(-0.6, 0, 0.6),
    labels = c("-0.6", "0", "0.6")
  )
)

# Abundance
cor <- read.csv("E:/2025.8-2026.7/iScience_revision/14_pheno_analysis/analysis_abundance/cor_all_final.csv", header = TRUE, row.names = 1)
pval <- read.csv("E:/2025.8-2026.7/iScience_revision/14_pheno_analysis/analysis_abundance/qval_all_final.csv", header = TRUE, row.names = 1)
type <- read.csv("E:/2025.8-2026.7/iScience_revision/14_pheno_analysis/pheno_final_cate_46.csv", header = TRUE)
unique(type$Type)

rownames(cor) <- gsub("_", " ", rownames(cor))
rownames(pval) <- gsub("_", " ", rownames(pval))
species_order <- gsub("_", " ", c(
  "Pseudoalteromonas_sp._3J6",
  "Prescottella_equi",
  "Bdellovibrio_bacteriovorus",
  "Parabacteroides_distasonis",
  "Enterobacter_hormaechei",
  "Xylella_fastidiosa",
  "Wolbachia_pipientis",
  "Yersinia_enterocolitica",
  "Elizabethkingia_anophelis",
  "Treponema_pallidum",
  "Klebsiella_michiganensis",
  "Ochrobactrum_quorumnocens"
))

pval <- pval[, match(colnames(cor), colnames(pval))]

colnames(cor) <- type$Pheno[match(colnames(cor), type$pheno)]
colnames(pval) <- type$Pheno[match(colnames(pval), type$pheno)]

if (!all(colnames(cor) == colnames(pval))) {
  stop("Column names of cor and pval are not consistent after renaming.")
}

cor <- cor[, order(match(colnames(cor), type$Pheno))]
pval <- pval[, order(match(colnames(pval), type$Pheno))]

annotation_col <- data.frame(Type = type$Type[match(colnames(cor), type$Pheno)])
rownames(annotation_col) <- colnames(cor)

cor <- cor[species_order, , drop = FALSE]
pval <- pval[species_order, , drop = FALSE]

significance_matrix <- ifelse(is.na(pval), "", 
                              ifelse(pval < 0.001, "***",
                                     ifelse(pval < 0.01, "**",
                                            ifelse(pval < 0.05, "*", ""))))
rownames(significance_matrix) <- rownames(pval)
colnames(significance_matrix) <- colnames(pval)

phenotype_order <- c("Host", "Pregnancy", "Diet", "Disease/Medication",
                     "Adaptive", "Cognitive", "Behavior", "Social/Communication")
annotation_col$Type <- factor(annotation_col$Type, levels = phenotype_order)

if (!all(colnames(cor) %in% rownames(annotation_col))) {
  stop("Error: Some columns in 'cor' do not have matching rows in 'annotation_col'.")
}

cor2 <- cor

quantile(as.matrix(cor2), 
         probs = c(0, 0.01, 0.05, 0.25, 0.5, 0.75, 0.95, 0.99, 1),
         na.rm = TRUE)
         # 0%          1%          5%         25%         50%         75%         95%         99%        100% 
# -3.05962942 -1.94168714 -1.11145045 -0.10052026 -0.01568966  0.02486201  0.15437841  0.42579165  0.77143278 

col_fun <- colorRamp2(c(-0.1, 0, 0.1), c("#0f86a9", "white", "#FC8452"))

ht2 <- Heatmap(
  cor2,
  name = "Correlation",
  col = col_fun,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_names = TRUE,
  show_column_names = TRUE,
  column_split = annotation_col$Type,
  border = TRUE,
  rect_gp = gpar(col = "white", lwd = 1.5),
  border_gp = gpar(col = "#0f86a9", lty = 2, lwd = 1.2),
  cell_fun = function(j, i, x, y, width, height, fill) {
    grid.text(significance_matrix[i, j], x, y, gp = gpar(fontsize = 10))
  },
  row_names_gp = gpar(fontsize = 13, fontface = "italic"),
  column_title = NULL,
  column_names_rot = 45,
  column_names_gp = gpar(fontsize = 13),
  column_title_gp = gpar(fontsize = 14, fontface = "bold"),
  heatmap_legend_param = list(
    title = "Correlation",
    title_gp = gpar(fontsize = 14, fontface = "bold"),
    labels_gp = gpar(fontsize = 12),
    at = c(-0.1, 0, 0.1),
    labels = c("-0.1", "0", "0.1")
  )
)

ht3 <- ht1 %v% ht2

pdf("E:/2025.8-2026.7/iScience_revision/Figures/Figure S5.pdf", width = 18, height = 10)
draw(ht3)
dev.off()