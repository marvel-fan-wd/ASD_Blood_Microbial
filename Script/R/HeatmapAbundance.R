library(ComplexHeatmap);library(circlize);library(grid)

cor <- read.csv("E:/2025.8-2026.7/iScience_revision/Abundance_phenotype_association_beta.csv", header = TRUE, row.names = 1)
pval <- read.csv("E:/2025.8-2026.7/iScience_revision/Abundance_phenotype_association_pval.csv", header = TRUE, row.names = 1)
type <- read.csv("E:/2025.8-2026.7/iScience_revision/pheno_analysis/pheno.csv", header = TRUE)
unique(type$Type)

rownames(cor) <- gsub("_", " ", rownames(cor))
rownames(pval) <- gsub("_", " ", rownames(pval))

pval <- pval[, match(colnames(cor), colnames(pval))]
colnames(cor) <- type$Pheno[match(colnames(cor), type$pheno)]
colnames(pval) <- type$Pheno[match(colnames(pval), type$pheno)]

cor <- cor[, order(match(colnames(cor), type$Pheno))]
pval <- pval[, order(match(colnames(pval), type$Pheno))]

annotation_col <- data.frame(Type = type$Type[match(colnames(cor), type$Pheno)])
rownames(annotation_col) <- colnames(cor)

significance_matrix <- ifelse(is.na(pval), "", 
                              ifelse(pval < 0.001, "***",
                                     ifelse(pval < 0.01, "**",
                                            ifelse(pval < 0.05, "*", ""))))
rownames(significance_matrix) <- rownames(pval)
colnames(significance_matrix) <- colnames(pval)

phenotype_order <- c("Host", "Pregnancy", "Diet", "Disease/Medication",
                     "Adaptive", "Cognitive", "Behavior", "Social/Communication")
annotation_col$Type <- factor(annotation_col$Type, levels = phenotype_order)

col_fun <- colorRamp2(c(-0.1, 0, 0.1), c("#0f86a9", "white", "#FC8452"))

Heatmap(
  cor,
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
