####------------------------------------------------------- (14).phenotypic association analysis -------------------------------------------------------
rm(list=ls())
library(cowplot)
library(logistf)
library(ComplexHeatmap)
library(circlize)
library(grid)

# Relative abundance--------------------------------------------------------------------
## 1.Continuous traits
# Prepare abandunce matrix
ma=read.table("E:/2025.8-2026.7/iScience_revision/13_differential_analysis/3_DA/Maaslin2/all_results.tsv", header = TRUE, sep = "\t")
ma_sig=ma[ma$qval<0.05 & ma$metadata == "Role",]

preva=read.csv("E:/2025.8-2026.7/iScience_revision/13_differential_analysis/3_DA/Prevalence/prevalence_diff_result_CLR.csv",header=T,row.names = 1)
colnames(preva)
preva_sig=preva[preva$q_poss<0.05,]

diff <- union(ma_sig$feature, preva_sig$Species)
diff # 12
# [1] "Pseudoalteromonas_sp._3J6"  "Prescottella_equi"          "Parabacteroides_distasonis" "Klebsiella_michiganensis"  
# [5] "Xylella_fastidiosa"         "Enterobacter_hormaechei"    "Wolbachia_pipientis"        "Treponema_pallidum"        
# [9] "Yersinia_enterocolitica"    "Elizabethkingia_anophelis"  "Bdellovibrio_bacteriovorus" "Ochrobactrum_quorumnocens" 

otu_raw <- read.csv("E:/2025.8-2026.7/iScience_revision/11_depth_normalize/Species_relative_abundance_100.csv")
cate=read.csv("E:/2025.8-2026.7/iScience_revision/1_rawdata/SSC_kraken_cate_16344.csv",row.names = 1)
colnames(cate)[colnames(cate) == "V1"] <- "X"
cate=cate[,c("Species","X")]
otu=merge(otu_raw,cate,by="X")
row.names(otu)=otu$Species
otu <- otu[, !(colnames(otu) %in% c("Species", "X"))]

diff_df=otu[rownames(otu) %in% diff, ]
t_diff_df=as.data.frame(t(diff_df))
t_diff_df$SPID=row.names(t_diff_df)

# Prepare phenotype matrix
sa_pheno=read.csv("E:/2025.8-2026.7/iScience_revision/14_pheno_analysis/pheno_conti_38.csv",header = T)
all_pheno=read.csv("E:/2025.8-2026.7/iScience_revision/14_pheno_analysis/pheno_46.csv",header = T)
# Covariates
covariates <- all_pheno[,c(1:4)]
covariates[covariates == ""] <- NA
covariates$Sex <- ifelse(covariates$Sex == "female", 0, 1)
covariates$ethnicity <- ifelse(covariates$ethnicity == "hispanic", 0, 1)
covariates_data <- apply(covariates[, 2:4], 2, as.numeric)
colnames(covariates_data)

colnames(sa_pheno)
aim0=merge(t_diff_df,sa_pheno,by="SPID")
aim=aim0

species_data <- apply(aim[, 2:13], 2, as.numeric)
phenotype_data <- apply(aim[, 15:51], 2, as.numeric);dim(phenotype_data) # 37

# Initialization
cor_matrix <- matrix(NA, nrow = ncol(species_data), ncol = ncol(phenotype_data))
p_matrix <- matrix(NA, nrow = ncol(species_data), ncol = ncol(phenotype_data))
q_matrix <- matrix(NA, nrow = ncol(species_data), ncol = ncol(phenotype_data))
err_matrix <- matrix(NA, nrow = ncol(species_data), ncol = ncol(phenotype_data))

rownames(cor_matrix) <- colnames(species_data)
colnames(cor_matrix) <- colnames(phenotype_data)
rownames(p_matrix) <- colnames(species_data)
colnames(p_matrix) <- colnames(phenotype_data)
rownames(q_matrix) <- colnames(species_data)
colnames(q_matrix) <- colnames(phenotype_data)
rownames(err_matrix) <- colnames(species_data)
colnames(err_matrix) <- colnames(phenotype_data)

# Linear regression
for (i in 1:ncol(species_data)) {
  for (j in 1:ncol(phenotype_data)) {
	df <- as.data.frame(cbind(pheno = phenotype_data[, j], microbe = log2(species_data[, i] + 1), covariates_data)) # log2 transformation
    fit <- lm(pheno ~ microbe + Age + Sex + ethnicity, data = df)
    s <- summary(fit)$coefficients
    cor_matrix[i, j] <- s["microbe", "Estimate"]
    p_matrix[i, j] <- s["microbe", "Pr(>|t|)"]
	err_matrix[i, j] <- s["microbe", "Std. Error"]
  }
}

q_matrix <- apply(p_matrix, 2, function(p) p.adjust(p, method = "fdr"))

p_matrix[p_matrix > 0.05] <- NA
q_matrix[q_matrix > 0.05] <- NA

fwrite(as.data.frame(cor_matrix), "E:/2025.8-2026.7/iScience_revision/14_pheno_analysis/analysis_abundance/continuous_cor_log2.csv", row.names = T) 
fwrite(as.data.frame(p_matrix), "E:/2025.8-2026.7/iScience_revision/14_pheno_analysis/analysis_abundance/continuous_p_log2.csv", row.names = T) 
fwrite(as.data.frame(q_matrix), "E:/2025.8-2026.7/iScience_revision/14_pheno_analysis/analysis_abundance/continuous_q_fdr_log2.csv", row.names = T)
fwrite(as.data.frame(err_matrix), "E:/2025.8-2026.7/iScience_revision/14_pheno_analysis/analysis_abundance/continuous_err_log2.csv", row.names = T)

## 2.Binary traits
# Prepare phenotype matrix
sa_pheno=read.csv("E:/2025.8-2026.7/iScience_revision/14_pheno_analysis/pheno_logic_8.csv",header = T) 
colnames(sa_pheno)
sa_pheno1 <- sa_pheno
apply(sa_pheno1, 2, table)
sa_pheno1[sa_pheno1 == ""] <- NA

aim=merge(t_diff_df,sa_pheno1,by="SPID")

all_pheno=read.csv("E:/2025.8-2026.7/iScience_revision/14_pheno_analysis/pheno_46.csv",header = T)
# Covariates
covariates <- all_pheno[,c(1:4)]
covariates[covariates == ""] <- NA
covariates$Sex <- ifelse(covariates$Sex == "female", 0, 1)
covariates$ethnicity <- ifelse(covariates$ethnicity == "hispanic", 0, 1)
covariates_data <- apply(covariates[, 2:4], 2, as.numeric)
colnames(covariates_data)

species_data <- apply(aim[, 2:13], 2, as.numeric)
phenotype_data <- apply(aim[, 16:21], 2, as.numeric)

# Initialization
cor_matrix <- matrix(NA, nrow = ncol(species_data), ncol = ncol(phenotype_data))
p_matrix <- matrix(NA, nrow = ncol(species_data), ncol = ncol(phenotype_data))
q_matrix <- matrix(NA, nrow = ncol(species_data), ncol = ncol(phenotype_data))
err_matrix <- matrix(NA, nrow = ncol(species_data), ncol = ncol(phenotype_data))

rownames(cor_matrix) <- colnames(species_data)
colnames(cor_matrix) <- colnames(phenotype_data)
rownames(p_matrix) <- colnames(species_data)
colnames(p_matrix) <- colnames(phenotype_data)
rownames(q_matrix) <- colnames(species_data)
colnames(q_matrix) <- colnames(phenotype_data)
rownames(err_matrix) <- colnames(species_data)
colnames(err_matrix) <- colnames(phenotype_data)

# Logistic regression
for (i in 1:ncol(species_data)) {
  for (j in 1:ncol(phenotype_data)) {
	df <- as.data.frame(cbind(pheno = phenotype_data[, j], microbe = log2(species_data[, i] + 1), covariates_data)) # log2 transformation
	df$pheno <- ifelse(df$pheno == 0, FALSE, TRUE)
    fit <- glm(pheno ~ microbe + Age + Sex + ethnicity, data = df, family = binomial)
    s <- summary(fit)$coefficients
    cor_matrix[i, j] <- s["microbe", "Estimate"]
    p_matrix[i, j] <- s["microbe", "Pr(>|z|)"]
	err_matrix[i, j] <- s["microbe", "Std. Error"]
  }
}

q_matrix <- apply(p_matrix, 2, function(p) p.adjust(p, method = "fdr"))

p_matrix[p_matrix > 0.05] <- NA
q_matrix[q_matrix > 0.05] <- NA

fwrite(as.data.frame(cor_matrix), "E:/2025.8-2026.7/iScience_revision/14_pheno_analysis/analysis_abundance/logic_cor_log2.csv", row.names = T) 
fwrite(as.data.frame(p_matrix), "E:/2025.8-2026.7/iScience_revision/14_pheno_analysis/analysis_abundance/logic_p_log2.csv", row.names = T) 
fwrite(as.data.frame(q_matrix), "E:/2025.8-2026.7/iScience_revision/14_pheno_analysis/analysis_abundance/logic_q_fdr_log2.csv", row.names = T)
fwrite(as.data.frame(err_matrix), "E:/2025.8-2026.7/iScience_revision/14_pheno_analysis/analysis_abundance/logic_err_log2.csv", row.names = T)

# Initialization
cor_matrix <- matrix(NA, nrow = ncol(species_data), ncol = ncol(phenotype_data))
p_matrix <- matrix(NA, nrow = ncol(species_data), ncol = ncol(phenotype_data))
q_matrix <- matrix(NA, nrow = ncol(species_data), ncol = ncol(phenotype_data))
updown_matrix <- matrix(NA, nrow = ncol(species_data), ncol = ncol(phenotype_data))

rownames(cor_matrix) <- colnames(species_data)
colnames(cor_matrix) <- colnames(phenotype_data)
rownames(p_matrix) <- colnames(species_data)
colnames(p_matrix) <- colnames(phenotype_data)
rownames(q_matrix) <- colnames(species_data)
colnames(q_matrix) <- colnames(phenotype_data)
rownames(updown_matrix) <- colnames(species_data)
colnames(updown_matrix) <- colnames(phenotype_data)

# Fitrh logistic regression
for (i in 1:ncol(species_data)) {
  for (j in 1:ncol(phenotype_data)) {
	df <- as.data.frame(cbind(pheno = phenotype_data[, j], microbe = log2(species_data[, i] + 1), covariates_data)) # log2 transformation
	df$pheno <- ifelse(df$pheno == 0, FALSE, TRUE)
	fit <- logistf(pheno ~ microbe + Age + Sex + ethnicity, data = df)
	s <- fit$coefficients
	cor_matrix[i, j] <- s["microbe"]
	p_matrix[i, j] <- fit$prob["microbe"]
	updown_matrix[i, j] <- paste0("(", confint(fit)["microbe",][1], ", ", confint(fit)["microbe",][2], ")")
	}
}

q_matrix <- apply(p_matrix, 2, function(p) p.adjust(p, method = "fdr"))

p_matrix[p_matrix > 0.05] <- NA
q_matrix[q_matrix > 0.05] <- NA

fwrite(as.data.frame(cor_matrix), "E:/2025.8-2026.7/iScience_revision/14_pheno_analysis/analysis_abundance/logicf_cor_log2.csv", row.names = T) 
fwrite(as.data.frame(p_matrix), "E:/2025.8-2026.7/iScience_revision/14_pheno_analysis/analysis_abundance/logicf_p_log2.csv", row.names = T) 
fwrite(as.data.frame(q_matrix), "E:/2025.8-2026.7/iScience_revision/14_pheno_analysis/analysis_abundance/logicf_q_fdr_log2.csv", row.names = T)
fwrite(as.data.frame(updown_matrix), "E:/2025.8-2026.7/iScience_revision/14_pheno_analysis/analysis_abundance/logicf_updown_log2.csv", row.names = T)

# Result combination
continuous_cor=read.csv("E:/2025.8-2026.7/iScience_revision/14_pheno_analysis/analysis_abundance/continuous_cor_log2.csv", header=T)
continuous_p=read.csv("E:/2025.8-2026.7/iScience_revision/14_pheno_analysis/analysis_abundance/continuous_p_log2.csv", header=T)
continuous_q=read.csv("E:/2025.8-2026.7/iScience_revision/14_pheno_analysis/analysis_abundance/continuous_q_fdr_log2.csv",header=T)

logic_cor=read.csv("E:/2025.8-2026.7/iScience_revision/14_pheno_analysis/analysis_abundance/logic_cor_log2.csv", header=T)
logic_p=read.csv("E:/2025.8-2026.7/iScience_revision/14_pheno_analysis/analysis_abundance/logic_p_log2.csv", header=T) # 3 associations significant (p < 0.05 firth correction)
logic_q=read.csv("E:/2025.8-2026.7/iScience_revision/14_pheno_analysis/analysis_abundance/logic_q_fdr_log2.csv", header=T) # 0 associations significant

logicf_cor=read.csv("E:/2025.8-2026.7/iScience_revision/14_pheno_analysis/analysis_abundance/logicf_cor_log2.csv", header=T)
logicf_p=read.csv("E:/2025.8-2026.7/iScience_revision/14_pheno_analysis/analysis_abundance/logicf_p_log2.csv", header=T) # 3 associations significant
logicf_q=read.csv("E:/2025.8-2026.7/iScience_revision/14_pheno_analysis/analysis_abundance/logicf_q_fdr_log2.csv", header=T) # 0 associations significant

# logistic p < 0.05 => firth logistic p
logic_cor_final <- logic_cor
logic_cor_final[9,3] <- logicf_cor[9,3]
logic_cor_final[10,2] <- logicf_cor[10,2]
logic_cor_final[10,3] <- logicf_cor[10,3]

corval=merge(continuous_cor,logic_cor_final,by="X")
pval=merge(continuous_p,logicf_p,by="X")
qval=merge(continuous_q,logic_q,by="X")

fwrite(as.data.frame(corval), "E:/2025.8-2026.7/iScience_revision/14_pheno_analysis/analysis_abundance/cor_all_final.csv", row.names = F)
fwrite(as.data.frame(pval), "E:/2025.8-2026.7/iScience_revision/14_pheno_analysis/analysis_abundance/pval_all_final.csv", row.names = F)
fwrite(as.data.frame(qval), "E:/2025.8-2026.7/iScience_revision/14_pheno_analysis/analysis_abundance/qval_all_final.csv", row.names = F)

# Visualization
corval <- read.csv("E:/2025.8-2026.7/iScience_revision/14_pheno_analysis/analysis_abundance/cor_all_final.csv", header = TRUE, row.names = 1)
pval <- read.csv("E:/2025.8-2026.7/iScience_revision/14_pheno_analysis/analysis_abundance/pval_all_final.csv", header = TRUE, row.names = 1)
type <- read.csv("E:/2025.8-2026.7/iScience_revision/14_pheno_analysis/pheno_final_cate_46.csv", header = TRUE)
unique(type$Type)

rownames(corval) <- gsub("_", " ", rownames(corval))
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

pval <- pval[, match(colnames(corval), colnames(pval))]

colnames(corval) <- type$Pheno[match(colnames(corval), type$pheno)]
colnames(pval) <- type$Pheno[match(colnames(pval), type$pheno)]

if (!all(colnames(corval) == colnames(pval))) {
  stop("Column names of corval and pval are not consistent after renaming.")
}

corval <- corval[, order(match(colnames(corval), type$Pheno))]
pval <- pval[, order(match(colnames(pval), type$Pheno))]

annotation_col <- data.frame(Type = type$Type[match(colnames(corval), type$Pheno)])
rownames(annotation_col) <- colnames(corval)

corval <- corval[species_order, , drop = FALSE]
pval <- pval[species_order, , drop = FALSE]

significance_matrix <- ifelse(is.na(pval), "", 
                              ifelse(pval < 0.001, "***",
                                     ifelse(pval < 0.01, "**",
                                            ifelse(pval < 0.05, "*", ""))))
rownames(significance_matrix) <- rownames(pval)
colnames(significance_matrix) <- colnames(pval)

col_fun <- colorRamp2(c(-0.1, 0, 0.1), c("#0f86a9", "white", "#FC8452"))

phenotype_order <- c("Host", "Pregnancy", "Diet", "Disease/Medication",
                     "Adaptive", "Cognitive", "Behavior", "Social/Communication")
annotation_col$Type <- factor(annotation_col$Type, levels = phenotype_order)

if (!all(colnames(corval) %in% rownames(annotation_col))) {
  stop("Error: Some columns in 'corval' do not have matching rows in 'annotation_col'.")
}

ht <- Heatmap(
  corval,
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
  row_names_gp = gpar(fontsize = 10, fontface = "italic"),
  column_names_rot = 45,
  column_names_gp = gpar(fontsize = 10),
  column_title_gp = gpar(fontsize = 12, fontface = "bold"),
  heatmap_legend_param = list(
    title = "Correlation",
    title_gp = gpar(fontsize = 12, fontface = "bold"),
    labels_gp = gpar(fontsize = 10),
    at = c(-0.1, 0, 0.1),
    labels = c("-0.1", "0", "0.1")
  )
)

pdf("E:/2025.8-2026.7/iScience_revision/14_pheno_analysis/analysis_abundance/heatmap_output_log2.pdf", width = 15, height = 6.6)
draw(ht)
dev.off()