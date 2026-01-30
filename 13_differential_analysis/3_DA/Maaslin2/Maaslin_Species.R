####------------------------------------------------------- (13).differential analysis -------------------------------------------------------
#################################################### 13.3 relative abundance-MaAsLin2 ####################################################
rm(list=ls())
library(Maaslin2)
library(dplyr)
library(ggplot2)
library(ggrepel)

# Species--------------------------------------------------------------------
metadisc <- read.csv("E:/2025.8-2026.7/iScience_revision/1_rawdata/SSC_sample_phenotype_age&bmi_impute.csv",header = T,row.names = 1)
metadata <- as.data.frame(metadisc[, c("SFID", "Role","Sex","age_at_ados_pmm","bmi_pmm","Predicted.ancestry"), drop = FALSE])
table(metadata$Role)
class(metadata$Role)

metadata$Role <- ifelse(metadata$Role == "s1", 0, 1)
metadata$Role <- as.factor(metadata$Role)
metadata$SFID <- as.factor(metadata$SFID)

otu_raw <- read.csv("E:/2025.8-2026.7/iScience_revision/11_depth_normalize/Species_relative_abundance_100.csv")
cate=read.csv("E:/2025.8-2026.7/iScience_revision/1_rawdata/SSC_kraken_cate_16344.csv",row.names = 1)

colnames(cate)[colnames(cate) == "V1"] <- "X"
cate <- cate[, c("Species", "X")]
otu <- merge(otu_raw, cate, by = "X")
row.names(otu) <- otu$Species
otu <- otu[, !(colnames(otu) %in% c("Species", "X"))]

df_read0 <- as.data.frame(t(otu), stringsAsFactors = FALSE)
df_read <- df_read0[rownames(metadata), ]

df_input_read <- df_read
df_input_meta <- metadata

# Binary transformation
df_input_meta$Ancestry_bin <- ifelse(df_input_meta$Predicted.ancestry == "EUR", "EUR", "nonEUR")
df_input_meta$Ancestry_bin <- relevel(factor(df_input_meta$Ancestry_bin), ref = "EUR")

# Run the Maaslin2 analysis with Role (ASD as the reference level), Sex, age_at_ados_pmm, bmi_pmm, and Predicted ancestry (EUR as the reference level) as fixed effect
fit_data1 <- Maaslin2(input_data = df_input_read, 
                      input_metadata = df_input_meta,
                      min_abundance = 0,
                      min_prevalence = 0.0001,
                      max_significance = 0.05,
                      normalization = "TSS",
                      transform = "LOG",
                      analysis_method = "LM",
                      correction = "BH",
                      standardize = TRUE,
                      plot_heatmap = TRUE,
                      plot_scatter = TRUE,
                      fixed_effects = c("Role","Sex","age_at_ados_pmm","bmi_pmm","Ancestry_bin"),
                      random_effects = c("SFID"),
                      output = "E:/2025.8-2026.7/iScience_revision/13_differential_analysis/3_DA/Maaslin2/")

# Take the intersection with the prevalence difference species
rm(list=ls())
ma <- read.delim("E:/2025.8-2026.7/iScience_revision/13_differential_analysis/3_DA/Maaslin2/all_results.tsv", header = TRUE)
ma_sig=ma[ma$qval<0.05,] # 12: 11 Role; 1 BMI
ma12=ma_sig
colnames(ma12)[1]=c("Species")
ma12[12,1] <- "Thalassobaculum_sp._OXR-137"

cate=read.csv("E:/2025.8-2026.7/iScience_revision/1_rawdata/SSC_kraken_cate_16344.csv",row.names = 1)
ma12_cate=merge(ma12, cate, by="Species") #12
fwrite(ma12_cate, "E:/2025.8-2026.7/iScience_revision/13_differential_analysis/3_DA/Maaslin2/abundance_diff_result_sig.csv", row.names = F)
fwrite(ma12_cate[ma12_cate$metadata == "Role", ], "E:/2025.8-2026.7/iScience_revision/13_differential_analysis/3_DA/Maaslin2/abundance_diff_result_sig_Role.csv", row.names = F)

pre=read.csv("E:/2025.8-2026.7/iScience_revision/13_differential_analysis/3_DA/Prevalence/prevalence_diff_result_CLR.csv",row.names = 1)
pre_sig <- pre[pre$q_poss < 0.05, ]
intersect(ma_sig$feature, pre_sig$Species)#11
 # [1] "Pseudoalteromonas_sp._3J6"  "Prescottella_equi"          "Parabacteroides_distasonis" "Klebsiella_michiganensis"  
 # [5] "Xylella_fastidiosa"         "Wolbachia_pipientis"        "Enterobacter_hormaechei"    "Elizabethkingia_anophelis" 
 # [9] "Yersinia_enterocolitica"    "Treponema_pallidum"         "Bdellovibrio_bacteriovorus"

# Visualization
ma$feature <- gsub("_", " ", ma$feature)
ma$log_qval <- -log10(ma$qval)
ma$color <- "grey"
ma$color[ma$qval < 0.05 & ma$coef > 0] <- "red"
ma$color[ma$qval < 0.05 & ma$coef < 0] <- "blue"
ma$label <- ifelse(ma$qval < 0.05, ma$feature, NA)
ma1 <- ma[ma$metadata == "Role", ]

ggplot(ma1, aes(x = coef, y = log_qval, color = color)) +
  geom_point(alpha = 0.85, size = 3) +
  scale_color_manual(values = c("grey" = "grey", "red" = "red", "blue" = "blue")) +
  ylim(c(-1, 24)) +
  geom_hline(yintercept = -log10(0.05), col = "black", lwd = 0.8, linetype = "dashed") +
  labs(
    title = "Maaslin2",
    x = "Coefficient",
    y = "-log10(qval)"
  ) +
  geom_text_repel(
    aes(label = label),
    size = 3,
    box.padding = unit(0.3, "lines"),
    point.padding = unit(0.3, "lines"),
    segment.color = "grey50",
    max.overlaps = 10000
  ) +
  theme(
    panel.background = element_blank(),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    axis.ticks.length = unit(0.2, "cm"),
    axis.ticks = element_line(size = 0.8),
    axis.text = element_text(size = 12, colour = "black"),
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    legend.position = "none"
  )

ggsave("prevalence_diff_Maaslin2.pdf", device = "pdf", path = "E:/2025.8-2026.7/iScience_revision/13_differential_analysis/3_DA/Maaslin2/", width=3.8, height=4.2)