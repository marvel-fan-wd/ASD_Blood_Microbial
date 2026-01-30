####------------------------------------------------------- (13).differential analysis -------------------------------------------------------
#################################################### 13.3 relative abundance-MaAsLin2 ####################################################
rm(list=ls())
library(Maaslin2)
library(dplyr)
library(ggplot2)
library(ggrepel)

# Genus--------------------------------------------------------------------
cate=read.csv("E:/2025.8-2026.7/iScience_revision/1_rawdata/SSC_kraken_cate_16344.csv",row.names = 1)
colnames(cate)[colnames(cate) == "V1"] <- "X"
otu_raw <- read.csv("E:/2025.8-2026.7/iScience_revision/11_depth_normalize/Species_relative_abundance_100.csv")
otu=merge(otu_raw,cate,by="X") 
otu=otu[,c(2:3893,3899)];length(unique(otu$Genus)) # 57
otu <- otu %>%
  group_by(Genus) %>%
  summarise(across(everything(), \(x) sum(x, na.rm = TRUE))) %>%
  ungroup()
otu=as.data.frame(otu)
row.names(otu)=otu$Genus
otu <- otu[, !(colnames(otu) %in% c("Genus"))]

metadisc <- read.csv("E:/2025.8-2026.7/iScience_revision/1_rawdata/SSC_sample_phenotype_age&bmi_impute.csv",header = T,row.names = 1)
metadata <- as.data.frame(metadisc[, c("SFID","Role","Sex","age_at_ados_pmm","bmi_pmm","Predicted.ancestry"), drop = FALSE])
table(metadata$Role)
class(metadata$Role)

metadata$Role <- ifelse(metadata$Role == "s1", 0, 1)
metadata$Role <- as.factor(metadata$Role)

df_read0=as.data.frame(t(otu),stringsAsFactors = FALSE)
df_read=df_read0[rownames(metadata),]

df_input_read = df_read
df_input_meta = metadata

# Binary transformation
df_input_meta$Ancestry_bin <- ifelse(df_input_meta$Predicted.ancestry == "EUR", "EUR", "nonEUR")
df_input_meta$Ancestry_bin <- relevel(factor(df_input_meta$Ancestry_bin), ref = "EUR")

fit_data1 = Maaslin2(input_data=df_input_read, input_metadata=df_input_meta,
                     min_abundance=0, min_prevalence=0.0001,max_significance=0.05,
                     normalization="TSS",transform="LOG",analysis_method="LM",
                     correction="BH",standardize=T,plot_heatmap=T,plot_scatter=T,
                     fixed_effects=c("Role","Sex","age_at_ados_pmm","bmi_pmm","Ancestry_bin"),
					 random_effects = c("SFID"),
                     output="E:/2025.8-2026.7/iScience_revision/13_differential_analysis/3_DA/Maaslin2_Genus/")

# Take the intersection with the prevalence difference genera
rm(list=ls())
ma <- read.delim("E:/2025.8-2026.7/iScience_revision/13_differential_analysis/3_DA/Maaslin2_Genus/all_results.tsv", header = TRUE)
ma_sig=ma[ma$qval<0.05,] # 13: 11 Role; 1 BMI; 1 Age (Thalassobaculum: 2)

pre=read.csv("E:/2025.8-2026.7/iScience_revision/13_differential_analysis/3_DA/Prevalence_Genus/prevalence_diff_result_Genus_CLR.csv",row.names = 1)
pre_sig <- pre[pre$poss_q < 0.05, ]
intersect(ma_sig$feature, rownames(pre_sig)) #11
 # [1] "Pseudoalteromonas" "Prescottella"      "Parabacteroides"   "Klebsiella"        "Xylella"          
 # [6] "Enterobacter"      "Wolbachia"         "Treponema"         "Yersinia"          "Elizabethkingia"  
# [11] "Bdellovibrio"  

ma11_cate <- read.csv("E:/2025.8-2026.7/iScience_revision/13_differential_analysis/3_DA/Maaslin2/abundance_diff_result_sig_Role.csv")
length(intersect(ma11_cate$Genus,ma_sig$feature)) # 11

# Visualization
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
  geom_hline(yintercept = -log10(0.05), col = "black", lwd = 0.8) +
  geom_text_repel(
    aes(label = label),
    size = 3,
    box.padding = unit(0.3, "lines"),
    point.padding = unit(0.3, "lines"),
    segment.color = "grey50",
    max.overlaps = 10000
  ) +
  labs(
    title = "Maaslin2 Genus",
    x = "Coefficient",
    y = "-log10(qval)"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    legend.position = "none"
  )

ggsave("Genus_prevalence_diff_Maaslin2.pdf", device = "pdf", path = "E:/2025.8-2026.7/iScience_revision/13_differential_analysis/3_DA/Maaslin2_Genus/", width = 3.8, height = 4)