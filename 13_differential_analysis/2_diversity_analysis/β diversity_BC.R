####------------------------------------------------------- (13).differential analysis -------------------------------------------------------
#################################################### 13.2 diversity analysis-beta diversity ####################################################
rm(list=ls())
library(vegan)
library(ggplot2)

otu_raw <- read.csv("E:/2025.8-2026.7/iScience_revision/11_depth_normalize/Species_relative_abundance_100.csv", row.names = 1)
# Calculate the average number of species for each sample
mean(apply(otu_raw > 0, 2, sum)) # 5.484327

otu_raw <- otu_raw[, colSums(otu_raw) != 0]
otu <- t(otu_raw)
otu.distance <- vegdist(otu)
otu.distance <- as.matrix(round(otu.distance,digits = 3))

meta <- read.csv("E:/2025.8-2026.7/iScienceÎÄÕÂ·µÐÞ/Scripts_final/1_rawdata/SSC_sample_phenotype_age&bmi_impute.csv",header = T)
rownames(meta)=meta$SPID
meta2 <- meta[colnames(otu_raw), , drop = FALSE]

# PCOA
pcoa <- cmdscale (otu.distance,eig=TRUE)
pc12 <- pcoa$points[,1:2]
pc <- round(pcoa$eig/sum(pcoa$eig)*100,digits=2)#½âÊÍ¶È

pc12 <- as.data.frame(pc12)
pc12$samples <- row.names(pc12)
head(pc12)
                    # V1         V2   samples
# SS0012979 -0.037456222 0.07790324 SS0012979
# SSC00081  -0.008595475 0.05918308  SSC00081
# SSC00082  -0.015932929 0.06032064  SSC00082
# SSC00115  -0.006598969 0.06103132  SSC00115
# SSC00119  -0.010249951 0.04588307  SSC00119
# SSC00137  -0.008605407 0.04695205  SSC00137

group <- as.data.frame(meta2[, c("SPID","SFID","Role","Sex","age_at_ados_pmm","bmi_pmm","Predicted.ancestry"), drop = FALSE])
colnames(group) <- c("samples","SFID","Role","Sex","age_at_ados_pmm","bmi_pmm","Predicted.ancestry")

df <- merge(pc12, group, by="samples")
head(df)
    # samples          V1          V2  SFID Role    Sex age_at_ados_pmm bmi_pmm Predicted.ancestry
# 1 SS0012979 -0.03745622  0.07790324 13750   s1   male           50.04    14.4                OTH
# 2 SS0012980  0.50030077 -0.15296700 12160   s1   male           63.96    23.1                EUR
# 3 SS0012982  0.37617817 -0.36191276 11266   s1 female           69.96    14.6                EUR
# 4 SS0012986 -0.13883726 -0.07574521 13867   s1   male          150.96    22.2                EUR
# 5 SS0012989  0.33561379 -0.08612189 13521   s1 female           48.96    16.4                EUR
# 6 SS0013009  0.21093558 -0.01300054 13348   p1   male           54.00    15.4                EUR

ggplot(data=df,aes(x=V1,y=V2,color=Role))+
  geom_point(size=1.5,alpha=0.6)+
  theme(panel.grid = element_blank())+
  labs(x=paste0("PC1 ",pc[1],"%"),
       y=paste0("PC2 ",pc[2],"%"),
       title = "PCOA:BC Distance")+
  scale_color_manual(values = c("#F6C63C","lightblue")) +
  theme(plot.title = element_text(hjust = 0.5,size = 15),
        panel.background = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, size = 1),
        axis.ticks = element_line(size=0.6, colour = "black"),
        axis.text.x = element_text(size = 10, colour = "black"),
        axis.text.y = element_text(size = 10, colour = "black"),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12))
ggsave("beta_BC_PCOA_100_ASD.pdf",  device = "pdf", path = "E:/2025.8-2026.7/iScience_revision/13_differential_analysis/2_diversity_analysis/", width=5.5, height=5)

# Adonis test
Adonis <- adonis2(otu.distance~Role, data=df, distance = "bray", permutations = 999, strata = df$SFID)
Adonis
# Permutation test for adonis under reduced model
# Blocks:  strata 
# Permutation: free
# Number of permutations: 999

# adonis2(formula = otu.distance ~ Role, data = df, permutations = 999, strata = df$SFID, distance = "bray")
           # Df SumOfSqs      R2      F Pr(>F)
# Model       1     0.41 0.00033 0.9196  0.238
# Residual 2811  1250.07 0.99967              
# Total    2812  1250.47 1.00000      
save(Adonis, file="E:/2025.8-2026.7/iScience_revision/13_differential_analysis/2_diversity_analysis/beta_BC_adonis.RData")
load("E:/2025.8-2026.7/iScience_revision/13_differential_analysis/2_diversity_analysis/beta_BC_adonis.RData")