####------------------------------------------------------- Figure S4 -------------------------------------------------------
rm(list=ls())
library(dplyr)
library(tidyr)
library(data.table)
library(ggplot2)
library(ggrepel)

metadisc <- read.csv("E:/2025.8-2026.7/iScience_revision/1_rawdata/SSC_sample_phenotype_age&bmi_impute.csv",header = T,row.names = 1)
metadata <- as.data.frame(metadisc[, c("Role","Sex","age_at_ados_pmm","bmi_pmm","Predicted.ancestry"), drop = FALSE])
table(metadata$Role)
class(metadata$Role)

metadata$Role <- ifelse(metadata$Role == "s1", 0, 1)
metadata$Role <- as.numeric(metadata$Role)

otu_raw <- read.csv("E:/2025.8-2026.7/iScience_revision/11_depth_normalize/Species_relative_abundance_100.csv")
cate=read.csv("E:/2025.8-2026.7/iScience_revision/1_rawdata/SSC_kraken_cate_16344.csv",row.names = 1)
colnames(cate)[colnames(cate) == "V1"] <- "X"
cate=cate[,c("Species","X")]
otu=merge(otu_raw,cate,by="X")
row.names(otu)=otu$Species
otu <- otu[, !(colnames(otu) %in% c("Species", "X"))]
df_read0=as.data.frame(t(otu),stringsAsFactors = FALSE)
df_read=df_read0[rownames(metadata),]
df_read[df_read> 0] <- 1

aim=merge(df_read, metadata, by.x = "row.names", by.y = "row.names")
colnames(aim)[1]="SPID"
tt1=aim
pheno_group = as.data.frame(metadisc[, c("Role"), drop = FALSE])
pheno_group$Role <- gsub("p1", "ASD", gsub("s1", "NT", pheno_group$Role))

df_combined <- merge(df_read, pheno_group, by = "row.names", all.x = TRUE)  
df_combined <- df_combined[, -1]  

results <- data.frame(Species = colnames(df_read), 
                             ASD_N = rep(NA, ncol(df_read)), 
                             NT_N = rep(NA, ncol(df_read)),
                             stringsAsFactors = FALSE)

for (i in 1:ncol(df_read)) {
  species_data <- df_combined[, c(i, ncol(df_combined))]  
  colnames(species_data) <- c("Species", "Role")

  results$ASD_N[i] <- sum(species_data$Species[species_data$Role == "ASD"] == 1)  
  results$NT_N[i] <- sum(species_data$Species[species_data$Role == "NT"] == 1)   
}

colnames(results)
results <- results[, c("Species", "ASD_N", "NT_N")]

poss_p1 <- read.csv("E:/2025.8-2026.7/iScience_revision/13_differential_analysis/3_DA/Prevalence/prevalence_diff_result_CLR.csv", row.names = 1)
maaslin <- read.delim("E:/2025.8-2026.7/iScience_revision/13_differential_analysis/3_DA/Maaslin2/all_results.tsv")
colnames(maaslin)[1] <- "Species"

diff_all <- union(poss_p1[poss_p1$q_poss < 0.05,]$Species, maaslin[maaslin$qval < 0.05 & maaslin$metadata == "Role",]$Species)
diff_all
 # [1] "Prescottella_equi"          "Parabacteroides_distasonis" "Elizabethkingia_anophelis"  "Bdellovibrio_bacteriovorus"
 # [5] "Ochrobactrum_quorumnocens"  "Wolbachia_pipientis"        "Pseudoalteromonas_sp._3J6"  "Enterobacter_hormaechei"   
 # [9] "Klebsiella_michiganensis"   "Yersinia_enterocolitica"    "Xylella_fastidiosa"         "Treponema_pallidum" 
 
aim=merge(poss_p1,results,by="Species")

final_result=aim[aim$q_poss<0.05,] # 12

final_result$ASD_N <- final_result$ASD_N / 1946 * 100
final_result$NT_N <- final_result$NT_N / 1946 * 100

final_result <- final_result %>%
  arrange(desc(q_poss))

final_result$Species <- gsub("_", " ", final_result$Species)
final_result$Species <- factor(final_result$Species, levels = final_result$Species)

long_data <- final_result %>%
  pivot_longer(cols = c("ASD_N", "NT_N"), names_to = "Group", values_to = "Value")

p <- ggplot(long_data, aes(x = Species, y = Value, fill = Group)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), 
           width = 0.7, color = "black") +  
  scale_fill_manual(values = c("ASD_N" = "#f1761d", "NT_N" = "#3273ae"), 
                    labels = c("ASD", "NT")) +  
  labs(x = "", y = "Prevalence(%)", fill = "Group") +  
  scale_y_continuous(limits = c(0, 20)) + 
  coord_flip() +
  theme(panel.background = element_blank(),
        
        panel.border = element_rect(color = "black", fill = NA, size = 1),
        axis.ticks.length = unit(0.2, "cm"),        
        axis.ticks = element_line(size = 0.8),      
        axis.text.x = element_text(angle = 0, hjust = 1, vjust = 1, colour = "black", size = 12),
        axis.text.y  = element_text(size = 12, colour = "black", face = "italic"),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        legend.position = c(0.85, 0.132))

ggsave("Figure S4.pdf", plot = p, device = "pdf", path = "E:/2025.8-2026.7/iScience_revision/Figures/", width=6, height=4)
