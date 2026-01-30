####------------------------------------------------------- (13).differential analysis -------------------------------------------------------
#################################################### 13.3 prevalence-CLR ####################################################
rm(list=ls())
library(data.table)
library(dplyr)
library(reshape2)
library(tidyr)
library(ggplot2)
library(survival)
library(ggrepel)

# Species--------------------------------------------------------------------
## 1.Data prepare
metadisc <- read.csv("E:/2025.8-2026.7/iScienceÎÄÕÂ·µÐÞ/Newflow/1_rawdata/SSC_sample_phenotype_age&bmi_impute.csv",header = T,row.names = 1)
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

## 2.Conditional logistic regression
metadata1 <- as.data.frame(metadisc[, c("SFID", "Role","Sex","age_at_ados_pmm","bmi_pmm","Predicted.ancestry"), drop = FALSE])
metadata1$Role <- ifelse(metadata1$Role == "s1", 0, 1)
metadata1$Role <- as.numeric(metadata1$Role)
aim1=merge(df_read, metadata1, by.x = "row.names", by.y = "row.names")
colnames(aim1)[1]="SPID"
tt11=aim1

poss_p1 <- data.frame(Species = character(),
                     Coefficient = numeric(),
                     P_poss = numeric())

for (i in 2:101) {
  species_data <- tt11[, c(1, i, 102:107)] # FamilyID
  colnames(species_data)[2] <- "Species"

  species_data$Role <- as.factor(species_data$Role)
  species_data$Role <- relevel(species_data$Role, ref = "0")

  model <- clogit(Species ~ Role + Sex + age_at_ados_pmm + bmi_pmm + Predicted.ancestry + strata(SFID),
                  data = species_data)

  coef_val <- summary(model)$coefficients["Role1", "coef"]
  p_val    <- summary(model)$coefficients["Role1", "Pr(>|z|)"]

  poss_p1 <- rbind(poss_p1,
                  data.frame(Species = colnames(tt11)[i],
                             Coefficient = coef_val,
                             P_poss = p_val))
}

poss_p1$q_poss <- p.adjust(poss_p1$P_poss, method = "BH")

fwrite(poss_p1, "E:/2025.8-2026.7/iScience_revision/13_differential_analysis/3_DA/Prevalence/prevalence_diff_result_CLR.csv", row.names = TRUE) 

## 3.Visualization of results of conditional logistic regression
colnames(poss_p1)[4] <- "poss_q"
poss_p1$log_poss_q <- -log10(poss_p1$poss_q)
poss_p1$color <- "grey"
poss_p1$color[poss_p1$poss_q < 0.05 & poss_p1$Coefficient > 0] <- "red"
poss_p1$color[poss_p1$poss_q < 0.05 & poss_p1$Coefficient < 0] <- "blue"
poss_p1$label <- ifelse(poss_p1$poss_q < 0.05, poss_p1$Species, NA)

ggplot(poss_p1, aes(x = Coefficient, y = log_poss_q, color = color)) +
  geom_point(alpha = 0.85, size = 3) +
  scale_color_manual(values = c("grey" = "grey", "red" = "red", "blue" = "blue")) +
  ylim(c(-1, 20.5)) +
  geom_hline(yintercept = -log10(0.05), col = "black", lwd = 0.8, linetype = "dashed") +
  geom_text_repel(
    aes(label = label),
    size = 3,
    box.padding = unit(0.3, "lines"),
    point.padding = unit(0.3, "lines"),
    segment.color = "grey50",
    max.overlaps = 10000
  ) +
  labs(title = "Conditional logistic regression",x = "Coefficient",y = "-log10 q value") +
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

ggsave("prevalence_diff_CLR.pdf",  device = "pdf", path = "E:/2025.8-2026.7/iScience_revision/13_differential_analysis/3_DA/Prevalence/", width=3, height=4.2)

## 4.Visualization of the intergroup prevalence of significant species
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

# Significant species
poss_p1 <- read.csv("E:/2025.8-2026.7/iScience_revision/13_differential_analysis/3_DA/Prevalence/prevalence_diff_result_CLR.csv", row.names = 1)
aim=merge(poss_p1, results, by="Species")

final_result=aim[aim$q_poss<0.05,] # 12
final_result$ASD_N <- final_result$ASD_N / 1946 * 100
final_result$NT_N <- final_result$NT_N / 1946 * 100

final_result <- final_result %>%
  arrange(q_poss)

final_result$Species <- factor(final_result$Species, levels = final_result$Species)

long_data <- final_result %>%
  pivot_longer(cols = c("ASD_N", "NT_N"), names_to = "Group", values_to = "Value")

ggplot(long_data, aes(x = Species, y = Value, fill = Group)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), 
           width = 0.7, color = "black") +
  scale_fill_manual(values = c("ASD_N" = "#f1761d", "NT_N" = "#3273ae"), 
                    labels = c("ASD", "NT")) +
  labs(x = "", y = "Prevalence(%)", fill = "Group") +
  theme(panel.background = element_blank(),
        plot.title = element_text(hjust = 0.5),
        panel.border = element_rect(color = "black", fill = NA, size = 1),
        axis.ticks.length = unit(0.2, "cm"),
        axis.ticks = element_line(size = 0.8),
        axis.text.x = element_text(angle = 55, hjust = 1, vjust = 1, size = 12),
        axis.text.y  = element_text(size = 12, colour = "black"),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        legend.position = "top",
        plot.margin = unit(c(1, 1, 2, 1.5), "cm"))

ggsave("E:/2025.8-2026.7/iScience_revision/13_differential_analysis/3_DA/Prevalence/prevalence_diffspecies12.pdf", device = "pdf", width = 4, height = 6.5)