####------------------------------------------------------- (13).differential analysis -------------------------------------------------------
#################################################### 13.3 prevalence-CLR ####################################################
rm(list=ls())
library(data.table)
library(dplyr)
library(reshape2)
library(survival)
library(ggplot2)
library(ggrepel)
library(tidyr)

# Genus--------------------------------------------------------------------
## 1.Data prepare
metadisc <- read.csv("E:/2025.8-2026.7/iScience_revision/1_rawdata/SSC_sample_phenotype_age&bmi_impute.csv",header = T,row.names = 1)
metadata <- as.data.frame(metadisc[, c("SFID", "Role","Sex","age_at_ados_pmm","bmi_pmm","Predicted.ancestry"), drop = FALSE])
table(metadata$Role)
class(metadata$Role)

metadata$Role <- ifelse(metadata$Role == "s1", 0, 1)
metadata$Role <- as.numeric(metadata$Role)

otu_raw <- read.csv("E:/2025.8-2026.7/iScience_revision/11_depth_normalize/Species_relative_abundance_100.csv")
cate=read.csv("E:/2025.8-2026.7/iScience_revision/1_rawdata/SSC_kraken_cate_16344.csv",row.names = 1)
colnames(cate)[colnames(cate) == "V1"] <- "X"
otu=merge(otu_raw,cate,by="X")
otu=otu[,c(2:3893,3899)];length(unique(otu$Genus)) # 57
otu <- otu %>%
  group_by(Genus) %>%
  summarise(across(everything(), \(x) sum(x, na.rm = TRUE))) %>%
  ungroup()
otu=as.data.frame(otu)
row.names(otu)=otu$Genus
otu <- otu[, !(colnames(otu) %in% c("Genus"))]

df_read0=as.data.frame(t(otu),stringsAsFactors = FALSE)
df_read=df_read0[rownames(metadata),]
df_read[df_read> 0] <- 1 

aim=merge(df_read, metadata, by.x = "row.names", by.y = "row.names")
colnames(aim)[1]="SPID"
tt1=aim

## 2.Conditional logistic regression
poss_p <- data.frame(Genus = character(), Coefficient = numeric(), P_poss = numeric())

for (i in 2:58) {
  species_data <- tt1[, c(1, i, 59:64)]
  colnames(species_data)[2] <- "Genus"
  species_data$Role=as.factor(species_data$Role)
  species_data$Role<- relevel(species_data$Role, ref = "0")
  
  model <- clogit(Genus ~ Role + Sex + age_at_ados_pmm + bmi_pmm + Predicted.ancestry + strata(SFID),
                  data = species_data)

  coef_val <- summary(model)$coefficients["Role1", "coef"]
  p_val    <- summary(model)$coefficients["Role1", "Pr(>|z|)"]

  poss_p <- rbind(poss_p, data.frame(Genus = colnames(tt1)[i], Coefficient = coef_val, P_poss = p_val))
}

poss_p$q_poss <- p.adjust(poss_p$P_poss, method = "BH")
colnames(poss_p)
colnames(poss_p)=c("Genus","Coefficient", "poss_p","poss_q")
sum(poss_p$poss_p < 0.05) # 14
sum(poss_p$poss_q < 0.05) # 13

pheno_group = as.data.frame(metadisc[, c("Role"), drop = FALSE])
pheno_group$Role <- gsub("p1", "ASD", gsub("s1", "NT", pheno_group$Role))

df_combined <- merge(df_read, pheno_group, by = "row.names", all.x = TRUE)
df_combined <- df_combined[, -1]

results <- data.frame(Genus = colnames(df_read), 
                             ASD_N = rep(NA, ncol(df_read)), 
                             NT_N = rep(NA, ncol(df_read)),
                             stringsAsFactors = FALSE)

for (i in 1:ncol(df_read)) {
  species_data <- df_combined[, c(i, ncol(df_combined))]
  colnames(species_data) <- c("Genus", "Role")

  results$ASD_N[i] <- sum(species_data$Genus[species_data$Role == "ASD"] == 1)
  results$NT_N[i] <- sum(species_data$Genus[species_data$Role == "NT"] == 1)
}

colnames(results)
results <- results[, c("Genus", "ASD_N", "NT_N")]

aim=merge(poss_p, results, by="Genus")
fwrite(aim, "E:/2025.8-2026.7/iScience_revision/13_differential_analysis/3_DA/Prevalence_Genus/prevalence_diff_result_Genus.csv", row.names = F) 
fwrite(poss_p, "E:/2025.8-2026.7/iScience_revision/13_differential_analysis/3_DA/Prevalence_Genus/prevalence_diff_result_Genus_CLR.csv", row.names = F)

## 3.Visualization of results of conditional logistic regression
aim$log_poss_q <- -log10(aim$poss_q)
aim$color <- "grey"
aim$color[aim$poss_q < 0.05 & aim$Coefficient > 0] <- "red"
aim$color[aim$poss_q < 0.05 & aim$Coefficient < 0] <- "blue"
aim$label <- ifelse(aim$poss_q < 0.05, aim$Genus, NA)

ggplot(aim, aes(x = Coefficient, y = log_poss_q, color = color)) +
  geom_point(alpha = 0.85, size = 3) +
  scale_color_manual(values = c("grey" = "grey", "red" = "red", "blue" = "blue")) +
  ylim(c(-1, 21)) +
  geom_hline(yintercept = -log10(0.05), col = "black", lwd = 0.8) +
  geom_text_repel(
    aes(label = label),
    size = 3,
    box.padding = unit(0.3, "lines"),
    point.padding = unit(0.3, "lines"),
    segment.color = "grey50",
    max.overlaps = 10000
  ) +
  labs(title = "Conditional Logistic Regression Genus",x = "Coefficient",y = "-log10 q value") +
  theme_bw() +
  theme(
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    legend.position = "none"
  )

ggsave("Genus_prevalence_diff_CLR.pdf",  device = "pdf", path = "E:/2025.8-2026.7/iScience_revision/13_differential_analysis/3_DA/Prevalence_Genus/",width=4,height=4)

## 4.Visualization of the intergroup prevalence of significant genera
final_result= aim[aim$poss_q < 0.05, ] # 13

final_result$ASD_N <- final_result$ASD_N / 1946 * 100
final_result$NT_N <- final_result$NT_N / 1946 * 100

final_result <- final_result %>%
  arrange(poss_q)

final_result$Genus <- factor(final_result$Genus, levels = final_result$Genus)

long_data <- final_result %>%
  pivot_longer(cols = c("ASD_N", "NT_N"), names_to = "Group", values_to = "Value")

ggplot(long_data, aes(x = Genus, y = Value, fill = Group)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), 
           width = 0.7, color = "black") +
  scale_fill_manual(values = c("ASD_N" = "#f1761d", "NT_N" = "#3273ae"), 
                    labels = c("ASD", "NT")) +
  labs(x = "", y = "Prevalence(%)", fill = "Group") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 10),
    plot.title = element_text(hjust = 0.5),
    plot.margin = unit(c(1, 1, 2, 1.5), "cm"),
    legend.position = "top"
  )
ggsave("E:/2025.8-2026.7/iScience_revision/13_differential_analysis/3_DA/Prevalence_Genus/Genus_prevalence_diff13.pdf", device = "pdf", width = 6, height = 6)