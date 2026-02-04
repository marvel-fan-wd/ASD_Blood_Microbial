rm(list=ls())
cor_tt <- read.delim("E:/2025.8-2026.7/iScience_revision/Validation.tsv")

## ------------------------------------------------------------- prevalence ------------------------------------------------------------- 
cor.test(cor_tt$prevalence_discovery, cor_tt$prevalence_validation, method = "spearman")
# intersection
intersect(cor_tt[order(cor_tt$prevalence_discovery, decreasing = T),]$Genus[1:30], cor_tt[order(cor_tt$prevalence_validation, decreasing = T),]$Genus[1:30])

## ------------------------------------------------------------- abundance ------------------------------------------------------------- 
cor.test(cor_tt$mean_abundance_discovery, cor_tt$mean_abundance_validation, method = "spearman")
# intersection

intersect(cor_tt[order(cor_tt$mean_abundance_discovery, decreasing = T),]$Genus[1:30], cor_tt[order(cor_tt$mean_abundance_validation, decreasing = T),]$Genus[1:30])
