####------------------------------------------------------- (13).differential analysis -------------------------------------------------------
#################################################### 13.2 diversity analysis-alpha diversity ####################################################
rm(list=ls())
library(vegan)
library(picante)      
library(dplyr)
library(ggpubr)

metadisc <- read.csv("E:/2025.8-2026.7/iScience_revision/1_rawdata/SSC_sample_phenotype_age&bmi_impute.csv",header = T)
meta_part <- as.data.frame(metadisc[, c("SPID","SFID","Role"), drop = FALSE])

df <- read.csv("E:/2025.8-2026.7/iScience_revision/11_depth_normalize/Species_relative_abundance_100.csv",row.names = 1,header = T)
cate=read.csv("E:/2025.8-2026.7/iScience_revision/1_rawdata/SSC_kraken_cate_16344.csv",row.names = 1)

Shannon <- diversity(df, index = "shannon", MARGIN = 2)
Simpson <- diversity(df, index = "simpson", MARGIN = 2)
Richness <- specnumber(df, MARGIN = 2)
alpha <- as.data.frame(cbind(Shannon, Simpson, Richness))
alpha$SPID=rownames(alpha)
alpha2 <- merge(alpha,meta_part,by="SPID")
alpha2$Role <- gsub("p1", "ASD", gsub("s1", "NT", alpha2$Role))
colnames(alpha2)

alpha2 %>%
  group_by(Role) %>%
  summarise(mean_shannon = mean(Shannon, na.rm = TRUE),
            mean_simpson = mean(Simpson, na.rm = TRUE),
            mean_richness = mean(Richness, na.rm = TRUE))
  # Role  mean_shannon mean_simpson mean_richness
  # <chr>        <dbl>        <dbl>         <dbl>
# 1 ASD          0.776        0.624          5.31
# 2 NT           0.826        0.617          5.66

# Wilcoxon signed-rank test
library(dplyr)
alpha2_asd <- alpha2 %>% filter(Role == "ASD")
alpha2_nt <- alpha2 %>% filter(Role == "NT")
paired_data <- merge(alpha2_asd, alpha2_nt, by = "SFID", suffixes = c("_ASD", "_NT"))
# Shannon
wilcox.test(paired_data$Shannon_ASD, paired_data$Shannon_NT, paired = TRUE)
# Simpson
wilcox.test(paired_data$Simpson_ASD, paired_data$Simpson_NT, paired = TRUE)
# Richness
wilcox.test(paired_data$Richness_ASD, paired_data$Richness_NT, paired = TRUE)

# Paired Wilcoxon rank test
# Shannon p-value = 0.0004121
# Simpson p-value = 0.2643
# Richness p-value = 1.163e-06

# Visualization
df <- alpha2
df$Role <- factor(df$Role, levels = c("ASD", "NT"))

# Shannon
p <- ggbarplot(
  df, 
  x = "Role", 
  y = "Shannon",
  alpha = 0.8,
  color = "Role",
  fill = "Role",
  palette = c("ASD" = "#f1761d", "NT" = "#3273ae"),
  add = c("mean_se"),
  add.params = list(size = 1),
  xlab = "",
  ylab = "Shannon",
  position = position_dodge()
)

p +
  coord_cartesian(ylim = c(0, 1)) +
  theme(
    panel.background = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    axis.ticks.length = unit(0.2, "cm"),
    axis.ticks = element_line(size = 0.8),
    axis.text = element_text(size = 14, colour = "black"),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
	legend.position = "none"
  )

ggsave("alpha_Shannon_bar.pdf", device = "pdf", path = "E:/2025.8-2026.7/iScience_revision/13_differential_analysis/2_diversity_analysis/",width=3.4,height=4.2)

# Simpson
p1 <- ggbarplot(
  df, 
  x = "Role",
  y = "Simpson",
  alpha = 0.8,
  color = "Role",
  fill = "Role",
  palette = c("ASD" = "#f1761d", "NT" = "#3273ae"),
  add = c("mean_se"),
  add.params = list(size = 1),
  xlab = "",
  ylab = "Simpson",
  position = position_dodge()
)

p1 +
  coord_cartesian(ylim = c(0, 0.75)) +
  theme(
    panel.background = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    axis.ticks.length = unit(0.2, "cm"),
    axis.ticks = element_line(size = 0.8),
    axis.text = element_text(size = 12, colour = "black"),
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12),
	legend.position = "none"
  )

ggsave("alpha_Simpson_bar.pdf",  device = "pdf", path = "E:/2025.8-2026.7/iScience_revision/13_differential_analysis/2_diversity_analysis/",width=3.4,height=4.2)
# Richness
p2 <- ggbarplot(
  df, 
  x = "Role",
  y = "Richness",
  alpha = 0.8,
  color = "Role",
  fill = "Role",
  palette = c("ASD" = "#f1761d", "NT" = "#3273ae"),
  add = c("mean_se"),
  add.params = list(size = 1),
  xlab = "",
  ylab = "Richness",
  position = position_dodge()
)

p2 +
  coord_cartesian(ylim = c(0, 6.5)) +
  theme(
    panel.background = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    axis.ticks.length = unit(0.2, "cm"),
    axis.ticks = element_line(size = 0.8),
    axis.text = element_text(size = 12, colour = "black"),
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12),
	legend.position = "none"
  )

ggsave("alpha_Richness_bar.pdf",  device = "pdf", path = "E:/2025.8-2026.7/iScience_revision/13_differential_analysis/2_diversity_analysis/",width=3.4,height=4.2)

# Load reads
df0 <- read.csv("E:/2025.8-2026.7/iScience_revision/10_prevalence_filter/SSC_specieQC_prevalence_filter_0.2_100.csv",row.names = 1,header = T)

df1 <- df0
dft <- as.data.frame(t(df1))
dft$sum <- rowSums(dft)
dft$SPID <- rownames(dft)
dftload <- dft[,c("SPID", "sum")]
alpha3 <- merge(dftload, meta_part,by="SPID")
alpha3$Role <- gsub("p1", "ASD", gsub("s1", "NT", alpha3$Role))

alpha3 %>%
  group_by(Role) %>%
  summarise(mean_sum = mean(sum, na.rm = TRUE))
  # Role  mean_sum
  # <chr>    <dbl>
# 1 ASD       712.
# 2 NT        678.

alpha3_asd <- alpha3 %>% filter(Role == "ASD")
alpha3_nt <- alpha3 %>% filter(Role == "NT")
paired_data <- merge(alpha3_asd, alpha3_nt, by = "SFID", suffixes = c("_ASD", "_NT"))
# Load
wilcox.test(paired_data$sum_ASD, paired_data$sum_NT, paired = TRUE) # reads p-value = 0.01039

df2 <- alpha3
df2$Role <- factor(df2$Role, levels = c("ASD", "NT"))

p3 <- ggbarplot(
  df2, 
  x = "Role",
  y = "sum",
  alpha = 0.8,
  color = "Role",
  fill = "Role",
  palette = c("ASD" = "#f1761d", "NT" = "#3273ae"),
  add = c("mean_se"),
  add.params = list(size = 1),
  xlab = "",
  ylab = "Reads",
  position = position_dodge()
)

p3 + coord_cartesian(ylim = c(0, 900)) +
  theme(
    panel.background = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    axis.ticks.length = unit(0.2, "cm"),
    axis.ticks = element_line(size = 0.8),
    axis.text = element_text(size = 12, colour = "black"),
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12),
	legend.position = "none"
  )

ggsave("alpha_load_reads_bar.pdf", device = "pdf", path = "E:/2025.8-2026.7/iScience_revision/13_differential_analysis/2_diversity_analysis/",width=3.4,height=4.2)