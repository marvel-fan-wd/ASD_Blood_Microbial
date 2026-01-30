####------------------------------------------------------- Figure S3 -------------------------------------------------------
rm(list=ls())
library(vegan)
library(picante)      
library(dplyr)
library(ggpubr)
library(cowplot)

metadisc <- read.csv("E:/2025.8-2026.7/iScience_revision/1_rawdata/SSC_sample_phenotype_age&bmi_impute.csv",header = T)
meta_part <- as.data.frame(metadisc[, c("SPID","SFID","Role"), drop = FALSE])

df <- read.csv("E:/2025.8-2026.7/iScience_revision/11_depth_normalize/Species_relative_abundance_100.csv",row.names = 1,header = T)
cate=read.csv("E:/2025.8-2026.7/iScience_revision/1_rawdata/SSC_kraken_cate_16344.csv",row.names = 1)

Simpson <- diversity(df, index = "simpson", MARGIN = 2)
Richness <- specnumber(df, MARGIN = 2)
alpha <- as.data.frame(cbind(Simpson, Richness))
alpha$SPID=rownames(alpha)
alpha2 <- merge(alpha,meta_part,by="SPID")
alpha2$Role <- gsub("p1", "ASD", gsub("s1", "NT", alpha2$Role))
colnames(alpha2)

alpha2 %>%
  group_by(Role) %>%
  summarise(mean_simpson = mean(Simpson, na.rm = TRUE),
            mean_richness = mean(Richness, na.rm = TRUE))
  # Role  mean_simpson mean_richness
  # <chr>        <dbl>         <dbl>
# 1 ASD          0.624          5.31
# 2 NT           0.617          5.66

df <- alpha2
df$Role <- factor(df$Role, levels = c("ASD", "NT"))  

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

df2 <- alpha3
df2$Role <- factor(df2$Role, levels = c("ASD", "NT"))  

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

pb <- p1 +
  coord_cartesian(ylim = c(0, 0.79)) +
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

pc <- p2 +
  coord_cartesian(ylim = c(0, 7.2)) +
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
  ylab = "Microbial DNA load",  
  position = position_dodge()  
)

pa <- p3 + coord_cartesian(ylim = c(0, 910)) +  
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
  

plotall <- plot_grid(pa, pb, pc, labels = c("A", "B", "C"), label_size = 18, ncol = 3, rel_widths = c(1.07, 1.01, 1))

ggsave("Figure S3.pdf", plot = plotall, device = "pdf", path = "E:/2025.8-2026.7/iScience_revision/Figures/", width=10, height=4)
