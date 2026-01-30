####------------------------------------------------------- Figure 3 -------------------------------------------------------
rm(list=ls())
library(vegan)
library(picante)
library(dplyr)
library(cowplot)
library(tidyr)
library(data.table)
library(ggplot2)
library(reshape2)
library(ggpubr)
library(ggrepel)

#------------------------------------------------------------------ A ------------------------------------------------------------------
metadisc <- read.csv("E:/2025.8-2026.7/iScience_revision/1_rawdata/SSC_sample_phenotype_age&bmi_impute.csv",header = T)
meta_part <- as.data.frame(metadisc[, c("SPID","SFID","Role"), drop = FALSE])

df <- read.csv("E:/2025.8-2026.7/iScience_revision/11_depth_normalize/Species_relative_abundance_100.csv",row.names = 1,header = T)
cate=read.csv("E:/2025.8-2026.7/iScience_revision/1_rawdata/SSC_kraken_cate_16344.csv",row.names = 1)

Shannon <- diversity(df, index = "shannon", MARGIN = 2)
alpha <- as.data.frame(Shannon)
alpha$SPID=rownames(alpha)
alpha2 <- merge(alpha,meta_part,by="SPID")
alpha2$Role <- gsub("p1", "ASD", gsub("s1", "NT", alpha2$Role))
colnames(alpha2)

alpha2 %>%
  group_by(Role) %>%
  summarise(mean_shannon = mean(Shannon, na.rm = TRUE))
  # Role  mean_shannon
  # <chr>        <dbl>
# 1 ASD          0.776
# 2 NT           0.826

df <- alpha2
df$Role <- factor(df$Role, levels = c("ASD", "NT"))  

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

pa <- p +
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

#------------------------------------------------------------------ B ------------------------------------------------------------------
poss_p1 <- read.csv("E:/2025.8-2026.7/iScience_revision/13_differential_analysis/3_DA/Prevalence/prevalence_diff_result_CLR.csv", row.names = 1)

colnames(poss_p1)[4] <- "poss_q"
poss_p1$log_poss_q <- -log10(poss_p1$poss_q)
poss_p1$color <- "grey"  
poss_p1$color[poss_p1$poss_q < 0.05 & poss_p1$Coefficient > 0] <- "red"  
poss_p1$color[poss_p1$poss_q < 0.05 & poss_p1$Coefficient < 0] <- "blue"  
poss_p1$label <- ifelse(poss_p1$poss_q < 0.05, poss_p1$Species, NA)  
poss_p1$label <- gsub("_", " ", poss_p1$label)

poss_p2 <- poss_p1[poss_p1$Coefficient < 10 & poss_p1$Coefficient > -105, ]

pb <- ggplot(poss_p2, aes(x = Coefficient, y = log_poss_q, color = color)) +
  geom_point(alpha = 0.7, size = 3) +  
  scale_color_manual(values = c("grey" = "grey", "red" = "#f1761d", "blue" = "#3273ae")) +  
  ylim(c(-1, 21)) +
  geom_hline(yintercept = -log10(0.05), col = "black", lwd = 0.8, linetype = "dashed") +  
  geom_text_repel(
    aes(label = label),
	size = 4,
    box.padding = unit(0.3, "lines"),  
    point.padding = unit(0.3, "lines"),  
    segment.color = "grey50",  
    max.overlaps = 10000  
  ) +
  labs(
  x = "Coefficient",y = "-log10(qval)") +
  theme(
    panel.background = element_blank(),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    axis.ticks.length = unit(0.2, "cm"),        
    axis.ticks = element_line(size = 0.8),      
    axis.text = element_text(size = 14, colour = "black"),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    legend.position = "none"
  )

#------------------------------------------------------------------ C ------------------------------------------------------------------
ma <- read.delim("E:/2025.8-2026.7/iScience_revision/13_differential_analysis/3_DA/Maaslin2/all_results.tsv", header = TRUE)  

ma$feature <- gsub("_", " ", ma$feature)
ma$log_qval <- -log10(ma$qval)

ma$color <- "grey"  
ma$color[ma$qval < 0.05 & ma$coef > 0] <- "red"  
ma$color[ma$qval < 0.05 & ma$coef < 0] <- "blue"  

ma$label <- ifelse(ma$qval < 0.05, ma$feature, NA)  

ma_role <- ma[ma$metadata == "Role",]

pc <- ggplot(ma_role, aes(x = coef, y = log_qval, color = color)) +
  geom_point(alpha = 0.7, size = 3) +  
  scale_color_manual(values = c("grey" = "grey", "red" = "#f1761d", "blue" = "#3273ae")) +  
  scale_y_continuous(limits = c(-1, 24)
  ) +
  geom_hline(yintercept = -log10(0.05), col = "black", lwd = 0.8, linetype = "dashed") +  
  labs(
    x = "Coefficient",
    y = "-log10(qval)"
  ) +
  geom_text_repel(
    aes(label = label),
	size = 4,
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
    axis.text = element_text(size = 14, colour = "black"),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    legend.position = "none"
  )

pfinal <- plot_grid(pa, pb, pc, labels = c("A", "B", "C"), label_size = 18, ncol = 3, rel_widths = c(0.8, 1, 1))

ggsave("Figure 3.pdf", plot = pfinal, device = "pdf", path = "E:/2025.8-2026.7/iScience_revision/Figures", width=10, height=4)