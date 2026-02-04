library(vegan);library(picante);library(dplyr);library(ggpubr)

df <- read.csv("E:/2025.8-2026.7/iScience_revision/Species_abundance.csv",row.names = 1,header = T)
Shannon <- diversity(df, index = "shannon", MARGIN = 2)
Simpson <- diversity(df, index = "simpson", MARGIN = 2)
Richness <- specnumber(df, MARGIN = 2)
alpha <- as.data.frame(cbind(Shannon, Simpson, Richness))
meta <- read.csv("E:/2025.8-2026.7/iScience_revision/SSC_sample_meta.csv",header = T)
alpha2 <- merge(alpha,meta,by="SPID")
alpha2$Role <- gsub("p1", "ASD", gsub("s1", "NT", alpha2$Role))

# Wilcoxon signed-rank test
library(dplyr)
alpha2_asd <- alpha2 %>% filter(Role == "ASD")
alpha2_nt <- alpha2 %>% filter(Role == "NT")
paired <- merge(alpha2_asd, alpha2_nt, by = "SFID", suffixes = c("_ASD", "_NT"))
wilcox.test(paired$Shannon_ASD, paired$Shannon_NT, paired = TRUE)

ggbarplot(
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
)+
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