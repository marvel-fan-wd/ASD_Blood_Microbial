library(vegan);library(ggplot2)

df <- read.csv("E:/2025.8-2026.7/iScience_revision/Species_abundance.csv",row.names = 1,header = T)
df <- t(df);otu.distance <- vegdist(df)
otu.distance <- as.matrix(round(otu.distance,digits = 3))
adonis2(otu.distance~Role, data=df, distance = "bray", permutations = 999, strata = df$SFID)

