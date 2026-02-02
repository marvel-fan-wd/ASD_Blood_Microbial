library(ggplot2);library(data.table);library(survival)
meta <- read.csv("E:/2025.8-2026.7/iScience_revision/SSC_sample_meta.csv",header = T,row.names = 1)
metapart <- as.data.frame(meta[, c("SFID", "group","sex","age","bmi","ancestry"), drop = FALSE])
df <- read.csv("E:/2025.8-2026.7/iScience_revision/Species_abundance.csv")
aim=merge(df, metapart, by.x = "row.names", by.y = "row.names")
aa <- data.frame(Species = character(), Coefficient = numeric(), p = numeric(),
                 stringsAsFactors = FALSE)

for (i in 2:101) {
  species_data <- aim[, c(1, i, 102:107)]
  colnames(species_data)[2] <- "Species"

  species_data$group <- as.factor(species_data$group)
  species_data$group <- relevel(species_data$group, ref = "0")

  model <- clogit(Species ~ group + sex + age + bmi + ancestry + strata(SFID),
                  data = species_data)
  sm <- summary(model)$coefficients
  grp_row <- grep("^group", rownames(sm), value = TRUE)[1]

  coef_val <- sm[grp_row, "coef"]
  p_val    <- sm[grp_row, "Pr(>|z|)"]

  aa <- rbind(aa,
              data.frame(Species = names(df)[i],
                         Coefficient = coef_val,
                         p = p_val,
                         stringsAsFactors = FALSE))
}

aa$qvalue <- p.adjust(aa$p, method = "BH")