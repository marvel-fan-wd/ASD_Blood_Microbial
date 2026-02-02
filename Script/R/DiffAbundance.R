library(Maaslin2);library(dplyr);library(ggplot2);library(ggrepel)

meta <- read.csv("E:/2025.8-2026.7/iScience_revision/SSC_sample_meta.csv",header = T,row.names = 1)
metadata <- as.data.frame(meta[, c("SFID", "group","sex","age","bmi","ancestry"), drop = FALSE])
metadata$Role <- ifelse(metadata$Role == "s1", 0, 1)
metadata$Role <- as.factor(metadata$Role);metadata$SFID <- as.factor(metadata$SFID)

df <- read.csv("E:/2025.8-2026.7/iScience_revision/Species_abundance.csv")
df_read0 <- as.data.frame(t(df), stringsAsFactors = FALSE)
df_read <- df_read0[rownames(metadata), ]

df_input_read <- df_read;df_input_meta <- metadata
Maaslin2(input_data = df_input_read, 
                      input_metadata = df_input_meta,
                      min_abundance = 0,
                      min_prevalence = 0.0001,
                      max_significance = 0.05,
                      normalization = "TSS",
                      transform = "LOG",
                      analysis_method = "LM",
                      correction = "BH",
                      standardize = TRUE,
                      plot_heatmap = TRUE,
                      plot_scatter = TRUE,
                      fixed_effects = c("group","sex","age","bmi","ancestry"),
                      random_effects = c("SFID"),
                      output = "E:/2025.8-2026.7/iScience_revision/differential_analysis/Maaslin2/")