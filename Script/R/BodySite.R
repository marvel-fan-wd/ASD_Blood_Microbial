library(ggplot2);library(dplyr)

site_counts <- read.delim("E:/2025.8-2026.7/iScience_revision/Body_sites.tsv")
ggplot(site_counts, aes(x = Body_site, y = No_of_species, fill = Body_site)) +
  geom_bar(stat = "identity") +
  labs(x = "Body site", y = "Number of species", title = "") +
  scale_fill_brewer(palette = "Set3") +
  theme_minimal() +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1))