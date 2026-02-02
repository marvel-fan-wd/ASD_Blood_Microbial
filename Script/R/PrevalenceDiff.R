library(ggplot2);library(data.table)
meta <- read.csv("E:/2025.8-2026.7/iScience_revision/SSC_sample_meta.csv",header = T,row.names = 1)
df <- read.csv("E:/2025.8-2026.7/iScience_revision/Species_abundance.csv")
aim=merge(df, metadata, by.x = "row.names", by.y = "row.names")

species_data <- aim[, 2:101];role_data <- aim[, 101]
species_count <- apply(species_data, 1, function(x) sum(x > 0))
count <- data.frame(SPID = aim$SPID,Role = role_data,
                    Species_Count = species_count)
count %>%
  group_by(Role) %>%
  summarise(Average_Species_Count = mean(Species_Count, na.rm = TRUE))
aim$contains_species <- rowSums(aim[, 2:101]) > 0 
contingency_table <- table(aim$Role,aim$contains_species)
fisher_result <- fisher.test(contingency_table)

data_summary <- aim %>%
  group_by(Role, contains_species) %>%
  summarise(Count = n(), .groups = "drop") %>%
  mutate(contains_species = ifelse(contains_species == 1, "Yes", "No")) 

ggplot(data_summary, aes(x = Role, y = Count, fill = contains_species)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = c("Yes" = "#f1a89a", "No" = "#D3D3D3")) +
  labs(x = "", y = "Number of samples ", fill = "", title = "") +
  geom_text(aes(label = Count), position = position_stack(vjust = 0.5), size = 5) +
  theme(panel.background = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, size = 1),
        axis.ticks.length = unit(0.2, "cm"),
        axis.ticks = element_line(size = 0.8),
        axis.text = element_text(size = 12, colour = "black"),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12))