library(data.table);library(ggplot2)

oxygen <- read.delim("E:/2025.8-2026.7/iScience_revision/oxygen_type.tsv")

ggplot(oxygen, aes(x = 2, y = n, fill = reorder(Oxygen, -percentage))) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar("y", start = 0) +
  geom_text(aes(label = paste0(round(percentage, 1), "%")), 
            position = position_stack(vjust = 0.5), size = 4) +
  scale_fill_brewer(palette = "Set3") +
  theme_void() +
  theme(
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 10)
  ) +
  xlim(0.5, 2.5) +
  labs(title = "", fill = "Oxygen Type")
