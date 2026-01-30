####------------------------------------------------------- (12).micro characterization -------------------------------------------------------
## 12.2 Oxygen type
rm(list=ls())
library(readxl)
library(ggplot2)
library(dplyr)
library(scales)

oxygen <- read_xlsx("E:/2025.8-2026.7/iScience_revision/12_micro_characterization/Oxygen_type/Specie_oxygen_type_100.xlsx")
unique(oxygen$Oxygen)

oxygen_counts <- oxygen %>%
  count(Oxygen) %>%
  mutate(percentage = n / sum(n) * 100, # Calculate percentage
         label = paste(n, "\n", round(percentage, 1), "%", sep = "")) # Only display numbers and percentages
colnames(oxygen_counts)

df=oxygen_counts %>% 
  arrange(desc(percentage)) # Sort in descending order by percentage

ggplot(df, aes(x = 2, y = n, fill = reorder(Oxygen, -percentage))) +
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

ggsave("E:/2025.8-2026.7/iScience_revision/12_micro_characterization/Oxygen_type/Oxygen_Type_donut.pdf", device = "pdf", width=6.2,height=4.5)