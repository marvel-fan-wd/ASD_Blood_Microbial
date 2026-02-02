library(ggplot2)
pre_top30 <- read.delim("E:/2025.8-2026.7/iScience_revision/PrevalenceTop.tsv")
ggplot(pre_top30, aes(x = reorder(Species, -all), y = Prevalence, fill = Group)) +
  geom_bar(stat = "identity", width = 0.7) +
  scale_fill_manual(values = c("NT" = "#3273ae","ASD" = "#f1761d")) +
  labs(x = "", y = "Prevalence(%)", title = "") +
  scale_y_continuous(limits = c(0, 20.5)) + 
  theme(panel.background = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, size = 1),
        axis.ticks.length = unit(0.2, "cm"),
        axis.ticks = element_line(size = 0.8),
        axis.text.x = element_text(angle = 55, hjust = 1, vjust = 1, size = 12),
        axis.text.y  = element_text(size = 12, colour = "black"),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        plot.margin = unit(c(1, 1, 2, 1), "cm"))