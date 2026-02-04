library(data.table);library(ggplot2)

#Phylum
data <- read.delim("E:/2025.8-2026.7/iScience_revision/Abundance_percent_phylum.tsv")
custom_colors <- c(
  "Pseudomonadota" = "#87BEFA",
  "Actinomycetota" = "#FF9999",
  "Spirochaetota" = "#7dc43a",
  "Bacteroidota" = "#DE9AEA",
  "Chlamydiota" = "#FF7F00",
  "Bacillota" = "#FFE7A3",
  "Residuals" = "#A6A6A6"
)

ggplot(data, aes(x = Role, y = Relatice_Abundance, fill = Phylum)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(x = "", y = "Relative abundance (%)") +
  scale_fill_manual(values = custom_colors)+
  scale_y_continuous(limits = c(0, 100), 
                       breaks = seq(0, 100, 20), labels = c("0", "20", "40", "60", "80", "100")) + 
  theme(plot.title = element_text(hjust = 0.5, vjust = 0.5,
								   face = "bold",
                                   size = 14, color = "black"),
		panel.background = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, size = 0.8),
        axis.ticks.length = unit(0.2, "cm"),
        axis.ticks = element_line(size = 0.8, colour = "black"),
        axis.text = element_text(size = 12, colour = "black"),
        axis.title.x = element_text(size = 12, colour = "black"),
        axis.title.y = element_text(size = 12, colour = "black"),
		legend.position = "none",
		plot.margin = unit(c(1, 4, 0, 8), "mm"))

# Class
data <- read.delim("E:/2025.8-2026.7/iScience_revision/Abundance_percent_class.tsv")
custom_colors1 <- c(
  "Residuals" = "#A6A6A6",
  "Negativicutes" = "#fae732",
  "Bacilli" = "#FFE7A3",
  "Chlamydiia" = "#FF7F00",
  "Sphingobacteriia" = "#984EA3",
  "Flavobacteriia" = "#DE9AEA",
  "Bacteroidia" = "#F8D4FF",
  "Coriobacteriia" = "#e6644b",
  "Actinomycetes" = "#FF9999",
  "Spirochaetia" = "#7dc43a",
  "Gammaproteobacteria" = "#377EB8",
  "Betaproteobacteria" = "#87BEFA",
  "Alphaproteobacteria" = "#D2E5FF"
)

ggplot(data, aes(x = Role, y = Relatice_Abundance, fill = Class)) +
  geom_bar(stat = "identity") +
  labs(x = "", y = "Relative abundance (%)", fill = "Class", title = "Class") +
  scale_fill_manual(values = custom_colors1)+
  scale_y_continuous(limits = c(0, 100), 
                       breaks = seq(0, 100, 20), labels = c("0", "20", "40", "60", "80", "100")) + 
  theme(plot.title = element_text(hjust = 0.5, vjust = 0.5,
								   face = "bold",
                                   size = 14, color = "black"),
		panel.background = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, size = 0.8),
        axis.ticks.length = unit(0.2, "cm"),
        axis.ticks = element_line(size = 0.8, colour = "black"),
        axis.text = element_text(size = 12, colour = "black"),
        axis.title.x = element_text(size = 12, colour = "black"),
        axis.title.y = element_text(size = 12, colour = "black"),
		legend.position = "none",
		plot.margin = unit(c(1, 4, 0, 8), "mm"))