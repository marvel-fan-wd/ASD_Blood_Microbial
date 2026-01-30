####------------------------------------------------------- (13).differential analysis -------------------------------------------------------
#################################################### 13.3 prevalence-summary-species ####################################################
rm(list=ls())
library(data.table)
library(dplyr)
library(reshape2)
library(tidyr)
library(ggplot2)

## 1.Histogram of top 30 prevelance
metadisc <- read.csv("E:/2025.8-2026.7/iScienceÎÄÕÂ·µÐÞ/Newflow/1_rawdata/SSC_sample_phenotype_age&bmi_impute.csv",header = T,row.names = 1)
metadata <- as.data.frame(metadisc[, c("Role"), drop = FALSE])
table(metadata$Role)
class(metadata$Role)
metadata$Role <- gsub("p1", "ASD", gsub("s1", "NT", metadata$Role))

otu_raw <- read.csv("E:/2025.8-2026.7/iScience_revision/11_depth_normalize/Species_relative_abundance_100.csv")
cate=read.csv("E:/2025.8-2026.7/iScience_revision/1_rawdata/SSC_kraken_cate_16344.csv",row.names = 1)
colnames(cate)[colnames(cate) == "V1"] <- "X"
cate=cate[,c("Species","X")]
otu=merge(otu_raw,cate,by="X")
row.names(otu)=otu$Species
otu <- otu[, !(colnames(otu) %in% c("Species", "X"))]
df_read0=as.data.frame(t(otu),stringsAsFactors = FALSE)
df_read=df_read0[rownames(metadata),]
df_read[df_read> 0] <- 1 

aim=merge(df_read, metadata, by.x = "row.names", by.y = "row.names")
colnames(aim)[1]="SPID"
row.names(aim)=aim$SPID

aim$Role <- as.character(aim$Role)
long_data <- aim %>%
  pivot_longer(cols = 2:101,
               names_to = "Species",
               values_to = "Presence")

# Calculate the number of samples of each species in the ASD and NT groups
result <- long_data %>%
  filter(Presence == 1) %>%
  group_by(Species, Role) %>%
  summarise(Count = n()) %>%
  spread(key = Role, value = Count, fill = 0)

final_result <- result
final_result$ASD=final_result$ASD/3892*100
final_result$NT=final_result$NT/3892*100
final_result <- final_result %>%
  mutate(all = (ASD + NT) / 2) 

top30_species <- final_result %>%
  arrange(desc(all)) %>% 
  head(30)

top30_long <- top30_species %>%
  pivot_longer(cols = c(ASD, NT),
               names_to = "Group",
               values_to = "Prevalence")

top30_long$Species <- gsub("_", " ", top30_long$Species)

ggplot(top30_long, aes(x = reorder(Species, -all), y = Prevalence, fill = Group)) +
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

ggsave("E:/2025.8-2026.7/iScience_revision/13_differential_analysis/3_DA/Prevalence/prevalence_top30.pdf", device = "pdf", width = 10, height = 8)

## 2.Histogram of sample species quantity
aim=merge(df_read, metadata, by.x = "row.names", by.y = "row.names")
colnames(aim)[1]="SPID"
row.names(aim)=aim$SPID
species_counts <- rowSums(aim[, 2:101])
data_for_plot <- data.frame(species_count = species_counts, Role = aim$Role)
plot_data <- data_for_plot %>%
  group_by(species_count, Role) %>%
  summarise(count = n()) %>%
  ungroup()
plot_data$Role <- factor(plot_data$Role, levels = c("ASD", "NT"))
filtered_data <- plot_data[plot_data$species_count != 0, ]
# Calculate the average
weighted_sum <- sum(filtered_data$species_count * filtered_data$count)
count_sum <- sum(filtered_data$count)
weighted_sum / count_sum # 7.587984
# ASD
asd_data <- filtered_data[filtered_data$Role == "ASD", ]
weighted_sum <- sum(asd_data $species_count * asd_data $count)
count_sum <- sum(asd_data $count)
weighted_sum / count_sum # 7.518923
# NT
nt_data <- filtered_data[filtered_data$Role == "NT", ]
weighted_sum <- sum(nt_data $species_count * nt_data $count)
count_sum <- sum(nt_data $count)
weighted_sum / count_sum # 7.653926

# Calculate the median
expanded_species_counts <- rep(filtered_data$species_count, filtered_data$count)
median(expanded_species_counts) # 4
# ASD
asd_data <- filtered_data[filtered_data$Role == "ASD", ]
median(rep(asd_data$species_count, asd_data$count)) # 4
# NT
nt_data <- filtered_data[filtered_data$Role == "NT", ]
median(rep(nt_data$species_count, nt_data$count)) # 4

# Calculate the number of species for each sample
species_counts <- rowSums(aim[, 2:101])
# Classify those greater than 20 as category 20+
species_counts[species_counts > 20] <- "20+"
data_for_plot <- data.frame(species_count = species_counts, Role = aim$Role)
plot_data <- data_for_plot %>%
  group_by(species_count, Role) %>%
  summarise(count = n(), .groups = "drop")

# Make sure that "Role" is listed as a factor and specify its display order
plot_data$Role <- factor(plot_data$Role, levels = c("ASD", "NT"))

# Convert the species number column into factors and manually set their display order
plot_data$species_count <- factor(plot_data$species_count, levels = c(0:20, "20+"))

ggplot(plot_data, aes(x = species_count, y = count, fill = Role)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = c("ASD" = "#f1761d", "NT" = "#3273ae")) +
  labs(x = "Number of Species", y = "Number of Samples", title = "") +
  theme(panel.background = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, size = 1),
        axis.ticks.length = unit(0.2, "cm"),
        axis.ticks = element_line(size = 0.8),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 12),
        axis.text.y  = element_text(size = 12, colour = "black"),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        plot.margin = unit(c(1, 1, 2, 1), "cm") )

ggsave("E:/2025.8-2026.7/iScience_revision/13_differential_analysis/3_DA/Prevalence/sample_specie_number_100.pdf", device = "pdf", width=7, height=5)

# Detect status difference between groups
aim=merge(df_read, metadata, by.x = "row.names", by.y = "row.names")
colnames(aim)[1]="SPID"

species_data <- aim[, 2:101]
role_data <- aim[, 101]
species_count <- apply(species_data, 1, function(x) sum(x > 0))
count <- data.frame(SPID = aim$SPID,
                    Role = role_data,
                    Species_Count = species_count)
count %>%
  group_by(Role) %>%
  summarise(Average_Species_Count = mean(Species_Count, na.rm = TRUE))
   # Role Average_Species_Count
  # <dbl>                 <dbl>
# 1     0                  5.22
# 2     1                  6.73

# Fisher test
aim$contains_species <- rowSums(aim[, 2:101]) > 0  # Contain: TRUE; not contain: FALSE
contingency_table <- table(aim$Role,aim$contains_species)
      # FALSE TRUE
  # ASD   572 1374
  # NT    507 1439
fisher_result <- fisher.test(contingency_table)
#p-value = 0.02189 odds ratio = 1.181523  

data_summary <- aim %>%
  group_by(Role, contains_species) %>%
  summarise(Count = n(), .groups = "drop") %>%
  mutate(contains_species = ifelse(contains_species == 1, "Yes", "No"))  # Convert to classification tags

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

ggsave("specie_present_ratio.pdf",  device = "pdf", path = "E:/2025.8-2026.7/iScience_revision/13_differential_analysis/3_DA/Prevalence", width=3, height=4)