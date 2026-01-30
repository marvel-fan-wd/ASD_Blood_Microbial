####------------------------------------------------------- Figure 2 -------------------------------------------------------
rm(list=ls())
library(data.table)
library(dplyr)
library(reshape2)
library(cowplot)
library(ggplot2)
library(tidyr)
library(RColorBrewer)
library(readxl)
library(scales) 

# ----------------------------------------------------------- A -----------------------------------------------------------
# Phylum ------------------------------------------------------------------
metadisc <- read.csv("E:/2025.8-2026.7/iScience_revision/1_rawdata/SSC_sample_phenotype_age&bmi_impute.csv",header = T,row.names = 1)
metadata <- as.data.frame(metadisc[, c("Role", "SFID"), drop = FALSE])
table(metadata$Role)
class(metadata$Role)

metadata$Role <- gsub("p1", "ASD", gsub("s1", "NT", metadata$Role))
metadata$SPID=rownames(metadata)

otu_raw <- read.csv("E:/2025.8-2026.7/iScience_revision/11_depth_normalize/Species_relative_abundance_100.csv")
cate=read.csv("E:/2025.8-2026.7/iScience_revision/1_rawdata/SSC_kraken_cate_16344.csv",row.names = 1)
colnames(cate)[colnames(cate) == "V1"] <- "X"
otu=merge(otu_raw,cate,by="X")
apply(otu[, 3894:3900], 2, function(x) length(unique(x)))

# Kingdom  Phylum   Class   Order  Family   Genus Species 
      # 1       9      15      29      41      57     100 
otu_Phylum=otu[,c(2:3893,3895)]
length(unique(otu_Phylum$Phylum)) # 9

otu_Phylum <- otu_Phylum %>%
  group_by(Phylum) %>%
  summarise(across(everything(), \(x) sum(x, na.rm = TRUE))) %>%
  ungroup()
otu_Phylum=as.data.frame(otu_Phylum)
row.names(otu_Phylum)=otu_Phylum$Phylum
otu_Phylum <- otu_Phylum[, !(colnames(otu_Phylum) %in% c("Phylum"))]
PP=as.data.frame(t(otu_Phylum))
PP$SPID=row.names(PP)
PP2=merge(PP,metadata,by="SPID")
row.names(PP2)=PP2$SPID;PP2 <- PP2[, !(colnames(PP2) %in% c("SPID"))]
apply(PP2[, 1:9], 2, function(x) sum(x != 0))
             # Actinomycetota                   Bacillota                Bacteroidota            Bdellovibrionota 
                       # 1008                         140                         462                         107 
           # Campylobacterota Candidatus_Saccharibacteria                 Chlamydiota              Pseudomonadota 
                         # 84                          17                         292                        2260 
              # Spirochaetota 
                        # 685 

PP2 <- PP2[,c(1:10)]

custom_colors <- c(
  "Pseudomonadota" = "#87BEFA",
  "Actinomycetota" = "#FF9999",
  "Spirochaetota" = "#7dc43a",
  "Bacteroidota" = "#DE9AEA",
  "Chlamydiota" = "#FF7F00",
  "Bacillota" = "#FFE7A3",
  "Residuals" = "#A6A6A6"
)

# Residuals = Candidatus_Saccharibacteria + Bdellovibrionota + Campylobacterota

PP2_long <- PP2 %>%
  pivot_longer(-Role, names_to = "Species", values_to = "Abundance") %>%
  mutate(Species = ifelse(Species %in% names(custom_colors), Species, "Residuals")) %>%
  group_by(Role, Species) %>%
  summarise(Abundance = mean(Abundance, na.rm = TRUE), .groups = "drop") %>%
  group_by(Role) %>%
  mutate(Relative_Abundance = Abundance / sum(Abundance) * 100) #

sum(PP2_long$Relative_Abundance) # 200

PP2_long

   # Role  Species        Abundance Relative_Abundance
   # <chr> <chr>              <dbl>              <dbl>
 # 1 ASD   Actinomycetota     71.5               7.28 
 # 2 ASD   Bacillota          37.2               3.79 
 # 3 ASD   Bacteroidota       71.3               7.26 
 # 4 ASD   Chlamydiota        50.8               5.18 
 # 5 ASD   Pseudomonadota    665.               67.8  
 # 6 ASD   Residuals           4.98              0.507
 # 7 ASD   Spirochaetota      80.4               8.19 
 # 8 NT    Actinomycetota     41.1               4.05 
 # 9 NT    Bacillota          11.3               1.11 
# 10 NT    Bacteroidota       48.0               4.73 
# 11 NT    Chlamydiota        49.0               4.83 
# 12 NT    Pseudomonadota    657.               64.8  
# 13 NT    Residuals           3.77              0.372
# 14 NT    Spirochaetota     204.               20.1  

species_order <- PP2_long %>%
  filter(Role == "ASD") %>%
  group_by(Species) %>%
  summarise(Total_Abundance = sum(Relative_Abundance, na.rm = TRUE)) %>%
  arrange(Total_Abundance) %>%
  pull(Species)

PP2_long$Species <- factor(PP2_long$Species, levels = species_order)

PP2_long0 <- PP2_long

p1 <- ggplot(PP2_long0, aes(x = Role, y = Relative_Abundance, fill = Species)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(x = "", y = "Abundance (%)", fill = "Phylum", title = "Phylum") +
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

# Class --------------------------------------------------------------------
otu_Class=otu[,c(2:3893,3896)]
length(unique(otu_Class$Class)) # 15

otu_Class <- otu_Class %>%
  group_by(Class) %>%
  summarise(across(everything(), \(x) sum(x, na.rm = TRUE))) %>%
  ungroup()
otu_Class=as.data.frame(otu_Class)
row.names(otu_Class)=otu_Class$Class
otu_Class <- otu_Class[, !(colnames(otu_Class) %in% c("Class"))]
PP=as.data.frame(t(otu_Class))
PP$SPID=row.names(PP)
PP2=merge(PP,metadata,by="SPID")
row.names(PP2)=PP2$SPID;PP2 <- PP2[, !(colnames(PP2) %in% c("SPID"))]
apply(PP2[, 1:15], 2, function(x) sum(x != 0))
             # Actinomycetes        Alphaproteobacteria                    Bacilli                Bacteroidia 
                      # 1008                        791                        121                        247 
           # Bdellovibrionia         Betaproteobacteria Candidatus_Saccharimonadia                 Chlamydiia 
                       # 107                       1646                         17                        292 
            # Coriobacteriia      Epsilonproteobacteria             Flavobacteriia        Gammaproteobacteria 
                        # 15                         84                        301                        878 
             # Negativicutes           Sphingobacteriia               Spirochaetia 
                        # 20                         82                        685 

PP2 <- PP2[,c(1:16)]
PP2_long <- PP2 %>%
  pivot_longer(-Role, names_to = "Species", values_to = "Abundance") %>%
  group_by(Role, Species) %>%
  summarise(Abundance = mean(Abundance, na.rm = TRUE), .groups = "drop") %>%
  group_by(Role) %>%
  mutate(Relative_Abundance = Abundance / sum(Abundance) * 100)

cate=unique(otu[,c(3895:3896)])
colnames(cate)[2] <- "Species"
PP2_long2=merge(PP2_long,cate,by="Species")
colnames(PP2_long2)
PP2_long2
                      # Species Role   Abundance Relative_Abundance                      Phylum
# 1               Actinomycetes  ASD  68.1185644         6.87005201              Actinomycetota
# 2               Actinomycetes   NT  40.0771383         3.92218576              Actinomycetota
# 3         Alphaproteobacteria  ASD 386.3717713        38.96726520              Pseudomonadota
# 4         Alphaproteobacteria   NT 398.6026325        39.00961089              Pseudomonadota
# 5                     Bacilli  ASD   1.7138184         0.17284600                   Bacillota
# 6                     Bacilli   NT   2.6602530         0.26034809                   Bacillota
# 7                 Bacteroidia   NT  39.8215708         3.89717441                Bacteroidota
# 8                 Bacteroidia  ASD  63.8078367         6.43529646                Bacteroidota
# 9             Bdellovibrionia  ASD   2.8344295         0.28586448            Bdellovibrionota
# 10            Bdellovibrionia   NT   1.5098779         0.14776558            Bdellovibrionota
# 11         Betaproteobacteria  ASD 165.3820481        16.67949526              Pseudomonadota
# 12         Betaproteobacteria   NT 173.0031013        16.93110661              Pseudomonadota
# 13 Candidatus_Saccharimonadia  ASD   4.9321835         0.49743205 Candidatus_Saccharibacteria
# 14 Candidatus_Saccharimonadia   NT   4.6867497         0.45867304 Candidatus_Saccharibacteria
# 15                 Chlamydiia   NT  48.9956511         4.79500417                 Chlamydiota
# 16                 Chlamydiia  ASD  50.8109632         5.12450553                 Chlamydiota
# 17             Coriobacteriia  ASD   3.3513477         0.33799792              Actinomycetota
# 18             Coriobacteriia   NT   1.0267977         0.10048850              Actinomycetota
# 19      Epsilonproteobacteria  ASD   7.1769225         0.72382369            Campylobacterota
# 20      Epsilonproteobacteria   NT   5.1101604         0.50011052            Campylobacterota
# 21             Flavobacteriia  ASD   6.6673612         0.67243223                Bacteroidota
# 22             Flavobacteriia   NT   7.2700119         0.71148636                Bacteroidota
# 23        Gammaproteobacteria   NT  85.3136839         8.34930164              Pseudomonadota
# 24        Gammaproteobacteria  ASD 113.6860857        11.46573373              Pseudomonadota
# 25              Negativicutes  ASD  35.4611818         3.57641364                   Bacillota
# 26              Negativicutes   NT   8.5941189         0.84107130                   Bacillota
# 27           Sphingobacteriia  ASD   0.7929756         0.07997502                Bacteroidota
# 28           Sphingobacteriia   NT   0.8612581         0.08428781                Bacteroidota
# 29               Spirochaetia  ASD  80.4216038         8.11086678               Spirochaetota
# 30               Spirochaetia   NT 204.2732196        19.99138531               Spirochaetota

PP2_long2 <- PP2_long2 %>%
  mutate(
    Species = ifelse(
      Phylum %in% c("Actinomycetota", "Pseudomonadota", "Spirochaetota", "Bacteroidota", "Chlamydiota", "Bacillota"),
      Species,
      "Residuals"
    )
  )

PP2_long3 <- PP2_long2 %>%
  group_by(Species, Role) %>%
  summarise(
    Relative_Abundance = sum(Relative_Abundance, na.rm = TRUE),
    .groups = "drop"
  )
PP2_long4=merge(PP2_long3,cate,by="Species",all.x = T)
length(unique(PP2_long4$Species)) # 
class2phylum <- unique(PP2_long4[,c("Species", "Phylum")])
class2phylum <- class2phylum[order(class2phylum$Phylum), ]
class2phylum 
               # Species         Phylum
# 1        Actinomycetes Actinomycetota
# 13      Coriobacteriia Actinomycetota
# 5              Bacilli      Bacillota
# 19       Negativicutes      Bacillota
# 7          Bacteroidia   Bacteroidota
# 15      Flavobacteriia   Bacteroidota
# 23    Sphingobacteriia   Bacteroidota
# 11          Chlamydiia    Chlamydiota
# 3  Alphaproteobacteria Pseudomonadota
# 9   Betaproteobacteria Pseudomonadota
# 17 Gammaproteobacteria Pseudomonadota
# 25        Spirochaetia  Spirochaetota
# 21           Residuals           <NA>

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

species_order <- names(custom_colors1)

PP2_long4$Species <- factor(PP2_long4$Species, levels = species_order)

               # Species Role Relative_Abundance         Phylum
# 1        Actinomycetes  ASD         6.87005201 Actinomycetota
# 2        Actinomycetes   NT         3.92218576 Actinomycetota
# 3  Alphaproteobacteria  ASD        38.96726520 Pseudomonadota
# 4  Alphaproteobacteria   NT        39.00961089 Pseudomonadota
# 5              Bacilli  ASD         0.17284600      Bacillota
# 6              Bacilli   NT         0.26034809      Bacillota
# 7          Bacteroidia  ASD         6.43529646   Bacteroidota
# 8          Bacteroidia   NT         3.89717441   Bacteroidota
# 9   Betaproteobacteria  ASD        16.67949526 Pseudomonadota
# 10  Betaproteobacteria   NT        16.93110661 Pseudomonadota
# 11          Chlamydiia  ASD         5.12450553    Chlamydiota
# 12          Chlamydiia   NT         4.79500417    Chlamydiota
# 13      Coriobacteriia  ASD         0.33799792 Actinomycetota
# 14      Coriobacteriia   NT         0.10048850 Actinomycetota
# 15      Flavobacteriia  ASD         0.67243223   Bacteroidota
# 16      Flavobacteriia   NT         0.71148636   Bacteroidota
# 17 Gammaproteobacteria  ASD        11.46573373 Pseudomonadota
# 18 Gammaproteobacteria   NT         8.34930164 Pseudomonadota
# 19       Negativicutes  ASD         3.57641364      Bacillota
# 20       Negativicutes   NT         0.84107130      Bacillota
# 21           Residuals  ASD         1.50712023           <NA>
# 22           Residuals   NT         1.10654915           <NA>
# 23    Sphingobacteriia  ASD         0.07997502   Bacteroidota
# 24    Sphingobacteriia   NT         0.08428781   Bacteroidota
# 25        Spirochaetia  ASD         8.11086678  Spirochaetota
# 26        Spirochaetia   NT        19.99138531  Spirochaetota

p2 <- ggplot(PP2_long4, aes(x = Role, y = Relative_Abundance, fill = Species)) +
  geom_bar(stat = "identity") +
  labs(x = "", y = "Abundance (%)", fill = "Class", title = "Class") +
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

pa1 <- plot_grid(p1, p2, labels = c("A", ""), label_size = 18, ncol = 2, rel_widths = c(1, 1))

empty_plot <- ggplot() + 
  theme_void()

pa <- plot_grid(pa1, empty_plot, labels = c("", ""), nrow = 2, rel_heights = c(1, 0.55))

# ----------------------------------------------------------- B -----------------------------------------------------------
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


result <- long_data %>%
  filter(Presence == 1) %>%
  group_by(Species, Role) %>%
  summarise(Count = n()) %>%
  spread(key = Role, value = Count, fill = 0)
colnames(result)

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

pb <- ggplot(top30_long, aes(x = reorder(Species, -all), y = Prevalence, fill = Group)) +
  geom_bar(stat = "identity", width = 0.7) +
  scale_fill_manual(values = c("NT" = "#3273ae","ASD" = "#f1761d")) +
  labs(x = "", y = "Prevalence(%)", title = "") +
  scale_y_continuous(limits = c(0, 20), expand = c(0.04, 0)) + 
  scale_x_discrete(expand = c(0.03, 0)) +
  theme(panel.background = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, size = 0.8),
        axis.ticks.length = unit(0.2, "cm"),
        axis.ticks = element_line(size = 0.8, color = "black"),
        axis.text.x = element_text(angle = 55, hjust = 1, vjust = 1,
                                   size = 12, color = "black", face = "italic"),
        axis.text.y  = element_text(size = 12, colour = "black"),
        axis.title.x = element_text(size = 12, color = "black"),
        axis.title.y = element_text(size = 12, color = "black"),
        legend.position = "top",
        plot.margin = unit(c(0, 5, 0, 10), "mm"),
		legend.title = element_text(size = 14, face = "bold", color = "black"),
		legend.text = element_text(size = 12, color = "black")
		)

# ----------------------------------------------------------- C -----------------------------------------------------------
species_counts <- rowSums(aim[, 2:101])
data_for_plot <- data.frame(species_count = species_counts, Role = aim$Role)

plot_data <- data_for_plot %>%
  group_by(species_count, Role) %>%
  summarise(count = n()) %>%
  ungroup()
plot_data$Role <- factor(plot_data$Role, levels = c("ASD", "NT"))

species_counts <- rowSums(aim[, 2:101])

species_counts[species_counts > 20] <- "20+"

data_for_plot <- data.frame(species_count = species_counts, Role = aim$Role)

plot_data <- data_for_plot %>%
  group_by(species_count, Role) %>%
  summarise(count = n(), .groups = "drop")

plot_data$Role <- factor(plot_data$Role, levels = c("ASD", "NT"))

plot_data$species_count <- factor(plot_data$species_count, levels = c(0:20, "20+"))

pc <- ggplot(plot_data, aes(x = species_count, y = count, fill = Role)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = c("ASD" = "#f1761d", "NT" = "#3273ae")) +
  labs(x = "Number of Species", y = "Number of Samples", title = "") +
  scale_y_continuous(limits = c(0, 1100), breaks = seq(0, 1000, 250), expand = c(0.04, 0)) + 
  scale_x_discrete(expand = c(0.035, 0)) +
  theme(panel.background = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, size = 1),
        axis.ticks.length = unit(0.2, "cm"),
        axis.ticks = element_line(size = 0.8),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1,
                                   size = 12, color = "black"),
        axis.text.y  = element_text(size = 12, colour = "black"),
        axis.title.x = element_text(size = 12, color = "black"),
        axis.title.y = element_text(size = 12, color = "black"),
        legend.position = "none")

pbc <- plot_grid(pb, pc, labels = c("B", "C"), label_size = 18, nrow = 2, rel_heights = c(1.2, 0.6), align = "v")

pabc <- plot_grid(pa, pbc, labels = c("", ""), ncol = 2, rel_widths = c(0.8, 1))

# ----------------------------------------------------------- D -----------------------------------------------------------
species_data <- aim[, 2:101]
role_data <- aim[, 102]

species_count <- apply(species_data, 1, function(x) sum(x > 0))

count <- data.frame(SPID = aim$SPID,
                    Role = role_data,
                    Species_Count = species_count)

count %>%
  group_by(Role) %>%
  summarise(Average_Species_Count = mean(Species_Count, na.rm = TRUE))
  # Role  Average_Species_Count
  # <chr>                 <dbl>
# 1 ASD                    5.31
# 2 NT                     5.66

aim$contains_species <- rowSums(aim[, 2:101]) > 0
contingency_table <- table(aim$Role,aim$contains_species )
      # FALSE TRUE
  # ASD   572 1374
  # NT    507 1439
fisher_result <- fisher.test(contingency_table)
#p-value = 0.02189 odds ratio = 1.181523

data_summary <- aim[,-1] %>%
  group_by(Role, contains_species) %>%
  summarise(Count = n(), .groups = "drop") %>%
  mutate(contains_species = ifelse(contains_species == 1, "Yes", "No"))

pd <- ggplot(data_summary, aes(x = Role, y = Count, fill = contains_species)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = c("Yes" = "#f1a89a", "No" = "#D3D3D3")) +
  labs(x = "", y = "Number of samples ", fill = "", title = "") +
  geom_text(aes(label = Count), position = position_stack(vjust = 0.5), size = 4.5) +
  scale_y_continuous(limits = c(0, 2150), breaks = seq(0, 2000, 500)) + 
  theme(panel.background = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, size = 0.8),
        axis.ticks.length = unit(0.2, "cm"),
        axis.ticks = element_line(size = 0.8),
        axis.text = element_text(size = 12, colour = "black"),
        axis.title.x = element_text(size = 12, color = "black"),
        axis.title.y = element_text(size = 12, color = "black"),
		legend.text = element_text(size = 12, color = "black"),
		plot.margin = unit(c(0, 0, 0, 6), "mm"))

# ----------------------------------------------------------- E -----------------------------------------------------------
cate1=read.csv("E:/2025.8-2026.7/iScience_revision/1_rawdata/SSC_kraken_cate_16344.csv",row.names = 1)
colnames(cate1)[colnames(cate1) == "V1"] <- "X"
otu1=merge(otu_raw,cate1,by="X")
row.names(otu1)=otu1$Species
otu1 <- otu1[, -c(2:3893)]
otu1$Species <- gsub("_", " ", otu1$Species)

body=read.csv("E:/2025.8-2026.7/iScience_revision/12_micro_characterization/Body_site/experiment human associate.csv")
colnames(body)
body_part=body[,c("organism_name","Site")]
aa1=merge(otu1,body_part,by.x = "Species",by.y = "organism_name",all.x = T)
length(unique(aa1$Species)) # 100

aa1_1 <- aa1 %>%
  filter(!is.na(Site))
length(unique(aa1_1$Species)) # 23
colnames(aa1_1)
aa1_1=aa1[,c("X","Kingdom","Phylum","Class","Order","Family","Genus","Species","Site")]

aa1_2 <- aa1 %>%
  filter(is.na(Site))
colnames(aa1_2)
aa1_2=aa1_2[,c("X","Kingdom","Phylum","Class","Order","Family","Genus","Species")]
aa2=merge(aa1_2,body_part,by.x = "Genus",by.y = "organism_name",all.x = T)
aa2=aa2[,c("X","Kingdom","Phylum","Class","Order","Family","Genus","Species","Site")]
length(unique(aa2$Species)) # 77

aa3=rbind(aa1_1,aa2)
length(unique(aa3$Species)) # 100

aa4=unique(aa3)
aa4$Site[aa4$Site == ""] <- NA
table(aa4$Site)
# Blood Genitourinary tract                 Gut                Oral   Respiratory tract                Skin 
# 18                  25                  33                  39                   7                  19 

site_counts <- aa4 %>%
  filter(!is.na(Site)) %>% 
  count(Site) %>%  
  arrange(desc(n)) 

site_counts$Site <- factor(site_counts$Site, levels = site_counts$Site)

pe <- ggplot(site_counts, aes(x = Site, y = n, fill = Site)) +
  geom_bar(stat = "identity") +
  labs(x = "Body site", y = "Number of species", title = "") +
  scale_fill_brewer(palette = "Set3") +
  scale_y_continuous(limits = c(0, 42), breaks = seq(0, 42, 10)) + 
  theme_bw() +
  theme(legend.position = "none",
		panel.border = element_rect(color = "black", fill = NA, size = 0.8),
		axis.title.x = element_text(size = 12, color = "black"),
        axis.title.y = element_text(size = 12, color = "black"),
        axis.text.y = element_text(size = 12, color = "black"),
		axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
		axis.ticks.length = unit(0.2, "cm"),
        axis.ticks.y = element_line(size = 0.8),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank())

# ----------------------------------------------------------- F -----------------------------------------------------------
oxygen <- read_xlsx("E:/2025.8-2026.7/iScience_revision/12_micro_characterization/Oxygen_type/Specie_oxygen_type_100.xlsx")
unique(oxygen$Oxygen)
oxygen_counts <- oxygen %>%
  count(Oxygen) %>%
  mutate(percentage = n / sum(n) * 100,
         label = paste(n, "\n", round(percentage, 1), "%", sep = ""))
colnames(oxygen_counts)

df=oxygen_counts %>% 
  arrange(desc(percentage))

pf <- ggplot(df, aes(x = 2, y = n, fill = reorder(Oxygen, -percentage))) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar("y", start = 0) +
  geom_text(aes(label = paste0(round(percentage, 1), "%")), 
            position = position_stack(vjust = 0.5), size = 4.5) +
  scale_fill_brewer(palette = "Set3") +
  theme_void() +
  theme(
    legend.title = element_text(size = 14, face = "bold", color = "black"),
    legend.text = element_text(size = 12, color = "black"),
	plot.margin = unit(c(0, 4, 0, 0), "mm")
  ) +
  xlim(0.5, 2.5) +
  labs(title = "", fill = "Oxygen Type")

pdef <- plot_grid(pd, pe, pf, labels = c("D", "E", "F"), label_size = 18, ncol = 3, rel_widths = c(1.1, 1.5, 1.8), align = "h")

pfinal <- plot_grid(pabc, pdef, labels = c("", ""), nrow = 2, rel_heights = c(2.2, 1))
ggsave("Figure 2.pdf", device = "pdf", path = "E:/2025.8-2026.7/iScience_revision/Figures", width=12, height=12)