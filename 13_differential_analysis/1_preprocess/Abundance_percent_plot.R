####------------------------------------------------------- (13).differential analysis -------------------------------------------------------
############################################################ 13.1 preprocess ############################################################
rm(list=ls())
library(data.table)
library(dplyr)
library(reshape2)
library(ggplot2)
library(RColorBrewer)
library(tidyr)

# Phylum --------------------------------------------------------------------
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

## 1.Differential analysis of Phyla
species_data <- PP2[, 1:9]
Phylum_results <- data.frame(Species = colnames(species_data), 
                             p_value = numeric(ncol(species_data)))
for (i in 1:ncol(species_data)) {
  temp <- PP2[,c(i,10,11)]
  temp_asd <- temp %>% filter(Role == "ASD")
  temp_nt <- temp %>% filter(Role == "NT")
  paired_data <- merge(temp_asd, temp_nt, by = "SFID", suffixes = c("_ASD", "_NT"))
  paired_test <- wilcox.test(paired_data[,2], paired_data[,4], paired = TRUE, exact = FALSE, na.action = na.omi)
  Phylum_results$p_value[i] <- paired_test$p.value
}

Phylum_results$q_value <- p.adjust(Phylum_results$p_value, method = "BH")
Phylum_results$q_value <- ifelse(Phylum_results$q_value > 0.05, NA, Phylum_results$q_value)
Phylum_results=Phylum_results[order(Phylum_results$q_value), ]
Phylum_results
                      # Species      p_value      q_value
# 3                Bacteroidota 1.283958e-10 1.155562e-09
# 9               Spirochaetota 6.344917e-09 2.855213e-08
# 4            Bdellovibrionota 2.935944e-03 8.807831e-03
# 1              Actinomycetota 1.479790e-02 3.329528e-02
# 2                   Bacillota 4.211746e-01           NA
# 5            Campylobacterota 7.952131e-01           NA
# 6 Candidatus_Saccharibacteria 9.321069e-01           NA
# 7                 Chlamydiota 4.385378e-02           NA
# 8              Pseudomonadota 1.559011e-01           NA

## 2.Check the relative distribution of all Phyla and sort them
PP2 <- PP2[,c(1:10)]
PP4 <- PP2
PP4$Role <- "All"
PP2_long_all1 <- PP4 %>%
  pivot_longer(-Role, names_to = "Species", values_to = "Abundance") %>%
  group_by(Role, Species) %>%
  summarise(Abundance = mean(Abundance, na.rm = TRUE), .groups = "drop") %>%
  group_by(Role) %>%
  mutate(Relative_Abundance = Abundance / sum(Abundance) * 100)

sum(PP2_long_all1$Relative_Abundance) # 100

PP2_long_all1[order(PP2_long_all1$Relative_Abundance),]
  # Role  Species                     Abundance Relative_Abundance
  # <chr> <chr>                           <dbl>              <dbl>
# 1 All   Bdellovibrionota                 2.17              0.216
# 2 All   Candidatus_Saccharibacteria      4.81              0.478
# 3 All   Campylobacterota                 6.14              0.610
# 4 All   Bacillota                       24.2               2.41 
# 5 All   Chlamydiota                     49.9               4.96 
# 6 All   Actinomycetota                  56.3               5.59 
# 7 All   Bacteroidota                    59.6               5.92 
# 8 All   Spirochaetota                  142.               14.1  
# 9 All   Pseudomonadota                 661.               65.7  
colnames(PP2_long_all1) <- c("All", "Species", "Abundance_all", "Relative_Abundance_all")

## 3.Check the relative distribution of all Phyla in different groups and sort them
PP2_long_group <- PP2 %>%
  pivot_longer(-Role, names_to = "Species", values_to = "Abundance") %>%
  group_by(Role, Species) %>%
  summarise(Abundance = mean(Abundance, na.rm = TRUE), .groups = "drop") %>%
  group_by(Role) %>%
  mutate(Relative_Abundance = Abundance / sum(Abundance) * 100)

sum(PP2_long_group$Relative_Abundance) # 200

PP2_long_group
   # Role  Species                     Abundance Relative_Abundance
   # <chr> <chr>                           <dbl>              <dbl>
 # 1 ASD   Actinomycetota                  71.5               7.21 
 # 2 ASD   Bacillota                       37.2               3.75 
 # 3 ASD   Bacteroidota                    71.3               7.19 
 # 4 ASD   Bdellovibrionota                 2.83              0.286
 # 5 ASD   Campylobacterota                 7.18              0.724
 # 6 ASD   Candidatus_Saccharibacteria      4.93              0.497
 # 7 ASD   Chlamydiota                     50.8               5.12 
 # 8 ASD   Pseudomonadota                 665.               67.1  
 # 9 ASD   Spirochaetota                   80.4               8.11 
# 10 NT    Actinomycetota                  41.1               4.02 
# 11 NT    Bacillota                       11.3               1.10 
# 12 NT    Bacteroidota                    48.0               4.69 
# 13 NT    Bdellovibrionota                 1.51              0.148
# 14 NT    Campylobacterota                 5.11              0.500
# 15 NT    Candidatus_Saccharibacteria      4.69              0.459
# 16 NT    Chlamydiota                     49.0               4.80 
# 17 NT    Pseudomonadota                 657.               64.3  
# 18 NT    Spirochaetota                  204.               20.0  

PP2_long_ASD <- PP2_long_group[PP2_long_group$Role == "ASD",]
colnames(PP2_long_ASD) <- c("ASD", "Species",  "Abundance_ASD", "Relative_Abundance_ASD")
PP2_long_NT <- PP2_long_group[PP2_long_group$Role == "NT",]
colnames(PP2_long_NT) <- c("NT", "Species",  "Abundance_NT", "Relative_Abundance_NT")
PP2_long_group_new <- merge(PP2_long_ASD, PP2_long_NT, by = "Species")
PP2_long_final <- merge(merge(PP2_long_group_new, Phylum_results, by = "Species"), PP2_long_all1, by = "Species")[,c(-2, -5, -10)]
PP2_long_final
                      # Species Abundance_ASD Relative_Abundance_ASD Abundance_NT Relative_Abundance_NT      p_value
# 1              Actinomycetota     71.469912              7.2080499    41.103936             4.0226743 1.479790e-02
# 2                   Bacillota     37.175000              3.7492596    11.254372             1.1014194 4.211746e-01
# 3                Bacteroidota     71.268173              7.1877037    47.952841             4.6929486 1.283958e-10
# 4            Bdellovibrionota      2.834429              0.2858645     1.509878             0.1477656 2.935944e-03
# 5            Campylobacterota      7.176923              0.7238237     5.110160             0.5001105 7.952131e-01
# 6 Candidatus_Saccharibacteria      4.932184              0.4974321     4.686750             0.4586730 9.321069e-01
# 7                 Chlamydiota     50.810963              5.1245055    48.995651             4.7950042 4.385378e-02
# 8              Pseudomonadota    665.439905             67.1124942   656.919418            64.2900191 1.559011e-01
# 9               Spirochaetota     80.421604              8.1108668   204.273220            19.9913853 6.344917e-09
       # q_value Abundance_all Relative_Abundance_all
# 1 3.329528e-02     56.286924              5.5914108
# 2           NA     24.214686              2.4054300
# 3 1.155562e-09     59.610507              5.9215677
# 4 8.807831e-03      2.172154              0.2157766
# 5           NA      6.143541              0.6102850
# 6           NA      4.809467              0.4777611
# 7           NA     49.903307              4.9572773
# 8           NA    661.179661             65.6800341
# 9 2.855213e-08    142.347412             14.1404574

write.table(PP2_long_final, "E:/2025.8-2026.7/iScience_revision/13_differential_analysis/1_preprocess/Phylum_results_all_paired.tsv", sep = "\t", quote = F, row.names = F)

# 4.Based on the statistical results, the three Phyla with a relative abundance less than 1% are combined and classified as Others
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

PP3 <- PP2
PP3$Role <- "All"
PP2_long_all <- PP3 %>%
  pivot_longer(-Role, names_to = "Species", values_to = "Abundance") %>%
  mutate(Species = ifelse(Species %in% names(custom_colors), Species, "Residuals")) %>%
  group_by(Role, Species) %>%
  summarise(Abundance = mean(Abundance, na.rm = TRUE), .groups = "drop") %>%
  group_by(Role) %>%
  mutate(Relative_Abundance = Abundance / sum(Abundance) * 100)

sum(PP2_long_all$Relative_Abundance) # 100

# The relative abundance distribution of all Phyla after merging the last three Phyla in terms of abundance ranking
PP2_long_all[order(PP2_long_all$Relative_Abundance),]
  # Role  Species        Abundance Relative_Abundance
  # <chr> <chr>              <dbl>              <dbl>
# 1 All   Residuals           4.38              0.438
# 2 All   Bacillota          24.2               2.43 
# 3 All   Chlamydiota        49.9               5.00 
# 4 All   Actinomycetota     56.3               5.64 
# 5 All   Bacteroidota       59.6               5.97 
# 6 All   Spirochaetota     142.               14.3  
# 7 All   Pseudomonadota    661.               66.3 

# 5.Visualization
PP2_long <- PP2 %>%
  pivot_longer(-Role, names_to = "Species", values_to = "Abundance") %>%
  mutate(Species = ifelse(Species %in% names(custom_colors), Species, "Residuals")) %>%
  group_by(Role, Species) %>%
  summarise(Abundance = mean(Abundance, na.rm = TRUE), .groups = "drop") %>%
  group_by(Role) %>%
  mutate(Relative_Abundance = Abundance / sum(Abundance) * 100)

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

# Calculate the species abundance ranking (including "Others") in the ASD group
species_order <- PP2_long %>%
  filter(Role == "ASD") %>%
  group_by(Species) %>%
  summarise(Total_Abundance = sum(Relative_Abundance, na.rm = TRUE)) %>%
  arrange(Total_Abundance) %>%
  pull(Species)

# Set the sequence of species factors according to the ASD abundance order
PP2_long$Species <- factor(PP2_long$Species, levels = species_order)

PP2_long0 <- PP2_long

ggplot(PP2_long, aes(x = Role, y = Relative_Abundance, fill = Species)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(x = "", y = "Abundance (%)", fill = "Phylum", title = "Phylum") +
  scale_fill_manual(values = custom_colors)+
  scale_y_continuous(limits = c(0, 100), 
                       breaks = seq(0, 100, 20), labels = c("0", "20", "40", "60", "80", "100")) + 
  theme(plot.title = element_text(hjust = 0.5, vjust = 0.5,
								   face = "bold", size = 14, color = "black"),
		panel.background = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, size = 0.8),
        axis.ticks.length = unit(0.2, "cm"),
        axis.ticks = element_line(size = 0.8, colour = "black"),
        axis.text = element_text(size = 12, colour = "black"),
        axis.title.x = element_text(size = 12, colour = "black"),
        axis.title.y = element_text(size = 12, colour = "black"),
		plot.margin = unit(c(1, 3, 0, 6), "mm"))

ggsave("Abundance_percent_Phylum.pdf", device = "pdf", path = "E:/2025.8-2026.7/iScience_revision/13_differential_analysis/1_preprocess/", width=4.2,height=4.8)

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

## 1.Differential analysis of Classes
species_data <- PP2[, 1:15]
Class_results <- data.frame(Class = colnames(species_data), 
                             p_value = numeric(ncol(species_data)))
for (i in 1:ncol(species_data)) {
  temp <- PP2[,c(i,16,17)]
  temp_asd <- temp %>% filter(Role == "ASD")
  temp_nt <- temp %>% filter(Role == "NT")
  paired_data <- merge(temp_asd, temp_nt, by = "SFID", suffixes = c("_ASD", "_NT"))
  paired_test <- wilcox.test(paired_data[,2], paired_data[,4], paired = TRUE, exact = FALSE, na.action = na.omi)
  Class_results$p_value[i] <- paired_test$p.value
}

Class_results$q_value <- p.adjust(Class_results$p_value, method = "BH")
Class_results$q_value <- ifelse(Class_results$q_value > 0.05, NA, Class_results$q_value)
Class_results=Class_results[order(Class_results$q_value), ]
Class_results
                        # Class      p_value      q_value
# 15               Spirochaetia 6.344917e-09 9.517375e-08
# 4                 Bacteroidia 1.349269e-06 1.011952e-05
# 2         Alphaproteobacteria 3.861197e-04 1.930599e-03
# 11             Flavobacteriia 1.251551e-03 4.693317e-03
# 5             Bdellovibrionia 2.935944e-03 8.807831e-03
# 1               Actinomycetes 1.503556e-02 3.758890e-02
# 3                     Bacilli 3.013984e-01           NA
# 6          Betaproteobacteria 8.786507e-01           NA
# 7  Candidatus_Saccharimonadia 9.321069e-01           NA
# 8                  Chlamydiia 4.385378e-02           NA
# 9              Coriobacteriia 4.016782e-01           NA
# 10      Epsilonproteobacteria 7.952131e-01           NA
# 12        Gammaproteobacteria 2.886679e-02           NA
# 13              Negativicutes 6.012497e-01           NA
# 14           Sphingobacteriia 5.648415e-01           NA
colnames(Class_results)[1] <- "Class"
cate=unique(otu[,c(3895:3896)])
Class_results2=merge(Class_results,cate,by="Class")
colnames(Class_results2)
Class_results2=Class_results2[,c("Phylum","Class","p_value","q_value")]
Class_results2=Class_results2[order(Class_results2$q_value), ]
Class_results2
                        # Phylum                      Class      p_value      q_value
# 15               Spirochaetota               Spirochaetia 6.344917e-09 9.517375e-08
# 4                 Bacteroidota                Bacteroidia 1.349269e-06 1.011952e-05
# 2               Pseudomonadota        Alphaproteobacteria 3.861197e-04 1.930599e-03
# 11                Bacteroidota             Flavobacteriia 1.251551e-03 4.693317e-03
# 5             Bdellovibrionota            Bdellovibrionia 2.935944e-03 8.807831e-03
# 1               Actinomycetota              Actinomycetes 1.503556e-02 3.758890e-02
# 3                    Bacillota                    Bacilli 3.013984e-01           NA
# 6               Pseudomonadota         Betaproteobacteria 8.786507e-01           NA
# 7  Candidatus_Saccharibacteria Candidatus_Saccharimonadia 9.321069e-01           NA
# 8                  Chlamydiota                 Chlamydiia 4.385378e-02           NA
# 9               Actinomycetota             Coriobacteriia 4.016782e-01           NA
# 10            Campylobacterota      Epsilonproteobacteria 7.952131e-01           NA
# 12              Pseudomonadota        Gammaproteobacteria 2.886679e-02           NA
# 13                   Bacillota              Negativicutes 6.012497e-01           NA
# 14                Bacteroidota           Sphingobacteriia 5.648415e-01           NA

## 2.Check the relative distribution of all Phyla and sort them
PP2 <- PP2[,c(1:16)]
PP4 <- PP2
PP4$Role <- "All"
PP2_long_all1 <- PP4 %>%
  pivot_longer(-Role, names_to = "Class", values_to = "Abundance") %>%
  group_by(Role, Class) %>%
  summarise(Abundance = mean(Abundance, na.rm = TRUE), .groups = "drop") %>%
  group_by(Role) %>%
  mutate(Relative_Abundance = Abundance / sum(Abundance) * 100)

sum(PP2_long_all1$Relative_Abundance) # 100

PP2_long_all1[order(PP2_long_all1$Relative_Abundance),]
   # Role  Class                      Abundance Relative_Abundance
   # <chr> <chr>                          <dbl>              <dbl>
 # 1 All   Sphingobacteriia               0.827             0.0822
 # 2 All   Bdellovibrionia                2.17              0.216 
 # 3 All   Bacilli                        2.19              0.217 
 # 4 All   Coriobacteriia                 2.19              0.217 
 # 5 All   Candidatus_Saccharimonadia     4.81              0.478 
 # 6 All   Epsilonproteobacteria          6.14              0.610 
 # 7 All   Flavobacteriia                 6.97              0.692 
 # 8 All   Negativicutes                 22.0               2.19  
 # 9 All   Chlamydiia                    49.9               4.96  
# 10 All   Bacteroidia                   51.8               5.15  
# 11 All   Actinomycetes                 54.1               5.37  
# 12 All   Gammaproteobacteria           99.5               9.88  
# 13 All   Spirochaetia                 142.               14.1   
# 14 All   Betaproteobacteria           169.               16.8   
# 15 All   Alphaproteobacteria          392.               39.0 

## 3.Check the relative distribution of all Classes in different groups
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

PP2_long2_ASD <- PP2_long2[PP2_long2$Role == "ASD",]
colnames(PP2_long2_ASD) <- c("Class", "ASD", "Abundance_ASD", "Relative_Abundance_ASD", "Phylum")
PP2_long2_NT <- PP2_long2[PP2_long2$Role == "NT",]
colnames(PP2_long2_NT) <- c("Class", "NT", "Abundance_NT", "Relative_Abundance_NT", "Phylum")
PP2_long2_group_new <- merge(PP2_long2_ASD[,-2], PP2_long2_NT[,c(-2, -5)], by = "Class")
PP2_long2_final <- merge(merge(PP2_long2_group_new, Class_results2[,-1], by = "Class"), PP2_long_all1[,-1], by = "Class")
PP2_long2_final
                        # Class Abundance_ASD Relative_Abundance_ASD                      Phylum Abundance_NT
# 1               Actinomycetes    68.1185644             6.87005201              Actinomycetota   40.0771383
# 2         Alphaproteobacteria   386.3717713            38.96726520              Pseudomonadota  398.6026325
# 3                     Bacilli     1.7138184             0.17284600                   Bacillota    2.6602530
# 4                 Bacteroidia    63.8078367             6.43529646                Bacteroidota   39.8215708
# 5             Bdellovibrionia     2.8344295             0.28586448            Bdellovibrionota    1.5098779
# 6          Betaproteobacteria   165.3820481            16.67949526              Pseudomonadota  173.0031013
# 7  Candidatus_Saccharimonadia     4.9321835             0.49743205 Candidatus_Saccharibacteria    4.6867497
# 8                  Chlamydiia    50.8109632             5.12450553                 Chlamydiota   48.9956511
# 9              Coriobacteriia     3.3513477             0.33799792              Actinomycetota    1.0267977
# 10      Epsilonproteobacteria     7.1769225             0.72382369            Campylobacterota    5.1101604
# 11             Flavobacteriia     6.6673612             0.67243223                Bacteroidota    7.2700119
# 12        Gammaproteobacteria   113.6860857            11.46573373              Pseudomonadota   85.3136839
# 13              Negativicutes    35.4611818             3.57641364                   Bacillota    8.5941189
# 14           Sphingobacteriia     0.7929756             0.07997502                Bacteroidota    0.8612581
# 15               Spirochaetia    80.4216038             8.11086678               Spirochaetota  204.2732196
   # Relative_Abundance_NT      p_value      q_value   Abundance Relative_Abundance
# 1             3.92218576 1.503556e-02 3.758890e-02  54.0978513         5.37395344
# 2            39.00961089 3.861197e-04 1.930599e-03 392.4872019        38.98875645
# 3             0.26034809 3.013984e-01           NA   2.1870357         0.21725499
# 4             3.89717441 1.349269e-06 1.011952e-05  51.8147037         5.14715092
# 5             0.14776558 2.935944e-03 8.807831e-03   2.1721537         0.21577664
# 6            16.93110661 8.786507e-01           NA 169.1925747        16.80719284
# 7             0.45867304 9.321069e-01           NA   4.8094666         0.47776111
# 8             4.79500417 4.385378e-02           NA  49.9033072         4.95727728
# 9             0.10048850 4.016782e-01           NA   2.1890727         0.21745734
# 10            0.50011052 7.952131e-01           NA   6.1435415         0.61028498
# 11            0.71148636 1.251551e-03 4.693317e-03   6.9686865         0.69225295
# 12            8.34930164 2.886679e-02           NA  99.4998848         9.88408477
# 13            0.84107130 6.012497e-01           NA  22.0276503         2.18817503
# 14            0.08428781 5.648415e-01           NA   0.8271168         0.08216384
# 15           19.99138531 6.344917e-09 9.517375e-08 142.3474117        14.14045742

write.table(PP2_long2_final, "E:/2025.8-2026.7/iScience_revision/13_differential_analysis/1_preprocess/Class_results_all_paired.tsv", sep = "\t", quote = F, row.names = F)

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

# Define Colors
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

# Create the factor order of the species in the order of custom_colors1
species_order <- names(custom_colors1)

# Set the species factor order in the order of custom_colors1
PP2_long4$Species <- factor(PP2_long4$Species, levels = species_order)
PP2_long4
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

ggplot(PP2_long4, aes(x = Role, y = Relative_Abundance, fill = Species)) +
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
		plot.margin = unit(c(1, 3, 0, 6), "mm"))

ggsave("Abundance_percent_Class.pdf", device = "pdf", path = "E:/2025.8-2026.7/iScience_revision/13_differential_analysis/1_preprocess/", width=4.2,height=4.8)