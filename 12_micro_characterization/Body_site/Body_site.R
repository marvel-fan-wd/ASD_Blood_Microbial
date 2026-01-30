####------------------------------------------------------- (12).micro characterization -------------------------------------------------------
## 12.1 Body site
rm(list=ls())
library(ggplot2)
library(dplyr)

otu_raw <- read.csv("E:/2025.8-2026.7/iScience_revision/10_prevalence_filter/SSC_specieQC_prevalence_filter_0.2_100.csv")
cate=read.csv("E:/2025.8-2026.7/iScience_revision/1_rawdata/SSC_kraken_cate_16344.csv",row.names = 1)
colnames(cate)[colnames(cate) == "V1"] <- "X"
otu=merge(otu_raw,cate,by="X")
row.names(otu)=otu$Species
otu <- otu[, -c(2:3893)]
otu$Species <- gsub("_", " ", otu$Species)

body=read.csv("E:/2025.8-2026.7/iScience_revision/12_micro_characterization/Body_site/experiment human associate.csv")
colnames(body)
body_part=body[,c("organism_name","Site")]
aa1=merge(otu,body_part,by.x = "Species",by.y = "organism_name",all.x = T)
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

# Count the number of each variable in the Site column and remove NA
site_counts <- aa4 %>%
  filter(!is.na(Site)) %>% 
  count(Site) %>%  
  arrange(desc(n)) 
site_counts$Site <- factor(site_counts$Site, levels = site_counts$Site)

# Visualization
ggplot(site_counts, aes(x = Site, y = n, fill = Site)) +
  geom_bar(stat = "identity") +
  labs(x = "Body site", y = "Number of species", title = "") +
  scale_fill_brewer(palette = "Set3") +
  theme_minimal() +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("E:/2025.8-2026.7/iScience_revision/12_micro_characterization/Body_site/Body_site_100.pdf", device = "pdf", width=4.5,height=4)