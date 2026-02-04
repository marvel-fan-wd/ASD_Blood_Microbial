rm(list=ls())
library(data.table)
library(stringr)

ff <- as.data.frame(fread("E:/2025.8-2026.7/iScience_revision/1_rawdata/SSC_rawdata.txt", header = TRUE))
colnames(ff)[1]=c("V1")
row_names <- ff$V1

extract_taxonomy_levels <- function(name) {
  levels <- strsplit(name, "\\|")[[1]]
  taxonomy <- list(Kingdom = NA, Phylum = NA, Class = NA, Order = NA, Family = NA, Genus = NA, Species = NA)
  for (level in levels) {
    level_type <- substring(level, 1, 1)
    level_name <- substring(level, 4)
    if (level_type == "k") {
      taxonomy$Kingdom <- ifelse(is.na(taxonomy$Kingdom), level_name, paste(taxonomy$Kingdom, level_name, sep = "|"))
    } else if (level_type == "p") {
      taxonomy$Phylum <- level_name
    } else if (level_type == "c") {
      taxonomy$Class <- level_name
    } else if (level_type == "o") {
      taxonomy$Order <- level_name
    } else if (level_type == "f") {
      taxonomy$Family <- level_name
    } else if (level_type == "g") {
      taxonomy$Genus <- level_name
    } else if (level_type == "s") {
      taxonomy$Species <- level_name
    }
  }
  return(unlist(taxonomy))
}

taxonomy_df <- do.call(rbind, lapply(row_names, extract_taxonomy_levels))
cate <- as.data.frame(taxonomy_df, stringsAsFactors = FALSE)
cate$V1 <- row_names
fwrite(cate, "E:/2025.8-2026.7/iScience_revision/1_rawdata/SSC_cate.csv", row.names=T)

Bact <- grep("k__Bacteria", ff$V1)
ff_Bacteria <- ff[Bact,]
Vir <- grep("k__Viruses", ff$V1)
ff_Viruses <- ff[Vir,]

Archae <- grep("k__Archaea", ff$V1)
ff_Archaea <- ff[Archae,]

Eukaryota_1 <- grep("k__Eukaryota", ff$V1)
Eukaryota <- ff[Eukaryota_1,]
Eukaryota_2 <- grep("k__Fungi", Eukaryota$V1)
Eukaryota_Fungi <- Eukaryota[Eukaryota_2,]

ff21 <- rbind(ff_Bacteria,ff_Viruses,ff_Archaea,Eukaryota_Fungi)
spp1 <- grep("s__", ff21$V1)
specie1 <- ff21[spp1,]
specie21=merge(cate,specie1,by="V1")
length(unique(specie21$Genus))
length(unique(specie21$Species))

ff2 <- rbind(ff_Bacteria)
spp <- grep("s__", ff2$V1)
specie <- ff2[spp,]
specie2=merge(cate,specie,by="V1")
length(unique(specie2$Genus))
length(unique(specie2$Species))

read <- specie2[, c(1, 9:3900)]
row.names(read) <- read$V1
read <- read[, -1]
fwrite(read,"E:/2025.8-2026.7/iScience_revision/SSC_specie1.csv",row.names=T)
