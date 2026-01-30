####------------------------------------------------------- (1).extract specie -------------------------------------------------------
rm(list=ls())
library(data.table)
library(stringr)

## 1.Decomposition of species level
ff <- as.data.frame(fread("E:/2025.8-2026.7/iScience_revision/1_rawdata/SSC_kraken_merged_count_3892_20260109.txt", header = TRUE))
colnames(ff)[1]=c("V1")
row_names <- ff$V1

# Define a function for extracting species information
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
fwrite(cate, "E:/2025.8-2026.7/iScience_revision/1_rawdata/SSC_kraken_cate_16344.csv", row.names=T)

## 2.Extract species at the species level
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
# Species information
ff21 <- rbind(ff_Bacteria,ff_Viruses,ff_Archaea,Eukaryota_Fungi)
spp1 <- grep("s__", ff21$V1)
specie1 <- ff21[spp1,]
specie21=merge(cate,specie1,by="V1")
length(unique(specie21$Genus)) # 2729
length(unique(specie21$Species)) # 12146

ff2 <- rbind(ff_Bacteria)
spp <- grep("s__", ff2$V1)
specie <- ff2[spp,]
specie2=merge(cate,specie,by="V1")
length(unique(specie2$Genus)) # 2110
length(unique(specie2$Species)) # 10841

read <- specie2[, c(1, 9:3900)]
row.names(read) <- read$V1
read <- read[, -1]
fwrite(read,"E:/2025.8-2026.7/iScience_revision/1_rawdata/SSC_specie_10841_count.csv",row.names=T)

## 3.Reads distribution of different kingdoms
ff_summary <- ff[ff$V1 %in% c("k__Bacteria", "k__Viruses", "k__Archaea", "k__Eukaryota", "k__Eukaryota|k__Metazoa"), ]
meta <- read.csv("E:/2025.8-2026.7/iScience_revision/1_rawdata/SSC_sample_phenotype_age&bmi_impute.csv",header = T)
# NT
sample_ids_s1 <- meta$SPID[meta$Role == "s1"]
ff_summary_s1 <- ff_summary[, c("V1", sample_ids_s1)]
ff_summary_s1[c(6,7,8),] <- 0
ff_summary_s1[c(6,7,8), 1] <- c("k__Fungi", "Bacteria_percent", "Viruses_Archaea_Fungi_percent")
for (i in 2:ncol(ff_summary_s1)){
  ff_summary_s1[6,i] <- ff_summary_s1[2,i] - ff_summary_s1[3,i]
  ff_summary_s1[7,i] <- ff_summary_s1[1,i] / (ff_summary_s1[1,i] + ff_summary_s1[2,i] + ff_summary_s1[4,i] + ff_summary_s1[5,i])
  ff_summary_s1[8,i] <- (ff_summary_s1[4,i] + ff_summary_s1[5,i] + ff_summary_s1[6,i]) / (ff_summary_s1[1,i] + ff_summary_s1[2,i] + ff_summary_s1[4,i] + ff_summary_s1[5,i])
}
ff_summary_s1$mean <- rowSums(ff_summary_s1[,-1])/(ncol(ff_summary_s1) - 1)
fwrite(ff_summary_s1,"E:/2025.8-2026.7/iScience_revision/1_rawdata/SSC_specie_sbling_kingdom.csv",row.names=T)

# ASD
sample_ids_p1 <- meta$SPID[meta$Role == "p1"]
ff_summary_p1 <- ff_summary[, c("V1", sample_ids_p1)]
ff_summary_p1[c(6,7,8),] <- 0
ff_summary_p1[c(6,7,8), 1] <- c("k__Fungi", "Bacteria_percent", "Viruses_Archaea_Fungi_percent")
for (i in 2:ncol(ff_summary_p1)){
  ff_summary_p1[6,i] <- ff_summary_p1[2,i] - ff_summary_p1[3,i]
  ff_summary_p1[7,i] <- ff_summary_p1[1,i] / (ff_summary_p1[1,i] + ff_summary_p1[2,i] + ff_summary_p1[4,i] + ff_summary_p1[5,i])
  ff_summary_p1[8,i] <- (ff_summary_p1[4,i] + ff_summary_p1[5,i] + ff_summary_p1[6,i]) / (ff_summary_p1[1,i] + ff_summary_p1[2,i] + ff_summary_p1[4,i] + ff_summary_p1[5,i])
}
ff_summary_p1$mean <- rowSums(ff_summary_p1[,-1])/(ncol(ff_summary_p1) - 1)
fwrite(ff_summary_p1,"E:/2025.8-2026.7/iScience_revision/1_rawdata/SSC_specie_proband_kingdom.csv",row.names=T)

####------------------------------------------------------- (2).abundance filter -------------------------------------------------------
rm(list=ls())
library(data.table)

## 1.Remove samples with microbial reads less than 100
read <- read.csv("E:/2025.8-2026.7/iScience_revision/1_rawdata/SSC_specie_10841_count.csv", row.names = 1)
read2 <- read[, colSums(read) > 100]

## 2.Set read pairs ¡Ü 10 as 0
read2[read2 <= 10] <- 0 

## 3.Calculate the relative abundance and remove the species with a relative abundance of ¡Ü 0.005
read2 <- sweep(read2, 2, colSums(read2), "/")
read2[read2 < 0.005] <- 0
read2 <- read2[rowSums(read2) > 0, ]

## 4.Summarize the number of species at specie level
cate=read.csv("E:/2025.8-2026.7/iScience_revision/1_rawdata/SSC_kraken_cate_16344.csv")
read2_cate=as.data.frame(rownames(read2))
colnames(read2_cate)=c("V1")
read2_cate1=merge(read2_cate,cate,by="V1")
length(unique(read2_cate1$Genus)) # 182
length(unique(read2_cate1$Species)) # 595

## 5.Restore read matrix
read3 <- read[rownames(read2), , drop = FALSE]
read3[read3 <= 10] <- 0
fwrite(read3, "E:/2025.8-2026.7/iScience_revision/2_abundance_filter/SSC_specie_595_count.csv",row.names=T)

####------------------------------------------------------- (3).prevalence prepare -------------------------------------------------------
rm(list=ls())
library(data.table)

read <- read.csv("E:/2025.8-2026.7/iScience_revision/2_abundance_filter/SSC_specie_595_count.csv", row.names = 1)
batch=read.csv("E:/2025.8-2026.7/iScience_revision/1_rawdata/SSC_sample_phenotype_age&bmi_impute.csv",header=T)
batch=batch[,c("SPID","Phase")];table(batch$Phase)
sample_ids <- batch$SPID;batch_info <- batch$Phase

## 1.Generate multiple matrices based on batch information and name them, while removing rows and those with a value of 0
temp_list <- list()
for (batch_id in unique(batch_info)) {
  batch_samples <- sample_ids[batch_info == batch_id]
  batch_matrix <- read[, colnames(read) %in% batch_samples, drop = FALSE]
  batch_matrix <- batch_matrix[rowSums(batch_matrix) > 0, ]
  temp_list[[as.character(batch_id)]] <- batch_matrix
}

# Convert each element in the list to a separate data.frame
list2env(temp_list, envir = .GlobalEnv)

fwrite(Phase1,"E:/2025.8-2026.7/iScience_revision/3_prevalence_prepare/SSC_specieQC_Phase_1.csv",row.names=T)
fwrite(Phase2,"E:/2025.8-2026.7/iScience_revision/3_prevalence_prepare/SSC_specieQC_Phase_2.csv",row.names=T)
fwrite(`Phase3-1`,"E:/2025.8-2026.7/iScience_revision/3_prevalence_prepare/SSC_specieQC_Phase_31.csv",row.names=T)
fwrite(`Phase3-2`,"E:/2025.8-2026.7/iScience_revision/3_prevalence_prepare/SSC_specieQC_Phase_32.csv",row.names=T)
fwrite(Phase4,"E:/2025.8-2026.7/iScience_revision/3_prevalence_prepare/SSC_specieQC_Phase_4.csv",row.names=T)
fwrite(Pilot,"E:/2025.8-2026.7/iScience_revision/3_prevalence_prepare/SSC_specieQC_Phase_p.csv",row.names=T)

## 2.Calculate the detection rate
output_dir2 <- "E:/2025.8-2026.7/iScience_revision/3_prevalence_prepare/"
output_dir3 <- "E:/2025.8-2026.7/iScience_revision/4_batch_filter_prepare/"

batch_matrices <- list(Phase1, Phase2, `Phase3-1`, `Phase3-2`, Phase4, Pilot)
batch_names <- c("Phase_1", "Phase_2", "Phase_31", "Phase_32", "Phase_4", "Phase_p")

for (i in seq_along(batch_matrices)) {
  ff2_spp <- batch_matrices[[i]]
  batch_name <- batch_names[i]
  # Check the number of samples. If there are less than 100 samples, skip this batch
  if (ncol(ff2_spp) < 100) {
    message(paste("Skipping", batch_name, "- fewer than 100 samples"))
    next
  }
  detect <- ff2_spp
  detect[detect > 0] <- 1
  detect["sum"] <- rowSums(detect)
  detect["detect_ratio"] <- detect$sum / ncol(detect)
  detect_ratio <- data.frame(
    species = row.names(detect),
    sum = detect$sum,
    all = ncol(detect),
    detect_ratio = detect$detect_ratio
  )
  
  output_detect_ratio <- paste0(output_dir3, "SSC_specieQC_prevalence_", batch_name, ".csv")
  fwrite(detect_ratio, output_detect_ratio, row.names = FALSE)
}

####------------------------------------------------------------ (4).batch filter prepare------------------------------------------------------------
## 1.Find species that occur in more than one batch
file4 <- list.files("E:/2025.8-2026.7/iScience_revision/4_batch_filter_prepare/", pattern = "SSC_specieQC_.*\\.csv")
BATCH <- NULL
for (ss in 1:length(file4)) {
  SS <- file4[ss]
  bat <- str_split(SS, "_")[[1]][5]
  bat <- gsub("\\.csv", "", bat)
  detect <- read.csv(paste0("E:/2025.8-2026.7/iScience_revision/4_batch_filter_prepare/SSC_specieQC_prevalence_Phase_", bat, ".csv"), row.names = 1, check.names = FALSE)
  detect[detect > 0] <- 1
  detect["sum"] <- rowSums(detect)
  detect["detect_ratio"] <- detect$sum / ncol(detect)
  detect_ratio <- data.frame(
    species = row.names(detect),
    sum = detect$sum,
    all = ncol(detect),
    detect_ratio = detect$detect_ratio
  )
  colnames(detect_ratio)[2:4] <- paste(colnames(detect_ratio)[2:4], bat, sep = "_")
  
  if (is.null(BATCH)) {
    BATCH <- detect_ratio
  } else {
    BATCH <- merge(BATCH, detect_ratio, by = "species", all = TRUE)
  }
}

row.names(BATCH) <- BATCH$species
BATCH <- BATCH[, which(colnames(BATCH) != "species")]

BATCH_preva2 <- BATCH[grep("detect_ratio_", colnames(BATCH))]
BATCH_preva2[is.na(BATCH_preva2)] <- 0

BATCH_preva2[BATCH_preva2 > 0] <- 1
appearance_counts <- rowSums(BATCH_preva2)
BATCH_preva3 <- BATCH_preva2[appearance_counts >= 2, ] # 453 species occur in multiple batches

####------------------------------------------------------------ (5).batch filter------------------------------------------------------------
## 1.Based on the batch filter species results, organize the original read file
file3 <- list.files("E:/2025.8-2026.7/iScience_revision/3_prevalence_prepare/")
for (ss in 1:length(file3)) {
  SS <- file3[ss]
  bat <- str_split(SS, "_")[[1]][4]
  bat <- gsub("\\.csv", "", bat)
  data_0 <- read.csv(paste0("E:/2025.8-2026.7/iScience_revision/3_prevalence_prepare/",SS),row.names=1)
  Species_batch <- data_0[which(row.names(data_0) %in% row.names(BATCH_preva3)), ]
  fwrite(Species_batch, paste0("E:/2025.8-2026.7/iScience_revision/5_batch_filter/SSC_specieQC_batch_filter_Phase_", bat, ".csv"), row.names = TRUE)
}

## 2.Summarize the number of species at specie level
# Count the remaining species of the batch filter
file_path <- "E:/2025.8-2026.7/iScience_revision/5_batch_filter/"
file5 <- list.files("E:/2025.8-2026.7/iScience_revision/5_batch_filter/")
col_union <- NULL

# Loop through the column names of each file and obtain the union
for (i in 1:length(file5)) {
  current_file <- fread(paste0(file_path, file5[i]), header = TRUE)
  col_union <- union(col_union, current_file$V1)
}
length(unique(col_union)) # 453

cate=read.csv("E:/2025.8-2026.7/iScience_revision/1_rawdata/SSC_kraken_cate_16344.csv")
read2_cate=as.data.frame(col_union)
colnames(read2_cate)=c("V1")
read2_cate1=merge(read2_cate,cate,by="V1")
length(unique(read2_cate1$Genus)) # 156
length(unique(read2_cate1$Species)) # 453

####------------------------------------------------------------ (6).read count filter ------------------------------------------------------------
rm(list = ls())

## 1.Combine results from multiple batches
file5 <- list.files("E:/2025.8-2026.7/iScience_revision/5_batch_filter/")
specie_all <- as.data.frame(fread(paste0("E:/2025.8-2026.7/iScience_revision/5_batch_filter/", file5[1]), header = TRUE))

for (i in 2:length(file5)) {
  current_file <- as.data.frame(fread(paste0("E:/2025.8-2026.7/iScience_revision/5_batch_filter/", file5[i]), header = TRUE))
  specie_all <- merge(specie_all, current_file, by = "V1", all = TRUE)
}

row.names(specie_all) <- specie_all$V1
specie_all <- specie_all[,which(colnames(specie_all) !="V1")]
setnafill(specie_all, fill = 0)

## 2.Screen out species with the number of reads in at least one sample >= 100
specie_all2 <- specie_all[apply(specie_all, 1, function(x) any(x >= 100)), ] # 412

## 3.Based on the read count filter species results, organize the original read file
file3 <- list.files("E:/2025.8-2026.7/iScience_revision/3_prevalence_prepare/")
for (ss in 1:length(file3)) {
  SS <- file3[ss]
  bat <- str_split(SS, "_")[[1]][4]
  bat <- gsub("\\.csv", "", bat)
  data_0 <- read.csv(paste0("E:/2025.8-2026.7/iScience_revision/3_prevalence_prepare/",SS),row.names=1)
  
  Species_batch <- data_0[which(row.names(data_0) %in% row.names(specie_all2)), ]
  fwrite(Species_batch, paste0("E:/2025.8-2026.7/iScience_revision/6_read_count_filter/SSC_specieQC_read_count_filter_Phase_", bat, ".csv"), row.names = TRUE)
}

## 4.Summarize the number of species at specie level
# Count the remaining species of the read count filter
cate=read.csv("E:/2025.8-2026.7/iScience_revision/1_rawdata/SSC_kraken_cate_16344.csv")
read2_cate=as.data.frame(rownames(specie_all2))
colnames(read2_cate)=c("V1")
read2_cate1=merge(read2_cate,cate,by="V1")
length(unique(read2_cate1$Genus)) # 107
length(unique(read2_cate1$Species)) # 336

####------------------------------------------------- (7).contaminate list filter -------------------------------------------------
rm(list = ls())

## 1.Screen out species with the number of reads in at least one sample >= 100
file6 <- list.files("E:/2025.8-2026.7/iScience_revision/6_read_count_filter/")
specie_all <- as.data.frame(fread(paste0("E:/2025.8-2026.7/iScience_revision/6_read_count_filter/", file6[1]), header = TRUE))

for (i in 2:length(file6)) {
  current_file <- as.data.frame(fread(paste0("E:/2025.8-2026.7/iScience_revision/6_read_count_filter/", file6[i]), header = TRUE))
  specie_all <- merge(specie_all, current_file, by = "V1", all = TRUE)
}

row.names(specie_all) <- specie_all$V1
specie_all <- specie_all[,which(colnames(specie_all) !="V1")]
setnafill(specie_all, fill = 0)
specie_all2 <- specie_all
specie_all2$V1 <- rownames(specie_all2)

## 2.Remove contaminant
cate=read.csv("E:/2025.8-2026.7/iScience_revision/1_rawdata/SSC_kraken_cate_16344.csv",row.names = 1)
specie_all4=merge(specie_all2,cate,by="V1")
contaminant_genus=c( "Mycoplasma", "Burkholderia","Bradyrhizobium","Mezorhizobium","Variovorax",
                     "Mycoplasmsa","Mycobacterium", "Staphylococcus","Streptomyces","Streptococcus",
                     "Pseudomonas","Achromobacter","Psuedomonas","Xanthomonas","Acidovorax","Mesorhizoium",
                     "Acinetobacter","Achromobacter","Pseudomonas","Stenotrophomonas","Ralstonia",
                     "Methylobacterium","Corynebacterium","Sphingomonas")
specie_all5 <- specie_all4[!specie_all4$Genus %in% contaminant_genus, ]
raw_species=read.csv("E:/2025.8-2026.7/iScience_revision/2_abundance_filter/SSC_specie_595_count.csv",header = T,row.names = 1);dim(raw_species) # 595
DEcontaminant <- raw_species[specie_all5$V1,];dim(DEcontaminant) # 197
Contaminant <- raw_species[!row.names(raw_species) %in% specie_all5$V1, ];dim(Contaminant) # 398

fwrite(Contaminant,"E:/2025.8-2026.7/iScience_revision/SSC_prevalence_filter_Contaminant_398_ratio.csv",row.names=T)
fwrite(DEcontaminant,"E:/2025.8-2026.7/iScience_revision/SSC_prevalence_filter_DEcontaminant_197_ratio.csv",row.names=T)

## 3.Based on the contaminate list filter species results, organize the original read file
file3 <- list.files("E:/2025.8-2026.7/iScience_revision/3_prevalence_prepare/")
for (ss in 1:length(file3)) {
  SS <- file3[ss]
  bat <- str_split(SS, "_")[[1]][4]
  bat <- gsub("\\.csv", "", bat)
  data_0 <- read.csv(paste0("E:/2025.8-2026.7/iScience_revision/3_prevalence_prepare/",SS),row.names=1)
  
  Species_batch <- data_0[which(row.names(data_0) %in% specie_all5$V1), ]
  fwrite(Species_batch, paste0("E:/2025.8-2026.7/iScience_revision/7_contam_filter/SSC_specieQC_contam_filter_Phase_", bat, ".csv"), row.names = TRUE)
}

## 4.Summarize the number of species at specie level
# Count the remaining species of the contaminate list filter
file_path <- "E:/2025.8-2026.7/iScience_revision/7_contam_filter/"
file5 <- list.files("E:/2025.8-2026.7/iScience_revision/7_contam_filter/")
col_union <- NULL

# Loop through the column names of each file and obtain the union
for (i in 1:length(file5)) {
  current_file <- fread(paste0(file_path, file5[i]), header = TRUE)
  col_union <- union(col_union, current_file$V1)
}
length(unique(col_union)) # 197

cate=read.csv("E:/2025.8-2026.7/iScience_revision/1_rawdata/SSC_kraken_cate_16344.csv")
read2_cate=as.data.frame(col_union)
colnames(read2_cate)=c("V1")
read2_cate1=merge(read2_cate,cate,by="V1")
length(unique(read2_cate1$Genus)) # 90
length(unique(read2_cate1$Species)) # 197

####------------------------------------------------------ (8).correlation filter prepare------------------------------------------------------
rm(list=ls())
library(compositions)
library(stringr)
library(data.table)

Contaminant=read.csv("E:/2025.8-2026.7/iScience_revision/SSC_prevalence_filter_Contaminant_398_ratio.csv", row.names = 1)
DEcontaminant=read.csv("E:/2025.8-2026.7/iScience_revision/SSC_prevalence_filter_DEcontaminant_197_ratio.csv", row.names = 1)

## 1.Generate correlation and p-value matrix
file3 <- list.files("E:/2025.8-2026.7/iScience_revision/3_prevalence_prepare/", pattern = "SSC_specieQC_.*\\.csv")
for (ss in 1:length(file3)) {
  SS <- file3[ss]
  
  bat <- str_split(SS, "_")[[1]][4]
  bat <- gsub("\\.csv", "", bat)
  
  data_0 <- read.csv(paste0("E:/2025.8-2026.7/iScience_revision/3_prevalence_prepare/", SS), row.names = 1)
  Species_con <- data_0[row.names(data_0) %in% row.names(Contaminant), ]
  Species_decon <- data_0[row.names(data_0) %in% row.names(DEcontaminant), ]
  
  # Calculate relative abundance
  total_microbial_reads <- colSums(data_0)
  data_relative_abundance <- sweep(data_0, 2, total_microbial_reads, FUN = "/")
  
  aim <- clr(data_relative_abundance)
  aim_data <- as.data.frame(aim)
  
  # Transposition
  Species_con <- data.frame(t(aim_data[row.names(Species_con), ]), check.names = FALSE)
  Species_decon <- data.frame(t(aim_data[row.names(Species_decon), ]), check.names = FALSE)
  
  # Initialization
  R <- data.frame(matrix(0, ncol(Species_con), ncol(Species_decon),
                         dimnames = list(colnames(Species_con), colnames(Species_decon))), check.names = FALSE)
  P <- data.frame(matrix(1, ncol(Species_con), ncol(Species_decon),
                         dimnames = list(colnames(Species_con), colnames(Species_decon))), check.names = FALSE)
  
  # Calculate the Spearman correlation and P-value
  for (i in 1:ncol(Species_con)) {
    for (j in 1:ncol(Species_decon)) {
      comp_data <- data.frame(con = as.numeric(Species_con[, i]),
                              decon = as.numeric(Species_decon[, j]))
      corr <- cor.test(comp_data$con, comp_data$decon, method = "spearman")
      R[i, j] <- corr$estimate
      P[i, j] <- corr$p.value
    }
  }
  
  # Output the correlation and p-value matrix
  fwrite(R, paste0("E:/2025.8-2026.7/iScience_revision/8_correlation_prepare/SSC_specieQC_cor_R_Phase_", bat, ".csv"),row.names = TRUE)
  fwrite(P, paste0("E:/2025.8-2026.7/iScience_revision/8_correlation_prepare/SSC_specieQC_cor_P_Phase_", bat, ".csv"),row.names = TRUE)
}

## 2.Verify the number of species
file_path <- "E:/2025.8-2026.7/iScience_revision/8_correlation_prepare/"
file5 <- list.files("E:/2025.8-2026.7/iScience_revision/8_correlation_prepare/")

# Contaminant
col_union <- NULL
for (i in 1:length(file5)) {
  current_file <- fread(paste0(file_path, file5[i]), header = TRUE)
  col_union <- union(col_union, current_file$V1)
}

length(unique(col_union)) # 398
# DEcontaminant
col_union <- NULL
for (i in 1:length(file5)) {
  current_file <- read.csv(paste0(file_path, file5[i]), row.names = 1)
  col_union <- union(col_union, colnames(current_file))
}

length(unique(col_union)) # 197

####------------------------------------------------- (9).correlation filter -------------------------------------------------
## 1.Remove contaminants according to the P and R matrices
rm(list=ls())

# Step1: Determine the total list of contaminants in all batches
removed_species_list <- list()
all_contaminants <- c()
file3 <- list.files("E:/2025.8-2026.7/iScience_revision/3_prevalence_prepare/", pattern = "SSC_specieQC_.*\\.csv")
for (ss in 1:length(file3)) {
  SS <- file3[ss]
  
  bat <- str_split(SS, "_")[[1]][4]
  bat <- gsub("\\.csv", "", bat)
  data_0 <- read.csv(paste0("E:/2025.8-2026.7/iScience_revision/3_prevalence_prepare/", SS), row.names = 1)
  
  R <- read.csv(paste0("E:/2025.8-2026.7/iScience_revision/8_correlation_prepare/SSC_specieQC_cor_R_Phase_", bat, ".csv"), row.names = 1, check.names = FALSE)
  P <- read.csv(paste0("E:/2025.8-2026.7/iScience_revision/8_correlation_prepare/SSC_specieQC_cor_P_Phase_", bat, ".csv"), row.names = 1, check.names = FALSE)
  R[is.na(R)] <- 0
  P[is.na(P)] <- 1
  S <- c()  
  for (i in 1:nrow(R)) {
    for (j in 1:ncol(R)) {
      r <- as.numeric(R[i, j])
      p <- as.numeric(P[i, j])
      if (r > 0.8 & p < 0.05) {  # Meet the significance condition
        S <- c(S, colnames(R)[j])  # Record the species to be removed
      }
    }
  }
  
  removed_species <- unique(S)
  RR <- R[, !(colnames(R) %in% removed_species)]
  Species_decontamination <- data_0[which(row.names(data_0) %in% colnames(RR)), ]
  all_contaminants <- c(all_contaminants, unique(S))
}

all_contaminants <- unique(all_contaminants) # 0.8sig 80

# Step2: Remove species from the total pollutant list in each batch
for (ss in 1:length(file3)) {
  SS <- file3[ss]
  
  bat <- str_split(SS, "_")[[1]][4]
  bat <- gsub("\\.csv", "", bat)
  data_0 <- read.csv(paste0("E:/2025.8-2026.7/iScience_revision/3_prevalence_prepare/", SS), row.names = 1)
  
  R <- read.csv(paste0("E:/2025.8-2026.7/iScience_revision/8_correlation_prepare/SSC_specieQC_cor_R_Phase_", bat, ".csv"), row.names = 1, check.names = FALSE)
  RR <- R[, !(colnames(R) %in% all_contaminants)]
  
  remaining_contaminants <- intersect(colnames(RR), all_contaminants)
  if (length(remaining_contaminants) > 0) {
    cat("Warning: Not all contaminants were removed in batch", bat, "\n")
    cat("Remaining contaminants:", remaining_contaminants, "\n")
  } else {
    cat("All contaminants removed successfully in batch", bat, "\n")
  }
  
  Species_decontamination <- data_0[which(row.names(data_0) %in% colnames(RR)), ]
  fwrite(Species_decontamination, paste0("E:/2025.8-2026.7/iScience_revision/9_correlation_filter/SSC_specieQC_cor_filter_R_Phase_", bat, ".csv"), row.names = TRUE)
}

## 2.Summarize the number of species at specie level
# Count the remaining species of the correlation filter
file_path <- "E:/2025.8-2026.7/iScience_revision/9_correlation_filter/"
file5 <- list.files("E:/2025.8-2026.7/iScience_revision/9_correlation_filter/")
col_union <- NULL

# Loop through the column names of each file and obtain the union
for (i in 1:length(file5)) {
  current_file <- fread(paste0(file_path, file5[i]), header = TRUE)
  col_union <- union(col_union, current_file$V1)
}
length(unique(col_union)) # 117

cate=read.csv("E:/2025.8-2026.7/iScience_revision/1_rawdata/SSC_kraken_cate_16344.csv")
read2_cate=as.data.frame(col_union)
colnames(read2_cate)=c("V1")
read2_cate1=merge(read2_cate,cate,by="V1")
length(unique(read2_cate1$Genus)) # 65
length(unique(read2_cate1$Species)) # 117

####------------------------------------------------------- (10).prevalence filter -------------------------------------------------------
rm(list = ls())
file5 <- list.files("E:/2025.8-2026.7/iScience_revision/9_correlation_filter/")
specie_all <- as.data.frame(fread(paste0("E:/2025.8-2026.7/iScience_revision/9_correlation_filter/", file5[1]), header = TRUE))

for (i in 2:length(file5)) {
  current_file <- as.data.frame(fread(paste0("E:/2025.8-2026.7/iScience_revision/9_correlation_filter/", file5[i]), header = TRUE))
  specie_all <- merge(specie_all, current_file, by = "V1", all = TRUE)
}

row.names(specie_all) <- specie_all$V1
specie_all <- specie_all[,which(colnames(specie_all) !="V1")]
setnafill(specie_all, fill = 0)

detect <- specie_all
detect[detect > 0] <- 1
detect["sum"] <- rowSums(detect)
detect["detect_ratio"] <- detect$sum / ncol(detect[,1:3892])
detect$V1=row.names(detect)

cate=read.csv("E:/2025.8-2026.7/iScience_revision/1_rawdata/SSC_kraken_cate_16344.csv",row.names = 1)
tt2=merge(detect,cate,by="V1")
# Screen species with a detection rate ¡Ý 0.2
tt3 <- tt2[tt2$detect_ratio <= 0.2, ]
aim=specie_all[tt3$V1,] # 100
fwrite(aim, "E:/2025.8-2026.7/iScience_revision/10_prevalence_filter/SSC_specieQC_prevalence_filter_0.2_100.csv", row.names = TRUE) 

# Screen the corresponding genus information of each species
da=read.csv("E:/2025.8-2026.7/iScience_revision/10_prevalence_filter/SSC_specieQC_prevalence_filter_0.2_100.csv",row.names = 1)
detect <- da
detect[detect > 0] <- 1
detect["sum"] <- rowSums(detect)
detect["detect_ratio"] <- detect$sum / ncol(detect[,1:3892])
detect$V1=row.names(detect)
tt=detect[,c("V1","detect_ratio")]
cate=read.csv("E:/2025.8-2026.7/iScience_revision/1_rawdata/SSC_kraken_cate_16344.csv",row.names = 1)
tt2=merge(tt,cate,by="V1")
length(unique(tt2$Genus)) # 57
fwrite(tt2, "E:/2025.8-2026.7/iScience_revision/10_prevalence_filter/SSC_specieQC_prevalence_filter_0.2_100_infor.csv", row.names = TRUE) 
