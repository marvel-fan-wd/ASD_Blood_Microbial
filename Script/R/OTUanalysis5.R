rm(list=ls())
library(compositions)
library(stringr)
library(data.table)

Contaminant=read.csv("E:/2025.8-2026.7/iScience_revision/SSC_prevalence_filter_Contaminant_398_ratio.csv", row.names = 1)
DEcontaminant=read.csv("E:/2025.8-2026.7/iScience_revision/SSC_prevalence_filter_DEcontaminant_197_ratio.csv", row.names = 1)
file3 <- list.files("E:/2025.8-2026.7/iScience_revision/prevalence_prepare/", pattern = "SSC_specieQC_.*\\.csv")
for (ss in 1:length(file3)) {
  SS <- file3[ss]
  
  bat <- str_split(SS, "_")[[1]][4]
  bat <- gsub("\\.csv", "", bat)
  
  data_0 <- read.csv(paste0("E:/2025.8-2026.7/iScience_revision/prevalence_prepare/", SS), row.names = 1)
  Species_con <- data_0[row.names(data_0) %in% row.names(Contaminant), ]
  Species_decon <- data_0[row.names(data_0) %in% row.names(DEcontaminant), ]
  
  total_microbial_reads <- colSums(data_0)
  data_relative_abundance <- sweep(data_0, 2, total_microbial_reads, FUN = "/")
  
  aim <- clr(data_relative_abundance)
  aim_data <- as.data.frame(aim)
  
  Species_con <- data.frame(t(aim_data[row.names(Species_con), ]), check.names = FALSE)
  Species_decon <- data.frame(t(aim_data[row.names(Species_decon), ]), check.names = FALSE)
  
  R <- data.frame(matrix(0, ncol(Species_con), ncol(Species_decon),
                         dimnames = list(colnames(Species_con), colnames(Species_decon))), check.names = FALSE)
  P <- data.frame(matrix(1, ncol(Species_con), ncol(Species_decon),
                         dimnames = list(colnames(Species_con), colnames(Species_decon))), check.names = FALSE)
  
  for (i in 1:ncol(Species_con)) {
    for (j in 1:ncol(Species_decon)) {
      comp_data <- data.frame(con = as.numeric(Species_con[, i]),
                              decon = as.numeric(Species_decon[, j]))
      corr <- cor.test(comp_data$con, comp_data$decon, method = "spearman")
      R[i, j] <- corr$estimate
      P[i, j] <- corr$p.value
    }
  }
  
  fwrite(R, paste0("E:/2025.8-2026.7/iScience_revision/correlation_prepare/SSC_specieQC_cor_R_Phase_", bat, ".csv"),row.names = TRUE)
  fwrite(P, paste0("E:/2025.8-2026.7/iScience_revision/correlation_prepare/SSC_specieQC_cor_P_Phase_", bat, ".csv"),row.names = TRUE)
}

removed_species_list <- list();all_contaminants <- c()
file3 <- list.files("E:/2025.8-2026.7/iScience_revision/prevalence_prepare/", pattern = "SSC_specieQC_.*\\.csv")
for (ss in 1:length(file3)) {
  SS <- file3[ss]
  
  bat <- str_split(SS, "_")[[1]][4]
  bat <- gsub("\\.csv", "", bat)
  data_0 <- read.csv(paste0("E:/2025.8-2026.7/iScience_revision/prevalence_prepare/", SS), row.names = 1)
  
  R <- read.csv(paste0("E:/2025.8-2026.7/iScience_revision/correlation_prepare/SSC_specieQC_cor_R_Phase_", bat, ".csv"), row.names = 1, check.names = FALSE)
  P <- read.csv(paste0("E:/2025.8-2026.7/iScience_revision/correlation_prepare/SSC_specieQC_cor_P_Phase_", bat, ".csv"), row.names = 1, check.names = FALSE)
  R[is.na(R)] <- 0
  P[is.na(P)] <- 1
  S <- c()  
  for (i in 1:nrow(R)) {
    for (j in 1:ncol(R)) {
      r <- as.numeric(R[i, j])
      p <- as.numeric(P[i, j])
      if (r > 0.8 & p < 0.05) {  
        S <- c(S, colnames(R)[j]) 
      }
    }
  }
  
  removed_species <- unique(S)
  RR <- R[, !(colnames(R) %in% removed_species)]
  Species_decontamination <- data_0[which(row.names(data_0) %in% colnames(RR)), ]
  all_contaminants <- c(all_contaminants, unique(S))
}

all_contaminants <- unique(all_contaminants)

for (ss in 1:length(file3)) {
  SS <- file3[ss]
  
  bat <- str_split(SS, "_")[[1]][4]
  bat <- gsub("\\.csv", "", bat)
  data_0 <- read.csv(paste0("E:/2025.8-2026.7/iScience_revision/prevalence_prepare/", SS), row.names = 1)
  
  R <- read.csv(paste0("E:/2025.8-2026.7/iScience_revision/correlation_prepare/SSC_specieQC_cor_R_Phase_", bat, ".csv"), row.names = 1, check.names = FALSE)
  RR <- R[, !(colnames(R) %in% all_contaminants)]
  
  remaining_contaminants <- intersect(colnames(RR), all_contaminants)
  if (length(remaining_contaminants) > 0) {
    cat("Warning: Not all contaminants were removed in batch", bat, "\n")
    cat("Remaining contaminants:", remaining_contaminants, "\n")
  } else {
    cat("All contaminants removed successfully in batch", bat, "\n")
  }
  
  Species_decontamination <- data_0[which(row.names(data_0) %in% colnames(RR)), ]
  fwrite(Species_decontamination, paste0("E:/2025.8-2026.7/iScience_revision/correlation_filter/SSC_specieQC_cor_filter_R_Phase_", bat, ".csv"), row.names = TRUE)
}

file_path <- "E:/2025.8-2026.7/iScience_revision/correlation_filter/"
file5 <- list.files("E:/2025.8-2026.7/iScience_revision/correlation_filter/")
specie_all <- as.data.frame(fread(paste0("E:/2025.8-2026.7/iScience_revision/correlation_filter/", file5[1]), header = TRUE))
for (i in 2:length(file5)) {
  current_file <- as.data.frame(fread(paste0("E:/2025.8-2026.7/iScience_revision/correlation_filter/", file5[i]), header = TRUE))
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

cate=read.csv("E:/2025.8-2026.7/iScience_revision/1_rawdata/SSC_cate.csv",row.names = 1)
tt2=merge(detect,cate,by="V1");tt3 <- tt2[tt2$detect_ratio <= 0.2, ]
aim=specie_all[tt3$V1,]
