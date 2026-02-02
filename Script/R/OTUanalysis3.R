rm(list=ls())
library(data.table)

read <- read.csv("E:/2025.8-2026.7/iScience_revision/SSC_specie2.csv", row.names = 1)
batch=read.csv("E:/2025.8-2026.7/iScience_revision/SSC_sample_meta.csv",header = T,row.names = 1)
batch=batch[,c("SPID","Phase")];table(batch$Phase)
sample_ids <- batch$SPID;batch_info <- batch$Phase
temp_list <- list()
for (batch_id in unique(batch_info)) {
  batch_samples <- sample_ids[batch_info == batch_id]
  batch_matrix <- read[, colnames(read) %in% batch_samples, drop = FALSE]
  batch_matrix <- batch_matrix[rowSums(batch_matrix) > 0, ]
  temp_list[[as.character(batch_id)]] <- batch_matrix
}
list2env(temp_list, envir = .GlobalEnv)

output_dir2 <- "E:/2025.8-2026.7/iScience_revision/prevalence_prepare/"
output_dir3 <- "E:/2025.8-2026.7/iScience_revision/batch_filter_prepare/"

batch_matrices <- list(Phase1, Phase2, `Phase3-1`, `Phase3-2`, Phase4, Pilot)
batch_names <- c("Phase_1", "Phase_2", "Phase_31", "Phase_32", "Phase_4", "Phase_p")

for (i in seq_along(batch_matrices)) {
  ff2_spp <- batch_matrices[[i]]
  batch_name <- batch_names[i]
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

file4 <- list.files("E:/2025.8-2026.7/iScience_revision/batch_filter_prepare/", pattern = "SSC_specieQC_.*\\.csv")
BATCH <- NULL
for (ss in 1:length(file4)) {
  SS <- file4[ss]
  bat <- str_split(SS, "_")[[1]][5]
  bat <- gsub("\\.csv", "", bat)
  detect <- read.csv(paste0("E:/2025.8-2026.7/iScience_revision/batch_filter_prepare/SSC_specieQC_prevalence_Phase_", bat, ".csv"), row.names = 1, check.names = FALSE)
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
BATCH_preva3 <- BATCH_preva2[appearance_counts >= 2, ] 

file3 <- list.files("E:/2025.8-2026.7/iScience_revision/prevalence_prepare/")
for (ss in 1:length(file3)) {
  SS <- file3[ss]
  bat <- str_split(SS, "_")[[1]][4]
  bat <- gsub("\\.csv", "", bat)
  data_0 <- read.csv(paste0("E:/2025.8-2026.7/iScience_revision/prevalence_prepare/",SS),row.names=1)
  Species_batch <- data_0[which(row.names(data_0) %in% row.names(BATCH_preva3)), ]
  fwrite(Species_batch, paste0("E:/2025.8-2026.7/iScience_revision/batch_filter/SSC_specieQC_batch_filter_Phase_", bat, ".csv"), row.names = TRUE)
}