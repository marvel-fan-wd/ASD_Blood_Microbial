rm(list = ls())
file5 <- list.files("E:/2025.8-2026.7/iScience_revision/batch_filter/")
specie_all <- as.data.frame(fread(paste0("E:/2025.8-2026.7/iScience_revision/batch_filter/", file5[1]), header = TRUE))

for (i in 2:length(file5)) {
  current_file <- as.data.frame(fread(paste0("E:/2025.8-2026.7/iScience_revision/batch_filter/", file5[i]), header = TRUE))
  specie_all <- merge(specie_all, current_file, by = "V1", all = TRUE)
}

row.names(specie_all) <- specie_all$V1
specie_all <- specie_all[,which(colnames(specie_all) !="V1")]
setnafill(specie_all, fill = 0)
specie_all2 <- specie_all[apply(specie_all, 1, function(x) any(x >= 100)), ]

file3 <- list.files("E:/2025.8-2026.7/iScience_revision/prevalence_prepare/")
for (ss in 1:length(file3)) {
  SS <- file3[ss]
  bat <- str_split(SS, "_")[[1]][4]
  bat <- gsub("\\.csv", "", bat)
  data_0 <- read.csv(paste0("E:/2025.8-2026.7/iScience_revision/prevalence_prepare/",SS),row.names=1)
  
  Species_batch <- data_0[which(row.names(data_0) %in% row.names(specie_all2)), ]
  fwrite(Species_batch, paste0("E:/2025.8-2026.7/iScience_revision/read_count_filter/SSC_specieQC_read_count_filter_Phase_", bat, ".csv"), row.names = TRUE)
}
