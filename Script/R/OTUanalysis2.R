rm(list=ls())
library(data.table)

read <- read.csv("E:/2025.8-2026.7/iScience_revision/SSC_specie1.csv", row.names = 1)
read2 <- read[, colSums(read) > 100];read2[read2 <= 10] <- 0 
read2 <- sweep(read2, 2, colSums(read2), "/")
read2[read2 < 0.005] <- 0
read2 <- read2[rowSums(read2) > 0, ]
cate=read.csv("E:/2025.8-2026.7/iScience_revision/1_rawdata/SSC_cate.csv")
read2_cate=as.data.frame(rownames(read2))
colnames(read2_cate)=c("V1")
read2_cate1=merge(read2_cate,cate,by="V1")
length(unique(read2_cate1$Genus))
length(unique(read2_cate1$Species))
read3 <- read[rownames(read2), , drop = FALSE]
read3[read3 <= 10] <- 0
fwrite(read3, "E:/2025.8-2026.7/iScience_revision/SSC_specie2.csv",row.names=T)