####------------------------------------------------------- (11).depth normalization -------------------------------------------------------
rm(list = ls())
library(data.table)
library(dplyr)

da=read.csv("E:/2025.8-2026.7/iScience_revision/10_prevalence_filter/SSC_specieQC_prevalence_filter_0.2_100.csv",row.names = 1)
da$class=row.names(da)

cate=read.csv("E:/2025.8-2026.7/iScience_revision/11_depth_normalize/SSC_specieQC_prevalence_filter_0.2_100_ref.csv", row.names = 1)
cate_ref=cate[,c(1,9)]
colnames(cate_ref)=c("class","ref")

sam=read.csv("E:/2025.8-2026.7/iScience_revision/11_depth_normalize/SSC_cram_mapped_read_samtools.csv")
colnames(sam)
da2=merge(da,cate_ref,by="class")
rownames(da2)=da2$class
da2=da2[,-1]

## 1.Align the order of da2 and sam according to the sample ID
aligned_sam <- sam %>% 
  filter(ID %in% colnames(da2)[1:3892]) %>% # sample filtering
  arrange(match(ID, colnames(da2)[1:3892])) # reordering

# check
if (!all(aligned_sam$ID == colnames(da2)[1:3892])) {
  stop("The sample ID alignment failed. Please check the input data")
}

## 2.Extract the required data
# Pay attention to unit matching
microbial_reads <- da2[, 1:3892] # microbe read count
microbial_genome_size <- da2[, 3893] # microbe genome size
human_reads <- aligned_sam$cram_mapped # human genome read count
human_genome_size <- 3000 # human genome size

## 3.Calculate relative abundance
absolute_abundance <- matrix(0, nrow = nrow(microbial_reads), ncol = ncol(microbial_reads))
for (i in 1:nrow(microbial_reads)) {
  microbial_row <- as.numeric(microbial_reads[i, ])
  current_genome_size <- as.numeric(microbial_genome_size[i])
  
  absolute_abundance[i, ] <- (2 * (microbial_row / current_genome_size)) / 
    (human_reads / human_genome_size) * 1e6
}

rownames(absolute_abundance) <- rownames(microbial_reads) # species
colnames(absolute_abundance) <- colnames(microbial_reads) # sample ID
class(absolute_abundance);absolute_abundance=as.data.frame(absolute_abundance);class(absolute_abundance)

fwrite(absolute_abundance,"E:/2025.8-2026.7/iScience_revision/11_depth_normalize/Species_relative_abundance_100.csv",row.names = T)
