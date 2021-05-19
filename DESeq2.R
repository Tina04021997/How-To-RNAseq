# Author: Tina Yang
# Date: May 19,2021
# RNA-seq data analization with DESeq2
# setwd("~/Desktop/Joey Lab/How-To-RNAseq")

##### Data preparation #####
# Load the counts txt file
GE <- read.table("RNAseq_GE.txt", header=T)
# Extract the wanted columns
GE_counts <- GE[,c(1,7:12)]
# Rename the columns
colnames(GE_counts) <- c("Gene","ZT02(R1)","ZT02(R2)","ZT02(R3)","ZT06(R1)","ZT06(R2)","ZT06(R3)")
# Save as csv file
write.csv(GE_counts, file = "GE_counts.csv", row.names = F)
