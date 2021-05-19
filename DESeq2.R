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

# Create the coldata.csv (sample information) with excel 
# Import coldata.csv and cts files
coldata <- read.csv("coldata.csv", row.names = 1)
cts <- as.matrix(read.csv("GE_counts.csv", row.names = 1))

# Correct the colnames of cts 
colnames(cts) <- row.names(coldata)

# Convert coldata from characters to factors
coldata$time.length <- factor(coldata$time.length)
coldata$type <- factor(coldata$type)

# Make sure that the columnm names of cts == row names of coldata
all((rownames(coldata)) %in% colnames(cts))
