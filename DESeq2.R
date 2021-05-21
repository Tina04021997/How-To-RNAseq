# Author: Tina Yang
# Date: May 21,2021
# RNA-seq data analization with DESeq2
# setwd("~/Desktop/Joey Lab/How-To-RNAseq")

##### Data preparation #####
# Load the counts txt file
GE <- read.table("RNAseq_GE.txt", header=T)

# Load packages
library("biomaRt")
library("dplyr")
library("tidyverse")

# Convert ENSMUSG number to Gene name
ENSMUSG <- GE$Geneid
require(biomaRt)

ensembl <- useMart('ensembl', dataset = 'mmusculus_gene_ensembl')

annot <- getBM(
  attributes = c(
    'mgi_symbol',
    'external_gene_name',
    'ensembl_gene_id',
    'entrezgene_id',
    'gene_biotype'),
  filters = 'ensembl_gene_id',
  values = ENSMUSG,
  mart = ensembl)

ge <- merge(
  x = as.data.frame(GE),
  y =  annot,
  by.y = 'ensembl_gene_id',
  all = F,
  by.x = 'Geneid',
  no.dups = T)

# Extract the wanted columns
GE_counts <- ge[,c(7:13)]
GE_counts <- GE_counts[,c(7,1:6)]

# Rename the columns
colnames(GE_counts) <- c("Gene","ZT02(R1)","ZT02(R2)","ZT02(R3)","ZT06(R1)","ZT06(R2)","ZT06(R3)")

# Remove rows with empty or NA values in Gene
GE_counts <- GE_counts[!(is.na(GE_counts$Gene) | GE_counts$Gene==""), ]

# Remove duplicate
GE_counts <- GE_counts %>% distinct(Gene, .keep_all = TRUE)

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

# Make sure that the columnm names of cts == row names of coldata
all(rownames(coldata) == colnames(cts))


##### Constructing DESeqDataSet #####
# Load packages
library("EnhancedVolcano")
library("BiocStyle")
library("knitr")
library("rmarkdown")
library("DESeq2")
library("magrittr")
library("tidyverse")
library("ggplot2")
library("pheatmap")
library("org.Mm.eg.db")
library("clusterProfiler")
library("enrichplot")
library("factoextra")

# Building the DESeqDataSet
dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = coldata,
                              design = ~ time.length)

# Pre-filtering (remove genes with less than 10 reads)
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

# Specify reference level
dds$time.length <- relevel(dds$time.length, ref = "2")


##### Differential expression analysis #####
# In res:log2 fold change (MLE): time.length 6 vs 2  means the estimates are of the logarithmic fold change log2(ZT06 vs ZT02)
dds <- DESeq(dds)
res <- results(dds)

# Under adjusted p value = 0.1, the number of up-regulated genes is 4914; number of down-regulated genes is 5317
summary(res)

## We can: reorder the results by the smallest p value
resOrdered <- res[order(res$pvalue),]
## We can: adjust the p value; LFC value)
res05 <- results(dds, alpha=0.05, lfcThreshold = 1)
## We can: see how many adjusted p value is < 0.1
sum(res$padj < 0.1, na.rm=TRUE)
## We can: subset DE genes
deg <- as.data.frame(subset(res, padj < 0.05, abs(log2FoldChange)>1))


# Exporting the results to csv files 
write.csv(resOrdered, file = "res_ZT02_ZT06.csv")
## We can: export results that pass a certain criterial 
write.csv(as.data.frame(subset(resOrdered, padj < 0.05)), 
          file="ZT02_ZT06_results.csv")


##### Volcano plot #####
keyvals <- ifelse(
  res$log2FoldChange < -2.5, 'royalblue',
  ifelse(res$log2FoldChange > 2.5, 'red',
         'grey'))
keyvals[is.na(keyvals)] <- 'black'
names(keyvals)[keyvals == 'red'] <- 'UP'
names(keyvals)[keyvals == 'grey'] <- 'INSIGNIFICANT'
names(keyvals)[keyvals == 'royalblue'] <- 'DOWN'

EnhancedVolcano(res,
                lab = rownames(res),
                x = 'log2FoldChange',
                y = 'pvalue',
                selectLab = rownames(res)[which(names(keyvals) %in% c('high', 'low'))],
                xlab = bquote(~Log[2]~ 'fold change'),
                ylab = bquote(~Log[10]~ adjusted~italic(P)),
                title = 'ZT02 vs ZT06',
                subtitle = 'CTRL',
                pCutoff = 10e-14,
                FCcutoff = 1.0,
                pointSize = 3.0,
                colCustom = keyvals,
                labSize = 5.0,
                labCol = 'black',
                labFace = 'bold',
                boxedLabels = TRUE,
                colAlpha = 0.5,
                legendPosition = 'right',
                legendLabSize = 11,
                legendIconSize = 4.0,
                drawConnectors = TRUE,
                widthConnectors = 0.5,
                gridlines.major = TRUE,
                gridlines.minor = FALSE,
                xlim = c(-10,10),
                border = 'full',
                colConnectors = 'grey')


##### PCA plots #####
de <- vsd[row.names(vsd) %in% row.names(deg),]
result <- de
pca_result <- prcomp(t(result))
pca_result_perc <- round(100*pca_result$sdev^2/sum(pca_result$sdev^2),1)
df_pca_result <- data.frame(PC1 = pca_result$x[,1], PC2 = pca_result$x[,2], sample = colnames(result), time.length = coldata$time.length)

ggplot(df_pca_result, aes(PC1, PC2, color = time.length)) +
  geom_point(size = 4) +
  ggtitle('DEGs') +
  labs(x = paste0('PC1: ',pca_result_perc[1],'% variance'), y = paste0('PC2: ',pca_result_perc[2],'% variance')) +
  theme(plot.title = element_text(hjust = 0.5, size = 18), axis.title = element_text(size = 14), axis.text = element_text(size = 12))

# Variance plot
# Take a look at the number of PCs 
pca_result[["x"]]
PCs <- factor(1:6, levels = 1:6)
df <- data.frame(PCs, pca_result_perc)
df <- df[1:6,]
ggplot(df, aes(PCs, pca_result_perc, group = 1)) +
  geom_bar(stat = 'identity', fill = 'royalblue') +
  ggtitle('Variances') +
  xlab('PCs') +
  ylab('Percentage of variances') +
  theme(plot.title = element_text(hjust = 0.5, size = 22), axis.title = element_text(size = 14), axis.text = element_text(size = 12)) +
  geom_text(aes(label = pca_result_perc), vjust = -0.3, size = 3.5) +
  geom_point() +
  geom_line()



##### Functional analysis (GO analysis)#####
UP <- deg %>% subset(log2FoldChange > 1) %>% tibble:: rownames_to_column('gene')
DOWN <- deg %>% subset(log2FoldChange < -1) %>% tibble:: rownames_to_column('gene')

# Check the keytypes in Mouse database
organism <- org.Mm.eg.db
keytypes(organism)

# Create enrichGO object
eGO_UP <- enrichGO(gene = UP$gene,
                   OrgDb = organism,
                   keyType = "SYMBOL",
                   ont = "ALL",
                   pAdjustMethod = "BH",
                   pvalueCutoff = 0.01,
                   qvalueCutoff = 0.05)

eGO_DOWN <- enrichGO(gene = DOWN$gene,
                     OrgDb = organism,
                     keyType = "SYMBOL",
                     ont = "ALL",
                     pAdjustMethod = "BH",
                     pvalueCutoff = 0.01,
                     qvalueCutoff = 0.05)

# Create dot plot
dotplot(eGO_UP, showCategory = 20) + ggtitle('GO_UPregulated')
dotplot(eGO_DOWN, showCategory = 20) + ggtitle('GO_DOWNregulated')

