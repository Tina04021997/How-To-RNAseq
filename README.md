# How-to-RNAseq
- This is the pipeline for RNA-seq DESeq2 analyzing, which mainly follows the workflow written by Michael I. Love et el. 


## Pipline workflow
Which includes:
1. QC test: FastQC + MultiQC
2. Alignment: STAR
3. Generate read summarization: featureCounts
4. Differential Expression (in R Studio): DESeq2

<img src="https://github.com/Tina04021997/How-To-RNAseq-exercise/blob/main/workflow.jpg" width="35%" height="35%">

## Environments
- STAR v2.7.8a
- Python v3.6.13
- Samtools v1.12
- Snakemake v5.7.0 
- Subread v2.0.1


## Input data
Download 8 fatsq PE files from GEO [GSE52778](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE52778).

We will use only "SRR1039508","SRR1039509", "SRR1039512","SRR1039513", "SRR1039516","SRR1039517","SRR1039520","SRR1039521" for comparison bwtween untreated and treated (DEX) groups.

Use **sra-tools** to download these files with **download.sh** script.

Download sra-tools by:
```
$ wget https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/2.11.0/sratoolkit.2.11.0-ubuntu64.tar.gz
$ tar -zxvf sratoolkit.2.11.0-ubuntu64.tar.gz
```

## Reference data
Ensemble v102 GRCm38.

## Run Snakefile
Run **Snakefile** by ```snakemake -p -j 20```

## References
Salmon Documentation
- https://salmon.readthedocs.io/en/latest/salmon.html

QC for snakemake
- https://snakemake-wrappers.readthedocs.io/en/stable/wrappers/multiqc.html 
