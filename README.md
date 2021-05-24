# How-to-RNAseq
- This is a practice package for RNA-seq DESeq2 analyzing, which mainly follows the workflow written by Michael I. Love et el. 
- Such script can be found in this [page](https://github.com/mikelove/rnaseqGene/blob/master/vignettes/rnaseqGene.Rmd), and the related paper can be found [here](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0099625).

## Pipline workflow
Which includes:
1. QC test: FastQC + MultiQC
2. Calculation of transcript abundance: Salmon
3. Differential Expression (in R Studio): DESeq2

<img src="https://github.com/Tina04021997/How-To-RNAseq-exercise/blob/main/workflow.jpg" width="35%" height="35%">

## Environments
- Python v3.6.13
- Samtools v1.12
- Snakemake v5.7.0
- FastQC v0.11.9
- MultiQC v1.10.1
- Salmon v1.4.0

## Input data
For this exercise, download 8 fatsq PE files from GEO [GSE52778](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE52778).

We will use only "SRR1039508","SRR1039509", "SRR1039512","SRR1039513", "SRR1039516","SRR1039517","SRR1039520","SRR1039521" for comparison bwtween untreated and treated (DEX) groups.

Use **sra-tools** to download these files with **download.sh** script.

Download sra-tools by:
```
$ wget https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/2.11.0/sratoolkit.2.11.0-ubuntu64.tar.gz
$ tar -zxvf sratoolkit.2.11.0-ubuntu64.tar.gz
```

## Reference data
Gencode.v37.transcripts.fa.gz
```
$ wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_37/gencode.v37.transcripts.fa.gz
```

## Run Snakefile
Run **Snakefile** by ```snakemake -p -j 20```

## References
Salmon Documentation
- https://salmon.readthedocs.io/en/latest/salmon.html

QC for snakemake
- https://snakemake-wrappers.readthedocs.io/en/stable/wrappers/multiqc.html 
