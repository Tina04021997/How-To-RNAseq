# How-to-RNAseq
This is the pipeline for RNA-seq DESeq2 analyzing, which uses the data from publication [A circadian clock regulates efflux by the
blood-brain barrier in mice and human cells](https://www.nature.com/articles/s41467-020-20795-9.pdf).

## Pipeline workflow
Which includes:
1. QC test: FastQC + MultiQC
2. Alignment: STAR
3. Generate read summarization: featureCounts
4. Differential Expression (in R Studio): DESeq2

<img src="https://github.com/Tina04021997/How-to-RNAseq/blob/main/flow.jpg" width="35%" height="35%">

## Environments
- STAR v2.7.8a
- Python v3.6.13
- Samtools v1.12
- Snakemake v5.7.0 
- Subread v2.0.1


## Input data
For this exercise, download ctrl ZT2, ZT6 each 3 pairs of fatsq files (total 12 fastq files) from [GEO](https://www.ncbi.nlm.nih.gov/sra?term=SRX6720701).

Now we'll have two conditions, three replicates:
- CTRL ZT2 replicate1 (SRR9973379) --> ZT02(R1)
- CTRL ZT2 replicate2 (SRR9973380) --> ZT02(R2)
- CTRL ZT2 replicate3 (SRR9973381) --> ZT02(R3)
- CTRL ZT6 replicate1 (SRR9973385) --> ZT06(R1)
- CTRL ZT6 replicate2 (SRR9973386) --> ZT06(R2)
- CTRL ZT6 replicate3 (SRR9973387) --> ZT06(R3)

Use **sra-tools** to download these files with **download.sh** script.

Download sra-tools by:
```
$ wget https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/2.11.0/sratoolkit.2.11.0-ubuntu64.tar.gz
$ tar -zxvf sratoolkit.2.11.0-ubuntu64.tar.gz
```

## Reference data
Ensemble v102 GRCm38.
```
$ wget -c http://ftp.ensembl.org/pub/release-102/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna.primary_assembly.fa.gz
$ wget -c http://ftp.ensembl.org/pub/release-102/gtf/mus_musculus/Mus_musculus.GRCm38.102.gtf.gz
$ gunzip Mus_musculus.GRCm38.dna.primary_assembly.fa.gz
$ gunzip Mus_musculus.GRCm38.102.gtf.gz
```

## Run Snakefile
Run **Snakefile** by ```snakemake -p -j 20```

## Download counts.txt file from Linux to desktop
Navigate to desktop's local terminal and enter the following scp line at command.
```
$scp LinuxUserName@avisIP:/LVM_data/tina/RNAseq/counts/counts.txt ~/Desktop/
```

## References
SnakeMake STAR script
- https://evodify.com/rna-seq-star-snakemake/ 

FeatureCounts Documentation
- http://bioinf.wehi.edu.au/featureCounts/  

DESeq2 Documentation
- http://www.bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html
