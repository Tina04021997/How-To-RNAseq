# Author: Tina Yang
# Date: May 18, 2021
# RNAseq analysis pipeline for PE reads:
    # Quality control: FastQC + MultiQC 
    # Alignment: STAR
    # Generate read summarization: featureCounts
# Reference: Ensemble v102 GRCm38

Genome_fa = "/LVM_data/tina/RNAseq/reference/Mus_musculus.GRCm38.dna.primary_assembly.fa"
Genome_GTF = "/LVM_data/tina/RNAseq/reference/Mus_musculus.GRCm38.102.gtf"

SAMPLES, = glob_wildcards("/LVM_data/tina/RNAseq/data/{sample}_1.fastq")

for sample in SAMPLES:
        print("Sample " + sample + " will be processed...STAY TUNED")
#SAMPLES = ["SRR9973379", "SRR9973380", "SRR9973381", "SRR9973385", "SRR9973386", "SRR9973387"]
SAMPLE_ZT02 = ["SRR9973379", "SRR9973380", "SRR9973381"]
SAMPLE_ZT06 = ["SRR9973385", "SRR9973386", "SRR9973387"]

rule all:
    input:
        "counts/counts.txt"

#Step1: Perform FastQC
rule fastqc:
     input:
         expand(["data/{sample}_1.fastq","data/{sample}_2.fastq"], sample = SAMPLES)
     output:
         "qc/fastqc/{sample}_1_fastqc.html",
         "qc/fastqc/{sample}_2_fastqc.html",
         "qc/fastqc/{sample}_1_fastqc.zip",
         "qc/fastqc/{sample}_2_fastqc.zip"
     log:
         "logs/fastqc/{sample}.log"
     threads: 1
     shell:
         "fastqc {input} -o qc/fastqc"

#Step2: Perform MultiQC
rule multiqc:
    input:
        expand(["qc/fastqc/{sample}_1_fastqc.html","qc/fastqc/{sample}_2_fastqc.html"], sample = SAMPLES),
        expand(["qc/fastqc/{sample}_1_fastqc.zip","qc/fastqc/{sample}_2_fastqc.zip"], sample = SAMPLES)
    output:
        "qc/multiqc/multiqc_report.html"
    log:
        "logs/multiqc/multiqc.log"
    threads: 1
    shell:
        "multiqc {input} -o qc/multiqc"
 
#Step3: Indexing
rule STAR_index:
    input:
        fa = "reference/Mus_musculus.GRCm38.dna.primary_assembly.fa",
        GTF = "reference/Mus_musculus.GRCm38.102.gtf"
    output:
        directory("index")
    threads: 20
    shell:
        "mkdir {output} "
        "STAR --runMode genomeGenerate "
        "--genomeDir {output} "
        "--genomeFastaFiles {input.fa} "
        "--sjdbGTFfile {input.GTF} "
        "--runThreadN {threads}"

#Step4: Mapping
# Params must be assigned since you want to specify the sample rather than randomly choosing them
# This step will generate STAR mapping files to the new dir "mapping" 
rule STAR_mapping:
    input:
        R1 = "/LVM_data/tina/RNAseq/data/{sample}_1.fastq",
        R2 = "/LVM_data/tina/RNAseq/data/{sample}_2.fastq",
        REF = "/LVM_data/tina/RNAseq/index"
    output:
        "mapped/{sample}_Aligned.out.bam"
    threads: 20
    params:
        prefix = "{sample}_",
        outdir = "mapped",
        ID = "{sample}"
    message:
        "___STAR mapping {params.prefix}__"
    shell:
        "cd {params.outdir} && "
        "STAR --genomeDir {input.REF} "
        "--readFilesIn {input.R1} {input.R2} "
        "--runThreadN {threads} "
        "--outSAMattrRGline ID:{params.ID} "
        "--outFileNamePrefix {params.prefix} "
        "--outSAMtype BAM Unsorted "
        "--outSAMunmapped Within KeepPairs "
        "--outFilterMismatchNmax 33 "
        "--seedSearchStartLmax 33 && cd .."
        
#Step5: Count reads
rule featureCounts:
    input:
       BAM = expand("/LVM_data/tina/RNAseq/mapped/{sample}_Aligned.out.bam", sample = SAMPLES),
       GTF = "reference/Mus_musculus.GRCm38.102.gtf"
    output:
       "counts/counts.txt"
    message:
        "___Creating counts.txt___"
    threads: 20
    shell:
        "featureCounts -p "
        "-t exon "
        "-g gene_id "
        "-a {input.GTF} "
        "-o {output} "
        "-T {threads} "
        "{input.BAM}"
