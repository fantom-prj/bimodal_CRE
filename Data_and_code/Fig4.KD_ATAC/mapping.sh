#!/bin/sh

cd /analysisdata/fantom6/Interactome/KD_ATAC/rawdata_2nd
#bwa index /analysisdata/fantom6/Interactome/resources/hg38.fa
bwa mem -t 8 /analysisdata/fantom6/Interactome/resources/hg38.fa scrONECUT2-1_S1_R1_001.fastq.gz scrONECUT2-1_S1_R2_001.fastq.gz | samtools view -h -b -o /analysisdata/fantom6/Interactome/KD_ATAC/BAM/Run2/Neuron_scr_rep1.bam
bwa mem -t 8 /analysisdata/fantom6/Interactome/resources/hg38.fa scrONECUT2-2_S2_R1_001.fastq.gz scrONECUT2-2_S2_R2_001.fastq.gz | samtools view -h -b -o /analysisdata/fantom6/Interactome/KD_ATAC/BAM/Run2/Neuron_scr_rep2.bam
bwa mem -t 8 /analysisdata/fantom6/Interactome/resources/hg38.fa siONECUT2-1-rep1_S3_R1_001.fastq.gz siONECUT2-1-rep1_S3_R2_001.fastq.gz | samtools view -h -b -o /analysisdata/fantom6/Interactome/KD_ATAC/BAM/Run2/Neuron_ONECUT2_1_rep1.bam
bwa mem -t 8 /analysisdata/fantom6/Interactome/resources/hg38.fa siONECUT2-1-rep2_S4_R1_001.fastq.gz siONECUT2-1-rep2_S4_R2_001.fastq.gz | samtools view -h -b -o /analysisdata/fantom6/Interactome/KD_ATAC/BAM/Run2/Neuron_ONECUT2_1_rep2.bam
bwa mem -t 8 /analysisdata/fantom6/Interactome/resources/hg38.fa siONECUT2-2-rep1_S5_R1_001.fastq.gz siONECUT2-2-rep1_S5_R2_001.fastq.gz | samtools view -h -b -o /analysisdata/fantom6/Interactome/KD_ATAC/BAM/Run2/Neuron_ONECUT2_2_rep1.bam
bwa mem -t 8 /analysisdata/fantom6/Interactome/resources/hg38.fa siONECUT2-2-rep2_S6_R1_001.fastq.gz siONECUT2-2-rep2_S6_R2_001.fastq.gz | samtools view -h -b -o /analysisdata/fantom6/Interactome/KD_ATAC/BAM/Run2/Neuron_ONECUT2_2_rep2.bam
bwa mem -t 8 /analysisdata/fantom6/Interactome/resources/hg38.fa scrTEAD24-rep1_S7_R1_001.fastq.gz scrTEAD24-rep1_S7_R2_001.fastq.gz | samtools view -h -b -o /analysisdata/fantom6/Interactome/KD_ATAC/BAM/Run2/iPSC_scr_rep1.bam
bwa mem -t 8 /analysisdata/fantom6/Interactome/resources/hg38.fa scrTEAD24-rep2_S8_R1_001.fastq.gz scrTEAD24-rep2_S8_R2_001.fastq.gz | samtools view -h -b -o /analysisdata/fantom6/Interactome/KD_ATAC/BAM/Run2/iPSC_scr_rep2.bam
bwa mem -t 8 /analysisdata/fantom6/Interactome/resources/hg38.fa siTEAD24-1-rep1_S9_R1_001.fastq.gz siTEAD24-1-rep1_S9_R2_001.fastq.gz | samtools view -h -b -o /analysisdata/fantom6/Interactome/KD_ATAC/BAM/Run2/iPSC_TEAD24_1_rep1.bam
bwa mem -t 8 /analysisdata/fantom6/Interactome/resources/hg38.fa siTEAD24-1-rep2_S10_R1_001.fastq.gz siTEAD24-1-rep2_S10_R2_001.fastq.gz | samtools view -h -b -o /analysisdata/fantom6/Interactome/KD_ATAC/BAM/Run2/iPSC_TEAD24_1_ep2.bam
bwa mem -t 8 /analysisdata/fantom6/Interactome/resources/hg38.fa siTEAD24-2-rep1_S11_R1_001.fastq.gz siTEAD24-2-rep1_S11_R2_001.fastq.gz | samtools view -h -b -o /analysisdata/fantom6/Interactome/KD_ATAC/BAM/Run2/iPSC_TEAD24_2_rep1.bam
bwa mem -t 8 /analysisdata/fantom6/Interactome/resources/hg38.fa siTEAD24-2-rep2_S12_R1_001.fastq.gz siTEAD24-2-rep2_S12_R2_001.fastq.gz | samtools view -h -b -o /analysisdata/fantom6/Interactome/KD_ATAC/BAM/Run2/iPSC_TEAD24_2_rep2.bam

