#!/bin/sh

cd /analysisdata/fantom6/Interactome/ChIP_seq/rawdata
#bwa index /analysisdata/fantom6/Interactome/resources/hg38.fa
bwa mem -t 8 /analysisdata/fantom6/Interactome/resources/hg38.fa TEAD4-rep1_1.fastq.gz TEAD4-rep1_2.fastq.gz | samtools view -h -b -o /analysisdata/fantom6/Interactome/ChIP_seq/BAM/TEAD4-rep1.bam
bwa mem -t 8 /analysisdata/fantom6/Interactome/resources/hg38.fa TEAD4-rep2_1.fastq.gz TEAD4-rep2_2.fastq.gz | samtools view -h -b -o /analysisdata/fantom6/Interactome/ChIP_seq/BAM/TEAD4-rep2.bam
bwa mem -t 8 /analysisdata/fantom6/Interactome/resources/hg38.fa IgG2a-rep1_1.fastq.gz IgG2a-rep1_2.fastq.gz | samtools view -h -b -o /analysisdata/fantom6/Interactome/ChIP_seq/BAM/IgG2a-rep1.bam
bwa mem -t 8 /analysisdata/fantom6/Interactome/resources/hg38.fa IgG2a-rep2_1.fastq.gz IgG2a-rep2_2.fastq.gz | samtools view -h -b -o /analysisdata/fantom6/Interactome/ChIP_seq/BAM/IgG2a-rep2.bam

cd /analysisdata/fantom6/Interactome/ChIP_seq/2nd_seq
#bwa index /analysisdata/fantom6/Interactome/resources/hg38.fa
bwa mem -t 8 /analysisdata/fantom6/Interactome/resources/hg38.fa ONECUT-rep1_1.fastq.gz ONECUT-rep1_2.fastq.gz | samtools view -h -b -o /analysisdata/fantom6/Interactome/ChIP_seq/BAM/ONECUT-rep1.bam
bwa mem -t 8 /analysisdata/fantom6/Interactome/resources/hg38.fa ONECUT-rep2_1.fastq.gz ONECUT-rep2_2.fastq.gz | samtools view -h -b -o /analysisdata/fantom6/Interactome/ChIP_seq/BAM/ONECUT-rep2.bam
bwa mem -t 8 /analysisdata/fantom6/Interactome/resources/hg38.fa IgG-rep1_1.fastq.gz IgG-rep1_2.fastq.gz | samtools view -h -b -o /analysisdata/fantom6/Interactome/ChIP_seq/BAM/IgG-rep1.bam
bwa mem -t 8 /analysisdata/fantom6/Interactome/resources/hg38.fa IgG-rep2_1.fastq.gz IgG-rep2_2.fastq.gz | samtools view -h -b -o /analysisdata/fantom6/Interactome/ChIP_seq/BAM/IgG-rep2.bam

cd /analysisdata/fantom6/Interactome/ChIP_seq/BAM
for file in *.bam; do samtools fixmate -m "$file" "${file%.bam}.fixmate.bam"; done
for file in *.fixmate.bam; do samtools sort "$file" > "${file%.bam}.sort.bam"; done
for file in *.fixmate.sort.bam; do samtools markdup -S "$file" "${file%.fixmate.sort.bam}.markdup.bam"; done
for file in *.markdup.bam; do samtools sort "$file" > "${file%.bam}.sort.bam"; done
for file in *.markdup.sort.bam; do samtools index "$file" ; done
for file in *.markdup.sort.bam; do echo "$file" ; samtools flagstat "$file"; done
for file in *.markdup.sort.bam; do samtools rmdup "$file" "${file%.bam}.nodup.bam"; done
