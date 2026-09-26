#!/bin/bash
#$ -cwd
#$ -q bigmem.q
#$ -pe def_slot 4

fq="/analysisdata/rawseq/fastq/SHARED/000012/fantom6atac/HN00163389/HN00163389_10X_RawData_Outs/F6-018-SA1-NRN/HHKCCCCX2"
id="NRN"
ref="/home/kouno/10x_ref/refdata-cellranger-arc-GRCh38-2020-A-2.0.0"

echo -e $fq $id $ref
ulimit -n 20000

export  PATH=/work/cellranger-atac-2.0.0/:$PATH
cellranger-atac count \
                 --fastqs=${fq}                             \
                 --id=${id}                                 \
                 --localcores=4                           \
                 --reference=${ref}

