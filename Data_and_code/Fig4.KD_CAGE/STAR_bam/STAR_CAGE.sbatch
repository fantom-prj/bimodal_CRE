#!/bin/bash
#SBATCH --job-name=STAR_CAGE
#SBATCH --cpus-per-task=12
#SBATCH --mem=48G
#SBATCH --time=12:00:00
#SBATCH --output=logs/STAR_%A_%a.out
#SBATCH --error=logs/STAR_%A_%a.err

set -euo pipefail

# Configuration
STAR=STAR
gdir=/analysisdata/fantom6/Interactome/resources/STAR_index_hg38_gencode39
sample_file=samples.tsv

thread_num=${SLURM_CPUS_PER_TASK}

mkdir -p map

# Read sample information
line=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$sample_file")
IFS=$'\t' read -r name R1 R2 <<< "$line"

echo "Processing sample: $name"

# Check input files
[[ -s "$R1" ]] || { echo "Missing R1: $R1"; exit 1; }
[[ -s "$R2" ]] || { echo "Missing R2: $R2"; exit 1; }

# Handle compressed FASTQ files
read_command=()

if [[ "$R1" == *.gz && "$R2" == *.gz ]]; then
    read_command=(--readFilesCommand zcat)
elif [[ "$R1" == *.gz || "$R2" == *.gz ]]; then
    echo "Error: R1 and R2 have different compression formats"
    exit 1
fi

# STAR alignment with coordinate sorting
"$STAR" \
    --twopassMode None \
    --runThreadN "$thread_num" \
    --outFilterMultimapNmax 20 \
    --outSAMtype BAM SortedByCoordinate \
    --outSAMunmapped None \
    --outReadsUnmapped Fastx \
    --genomeDir "$gdir" \
    --readFilesIn "$R1" "$R2" \
    "${read_command[@]}" \
    --outFileNamePrefix "map/${name}." \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --limitBAMsortRAM 30000000000

# Index sorted BAM
samtools index \
    -@ "$thread_num" \
    "map/${name}.Aligned.sortedByCoord.out.bam"

echo "Completed: $name"

