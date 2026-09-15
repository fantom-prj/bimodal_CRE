#!/bin/bash
#SBATCH --job-name=hic_to_bedpe
#SBATCH --output=HiC_logs/hic_to_bedpe_%A_%a.out
#SBATCH --error=HiC_logs/hic_to_bedpe_%A_%a.err
#SBATCH --array=0-2
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=02:00:00

set -euo pipefail

mkdir -p HiC_logs HiC_bedpe

# Input files (edit these)
ABS_FILES=(
  "/analysisdata/fantom6/Interactome/HiC_KI/wallace/iced/iPSC.5000_abs.bed"
  "/analysisdata/fantom6/Interactome/HiC_KI/wallace/iced/NSC.5000_abs.bed"
  "/analysisdata/fantom6/Interactome/HiC_KI/wallace/iced/NEU.5000_abs.bed"
)

MATRIX_FILES=(
  "/analysisdata/fantom6/Interactome/HiC_KI/wallace/iced/iPSC.5000_iced.matrix"
  "/analysisdata/fantom6/Interactome/HiC_KI/wallace/iced/NSC.5000_iced.matrix"
  "/analysisdata/fantom6/Interactome/HiC_KI/wallace/iced/NEU.5000_iced.matrix"
)

NAMES=(
  "iPSC"
  "NSC"
  "NEU"
)

ABS=${ABS_FILES[$SLURM_ARRAY_TASK_ID]}
MAT=${MATRIX_FILES[$SLURM_ARRAY_TASK_ID]}
BASE=${NAMES[$SLURM_ARRAY_TASK_ID]}

echo "Processing $BASE"
echo "ABS: $ABS"
echo "MATRIX: $MAT"

# Convert HiC-Pro to BEDPE
awk 'BEGIN{OFS="\t"}
FNR==NR {
    chr[$4]=$1
    start[$4]=$2
    end[$4]=$3
    next
}
{
    b1=$1
    b2=$2
    contact=$3

    if (b1 in chr && b2 in chr) {
        print chr[b1], start[b1], end[b1], chr[b2], start[b2], end[b2], ".", contact
    }
}' "$ABS" "$MAT" > "HiC_bedpe/${BASE}.bedpe"

# Sort
sort -k1,1 -k2,2n -k4,4 -k5,5n \
  "HiC_bedpe/${BASE}.bedpe" > "HiC_bedpe/${BASE}.sorted.bedpe"

echo "Done: HiC_bedpe/${BASE}.sorted.bedpe"
