#!/bin/bash
#SBATCH --job-name=bam_to_tagalign
#SBATCH --output=logs/bam_to_tagalign_%A_%a.out
#SBATCH --error=logs/bam_to_tagalign_%A_%a.err
#SBATCH --array=0-2
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=04:00:00

set -euo pipefail

mkdir -p logs tagalign bed

# List your 3 BAM files here
BAMS=(
  "/analysisdata/fantom6/Interactome/scATACcellranger_Kouno/iPS/outs/possorted_bam.bam"
  "/analysisdata/fantom6/Interactome/scATACcellranger_Kouno/NSC/outs/possorted_bam.bam"
  "/analysisdata/fantom6/Interactome/scATACcellranger_Kouno/NRN/outs/possorted_bam.bam"
)

NAMES=(
  "iPSC"
  "NSC"
  "NRN"
)

BAM=${BAMS[$SLURM_ARRAY_TASK_ID]}
BASE=${NAMES[$SLURM_ARRAY_TASK_ID]}


echo "Processing $BAM"

# 1. BAM -> sorted BED
bedtools bamtobed -i "$BAM" | \
  sort -k1,1 -k2,2n > "bed/${BASE}.sorted.bed"

# 2. bgzip compress
bgzip -c "bed/${BASE}.sorted.bed" > "bed/${BASE}.sorted.bed.gz"

# 3. tabix index (BED preset)
tabix -p bed "bed/${BASE}.sorted.bed.gz"

# 4. BED -> Tn5-shifted TagAlign
awk 'BEGIN{OFS="\t"}
{
  if ($6 == "+") {
    start = $2 + 4
    end = start + 1
    print $1, start, end, "N", 1000, "+"
  } else if ($6 == "-") {
    end = $3 - 5
    start = end - 1
    if (start >= 0) {
      print $1, start, end, "N", 1000, "-"
    }
  }
}' "bed/${BASE}.sorted.bed" > "tagalign/${BASE}.tn5.tagAlign"

sort -k1,1 -k2,2n "tagalign/${BASE}.tn5.tagAlign" | bgzip > "tagalign/${BASE}.tn5.tagAlign.gz"
tabix -p bed "tagalign/${BASE}.tn5.tagAlign.gz"

echo "Done: bed/${BASE}.sorted.bed.gz + index + tagalign"

