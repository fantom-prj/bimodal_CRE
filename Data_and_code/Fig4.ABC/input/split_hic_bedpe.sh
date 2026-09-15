#!/bin/bash
set -eo pipefail

BASE="/analysisdata/fantom6/Interactome/ABC/input/HiC_bedpe"

for sample in iPSC NSC NEU
do
  IN="${BASE}/${sample}.sorted.bedpe"
  OUT="${BASE}/${sample}"

  mkdir -p "${OUT}"

  for chr in $(seq 1 22) M X Y
  do
    mkdir -p "${OUT}/chr${chr}"

    awk -v c="chr${chr}" 'BEGIN{OFS="\t"} $1==c && $4==c' "${IN}" \
    | gzip -c > "${OUT}/chr${chr}/chr${chr}.bedpe.gz"
  done
done

