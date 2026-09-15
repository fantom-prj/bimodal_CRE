#!/bin/bash

cd /analysisdata/fantom6/Interactome/single_cell_wallace/ATAC/subset_BAM

for file in iPS*; do
~/subset-bam-1.1.0/subset-bam_linux \
--bam /analysisdata/fantom6/Interactome/scATACcellranger_Kouno/iPS/outs/possorted_bam.bam \
--bam-tag CB \
--cell-barcodes "$file" \
--cores 8 \
--out-bam "${file%.CB.csv}.bam"; done

for file in NSC*; do
~/subset-bam-1.1.0/subset-bam_linux \
--bam /analysisdata/fantom6/Interactome/scATACcellranger_Kouno/NSC/outs/possorted_bam.bam \
--bam-tag CB \
--cell-barcodes "$file" \
--cores 8 \
--out-bam "${file%.CB.csv}.bam"; done

for file in NRN*; do
~/subset-bam-1.1.0/subset-bam_linux \
--bam /analysisdata/fantom6/Interactome/scATACcellranger_Kouno/NRN/outs/possorted_bam.bam \
--bam-tag CB \
--cell-barcodes "$file" \
--cores 8 \
--out-bam "${file%.CB.csv}.bam"; done
