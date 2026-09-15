#!/bin/sh

cd /analysisdata/fantom6/Interactome/single_cell_wallace/ATAC/subset_BAM

for file in *CB.bam; do
    filename=$(basename "$file" _CB.bam)
    
    TOBIAS FootprintScores \
    --signal ./"$filename"/"$filename"_CB_corrected.bw \
    --regions /analysisdata/fantom6/Interactome/single_cell_wallace/ATAC/JC_merged/peaks.merged.all.bed \
    --output ./"$filename"/"$filename"_CB_footprints.bw \
    --cores 8 
done


