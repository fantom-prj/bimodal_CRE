#!/bin/sh

cd /analysisdata/fantom6/Interactome/single_cell_wallace/ATAC/subset_BAM

for file in *CB.bam; do
    # Extracting the filename without the extension
    filename=$(basename "$file" _CB.bam)
    
    # Create a directory with a specific name
    mkdir "$filename"

    TOBIAS ATACorrect \
    --bam "$file" \
    --genome /analysisdata/fantom6/Interactome/resources/minimap2_mapping/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta \
    --peaks /analysisdata/fantom6/Interactome/single_cell_wallace/ATAC/JC_merged/peaks.merged.all.bed \
    --blacklist ~/TOBIAS_snakemake/Blacklist/lists/hg38-blacklist.v2.bed \
    --outdir /analysisdata/fantom6/Interactome/single_cell_wallace/ATAC/subset_BAM/"$filename" \
    --cores 8
done

