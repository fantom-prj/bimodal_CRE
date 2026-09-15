#!/bin/bash
input_file="/osc-fs_home/manli/others/Wallace/F6_meme/input_bed_list.tsv"
output_dir="/osc-fs_home/manli/others/Wallace/F6_meme/input_fasta/"
genome_file="/osc-fs_home/manli/opt/genome/hg38/genome.fa"


mkdir -p "$output_dir"

while IFS= read -r url; do

    filename=$(basename "$url")
    
    output_file="$output_dir/${filename%.bed.gz}.fa"

    bedtools getfasta -fi "$genome_file" -bed "$url" -fo "$output_file"
    echo "Processed: $url"
done < "$input_file"
