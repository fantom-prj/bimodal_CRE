#!/bin/bash

export OMP_NUM_THREADS=64

input_dir="/osc-fs_home/manli/others/Wallace/F6_meme/input_fasta/"
#threshold=100000

#for fasta_file in "$input_dir"*.fa; do
#    filename=$(basename "$fasta_file" .fa)
#    output_dir="out/SEA_${filename}"
#
#    sea --p $fasta_file --m $motif_file --o $output_dir --thresh 100000000 --noseqs
#    echo $filename
#    echo $output_dir
#    echo "Processed: $fasta_file"
#done


# Function to process each FASTA file
process_fasta() {
    fasta_file="$1"
    filename=$(basename "$fasta_file" .fa)
    motif_file="/osc-fs_home/manli/others/Wallace/F6_meme/motif_files/JASPAR2024_CORE_vertebrates_non-redundant_pfms_meme.n879.txt"
    output_dir="out/SEA_${filename}"
    sea --p $fasta_file --m $motif_file --o $output_dir --thresh 10000000 --noseqs
    echo "Processed: $fasta_file"
}

export -f process_fasta

# Find all .fa files in the input directory and process them in parallel
find "$input_dir" -name "*.fa" | parallel -j 16 process_fasta
