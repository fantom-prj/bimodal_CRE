#!/bin/sh

cd /analysisdata/fantom6/Interactome/single_cell_wallace/Figure/Data_and_code/ChIP_seq/final_data/TEAD4
idr --samples TEAD4-rep1_peaks.narrowPeak TEAD4-rep2_peaks.narrowPeak \
    --input-file-type narrowPeak \
    --rank p.value \
    --output-file idr/TEAD4/TEAD4_2_output.txt \
    --plot

cd /analysisdata/fantom6/Interactome/single_cell_wallace/Figure/Data_and_code/ChIP_seq/final_data/ONECUT2
idr --samples ONECUT-rep1_peaks.narrowPeak ONECUT-rep2_peaks.narrowPeak \
    --input-file-type narrowPeak \
    --rank p.value \
    --output-file idr/ONECUT_2_output.txt \
    --plot
    