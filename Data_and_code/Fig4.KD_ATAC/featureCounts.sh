#!/bin/sh

cd /analysisdata/fantom6/Interactome/KD_ATAC/BAM
~/subread-2.0.6-Linux-x86_64/bin/featureCounts -p -F SAF -a /analysisdata/fantom6/Interactome/KD_ATAC/all_ATAC_peak.saf \
--fracOverlap 0.2 -o snATAC_peak_counts.tsv iPSC_scr_rep1.bam iPSC_scr_rep2.bam iPSC_TEAD24_1_rep1.bam iPSC_TEAD24_1_rep2.bam iPSC_TEAD24_2_rep1.bam iPSC_TEAD24_2_rep2.bam Neuron_scr_rep1.bam Neuron_scr_rep2.bam Neuron_ONECUT2_1_rep1.bam Neuron_ONECUT2_1_rep2.bam Neuron_ONECUT2_rep1.bam Neuron_ONECUT2_rep2.bam
