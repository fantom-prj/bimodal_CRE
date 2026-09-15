#!/bin/sh

cd /analysisdata/fantom6/Interactome/single_cell_wallace/ATAC/subset_BAM

TOBIAS BINDetect \
--motifs /analysisdata/fantom6/Interactome/single_cell_wallace/ATAC/find_summit/JASPAR2024_CORE_vertebrates_non-redundant_pfms_jaspar.txt \
--signals ./out_iPS_0_undifferentiaed/out_iPS_0_undifferentiaed_CB_footprints.bw \
./out_iPS_1_neuroepithelial/out_iPS_1_neuroepithelial_CB_footprints.bw \
./out_NSC_0_iPS_like/out_NSC_0_iPS_like_CB_footprints.bw \
./out_NSC_1_differentiating/out_NSC_1_differentiating_CB_footprints.bw \
./out_NSC_2_astrocyte_like/out_NSC_2_astrocyte_like_CB_footprints.bw \
./out_neuron_progenitor_n_schwann/out_neuron_progenitor_n_schwann_CB_footprints.bw \
./out_neuron_immature/out_neuron_immature_CB_footprints.bw \
./out_neuron_mature/out_neuron_mature_CB_footprints.bw \
./out_neuron_mature_2/out_neuron_mature_2_CB_footprints.bw \
./out_OPC_like/out_OPC_like_CB_footprints.bw \
--genome /analysisdata/fantom6/Interactome/resources/minimap2_mapping/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta \
--peaks /analysisdata/fantom6/Interactome/single_cell_wallace/ATAC/JC_merged/filtered_peaks.merged.all.bed \
--outdir ./BINDetect_compare_output \
--time-series \
--cores 8




