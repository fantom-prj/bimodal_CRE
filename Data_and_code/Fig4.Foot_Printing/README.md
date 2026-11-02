# Footprinting analyses on 879 TF motifs across 10 cell clusters

* footprint_process.R -> main code describing all the steps
* subset.bam.sh -> code for subsetting the iPSC/NSC/Neuron bam files into 10 cell clusters 
* merge.bam.sh -> code for merging the bam from 3 samples
* 20231229tobias.sh -> code for running TOBIAS ATAC correct
* 20231230tobias_step2.sh -> code for running TOBIAS footprint score
* 20240101tobias_step3.sh.sh -> code for running TOBIAS BINDetect
* Input data is not included here, please refer to https://fantom.gsc.riken.jp/6/suppl/Yip_et_al_2026_scTF
* For fastq and bam files, please refer to DDBJ (DRA026820 (DRR959950-DRR959961))

# Related Methods
* [Transcription Factor footprinting](#footprint)


# <a name="footprint"></a>Transcription Factor footprinting
Transcription factor footprinting was performed at the cell cluster level (n = 10) to resolve cell-type specific regulatory interactions. To generate cluster-specific signals, the global scATAC-seq BAM file was demultiplexed based on cell barcodes using subset-bam (v1.1.0). The resulting cluster-level BAM files were processed using the TOBIAS (v0.16.1) framework. Briefly, TN5 insertion biases were corrected using the ATACorrect function, utilizing the GRCh38 reference genome and masking ENCODE blacklist regions to generate bias-corrected accessibility signals. Footprint scores were subsequently calculated within the merged CRE regions using the FootprintScores function. To identify TF binding events, the BINDetect module was employed, integrating bias-corrected signals with the JASPAR 2024 CORE vertebrate non-redundant motif database (n = 879). A time-series analysis across the 10 differentiation clusters was implemented to identify differential TF occupancy and binding dynamics. Footprinting results were further validated by cross-referencing with ChIP-seq and knockdown-ATAC-seq data to ensure high-confidence regulatory associations.


```
[1] ~/subset-bam-1.1.0/subset-bam_linux \
--bam scATACcellranger_output.bam \
--bam-tag CB \
--cell-barcodes cell_barcode_file_for_each_cluster \
--cores 8 \
--out-bam cluster_output.bam

[2] TOBIAS ATACorrect \
--bam "$file" \
--genome GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta \
--peaks peaks.merged.all.bed \
--blacklist ~/TOBIAS_snakemake/Blacklist/lists/hg38-blacklist.v2.bed \
--outdir subset_BAM/"$filename" \
--cores 8

[3] TOBIAS FootprintScores \
--signal ./"$filename"/"$filename"_CB_corrected.bw \
--regions peaks.merged.all.bed \
--output ./"$filename"/"$filename"_CB_footprints.bw \
--cores 8

[4] TOBIAS BINDetect \
--motifs JASPAR2024_CORE_vertebrates_non-redundant_pfms_jaspar.txt \
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
--genome GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta \
--peaks filtered_peaks.merged.all.bed \
--outdir ./BINDetect_compare_output \
--time-series \
--cores 8
```
