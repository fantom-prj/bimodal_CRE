# Super enhancer identification from iPSC, NSC and Neuron using H3K27ac CUT&Tag data 

* SE.R -> main code describing all the steps
* *.json -> files contain the parameter and input files for ENCODE ATAC-seq pipeline 
* Results are presented in Fig.2
* Input data is not included here, please refer to https://fantom.gsc.riken.jp/6/suppl/Yip_et_al_2026_scTF
* For fastq and bam files, please refer to DDBJ DRA019568 (DRR614873-DRR614905)

# Related Methods
* [Super enhancer identification](#SE)


# <a name="SE"></a>Super enhancer identification
Super enhancers were identified independently from the 3 samples (iPSC, NSC, and neurons) as described previously. Briefly, histone mark H3K27ac was profiled from iPSC, NSC and Neuron by CUT&Tag. FASTQ files were subjected to the ENCODE ATAC-seq pipeline (v1.7.0) to identify H3K27ac peaks, which were used as the reference loci for counting the read number from the BAM files using ROSE (v1.3.1). H3K27ac peaks within a 10 kb distance were stitched together, and the total H3K27ac signal (read counts) was calculated for each stitched region, subtracting the background signal from the IgG control. The previously defined enhancer-like CREs were then cross-referenced with the super enhancer regions.
```
[1] ROSE_main.py -g hg38 -i [significant H3K27Ac peaks] -r [H3K27Ac bam file] -c [IgG bam file] -o [output folder] -s 10000 -t 2500
```
