# Differential accessibility based on bulk ATAC after siRNA KD of TEAD2/4 and ONECUT2 

* KD_ATAC.R -> main code describing all the steps
* mapping.sh -> code for BWA mapping fastq to bam file hg38
* featureCount.sh -> code for feature count to quantify each sample according to CRE defined from scATAC-seq
* *.pdf -> QC on the effect of GC content in the differential analyses
* Input data is not included here, please refer to https://fantom.gsc.riken.jp/6/suppl/Yip_et_al_2026_scTF
* For fastq and bam files, please refer to DDBJ (DRA026820 (DRR959950-DRR959961))

# Related Methods
* [Knockdown of transcription factors by siRNA](#KD_exp)
* [Real-time quantitative RT-PCR](#QPCR)
* [Bulk ATAC-seq for knockdown samples](#ATAC_exp)
* [Differentially accessible CREs identification from knockdown samples](#DE_ATAC_KD)


# <a name="KD_exp"></a>Knockdown of transcription factors by siRNA
Bulk ATAC-seq was performed after siRNA knockdown of specific TF. The cells were detached as described and one million cells were subjected to nuclei isolation according to 10x demonstrated protocol (Document CG000169 Rev D), with the concentration of BSA in PBS changed to 0.1%. The isolated nuclei (40,000~60,000 for each replicate) were subjected to library generation using ATAC-Seq Kit (Active motif #53150) according to the manufacturer’s instructions. The libraries were further size-selected using SPRIselect (Beckman Coulter, Cat.No. B23317) to remove large fragments as per the manufacturer’s instructions. The size distribution and quantification of the libraries were assessed by Agilent 2100 Bioanalyzer using High Sensitivity DNA kit (Agilent, Cat.No. 5067-4626) and Qubit 3.0 Fluorometer using Qubit dsDNA HS Assay Kit (Invitrogen, Cat.No. Q32851), respectively. Sequencing was done with NextSeq2000 in the following conditions: R1: 50 cycles; R2: 50 cycles; Index1: 8 cycles; Index2: 8 cycles.

# <a name="QPCR"></a>Real-time quantitative RT-PCR
The KD of TEAD2/4 in iPS cells were performed by transfecting siRNAs (Thermo Fisher Scientific, listed in Table S11) at a final concentration of 20 nM each with 6 l lipofectamine RNAiMAX (Thermo Fisher Scientific) in a 6-well where the cells were harvested 48 hr post-transfection for real-time quantitative RT-PCR and ATAC-seq. The KD of ONECUT2 in cortical neurons was performed on the 6th day after starting the differentiation from NSC. The medium was refreshed before transfection of siRNAs (20 nM with 6 ul lipofectamine RNAiMAX). The cells were harvested 48 hours post-transfection for real-time quantitative RT-PCR and ATAC-seq. For TEAD2/4, two siRNA combinations, including s16075(TEAD2) and s13965 (TEAD4), and s16076 (TEAD2) and s13966 (TEAD4) used with two technical replicates. For ONECUT2 KD, s18177 and s18178 were used with two technical replicates.

# <a name="ATAC_exp"></a>Bulk ATAC-seq for knockdown samples
RNA purification was performed with RNeasy Kit (Qiagen). The harvested cells were pelleted and lysed while the RNA was purified with RNeasy Kit (Qiagen) according to the manufacturer’s instructions. After quantifying the RNA by nanodrop (Thermo Fisher Scientific), real-time quantitative RT-PCR was performed by One Step PrimeScript RT-PCR Kit (Takara Bio Inc.). For each TF, 2 pairs of primer against the specific mRNA were used. The expression level was normalized by GAPDH while the knockdown efficiency was calculated from the fold-change between specific KD and KD by siRNA NC1 (Thermo Fisher Scientific). The primers’ sequences are listed in Table S11.

# <a name="DE_ATAC_KD"></a>Differentially accessible CREs identification from knockdown samples
The fastq files from bulk ATAC-seq after knockdown were mapped to the hg38 genome by BWA. The pair-end bam files were applied to featureCounts of subread to count the number of reads in the annotated CRE regions derived from our scATAC-seq data. The resulting count table was subjected to edgeR for differential accessibility. The knockdown effect of TEAD2/4 in iPSC was revealed by comparing siNC (2 replicates) with TEAD2/4 siRNA combination 1 (2 replicates) and combination 2 (2 replicates) in a design of “-1, 0.5, 0.5” using edgeR. The knockdown of ONECUT2 was revealed by comparing siNC (2 replicates) with ONECUT2 siRNA 1 (2 replicates) and 2 (2 replicates). In each comparison, the CREs were filtered by having at least 100 total reads across the 6 libraries and “filterByExpr()” in edgeR. These resulted in 382,880 and 235,151 evaluable CREs from TEAD2/4 knockdown and ONECUT2 knockdown respectively. The GC content of the CREs was measured and compared with ATAC count and log2 fold change to assess GC content bias. While bias was not observed across samples, data was only normalized by TMM without GC correction. The CREs with |log2FC| > 0.5 and FDR < 0.05 were considered significantly regulated upon KD.








