# Pre-processing and processing of 10x scATAC-seq data generated from three samples (iPSC, NSC and Neuron)

* scATAC_processing.R -> main code describing all the analyses
* cellranger/*.sh -> associated code running cellranger
* Results are presented in Fig.1
* Input data is not included here, please refer to https://fantom.gsc.riken.jp/6/suppl/Yip_et_al_2026_scTF

# Related Methods
* [Chromium for scATAC-seq](#10x_exp)
* [Single-cell ATAC-seq analytical pipeline](#scATAC)

# <a name="10x_exp"></a>Chromium for scATAC-seq
For the scATAC-seq, the 3 cell types were harvested at the same time where one million cells per each cell type were subjected to nuclei isolation according to the protocol developed by 10x Genomics (Document CG000169 Rev D), except that the concentration of BSA in PBS was changed to 0.1%. The isolated nuclei (~7,000 for each cell type) were loaded into the Chromium System in 3 reactions for cell barcode addition within GEMs using the Chromium Next GEM Single Cell ATAC Reagent Kits v1.1 (PN: 1000175) user guide (10x Genomics, CG000209 Rev F). Briefly, the transposition reaction was prepared by mixing the desired number of nuclei with ATAC Buffer (10x Genomics) and ATAC Enzyme (10x Genomics) and incubated for 60 min at 37 °C; nuclei were loaded with the goal of recovering 5,000 nuclei. DNA linear amplification was then performed by incubating the GEMs under the following thermal cycling conditions: 72 °C for 5 min, 98 °C for 30 s and 10 cycles of 98 °C for 10 s, 59 °C for 30 s and 72 °C for 1 min. After PCR indexing (iPS: SI-NA-A3; NSC: SI-NA-B3; Neuron: SI-NA-C3), sequencing libraries were subjected to a cleanup with SPRIselect reagent and verified by Agilent Bioanalyzer with a high-sensitivity chip (Agilent Technologies). The 3 libraries were pooled and subjected to HiseqX_Ten for sequencing with the following conditions: R1: 150 cycles; R2: 150 cycles; Index1: 8 cycles; Index2: 8 cycles.


# <a name="scATAC"></a>Single-cell ATAC-seq analytical pipeline
Raw data from scATAC-seq was processed through the Cell Ranger ATAC pipeline (v2.0.0) using pre-built reference from Cell Ranger refdata-cellranger-arc-GRCh38-2020-A-2.0.0. Signac was used to process the Cell Ranger output. Cells were then filtered for nCount_peaks > 1000 & < 20,000, FRiP > 0.5, is__cell_barcode equal to 1, nucleosome_signal < 4, and TSS.enrichment > 3. After quality control, 5,028 iPSCs, 5,040 NSCs, and 6,111 neurons were retained. These cells were then integrated by the lsi algorithm. The peaks called for each sample were then merged to an unified peak set (n = 463,966) which was then used to quantify reads for each cell in each sample to create an unified peak-cell matrix for downstream analysis. Gene activities were extracted by extending the gene regions 2000 bp upstream. Imputation were performed by Seurat using the sc5nRNA-seq data as reference. The cell clusters are then labeled manually according to the identity of cell clusters from the sc5nRNA-seq data.








