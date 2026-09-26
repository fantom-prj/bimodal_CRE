# Pre-processing and processing of 10x scATAC-seq data generated from three samples (iPSC, NSC and Neuron)

* scATAC_processing.R -> main code describing all the analyses
* cellranger/*.sh -> associated code running cellranger
* Results are presented in Fig.1
* Input data is not included here, please refer to https://fantom.gsc.riken.jp/6/suppl/Yip_et_al_2026_scTF

# Related Methods
* [Single-cell ATAC-seq analytical pipeline](#scATAC)

# <a name="scATAC"></a>Single-cell ATAC-seq analytical pipeline

Raw data from scATAC-seq was processed through the Cell Ranger ATAC pipeline (v2.0.0) using pre-built reference from Cell Ranger refdata-cellranger-arc-GRCh38-2020-A-2.0.0. Signac was used to process the Cell Ranger output. Cells were then filtered for nCount_peaks > 1000 & < 20,000, FRiP > 0.5, is__cell_barcode equal to 1, nucleosome_signal < 4, and TSS.enrichment > 3. After quality control, 5,028 iPSCs, 5,040 NSCs, and 6,111 neurons were retained. These cells were then integrated by the lsi algorithm. The peaks called for each sample were then merged to an unified peak set (n = 463,966) which was then used to quantify reads for each cell in each sample to create an unified peak-cell matrix for downstream analysis. Gene activities were extracted by extending the gene regions 2000 bp upstream. Imputation were performed by Seurat using the sc5nRNA-seq data as reference. The cell clusters are then labeled manually according to the identity of cell clusters from the sc5nRNA-seq data.








