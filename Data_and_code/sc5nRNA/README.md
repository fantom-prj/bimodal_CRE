# Pre-processing and processing of 10x sc5nRNA-seq data generated from three samples (iPSC, NSC and Neuron)

* sc5nRNA_seq_process.R -> main code describing all the analyses
* cellranger/*.sh -> associated code running cellranger
* scRNA_pre_processing.Rmd -> associated code for pre-process cellranger output into seurat object
* scRNA_QC_filtering.Rmd -> associated code for filtering cells
* scRNA_integration.Rmd -> associated code for integrating 3 samples into one seurat object
* scRNA_monocle.Rmd -> associated code for pseudotime
* skewC.R -> associated code for SkewC
* Results are presented in Fig.1
* Input data is not included here, please refer to https://fantom.gsc.riken.jp/6/suppl/Yip_et_al_2026_scTF

# Related Methods
* [Chromium for sc5nRNA-seq](#10xsep)
* [Single-cell RNA-seq analytical pipeline](#scRNA)

# <a name="10xsep"></a>Chromium for sc5nRNA-seq
For the sc5nRNA-seq, iPSC, NSC and cortical neuron were harvested at the same time as described above. The cells were resuspended in PBS with 0.1% BSA followed by filtration with 40 μm FLOWMI Tip strainers (BEL-ART). The cells were centrifuged at 150 × g for 10 min to remove the dead cell debris. Cell pellets were resuspended in PBS with 0.1% BSA at a stock concentration within the recommended range and loaded at a volume to target 4,500 cells. We generated two technical replicates per cell type by loading two channels of Chip G with the same cell mix. Libraries were prepared using the Chromium Next GEM Single Cell 5' Library and Gel Bead Kit v1.1 (PN: 1000165) according to the manufacturer’s instructions. The cDNA libraries were then indexed by PCR (iPS_rep1: SI-GA-E4, iPS_rep2: SI-GA-F4; NSC_rep1: SI-GA-C4, NSC_rep2: SI-GA-D4; Neuron_rep1: SI-GA-G4, Neuron_rep2: SI-GA-H4). After verification from Agilent Bioanalyzer with a high-sensitivity chip (Agilent Technologies), libraries were sequenced on a MGI sequencer using the following sequencing conditions: R1: 111 cycles (dark cycles 27-39); R2: 91 cycles; Index1: 8 cycles.

# <a name="scRNA"></a>Single-cell RNA-seq analytical pipeline
Raw data from the sc5nRNA-seq were applied to Cell Ranger (v6.0.0) using human genome reference refdata-gex-GRCh38-2020-A. Seurat (v4.3.0) was used to manage the Cell Ranger output. Scrublet4 was used to identify potential doublets. After removing the doublets, single cells were filtered to have minUMI ≥ 2000, minGenes ≥ 1000 and maxPctMt < 20. SkewC was applied to annotate the filtered cells as an additional quality control. SkewC identifies skewed cells based on gene body coverage skewness, which was observed in terms of coverage bias towards specific gene regions (5’ end, 3’ end, or middle of the gene body). After filtering, cells from the three cell types were integrated by Harmony and clustered by the Louvain algorithm with a resolution of 0.8. The cell clusters are then labeled manually according to the identity of cell clusters predicted by Azimuth. Normalized gene expression level per single-cell was obtained by Seurat with “LogNormalize”. <br>
Pseudotime trajectories were inferred from the sc5nRNA-seq dataset using Monocle. Harmony-corrected PCA and UMAP embeddings were used for trajectory construction, and a principal graph was learned across all cells without partitioning. The trajectory was rooted at the principal graph node containing the largest number of undifferentiated iPSC cells, and pseudotime values were assigned using order_cells.








