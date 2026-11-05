# TSS, TSS cluster and tCRE identification from sc5nRNA-seq data

* SCAFE_summary.R -> main code describing all the analyses [procedure in running SCAFE, SCAFE scafe.tool.sc.count function to get TSS count into CRE at single-cell level]
* sc5end/solo/2023070722.batch_run_scafe.F6_neuron_single_cell_solo.pl -> associated code getting CTSS from individual libraries
* sc5end/aggregate/run_full/20230711.batch_run_scafe_interactome.scafe_aggregate.sc5end.pl -> associated code for getting SCAFE tCRE across 3 cell types
* sc5end/count/20230711.batch_run_scafe_interactome.scafe_count.pl -> associated code for counting TSS to CRE and generate single cell matrix of TSS
* Input data is not included here, please refer to https://fantom.gsc.riken.jp/6/suppl/Yip_et_al_2026_scTF
* For bam files please refer to DDBJ DRA019607 (DRR618507-DRR618512)

# Related Methods
* [Identification of confident TSS cluster](#tCRE)
* [Establishment of CREs with bimodal signals](#count)



# <a name="tCRE"></a>Identification of confident TSS cluster

The output files generated from the Cell Ranger pipeline were applied to SCAFE (v1.0.1)12 to identify genuine transcription start sites (TSSs) and clusters. Briefly, TSSs were extracted to identify genuine TSS clusters from each replicate by running scafe.workflow.sc.solo. These clusters were then merged by running scafe.tool.cm.aggregate. Only the TSSs lying inside the genuine TSS clusters were used for quantification, which accounted for 68.4% of the total UMI. Transcribed CREs (tCREs) were indepednently identified from the three cell types using SCAFE. These tCRE regions were used as comparison (in Extended Fig. 2a) alone but not for final quantification.


# <a name="count"></a>Establishment of CREs with bimodal signals

The ATAC peak regions defined from the Cell Ranger ATAC pipeline were merged across the 3 cell-types. These merged peak regions represent the CRE regions (n = 463,866 from main chromosomes) of this study. The ATAC signal across all the single cells were recounted using the fragment files by CreateChromatinAssay of the Signac package as described above. To quantify the transcription activity per cell per CRE, CRE regions were applied to the scafe.tool.sc.count of the SCAFE package, counting the TSS located inside the genuine TSS clusters. Across the three cell-types, positions with the highest TSS count within CRE regions were considered as TSS summits. These steps generated ATAC and TSS quantification profiles across 463,866 CREs for each single-cell.



