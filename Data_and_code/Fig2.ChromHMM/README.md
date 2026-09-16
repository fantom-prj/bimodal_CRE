# chromatin state annotation on the genome based on histon marks and CTCF profiling from iPSC, NSC and Neuron 

* chromHMM.R -> main code describing all the analyses
* Results are presented in Fig.2
* Input data is not included here, please refer to https://fantom.gsc.riken.jp/6/suppl/Yip_et_al_2026_scTF
* bam files please refer to DDBJ DRA019568 (DRR614873-DRR614905)

# Related Methods
* [Chromatin state modeling](#chromatinstate)
* [TATA box prediction from non-transcribed CREs](#TATApredict)
* [Enrichment of chromatin states with transcription and regulatory elements ](#enrichment)
* [Genomic coordinate of CREs](#gcoordinate)

# <a name="chromatinstate"></a>Chromatin state modeling 

The bulk CUT&Tag data was processed by the ENCODE ATAC-seq pipeline (v1.7.0). Briefly, reads were mapped to hg38 human genome using Bowtie2 (v2.2.6). Following the removal of PCR duplicates, peak calling was performed on merged replicates for each cell type using MACS2 (v2.1.0) for peak calling. The significant peaks (p < 0.01) were merged and extended into a minimum window size of 150bp. 
Chromatin state was defined by ChromHMM (v1.24) at 200 nt resolution. The model integrated CUT&Tag data (with corresponding IgG controls) and scATAC-seq data, divided as three samples (iPSC, NSC, and neurons). Both CUT&Tag and scATAC-seq data were applied as mapped reads (bam files). A 16-state model was trained using the LearnModel function. States were manually annotated based on their enrichment for specific histone marks and external genomic features. Additionally, the summit coordination of the SCAFE TSS clusters derived from sc5nRNA-seq and the midpoint locations of SCREEN cCREs indicating promoter, enhancer and CTCF were applied as annotation to support curation.
```
[1] java -mx6000M -jar ChromHMM.jar BinarizeBed -gzip -b 200 hg38.txt bed_file_path input_summary.txt binary_folder
[2] java -mx6000M -jar ChromHMM.jar LearnModel -p 8 binary_folder output folder 16 hg38
```

