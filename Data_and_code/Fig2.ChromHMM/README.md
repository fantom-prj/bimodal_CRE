# Chromatin state annotation on the genome based on histone marks and CTCF profiling from iPSC, NSC and Neuron 

* chromHMM.R -> main code describing all the analyses
* bin.20231108.sh -> associated code for generating 200 nt bin files from bam files
* data_files/learnmodel.20240305.sh -> associated code for clustering the bins
* Results are presented in Fig.2
* Input data is not included here, please refer to https://fantom.gsc.riken.jp/6/suppl/Yip_et_al_2026_scTF
* bam files please refer to DDBJ DRA019568 (DRR614873-DRR614905)

# Related Methods
* [Chromatin state modeling](#chromatinstate)
* [Hierarchical CRE classification](#classification)

# <a name="chromatinstate"></a>Chromatin state modeling 

Bulk CUT&Tag data on the genomic profiles of CTCF, H3K4me3, H3K27me3, H3K4me1, and H3K27ac from the three samples (iPSC, NSC, and neuron) were used to define chromatin states genome wide. Briefly, reads were mapped to hg38 human genome using Bowtie2 (v2.2.6). Following the removal of PCR duplicates, the BAM files were applied to ChromHMM (v1.24)10 to define chromatin states at 200 nt resolution. Corresponding IgG controls were also subjected to ChromHMM. In addition, scATAC-seq data were divided as three samples (iPSC, NSC, and neurons) and applied as mapped reads (BAM files). After all the BAM files were binarized at 200 nt resoultion, a 16-state model was trained using the LearnModel function. To support curation, the summit coordination of the SCAFE TSS clusters derived from sc5nRNA-seq and the midpoint locations of SCREEN cCREs indicating promoter, enhancer and CTCF were applied as annotation. These 200nt-bins were annotated as 1) promoter, 2) promoter flanked, 3) bivalent promoter, 4) active enhancer, 5) repressed enhancer, 6) primed enhancer, 7) CTCF-alone, 8) genomic region and 9) repressed region.
```
[1] java -mx6000M -jar ChromHMM.jar BinarizeBed -gzip -b 200 hg38.txt bed_file_path input_summary.txt binary_folder
[2] java -mx6000M -jar ChromHMM.jar LearnModel -p 8 binary_folder output folder 16 hg38
```

# <a name="classification"></a>Hierarchical CRE classification

The final list of CREs (merged peaks derived from scATAC-seq) were assigned with chromatin states in a sample specific manner. For CREs overlapping multiple ChromHMM states, a representative state was selected based on a prioritized hierarchy: 1) promoter / promoter flanked, 2) bivalent promoter, 3) active enhancer, 4) repressed enhancer, 5) primed enhancer, and 6) CTCF-alone. By combining the chromatin states across the 3 samples, the CREs were ultimately categorized as, in the order of, “promoter-like”, "enhancer-like”, “CTCF-alone” or “unclassed”.

