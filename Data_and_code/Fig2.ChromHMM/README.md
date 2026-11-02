# Chromatin state annotation on the genome based on histone marks and CTCF profiling from iPSC, NSC and Neuron 

* chromHMM.R -> main code describing all the analyses
* independent_peak_call.R -> associated code for peak calling. These regions were not used in downstream analysis
* bin.20231110v4.sh -> associated code for generating 200 nt bin files from bed files
* data_files/learnmodel.20240305.sh -> associated code for clustering the bins
* Results are presented in Fig.2
* Input data is not included here, please refer to https://fantom.gsc.riken.jp/6/suppl/Yip_et_al_2026_scTF
* For bam files please refer to DDBJ DRA019568 (DRR614873-DRR614905)

# Related Methods
* [CUT&Tag](#CUT&Tag)
* [Chromatin state modeling](#chromatinstate)

# <a name="CUT&Tag"></a>CUT&Tag

The iPSC, NSC, and cortical neuron cells were harvested as described above. Then, the cells were subjected to nuclei isolation according to a protocol established by 10x Genomics (Document CG000169 Rev D), except the concentration of BSA in PBS was optimized to 0.1%. The isolated nuclei from two biological replicates (5 x 105 nuclei for each antibody) were used for library preparation by CUT&Tag-IT Assay Kit (Active motif, Cat.No. 53160), according to the manufacturer’s instructions with a slight modification. Briefly, the nuclei were immobilized to Concanavalin A coated beads for 10 min on an end-over-end rotator at room temperature, followed by magnetic separation to clear the liquid. The nuclei were resuspended with ice-cold antibody buffer supplemented with protease inhibitor cocktail and digitonin. Then, 1 μg specific antibody including CTCF, H3K4me3, H3K27me3, H3K4me1, and H3K27ac (Table S10) or 1 μg rabbit IgG isotype control antibody was added to each sample. After incubation overnight at 4°C, the nuclei were isolated by magnetic separation. The samples were incubated with 1:100 diluted guinea pig anti-rabbit secondary antibody in Dig-Wash buffer for 1 hr at room temperature, followed by washes with Dig-Wash buffer. The antibody-bound chromatin was incubated with the CUT&Tag-IT assembled pA-Tn5 transposomes in Dig-300 buffer for 1 hr at room temperature, followed by stringent washes with Dig-300 buffer. The antibody-bound chromatin was subsequently sheared, and sequence adaptors were added in a single step by tagmentation for 1 hr at 37°C. The samples were de-crosslinked with buffers containing 0.1% SDS, 16mM EDTA, and 1 μg Proteinase K for 1 hr at 55°C, followed by DNA extraction with the kit-supplied column purification. The libraries were amplified by PCR with indexing primers, and size selection was performed with SPRI beads to remove primer dimers. The size distribution and quantification of the libraries were assessed by the D1000 ScreenTape system (Agilent, Cat.No. 5067) and KAPA Library Quantification Kit (KAPA, Cat.No. KK4828), respectively. Sequencing was performed by HiSeqX_Ten (illumina) with the following conditions: R1: 150 cycles; R2: 150 cycles; Index1: 8 cycles; Index2: 8 cycles.


# <a name="chromatinstate"></a>Chromatin state modeling 

Bulk CUT&Tag data on the genomic profiles of CTCF, H3K4me3, H3K27me3, H3K4me1, and H3K27ac from the three samples (iPSC, NSC, and neuron) were used to define chromatin states genome wide. Briefly, reads were mapped to hg38 human genome using Bowtie2 (v2.2.6). Following the removal of PCR duplicates, the BAM files were converted into BED format and applied to ChromHMM (v1.24) to define chromatin states at a 200-nt resolution. Corresponding BED files of IgG controls were also subjected to ChromHMM. In addition, scATAC-seq data were divided as three samples (iPSC, NSC, and neurons) and applied as mapped reads (BED files) without negative control. After all the signals were binarized at a 200-nt resoultion, a 16-state model was trained using the LearnModel function. To support curation, the summit coordination of the SCAFE TSS clusters derived from sc5nRNA-seq and the midpoint locations of SCREEN cCREs indicating promoter, enhancer and CTCF were applied as annotation. These 200nt-bins were annotated as 1) promoter, 2) promoter flanked, 3) bivalent promoter, 4) active enhancer, 5) repressed enhancer, 6) primed enhancer, 7) CTCF-alone, 8) genomic region and 9) repressed region.
```
[1] java -mx6000M -jar ChromHMM.jar BinarizeBed -gzip -mixed hg38.txt [BED_path] Neuron.cellmarkfile.v4.txt [binary_folder]
[2] java -mx6000M -jar ChromHMM.jar LearnModel -p 8 [binary_folder] [output folder] 16 hg38
```



