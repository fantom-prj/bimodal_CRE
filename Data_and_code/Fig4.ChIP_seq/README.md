# ChIP-seq on TEAD4 in iPSC and ONECUT2 in  neuron 

* chipseq.process.R -> main code describing all the steps
* mapping.sh -> code for BWA mapping fastq to bam file hg38 and de-duplication
* macs2.sh -> code for macs2 peak-calling each replicate with applying IgG control as background
* idr.sh -> code for idr by applying the 2 replicates 
* Input data is not included here, please refer to https://fantom.gsc.riken.jp/6/suppl/Yip_et_al_2026_scTF
* For fastq and bam files, please refer to DDBJ (DRA026820 (DRR959942-DRR959949))

# Related Methods
* [ChIP-seq](#ChIP_exp)
* [Identification of TEAD4- and ONECUT2-bound CREs from ChIP-seq data](#ChIP_analysis)


# <a name="ChIP_exp"></a>ChIP-seq
The iPS, NSC, and cortical neuron cells were harvested as described above. The cell fixation and chromatin shearing were performed by truChIP Chromatin Shearing Kit (Covaris, Cat.No. 520127) according to the manufacturer’s instructions with some modifications. Briefly, 3 × 107 cells were crosslinked with 1% methanol-free formaldehyde (Electron Microscopy Sciences) in PBS for 5 min at room temperature, followed by quenching with 200 mM glycine for 5 min. The fixed cells were lysed, and the chromatin was sheared using the Covaris S220 ultrasonicator. To validate the chromatin shearing efficiency, the sheared chromatin was de-crosslinked in 50 mM NaCl for 20 min at 100°C, followed by ramping the temperature down to 50°C. A portion of 500 ng of this DNA was examined by 1.5% agarose gel electrophoresis. After confirming the chromatin shearing efficiency within the range of 200-1200 bp, the chromatin immunoprecipitation reaction was performed using the ChIP-IT High sensitivity (Active motif, Cat.No. 53040) with Spike-in Antibody (Active motif, Cat.No. 61686) and Spike-in chromatin (Active motif, Cat.No. 53083) according to the manufacturer’s instructions. Specific and IgG isotype control antibodies for TEAD4 and ONECUT2 used for ChIP-seq are listed in Table S10. Sequence libraries were prepared using the Next Gen DNA Library Kit (Active Motif, Cat.No. 53216). The resulting libraries were quantified by Qubit 3.0 Fluorometer using Qubit dsDNA HS Assay Kit (Invitrogen, Cat.No. Q32851) and were paired-end sequenced by NovaSeqX plus (illumina) with 150 bp read length and 8 bp single-indexed reads.

# <a name="ChIP_analysis"></a>Identification of TEAD4- and ONECUT2-bound CREs from ChIP-seq data
The fastq files from ChIP-seq were mapped to the hg38 genome by BWA. The BAM files were preceded with fixmate and markdup by samtools while the bam files free of duplicates were applied to peak calling. Peak calling was performed by MACS2 at FDR 0.01 for each replicate using the pair-end bam files from specific and isotypic antibodies. The peaks of the 2 replicates were merged by bedtools as a permissive output, and were subjected to IDR for reproducible peaks at p-value < 0.05. Peaks were intersected with CREs defined by scATAC-seq for downstream analyses.




