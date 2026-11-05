# Differential expression of CRE following TEAD2/4 KD in iPSC

* KD_CAGE.R -> main code describing all the analyses
* STAR_bam/STAR_CAGE.sbatch -> associated code for paired-end alignment
* SCAFE/batch_run_scafe_interactome_ssCAGE.bamtoctss.pl -> associated code for bam to ctss considering only uniquely mapped alignment
* Results are presented in Fig.4 & Ext.Fig.4
* Input data is not included here, please refer to https://fantom.gsc.riken.jp/6/suppl/Yip_et_al_2026_scTF
* For bam files please refer to DDBJ pending (pending)

# Related Methods
* [CAGE sequencing for TEAD2/4 knockdown samples](#CAGE)
* [Differentially expressed CREs identification from knockdown samples](#edgeR)



# <a name="CAGE"></a>CAGE sequencing for TEAD2/4 knockdown samples

The iPSCs sujbected to siRNA-mediated TEAD2/4 knockdown were harvested for CAGE sequencing. Six samples were collected, comprising two biological replicates each for each of three conditions including non-targeting control (NC1), and two independent TEAD2/4 siRNA combinations. Total RNA was purified as described above, and 5 μg RNA from each sample was used for library preparation according to the single-strand CAGE protocol.3 Libraries were sequenced on an Illumina NextSeq 2000 platform using paired-end 150-bp reads and a 8-bp index read.


# <a name="edgeR"></a>Differentially expressed CREs identification from knockdown samples

Paired-end CAGE reads from TEAD2/4 knockdown and control samples were aligned to the hg38 genome using STAR (v2.7.11b). The bam files were applied to SCAFE bam-to-ctss function in soft-clipping mode to identify CAGE transcription start sites (CTSS), retaining only uniquely mapped alignments. The resulting CTSS profiles were interescted with genuine TSS clusters previously identified from sc5nRNA-seq. These filtered CTSS conuts were applied to the annotated CRE regions (derived from scATAC-seq) to generate a CRE-level transcription count matrix.<br>
Differential expression was assessed using edgeR, comparing two replicates of the non-targeting (siNC) with two replicates each of two independent TEAD2/4 siRNA combinations. A generalized linear model was fitted with separate coefficients for the three experimental conditions, and the average effect of the two knockdown combinations relative to the control was tested using the contrast of “-1, 0.5, 0.5”. Lowly expressed CREs were filtered by “filterByExpr()” function in edgeR, retaining 21,599 evaluable CREs for differential analysis. The CREs with |log2FC| > 0.5 and FDR < 0.05 were considered significantly differentially expressed.

```
[1] STAR \
--twopassMode None \
--runThreadN 12 \
--outFilterMultimapNmax 20 \
--outSAMtype BAM SortedByCoordinate \
--outSAMunmapped None \
--outReadsUnmapped Fastx \
--genomeDir STAR_index_hg38_gencode39 \
--readFilesIn "$R1" "$R2" \
--readFilesCommand zcat \
--outFileNamePrefix "map/${name}." \
--alignIntronMin 20 \
--alignIntronMax 1000000 \
--limitBAMsortRAM 30000000000

[2] scafe.tool.bk.bam_to_ctss \
--TSS_mode=softclip \
--bamPath=$bamPath \
--unencoded_G_upstrm_nt=3 \
--max_thread=5 \
--min_MAPQ=225 \
--genome=hg38.gencode_v39 \
--max_softclip_length=3 \
--outputPrefix=$outputPrefix \
--outDir=$outDir
```




