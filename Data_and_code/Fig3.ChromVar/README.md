# Motif activity analyses at metacell level

* chromvar_metacell.R -> main code describing all the steps
* Input data is not included here, please refer to https://fantom.gsc.riken.jp/6/suppl/Yip_et_al_2026_scTF
* For fastq and bam files, please refer to DDBJ (DRA026820 (DRR959950-DRR959961))

# Related Methods
* [Motif activity analysis](#motif)
* [Cross-correlation between TF expression and motif activity](#corr)
* [Identification of evaluable TF Motifs](#evaluate)


# <a name="motif"></a>Motif activity analysis
Transcription factor motif activity was estimated across 393 metacells with the bimodal signals of CREs using chromVar (v1.20.2). These CREs were divided into promoter-like CRE linked to genes (n = 20,959) and enhancer-like (n = 20,951) according to the overall chromatin states. To account for the distinct regulatory properties of promoters and enhancers, chromVar was performed in a stratified manner across four conditions (promoter-ATAC, promoter-TSS, enhancer-ATAC, and enhancer-TSS). For each condition, background peak sets were independently generated to match the GC content and average accessibility specific to the included CRE subset. This approach prevented the high basal signal of promoters from masking the dynamic variability of enhancer-associated motifs. TF motifs were filtered for biological relevance based on two criteria: a bootstrap-adjusted p-value < 0.01 and a variability score exceeding the median for each specific analysis set. This rigorous filtering ensured that motifs with variability indistinguishable from the background null model were excluded. The resulting deviation Z-scores were extracted and utilized as a quantitative measure of TF motif activity for subsequent analyses.

# <a name="corr"></a>Cross-correlation between TF expression and motif activity
To quantify the relationship between TF mRNA levels (log normalized expression) and their corresponding motif activities (deviation Z-scores), we calculated cross-correlation coefficients across 304 metacells ordered by pseudotime from trajectory 1. Given the potential temporal delay between TF expression and subsequent changes in chromatin accessibility or transcription at CREs, we utilized the ccf function in R to account for pseudotime lags. A maximum lag of five metacells was allowed, and the peak absolute correlation coefficient across the trajectory was selected to represent the maximum regulatory coupling for each TF-motif pair.

# <a name="evaluate"></a>Identification of evaluable TF Motifs
Evaluable TF motifs were defined based on a stringent three-tier filtering pipeline applied independently across the four experimental conditions (Promoter-ATAC, Promoter-TSS, Enhancer-ATAC, and Enhancer-TSS). Only motifs satisfying all three criteria within a specific condition were considered evaluable and included in the final downstream regulatory analyses. These three criteria are:<br>
1.	Expression support: As described above, the corresponding TF(s) must reach a normalized expression level of 0.2 in at least three metacells. For dimeric motifs, both constituent TFs were required to meet this threshold.<br>
2.	Motif activity variability: Motifs were required to exhibit significant chromatin or transcriptional variability, defined by a chromVAR bootstrap-adjusted p-value (< 0.01) and a variability score exceeding the median for the respective analysis set.<br>
3.	Expression-activity correlation: The peak absolute cross-correlation coefficient between TF expression and motif activity (Z-score) must exceed 0.5.<br>





