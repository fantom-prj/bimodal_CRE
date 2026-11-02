# gathering analyses related to TEAD4 and ONECUT2 for downstream analyses

* TEAD_ONECUT_process.R -> main code describing all the steps
* Input data is not included here, please refer to https://fantom.gsc.riken.jp/6/suppl/Yip_et_al_2026_scTF
* For fastq and bam files, please refer to DDBJ (DRA026820 (DRR959950-DRR959961))

# Related Methods
* [Validation of TF footprinting using ChIP-seq](#AUC)
* [Chromatin accessibility of ChIP-seq peaks](#ChIP_peak)
* [Pathway enrichment analysis](#pathway)


# <a name="AUC"></a>Validation of TF footprinting using ChIP-seq
To evaluate the accuracy of TF footprinting, we compared TOBIAS-derived footprinting scores with experimentally determined TF-binding sites from ChIP-seq of TEAD4 and ONECUT2, towards footprints obtained from iPSC_Sphase and Immature_neuron cluster repectively. For each TF motif, bound and unbound motif instances identified by TOBIAS were combined, and the instance with the highest footprinting score was retained for each CRE. CREs overlapping ChIP-seq-supported binding sites were labeled as positive, while the remaining CREs were labeled as negative. Receiver operating characteristic (ROC) curves were generated using the footprinting scores to discriminate ChIP-seq-positive from ChIP-seq-negative CREs, and the area under the ROC curve (AUC) was calculated using the R package pROC.

# <a name="ChIP_peak"></a>Chromatin accessibility of ChIP-seq peaks
To assess the chromatin accessibility of TEAD4- and ONECUT2- binding regions identified from ChIP-seq. Peaks identified from the two replicates in ChIP-seq were merged and searched for motif of TEAD[1-4] and ONECUT[1-3] at a threshold of threshold of P < 1 × 10⁻⁴, according to JASPER 2024 database using motifmatchr package with hg38 in R. These peaks were considered as TEAD and ONECUT binding sites and intersected with accessible peaks defined from scATAC-seq of iPSC and Neuron respectively.

# <a name="pathway"></a>Pathway enrichment analysis
Genes directly associated with TEAD4- or ONECUT2-bound promoters, as identified by ChIP-seq, together with genes linked to TEAD4- or ONECUT2-bound enhancers through the ABC model, were used for pathway enrichment analyses. Gene Ontology (GO) enrichment was performed using topGO (v2.58.0) with the elim algorithm to account for the hierarchical structure of GO terms. KEGG and Reactome pathway enrichment analyses were performed using clusterProfiler (v4.14.6). Universe was set as the genes expressed with average log normalized expression > 0.1 from iPSC (For TEAD4) and Neuron (for ONECUT2). Unless otherwise specified, multiple-testing correction was performed using the Benjamini–Hochberg method.

