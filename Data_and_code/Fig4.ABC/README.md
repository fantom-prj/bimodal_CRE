# Enhancer regulating gene prediction from iPSC, NSC and Neuron 

* ABC_input_and_output.R -> main code describing all the steps
* ABC_sc.sh -> code for running ABD model workflow
* config.yaml, config_biosamples_NeuronSeries_sc.tsv and files inside input folder -> supportive parameter for running ABD model workflow
* Results are presented in Fig.4
* Input data is not included here, please refer to https://fantom.gsc.riken.jp/6/suppl/Yip_et_al_2026_scTF
* For fastq and bam files, please refer to DDBJ (H3K27ac CUT&Tag: DRA019568 (DRR614873-DRR614905); Hi-C: DRA019572 (DRR614942-DRR614953); scATAC-seq: DRA019608 (DRR618513-DRR618515))

# Related Methods
* [Activity-by-Contact (ABC) model](#ABC)


# <a name="ABC"></a>Activity-by-Contact (ABC) model
Enhancer-gene interactions were predicted using the Activity-by-Contact (ABC) model. Analyses were performed independently for iPSC, NSC and Neuron using scATAC-seq, H3K27ac CUT&Tag and Hi-C datasets. The ABC pipeline was executed using the Snakemake workflow under the hg38 genome assembly. Signals were normalized using the reference quantile-normalization distribution supplied with the ABC workflow. Gene annotations were based on the hg38 CollapsedGeneBounds reference provided by the ABC package.
We supplied a custom catalog of CREs identified from this study. These regions were formatted as BED intervals and substituted directly into the ABC workflow at the candidate-region generation stage. Predicted enhancer–gene interactions were extracted from the thresholded ABC output using the default score threshold of 0.02. Only interactions exceeding this threshold were considered high-confidence enhancer–gene links for downstream analyses. For comparative analyses across differentiation, ABC scores were calculated independently in iPSC, NSC and neuron samples using the same candidate regulatory element catalog.


