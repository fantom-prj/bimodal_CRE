# Methodology for "Bimodal Single-Cell Profiling Reveals the Uncoupled Regulation of Enhancers and Promoters in Human Neural Differentiation"


Transcription factors (TFs) establish cellular identity through coordinated regulation of chromatin accessibility and transcription, yet the extent to which these regulatory outputs are coupled across cis-regulatory elements (CREs) remains poorly understood. Here, we integrated single-cell 5’ RNA-seq and single-cell ATAC-seq across human neural differentiation to dynamically resolve chromatin accessibility and transcription. TF motif activities derived independently from these two modalities revealed widespread uncoupling, with only 30.6% of TF motifs identified through both accessibility- and transcription-derived activities. These bimodal TF activities further revealed distinct TF regulatory modes associated with chromatin accessibility or transcription, with enhancers shown as major sites of TF-mediated chromatin regulation. TEAD and ONECUT emerged as regulators of enhancer accessibility despite their occupancy at both promoters and enhancers. To resolves the regulatory relationship among TFs, enhancers and promoters, we reconstructed a gene regulatory network (GRN) integrating bimodal CREs accessibility and transcriptional activities. Notably, transcribed enhancers occupied highly connected positions within the and revealed regulatory routes that were not recovered from chromatin accessibility alone. Together, these findings establish chromatin accessibility and transcription as complementary but non-interchangeable dimensions of CRE activity and provide a framework for resolving TF regulatory modes and enhancer-mediated gene regulation during cell-state transitions.

## Method to code links

| Method | Code |
| --- | --- |
| Single-cell RNA-seq analytical pipeline | [code](https://github.com/fantom-prj/bimodal_CRE/tree/main/Data_and_code/sc5nRNA) |
| Single-cell ATAC-seq analytical pipeline | [code](https://github.com/fantom-prj/bimodal_CRE/tree/main/Data_and_code/scATAC) |
| Co-embedding and metacell grouping | [code](https://github.com/fantom-prj/bimodal_CRE/tree/main/Data_and_code/Fig1.coembed) |
| Chromatin state modeling | [code](https://github.com/fantom-prj/bimodal_CRE/tree/main/Data_and_code/Fig2.ChromHMM) |
| Identification of confident TSS clusters | [code](https://github.com/fantom-prj/bimodal_CRE/tree/main/Data_and_code/SCAFE) |
| Establishment of CREs with bimodal signals | [code](https://github.com/fantom-prj/bimodal_CRE/tree/main/Data_and_code/SCAFE) |
| Aggregation of single-cell signals into metacells | [code](https://github.com/fantom-prj/bimodal_CRE/tree/main/Data_and_code/Fig1.metacell) |
| Manual curation of TF function | [code](https://github.com/fantom-prj/bimodal_CRE/tree/main/Data_and_code/Fig1.metacell) |
| Hierarchical CRE classification | [code](https://github.com/fantom-prj/bimodal_CRE/tree/main/Data_and_code/Fig2.CRE_feature) |
| Regulatory element annotation for transcribed and non-transcribed CREs | [code](https://github.com/fantom-prj/bimodal_CRE/tree/main/Data_and_code/Fig2.CRE_feature) |
| TATA box prediction from non-transcribed CREs | [code](https://github.com/fantom-prj/bimodal_CRE/tree/main/Data_and_code/Fig2.CRE_feature) |
| Enrichment of chromatin states with transcription and regulatory elements | [code](https://github.com/fantom-prj/bimodal_CRE/tree/main/Data_and_code/Fig2.CRE_feature) |
| Genomic coordinates of CREs | [code](https://github.com/fantom-prj/bimodal_CRE/tree/main/Data_and_code/Fig2.CRE_feature) |
| Super-enhancer identification | [code](https://github.com/fantom-prj/bimodal_CRE/tree/main/Data_and_code/Fig2.SE_identification) |
| Connecting genes to CREs | [code](https://github.com/fantom-prj/bimodal_CRE/tree/main/Data_and_code/Fig2.CRE_feature) |
| Identification of TEAD4- and ONECUT2-bound CREs from ChIP-seq data | [code](https://github.com/fantom-prj/bimodal_CRE/tree/main/Data_and_code/Fig4.ChIP_seq) |
| Identification of differentially accessible CREs in knockdown samples | [code](https://github.com/fantom-prj/bimodal_CRE/tree/main/Data_and_code/Fig4.KD_ATAC) |
| Identification of differentially expressed CREs in knockdown samples | [code](https://github.com/fantom-prj/bimodal_CRE/tree/main/Data_and_code/Fig4.KD_CAGE) |
| Motif activity analysis | [code](https://github.com/fantom-prj/bimodal_CRE/tree/main/Data_and_code/Fig3.ChromVar) |
| Cross-correlation between TF expression and motif activity | [code](https://github.com/fantom-prj/bimodal_CRE/tree/main/Data_and_code/Fig3.ChromVar) |
| Identification of evaluable TF motifs | [code](https://github.com/fantom-prj/bimodal_CRE/tree/main/Data_and_code/Fig3.ChromVar) |
| Transcription factor footprinting | [code](https://github.com/fantom-prj/bimodal_CRE/tree/main/Data_and_code/Fig4.Foot_Printing) |
| Validation of TF footprinting using ChIP-seq | [code](https://github.com/fantom-prj/bimodal_CRE/tree/main/Data_and_code/Fig4.TEAD_ONECUT) |
| Chromatin accessibility of ChIP-seq peaks | [code](https://github.com/fantom-prj/bimodal_CRE/tree/main/Data_and_code/Fig4.TEAD_ONECUT) |
| Pathway enrichment analysis | [code](https://github.com/fantom-prj/bimodal_CRE/tree/main/Data_and_code/Fig4.TEAD_ONECUT) |
| Activity-by-Contact (ABC) model | [code](https://github.com/fantom-prj/bimodal_CRE/tree/main/Data_and_code/Fig4.ABC) |
| eQTL | code |
| Network reconstruction | code |
| Network variants and external validation | code |
| Filtering according to significance and goodness of fit | code |
| Global topology of TF-to-TF relationships | code |
| Cumulative TF effects and TF regulatory modes | code |
| Regulatory decoupling between accessibility and transcription | code |
| Enhancer roles and epigenomic feature enrichment | code |
| Concordance with Activity-by-Contact predictions | code |
| Contribution of enhancer transcriptional activity | code |
| Enhancer centrality and CRISPRi screen association | code |
| Multiple testing | code |


## Software and Tools

| Tool | Version | Source |
| --- | --- | --- |
| Cell Ranger | 6.0.0 | https://github.com/10XGenomics/cellranger |
| Seurat | 4.3.0 | https://github.com/satijalab/seurat |
| Scrublet | 0.2.1 | https://github.com/swolock/scrublet |
| SkewC | — | https://github.com/LSBDT/SkewC |
| Azimuth | — | https://github.com/satijalab/azimuth-references |
| Cell Ranger ATAC | 2.0.0 | https://www.10xgenomics.com/support/software/cell-ranger-atac/downloads |
| scDART | 0.1dev | https://github.com/PeterZZQ/scDART |
| SEACells | 0.3.3 | https://github.com/dpeerlab/SEACells |
| SCAFE | 1.0.1 | https://github.com/chung-lab/SCAFE |
| STAR | 2.7.11b | https://github.com/alexdobin/STAR |
| subset-bam | 1.1.0 | https://github.com/10XGenomics/subset-bam |
| TOBIAS | 0.16.1 | https://github.com/loosolab/TOBIAS |
| Homer | 4.11 | https://github.com/javrodriguez/HOMER |
| MEME-suite | 5.5.7 | https://meme-suite.org/meme/meme-software |
| Bowtie2 | 2.2.6 | https://github.com/BenLangmead/bowtie2 |
| MACS2 | 2.1.0 | https://github.com/jdavisturak/MACS2-2.1.1.20160309 |
| ROSE | 1.3.1 | https://github.com/stjude/ROSE |
| ChromHMM | 1.24 | https://ernstlab.github.io/ChromHMM/ |
| ABC-model | 1.1.2 | https://github.com/broadinstitute/ABC-Enhancer-Gene-Prediction |
| topGO | 2.58.0 | https://bioconductor.org/packages/release/bioc/html/topGO.html |
| clusterProfiler | 4.14.6 | https://bioconductor.org/packages/release/bioc/html/clusterProfiler.html |
| Subread | 2.0.6 | https://subread.sourceforge.net/ |
| BWA | 0.7.18-r1243-dirty | https://github.com/lh3/bwa |
| IDR | 2.0.2 | https://github.com/nboley/idr |
| Samtools | 1.11 | https://www.htslib.org/ |
| Bedtools | 2.30.0 | https://bedtools.readthedocs.io/en/latest/ |
| Tabix, Bgzip | 1.15.1 | https://www.htslib.org/ |
| bedGraphToBigWig | 2.8 | https://github.com/ENCODE-DCC/kentUtils/tree/master/src/utils |
| Bedparse | 0.2.3 | https://github.com/tleonardi/bedparse |
| R | 4.3.1 | https://cran.r-project.org/ |
| Perl | 5.26.2 | https://www.perl.org/get.html |


## Data Availability

### Core Data

| Dataset | Accession numbers |
| --- | --- |
| Single-cell 5′ RNA-seq | DRA019607 (DRR618507–DRR618512) |
| Single-cell ATAC-seq (iPSC, NSC, Neuron) | DRA019608 (DRR618513–DRR618515) |
| CUT&Tag (iPSC, NSC, Neuron; H3K27ac, H3K27me3, H3K4me1, H3K4me3, CTCF and control) | DRA019568 (DRR614873–DRR614905) |
| ChIP-seq (TEAD4 in iPSC; ONECUT2 in Neuron) | DRA026820 (DRR959942–DRR959949) |
| Bulk ATAC-seq after knockdown (TEAD2/4 in iPSC; ONECUT2 in Neuron) | DRA026820 (DRR959950–DRR959961) |
| Bulk CAGE after knockdown (TEAD2/4 in iPSC) | Pending |
| Hi-C (iPSC, NSC, Neuron) | DRA019572 (DRR614942–DRR614953) |


### Public Resources

| Resource | URL |
| --- | --- |
| eQTL Catalogue | https://www.ebi.ac.uk/eqtl/ |
| GWAS from CAUSALdb | http://www.mulinlab.org/causaldb/index.html |
| cCREs from SCREEN | https://screen.encodeproject.org/ |
| Disease genes from Open Targets | https://platform.opentargets.org/ |
| JASPAR 2024 database | https://jaspar2024.elixir.no/ |


## References
1.	Wang, C. et al. Scalable Production of iPSC-Derived Human Neurons to Identify Tau-Lowering Compounds by High-Content Screening. Stem Cell Reports 9, 1221–1233 (2017).
2.	Yip, C. W. et al. CFC-seq: identification of full-length capped RNAs unveil enhancer-derived transcription. Preprint at https://doi.org/10.1101/2024.10.31.620483 (2024).
3.	Takahashi, H., Nishiyori-Sueki, H., Ramilowski, J. A., Itoh, M. & Carninci, P. Low Quantity Single Strand CAGE (LQ-ssCAGE) Maps Regulatory Enhancers and Promoters. Methods Mol Biol 2351, 67–90 (2021).
4.	Hao, Y. et al. Integrated analysis of multimodal single-cell data. Cell 184, 3573-3587.e29 (2021).
5.	Wolock, S. L., Lopez, R. & Klein, A. M. Scrublet: Computational Identification of Cell Doublets in Single-Cell Transcriptomic Data. Cell Syst 8, 281-291.e9 (2019).
6.	Abugessaisa, I. et al. SkewC: Identifying cells with skewed gene body coverage in single-cell RNA sequencing data. iScience 25, 103777 (2022).
7.	Abugessaisa, I., Hasegawa, A., Katayama, S., Kere, J. & Kasukawa, T. Computational approach to evaluate scRNA-seq data quality and gene body coverage with SkewC. STAR Protoc 4, 102038 (2023).
8.	Stuart, T., Srivastava, A., Lareau, C. A. & Satija, R. Single-cell chromatin state analysis with Signac. Nat Methods 18, 1333–1341 (2021).
9.	Zhang, Z., Yang, C. & Zhang, X. scDART: integrating unmatched scRNA-seq and scATAC-seq data and learning cross-modality relationship simultaneously. Genome Biol 23, 139 (2022).
10.	Persad, S. et al. SEACells infers transcriptional and epigenomic cellular states from single-cell genomics data. Nat Biotechnol 41, 1746–1757 (2023).
11.	Ernst, J. & Kellis, M. ChromHMM: automating chromatin-state discovery and characterization. Nat Methods 9, 215–216 (2012).
12.	Moody, J. et al. SCAFE: a software suite for analysis of transcribed cis-regulatory elements in single cells. Bioinformatics 38, 5126–5128 (2022).
13.	Heinz, S. et al. Simple combinations of lineage-determining transcription factors prime cis-regulatory elements required for macrophage and B cell identities. Mol Cell 38, 576–589 (2010).
14.	Whyte, W. A. et al. Master transcription factors and mediator establish super-enhancers at key cell identity genes. Cell 153, 307–319 (2013).
15.	Schep, A. N., Wu, B., Buenrostro, J. D. & Greenleaf, W. J. chromVAR: inferring transcription-factor-associated accessibility from single-cell epigenomic data. Nat Methods 14, 975–978 (2017).
16.	Li, H. & Durbin, R. Fast and accurate long-read alignment with Burrows-Wheeler transform. Bioinformatics 26, 589–595 (2010).
17.	Danecek, P. et al. Twelve years of SAMtools and BCFtools. GigaScience 10, giab008 (2021).
18.	Zhang, Y. et al. Model-based analysis of ChIP-Seq (MACS). Genome Biol 9, R137 (2008).
19.	Li, Q., Brown, J. B., Huang, H. & Bickel, P. J. Measuring reproducibility of high-throughput experiments. https://doi.org/10.48550/ARXIV.1110.4705 (2011) doi:10.48550/ARXIV.1110.4705.
20.	Bentsen, M. et al. ATAC-seq footprinting unravels kinetics of transcription factor binding during zygotic genome activation. Nat Commun 11, 4267 (2020).
21.	Liao, Y., Smyth, G. K. & Shi, W. The Subread aligner: fast, accurate and scalable read mapping by seed-and-vote. Nucleic Acids Res 41, e108 (2013).
22.	Fulco, C. P. et al. Activity-by-contact model of enhancer–promoter regulation from thousands of CRISPR perturbations. Nat Genet 51, 1664–1669 (2019).
23.	Kamimoto, K. et al. Dissecting cell identity via network inference and in silico gene perturbation. Nature 614, 742–751 (2023).
24.	Bravo González-Blas, C. et al. SCENIC+: single-cell multiomic inference of enhancers and gene regulatory networks. Nat Methods 20, 1355–1367 (2023).
25.	Fleck, J. S. et al. Inferring and perturbing cell fate regulomes in human brain organoids. Nature 621, 365–372 (2023).
26.	Friedman, J., Hastie, T. & Tibshirani, R. Regularization Paths for Generalized Linear Models via Coordinate Descent. J. Stat. Soft. 33, (2010).
27.	Kerimov, N. et al. eQTL Catalogue: A Compendium of Uniformly Processed Human Gene Expression and Splicing QTLs. http://biorxiv.org/lookup/doi/10.1101/2020.01.29.924266 (2020) doi:10.1101/2020.01.29.924266.
28.	Wang, J. et al. CAUSALdb: a database for disease/trait causal variants identified using summary statistics of genome-wide association studies. Nucleic Acids Research gkz1026 (2019) doi:10.1093/nar/gkz1026.
29.	The ENCODE Project Consortium et al. Expanded encyclopaedias of DNA elements in the human and mouse genomes. Nature 583, 699–710 (2020).
30.	Koscielny, G. et al. Open Targets: a platform for therapeutic target identification and validation. Nucleic Acids Res 45, D985–D994 (2017).
31.	Rauluseviciute, I. et al. JASPAR 2024: 20th anniversary of the open-access database of transcription factor binding profiles. Nucleic Acids Res 52, D174–D182 (2024).
32.	Dobin, A. et al. STAR: ultrafast universal RNA-seq aligner. Bioinformatics 29, 15–21 (2013).
33.	Bailey, T. L. & Grant, C. E. SEA: Simple Enrichment Analysis of motifs. Preprint at https://doi.org/10.1101/2021.08.23.457422 (2021).
34.	Langmead, B. & Salzberg, S. L. Fast gapped-read alignment with Bowtie 2. Nat Methods 9, 357–359 (2012).
35.	Alexa, A., Rahnenführer, J. & Lengauer, T. Improved scoring of functional groups from gene expression data by decorrelating GO graph structure. Bioinformatics 22, 1600–1607 (2006).
36.	Xu, S. et al. Using clusterProfiler to characterize multiomics data. Nat Protoc 19, 3292–3320 (2024).
37.	Quinlan, A. R. & Hall, I. M. BEDTools: a flexible suite of utilities for comparing genomic features. Bioinformatics 26, 841–842 (2010).
38.	Li, H. Tabix: fast retrieval of sequence features from generic TAB-delimited files. Bioinformatics 27, 718–719 (2011).
39.	Kent, W. J. et al. The human genome browser at UCSC. Genome Res 12, 996–1006 (2002).
40.	Leonardi, T. Bedparse: feature extraction from BED files. JOSS 4, 1228 (2019).
41.	R Core Team. R: A language and environment for statistical computing. (2021).




