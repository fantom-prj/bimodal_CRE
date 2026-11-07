# Extraction of eQTL data

* eQTL.R -> main code describing all the analyses
* Input data is not included here, please refer to https://fantom.gsc.riken.jp/6/suppl/Yip_et_al_2026_scTF

# Related Methods
* [eQTL](#eQTL)



# <a name="eQTL"></a>eQTL

Processed eQTL data was obtained from the eQTL Catalog maintained by EMBL (https://www.ebi.ac.uk/eqtl/). The QTL of gene expression of a variant to a gene were filtered for the posterior inclusion probability greater than 0.3. Only the dataset related to iPSC and brain were used in this study. A full list of dataset incorporated is shown in Table S12. The variant locations were intersected with all the enhancer-like CRE regions and the Egenes were compared with the expressing genes connected with the promoter-like CREs.


