# Procedures of embedding sc5nRNA-seq and acATAC-seq single-cells & obtaining metacell with SEACells

* coembed.R -> main code describing input and output of scDART embedding
* seacell.R -> main code describing input and output of SEACell analysis
* *.scdart.py -> associated code for running scDART in each cell cluster
* *.SEACell.py -> associated code for running SEACell in each cell cluster
* Results are presented in Fig.1 and ext_Fig.1
* Input data is not included here, please refer to https://fantom.gsc.riken.jp/6/suppl/Yip_et_al_2026_scTF

# Related Methods
* [Co-embedding and metacell grouping](#coembed)


# <a name="coembed"></a>Co-embedding and metacell grouping

Single cells from sc5nRNA-seq and scATAC-seq were co-embedded using scDART at the cell cluster level. Briefly, gene annotations extended by 2 kb upstream of the 5’ ends were used to construct the gene-to-region relationship matrix of the top 2000 most variable genes from each cell cluster. Gene-expression and chromatin-accessibility signals corresponding to these 2000 variable genes were extracted from each single cell, as the input for scDART. The co-embeded latent coordinates of single-cells from sc5nRNA-seq and scATAC-seq from each cluster was applied to SEACells to group neighboring cells into metacells. <br>
By targeting approximately 100 single-cells per metacell, a total of 395 metacells were generated across the 10 clusters. Two  metacells conatining fewer than 3 single-cells from either the sc5nRNA-seq or scATAC-seq were excluded, resulting in 393 metacells for downstram analyses. Raw gene-expression, CRE accessibility and CRE transcription counts from constituent single-cells were summed within each metacell to generate metacell level profiles.










