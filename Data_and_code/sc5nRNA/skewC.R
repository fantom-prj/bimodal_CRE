
###############

#primary_folder=[primary_folder]
primary_folder="/analysisdata/fantom6/Interactome/single_cell_wallace/Figure/"
RNA=paste0(primary_folder,"Data_and_code/sc5nRNA/")

#===============================================================================
# identifiy skew cells from BAM files by skewC 
# 1. split BAM by 10x barcode
setwd(paste0(RNA,"output/sc5_SkewC/"))
cell_types = c("iPS1", "iPS2", "NSC1", "NSC2", "Neuron1", "Neuron2")
for (i in cell_types) {
  # 1. Split 10X by Barcode
  cmd1 = paste(
    "bash 0_split10XbyBarcode.sh",
    paste0("/analysisdata/fantom6/Interactome/sc5cellranger_Kouno/hg38/", i, "/outs/possorted_genome_bam.bam"),
    paste0("/analysisdata/fantom6/Interactome/sc5cellranger_Kouno/hg38/", i, "/outs/filtered_feature_bc_matrix/barcodes.tsv.gz"),
    paste0(i, "_raw"))
  system(cmd1)
  
  # 2. Calculate gene body coverage
  system(paste0("bash 1_geneBodyCoverage.sh hg38 ", i, "_intermediate ", i, "_raw"))
  
  # 3. Classify cells
  system(paste0("bash 2_SkewC.sh F6_Interactome_", i, " ", i, "_intermediate ", i))}

# _raw and _intermediate folders were not included

#===============================================================================
# parse the annotation
setwd(paste0(RNA,"output/sc5_SkewC/"))
files=list.files(pattern="SkewCAnnotation.rds", recursive=T)
files.names=sapply(strsplit(files,"\\/"),"[",1)
anno=data.frame()
for (i in 1: length(files)){
anno1=readRDS(files[i])
anno1$barcode=rownames(anno1)
anno1$orig.ident=files.names[i]
anno=rbind(anno, anno1)
}

write.table(anno[,c(2,3,1)],"sc5RNA_SkewC_QC_Annotation.tsv.gz", col.names=T, row.names=F, sep="\t",quote=F)

#===============================================================================


