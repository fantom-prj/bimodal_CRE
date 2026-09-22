
##########################
#primary_folder=[primary_folder]
primary_folder="/analysisdata/fantom6/Interactome/single_cell_wallace/Figure/"
SE_folder=paste0(primary_folder,"Data_and_code/SE_identification/")

#===============================================================================
# ENCODE ATAC-seq pipeline to identify H3K27ac peaks, FASTQ to peaks
system("caper run atac.wdl -i iPS_H3K27ac.json")
system("caper run atac.wdl -i NSC_H3K27ac.json")
system("caper run atac.wdl -i NRN_H3K27ac.json")

# output -> pooled.pval0.01.300K.bfilt.narrowPeak.gz

#======
files=list.files(pattern="pooled.pval0.01.300K.bfilt.narrowPeak.gz")
files.names=gsub("pooled.pval0.01.300K.bfilt.narrowPeak.gz","",files)
for (i in 1: length(files)){
data=read.delim(files[i], header=F, stringsAsFactors = F)
data$V2=data$V2+1
data$label="CUTnTag"
data1=data[,c(1,11,2,3,6,6,6,6)]
write.table(data1[order(data1$V1,data1$V2),],paste0(SE_folder,"input/",files.names[i],"27ac_narrowpeak.gff"),col.names=F, row.names=F, sep="¥t", quote=F)
}
# output -> 27ac_narrowpeak.gff

#======
# merged the filtered bam files from two replicates 
system("samtools merge NRN-K27ac-both_rmdup_sorted.bam NRN-K27ac-r1_rmdup_sorted.bam NRN-K27ac-r2_rmdup_sorted.bam")
system("samtools merge NSC-K27ac-both_rmdup_sorted.bam NSC-K27ac-r1_rmdup_sorted.bam NSC-K27ac-r2_rmdup_sorted.bam")
system("samtools merge iPS-K27ac-both_rmdup_sorted.bam iPS-K27ac-r1_rmdup_sorted.bam iPS-K27ac-r2_rmdup_sorted.bam")

system("for file in *both_rmdup_sorted.bam; do samtools index \"$file\"; done")

#===============================================================================
# run ROSE
# bam files please refer to DDBJ DRA019568 (DRR614873-DRR614905)
setwd(SE_folder)
system("ROSE_main.py -g hg38 -i input/ips_27ac_narrowpeak.gff -r iPS-K27ac-both_rmdup_sorted.bam –c iPS-IgG-r1_rmdup_sorted.bam -o iPS_k27ac5 -s 10000 -t 2500")
system("ROSE_main.py -g hg38 -i input/nsc_27ac_narrowpeak.gff -r NSC-K27ac-both_rmdup_sorted.bam –c NSC-IgG-r1_rmdup_sorted.bam -o NSC_k27ac5 -s 10000 -t 2500")
system("ROSE_main.py -g hg38 -i input/nrn_27ac_narrowpeak.gff -r NRN-K27ac-both_rmdup_sorted.bam –c NRN-IgG-r1_rmdup_sorted.bam -o NRN_k27ac5 -s 10000 -t 2500")





