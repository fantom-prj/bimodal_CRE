
##########################
#primary_folder=[primary_folder]
primary_folder="/analysisdata/fantom6/Interactome/single_cell_wallace/Figure/"
SE_folder=paste0(primary_folder,"Data_and_code/SE_identification/")

#===============================================================================
# bam files please refer to DDBJ DRA019568 (DRR614873-DRR614905)
setwd(SE_folder)
system("ROSE_main.py -g hg38 -i input/ips_27ac_narrowpeak.gff -r iPS-K27ac-both_rmdup_sorted.bam –c iPS-IgG-r1_rmdup_sorted.bam -o iPS_k27ac5 -s 10000 -t 2500")
system("ROSE_main.py -g hg38 -i input/nsc_27ac_narrowpeak.gff -r NSC-K27ac-both_rmdup_sorted.bam –c NSC-IgG-r1_rmdup_sorted.bam -o NSC_k27ac5 -s 10000 -t 2500")
system("ROSE_main.py -g hg38 -i input/nrn_27ac_narrowpeak.gff -r NRN-K27ac-both_rmdup_sorted.bam –c NRN-IgG-r1_rmdup_sorted.bam -o NRN_k27ac5 -s 10000 -t 2500")





