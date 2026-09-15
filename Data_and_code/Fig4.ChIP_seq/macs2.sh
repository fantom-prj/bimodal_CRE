#!/bin/bash

cd /analysisdata/fantom6/Interactome/single_cell_wallace/Figure/Data_and_code/ChIP_seq/final_data/TEAD4
#TEAD
macs2 callpeak \
-t TEAD4-rep1.markdup.sort.nodup.bam \
-c IgG2a-rep1.markdup.sort.nodup.bam \
-f BAMPE \
-g hs \
-n ./macs2/TEAD4-rep1 \
-B \
-q 0.01

macs2 callpeak \
-t TEAD4-rep2.markdup.sort.nodup.bam \
-c IgG2a-rep1.markdup.sort.nodup.bam \
-f BAMPE \
-g hs \
-n ./macs2/TEAD4-rep2 \
-B \
-q 0.01

#====
cd /analysisdata/fantom6/Interactome/single_cell_wallace/Figure/Data_and_code/ChIP_seq/final_data/ONECUT2
#ONECUT
macs2 callpeak \
-t ONECUT-rep1.markdup.sort.nodup.bam \
-c IgG-rep1.markdup.sort.nodup.bam \
-f BAMPE \
-g hs \
-n ./macs2/ONECUT-rep1 \
-B \
-q 0.01

macs2 callpeak \
-t ONECUT-rep2.markdup.sort.nodup.bam \
-c IgG-rep1.markdup.sort.nodup.bam \
-f BAMPE \
-g hs \
-n ./macs2/ONECUT-rep2 \
-B \
-q 0.01
