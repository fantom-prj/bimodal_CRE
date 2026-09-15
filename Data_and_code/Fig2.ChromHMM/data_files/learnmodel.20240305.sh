#!/bin/sh
cd /analysisdata/fantom6/Interactome/single_cell_wallace/Figure/Data_and_code/ChromHMM/data_files

java -mx6000M -jar ChromHMM.jar LearnModel -p 8 all_HG38 chromatin_state/k7 7 hg38
java -mx6000M -jar ChromHMM.jar LearnModel -p 8 all_HG38 chromatin_state/k8 8 hg38
java -mx6000M -jar ChromHMM.jar LearnModel -p 8 all_HG38 chromatin_state/k9 9 hg38
java -mx6000M -jar ChromHMM.jar LearnModel -p 8 all_HG38 chromatin_state/k10 10 hg38
java -mx6000M -jar ChromHMM.jar LearnModel -p 8 all_HG38 chromatin_state/k11 11 hg38
java -mx6000M -jar ChromHMM.jar LearnModel -p 8 all_HG38 chromatin_state/k12 12 hg38
java -mx6000M -jar ChromHMM.jar LearnModel -p 8 all_HG38 chromatin_state/k13 13 hg38
java -mx6000M -jar ChromHMM.jar LearnModel -p 8 all_HG38 chromatin_state/k14 14 hg38
java -mx6000M -jar ChromHMM.jar LearnModel -p 8 all_HG38 chromatin_state/k15 15 hg38
java -mx6000M -jar ChromHMM.jar LearnModel -p 8 all_HG38 chromatin_state/k16 16 hg38
