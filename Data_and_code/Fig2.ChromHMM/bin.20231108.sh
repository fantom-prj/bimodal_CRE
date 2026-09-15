#!/bin/sh
java -mx6000M -jar ChromHMM.jar BinarizeBam -gzip -mixed CHROMSIZES/hg38.txt /analysisdata/fantom6/Interactome/CUTnTag/BAM Neuron.cellmarkfile.txt ./data_files/all_HG38
