#!/bin/sh
java -mx12000M -jar ChromHMM.jar BinarizeBed -gzip CHROMSIZES/hg38.txt /analysisdata/fantom6/Interactome/CUTnTag/BED Neuron.cellmarkfile.v4.txt ./data_files/all_HG38