#!/bin/bash
#$ -cwd
#$ -q bigmem.q
#$ -pe def_slot 4

fq="/sequencedata/MGI/V300095456/Demultiplexed/L002/SI-GA-C4"
id="NSC1"
ref="/home/kouno/Visium/refdata-gex-GRCh38-2020-A"

echo -e $fq $id $ref
ulimit -n 20000

export  PATH=/work/cellranger-6.0.0/:$PATH
cellranger count \
                 --fastqs=${fq}                             \
                 --id=${id}                                 \
                 --localcores=4                            \
                 --transcriptome=${ref}

