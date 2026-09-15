#!/bin/bash
#SBATCH --job-name=ABC_3celltypes
#SBATCH --cpus-per-task=16
#SBATCH --mem=120G
#SBATCH --time=48:00:00
#SBATCH --output=logs_sc/ABC_%j.out
#SBATCH --error=logs_sc/ABC_%j.err

set -eo pipefail

mkdir -p logs_sc
cd /home/yip/tool2026/ABC-Enhancer-Gene-Prediction-main

source /home/yip/anaconda3/etc/profile.d/conda.sh
conda activate abc-env

snakemake \
  --cores ${SLURM_CPUS_PER_TASK} \
  --rerun-incomplete \
  --keep-going \
  /home/yip/tool2026/ABC-Enhancer-Gene-Prediction-main/result_sc/iPSC/Predictions/EnhancerPredictionsFull_threshold0.02_self_promoter.tsv \
  /home/yip/tool2026/ABC-Enhancer-Gene-Prediction-main/result_sc/NSC/Predictions/EnhancerPredictionsFull_threshold0.02_self_promoter.tsv \
  /home/yip/tool2026/ABC-Enhancer-Gene-Prediction-main/result_sc/Neuron/Predictions/EnhancerPredictionsFull_threshold0.02_self_promoter.tsv

