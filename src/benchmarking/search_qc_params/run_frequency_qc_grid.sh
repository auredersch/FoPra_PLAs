#!/bin/bash
#SBATCH --job-name=frequency_qc_grid
#SBATCH --output=/nfs/home/students/f.mathis/FoPra_PLAs/logs/frequency_qc_grid_%j.out
#SBATCH --error=/nfs/home/students/f.mathis/FoPra_PLAs/logs/frequency_qc_grid_%j.err
#SBATCH --chdir=/nfs/home/students/f.mathis/FoPra_PLAs
#SBATCH --time=02:00:00
#SBATCH --mem=128G
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1

set -euo pipefail
source ~/.bashrc
conda activate FoPra

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1

Rscript src/benchmarking/search_qc_params/Benchmarking_Frequency_QC_Grid.R
