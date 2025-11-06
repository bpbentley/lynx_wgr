#!/bin/bash
#SBATCH -J lynx_WGR_FEEMS
#SBATCH -o ./logs/%x_%a_%A.log
#SBATCH -e ./logs/%x_%a_%A.err
#SBATCH -p cpu
#SBATCH -t 48:00:00
#SBATCH -c 32  # Number of Cores per Task
#SBATCH --mem=240G  # Requested Memory

module load conda/latest
conda activate feems_e

python3 ./scripts/feems.sh

conda deactivate