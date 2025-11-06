#!/bin/bash
#SBATCH -J Ne
#SBATCH -o ./logs/Ne/%x_%a_%A.log
#SBATCH -e ./logs/Ne/%x_%a_%A.err
#SBATCH -p cpu
#SBATCH -t 24:00:00
#SBATCH -c 16  # Number of Cores per Task
#SBATCH --mem=180G  # Requested Memory

module load r/4.4.0

input=$1

echo "Running Ne estimation with input: $input"

Rscript /scratch3/workspace/bbentley_smith_edu-lynx/scripts/${input}_estimate_Ne.R
