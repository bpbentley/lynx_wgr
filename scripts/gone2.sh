#!/bin/bash
#SBATCH -J lynx_WGR_GONE2
#SBATCH -o ./logs/GONE2/%x_%a_%A.log 
#SBATCH -e ./logs/GONE2/%x_%a_%A.err
#SBATCH -p cpu
#SBATCH --nodes=1
#SBATCH -t 48:00:00 
#SBATCH -c 24 # Number of Cores per Task 
#SBATCH --mem=320G # Requested Memory
#SBATCH --array=1

POP=$(printf cluster${SLURM_ARRAY_TASK_ID})
INDIR=/scratch3/workspace/bbentley_smith_edu-lynx/gone/GONE_input

/home/bbentley_smith_edu/GONE2/gone2 -g 3 -r 1.9 $INDIR/${POP}_GONE.ped