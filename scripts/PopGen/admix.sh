#!/bin/bash
#SBATCH -J lynx_WGR_admixture
#SBATCH -o ./logs/admix/%x_%a_%A.log
#SBATCH -e ./logs/admix/%x_%a_%A.err
#SBATCH -p cpu
#SBATCH -t 12:00:00
#SBATCH -c 20  # Number of Cores per Task
#SBATCH --mem=200G  # Requested Memory
#SBATCH --array=1

K=${SLURM_ARRAY_TASK_ID}
INDIR=/scratch3/workspace/bbentley_smith_edu-lynx/outputs/popgen
OUTDIR=${INDIR}/admix

admixture --cv $INDIR/lynx_WGR_neutral_snps.bed $K | tee $OUTDIR/Uro_Adaptive_log${K}.out