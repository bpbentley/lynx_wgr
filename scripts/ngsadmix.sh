#!/bin/bash
#SBATCH -J lynx_WGR_NGSAdmix_Sep25
#SBATCH -o ./logs/NGSAdmix/Sep25/%x_%a_%A.log
#SBATCH -e ./logs/NGSAdmix/Sep25/%x_%a_%A.err
#SBATCH -p cpu
#SBATCH -t 08:00:00
#SBATCH -c 16  # Number of Cores per Task
#SBATCH --mem=200G  # Requested Memory
#SBATCH --array=1-50

module load conda/latest
conda activate bb_WGR

INDIR=/scratch3/workspace/bbentley_smith_edu-lynx/outputs/angsd
OUTDIR=/scratch3/workspace/bbentley_smith_edu-lynx/ngsadmix
K=$(sed -n ${SLURM_ARRAY_TASK_ID}p ./k_reps.txt | cut -f1)
rep=$(sed -n ${SLURM_ARRAY_TASK_ID}p ./k_reps.txt | cut -f2)

NGSadmix -likes $INDIR/all_samples_geno.beagle.gz -K $K -P 16 -o $OUTDIR/run_K${K}_rep${rep} -seed $RANDOM

conda deactivate
