#!/bin/bash
#SBATCH -J lynx_PCAngsd
#SBATCH -o ./logs/pcangsd/%x_%a_%A.log
#SBATCH -e ./logs/pcangsd/%x_%a_%A.err
#SBATCH -p cpu
#SBATCH -t 04:00:00
#SBATCH -c 16  # Number of Cores per Task
#SBATCH --mem=200G  # Requested Memory
#SBATCH --array=1-100

module load conda/latest
conda activate pcangsd

rep=${SLURM_ARRAY_TASK_ID}

mkdir -p pcangsd/output_logs

for k in {2..6}; do
mkdir -p pcangsd/K${k}/rep${rep}
pcangsd --beagle ./outputs/angsd/all_samples_geno.beagle.gz \
        --threads 16 \
        --out pcangsd/K${k}/rep${rep}/lynx_admix_K${k}_rep${rep} \
        --admix \
        --admix-K ${k} | tee pcangsd/output_logs/K${k}_rep${rep}.log
done

conda deactivate