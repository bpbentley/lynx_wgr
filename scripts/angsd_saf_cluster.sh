#!/bin/bash
#SBATCH -J lynx_WGR_ANGSD_Fst
#SBATCH -o ./logs/%x_%a_%A.log
#SBATCH -e ./logs/%x_%a_%A.err
#SBATCH -p cpu
#SBATCH -t 48:00:00
#SBATCH -c 32  # Number of Cores per Task
#SBATCH --mem=240G  # Requested Memory
#SBATCH --array=1-4

module load conda/latest
conda activate bb_WGR

POP=$(sed -n ${SLURM_ARRAY_TASK_ID}p ./clusters.txt)

angsd -b bamlists/${POP}_bams.txt  -anc ./ref/GCF_007474595.2_mLynCan4.pri.v2_genomic.fna -out angsd/${POP}/${POP} -dosaf 1 -gl 1