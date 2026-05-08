#!/bin/bash
#SBATCH -J lynx_WGR_ANGSD_Het
#SBATCH -o ./logs/ANGSD_Het/%x_%a_%A.log
#SBATCH -e ./logs/ANGSD_Het/%x_%a_%A.err
#SBATCH -p cpu
#SBATCH -t 24:00:00
#SBATCH -c 16  # Number of Cores per Task
#SBATCH --mem=120G  # Requested Memory
#SBATCH --array=2-57

module load conda/latest
conda activate bb_WGR

BAM=$(sed -n ${SLURM_ARRAY_TASK_ID}p ./all_samples_bam.txt)
IND=$(basename "$BAM" | sed 's/_dedup\.bam$//')

angsd -i ${BAM} -anc ./ref/GCF_007474595.2_mLynCan4.pri.v2_genomic.fna -C 50 -ref ./ref/GCF_007474595.2_mLynCan4.pri.v2_genomic.fna -minQ 20 -minmapq 30 -dosaf 1 -GL 2 -out ./angsd/het/${IND}

realSFS ./angsd/het/${IND}.saf.idx -anc ./ref/GCF_007474595.2_mLynCan4.pri.v2_genomic.fna -ref ./ref/GCF_007474595.2_mLynCan4.pri.v2_genomic.fna -fold 1 > ./angsd/het/${IND}.est.ml

