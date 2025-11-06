#!/bin/bash
#SBATCH -J lynx_cat_align
#SBATCH -o ./logs/%x_%a_%A.log
#SBATCH -e ./logs/%x_%a_%A.err
#SBATCH -p cpu
#SBATCH -t 12:00:00
#SBATCH -c 24  # Number of Cores per Task
#SBATCH --mem=200G  # Requested Memory

module load conda/latest
conda activate mummer

minimap2 -x asm5 ./ref/GCF_018350175.1_F.catus_clean.fna ./ref/GCF_007474595.2_mLynCan4.pri.v2_genomic.fna > ./ref/cat_lynx.paf

conda deactivate