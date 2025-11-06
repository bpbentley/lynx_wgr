#!/bin/bash
#SBATCH -J lynx_WGR_ANGSD_Fst
#SBATCH -o ./logs/%x_%a_%A.log
#SBATCH -e ./logs/%x_%a_%A.err
#SBATCH -p cpu
#SBATCH -t 48:00:00
#SBATCH -c 24  # Number of Cores per Task
#SBATCH --mem=240G  # Requested Memory

module load conda/latest
conda activate bb_WGR

realSFS fst angsd/NFLD/NFLD.saf.idx angsd/NSLR/NSLR.saf.idx angsd/ME/ME.saf.idx angsd/SSLR/SSLR.saf.idx -sfs angsd/NFLD.NSLR.ml angsd/NFLD.ME.ml angsd/NFLD.SSLR.ml angsd/NSLR.ME.ml angsd/NSLR.SSLR.ml angsd/ME.SSLR.ml  -fstout angsd/all_pops_fst

realSFS fst angsd/stats2 angsd/all_pops_fst.fst.idx -win 50000 -step 10000 > angsd/all_pops_fst_slidingwindow

conda deactivate