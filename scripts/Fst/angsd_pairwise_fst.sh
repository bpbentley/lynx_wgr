#!/bin/bash
#SBATCH -J lynx_WGR_ANGSD_PairwiseFst
#SBATCH -o ./logs/pw_fst/3rd/%x_%a_%A.log 
#SBATCH -e ./logs/pw_fst/3rd/%x_%a_%A.err
#SBATCH -p cpu
#SBATCH --nodes=1
#SBATCH -t 24:00:00 # Run for 3-4 days (add -q long)
#SBATCH -c 24 # Number of Cores per Task 
#SBATCH --mem=64G # >400GB for the first SFS calculation
#SBATCH --array=0-2   # 6 comparisons total

module load conda/latest
conda activate bb_WGR

# Define pairwise comparisons
pairs=(
  "NFLD NSLR"
  "NFLD ME"
  "NFLD SSLR"
  "NSLR ME"
  "NSLR SSLR"
  "ME SSLR"
)

# Select comparison for this array task
pair="${pairs[$SLURM_ARRAY_TASK_ID]}"
set -- $pair  # Split into two variables
pop1=$1
pop2=$2

echo "Running realSFS for $pop1 vs $pop2"

#realSFS angsd/${pop1}/${pop1}.saf.idx angsd/${pop2}/${pop2}.saf.idx -P 24 > angsd/${pop1}.${pop2}.ml

#realSFS fst index angsd/fst/${pop1}.saf.idx angsd/fst/${pop2}.saf.idx -sfs angsd/fst/${pop1}.${pop2}.ml -fstout angsd/fst/${pop1}_${pop2}
#realSFS fst stats angsd/fst/${pop1}_${pop2}.fst.idx

realSFS fst stats2 angsd/fst/${pop1}_${pop2}.fst.idx -win 50000 -step 10000 > angsd/fst/${pop1}_${pop2}_windows.txt

conda deactivate