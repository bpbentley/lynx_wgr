#!/bin/bash
#SBATCH -J lynx_WGR_bayescan
#SBATCH -o ./logs/bayescan/%x_%a_%A.log
#SBATCH -e ./logs/bayescan/%x_%a_%A.err
#SBATCH -p cpu
#SBATCH -t 08:00:00
#SBATCH -c 16  # Number of Cores per Task
#SBATCH --mem=250G  # Requested Memory

module load conda/latest
conda activate bb_WGR

OUTDIR=./outputs/popgen
OUTPREFIX="lynx_WGR"
THREADS=16

# === Convert VCF to BayeScan input format ===
#perl scripts/vcf2bayescan.pl -p popmap.txt -v outputs/popgen/lynx_WGR_pruned_snps.vcf
# Note that I ran this interactively first then moved the output file. Took ~2 mins

# === Run Bayescan ===
echo "Running BayeScan..."
bayescan2 $OUTDIR/${OUTPREFIX}_bayescan.txt \
  -od $OUTDIR \
  -o $OUTPREFIX \
  -threads $THREADS \
  -snp \
  -n 5000 \
  -thin 10 \
  -nbp 20 \
  -pilot 5000
 
conda deactivate

# === Output message ===
echo "BayeScan run complete. Results in $WORKDIR/${OUTPREFIX}_fst.txt and _sel.txt"