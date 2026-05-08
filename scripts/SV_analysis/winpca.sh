#!/bin/bash
#SBATCH -J lynx_winpca
#SBATCH -o ./logs/%x_%a_%A.log
#SBATCH -e ./logs/%x_%a_%A.err
#SBATCH -p cpu
#SBATCH -t 12:00:00
#SBATCH -c 8  # Number of Cores per Task
#SBATCH --mem=100G  # Requested Memory
#SBATCH --array=1-54

# === Load conda and the associated environment ===
module load conda/latest
conda activate bb_WGR

# === Specify the chromosome for the array and extract length ===
CHR=SUPER_${SLURM_ARRAY_TASK_ID}
LEN=$(sed -n ${SLURM_ARRAY_TASK_ID}p ./ref/GCF_007474595.2_mLynCan4.pri.v2_scaffold_lengths.txt | cut -f2)

# === Pull the variants from the chromosome from the full VCF ===
vcftools --gzvcf ./outputs/lynx_WGR_pruned_snps.vcf --chr ${CHR} --recode --out ./winpca/chr/${CHR}/${CHR}

# === Just take the biallelic SNPs and index the VCF ===
bcftools view --types snps -m 2 -M 2 ./winpca/chr/${CHR}/${CHR}.recode.vcf -Oz -o ./winpca/chr/${CHR}/${CHR}_biallelic.vcf.gz
bcftools index ./winpca/chr/${CHR}/${CHR}_biallelic.vcf.gz

# === Run the winpca PCA script and plot PC1 and heterozygosity for all samples ===
~/winpca/winpca pca --np ./outputs/winpca/out/${CHR} ./winpca/chr/${CHR}/${CHR}_biallelic.vcf.gz ${CHR}:1-${LEN}
~/winpca/winpca chromplot -i 1 ./outputs/winpca/out/${CHR} ${CHR}:1-${LEN} -m ./popmap -g Population
~/winpca/winpca chromplot -i 1 -p het ./outputs/${CHR} ${CHR}:1-${LEN} -m ./final_meta.txt -g Population

# === Deactivate conda environment ===
conda deactivate
echo "winpca analysis complete for" ${CHR}