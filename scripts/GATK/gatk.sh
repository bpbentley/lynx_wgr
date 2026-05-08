#!/bin/bash
#SBATCH -J lynx_WGR_GATK
#SBATCH -o ./logs/%x_%a_%A.log
#SBATCH -e ./logs/%x_%a_%A.err
#SBATCH -p cpu
#SBATCH -t 48:00:00
#SBATCH -c 16  # Number of Cores per Task
#SBATCH --mem=240G  # Requested Memory

#module load conda/latest
#conda activate bb_WGR

# Population Genetics and GEA Analysis Pipeline

# Configuration
REFERENCE="./ref/GCF_007474595.2_mLynCan4.pri.v2_genomic.fna"
BAM_DIR="/scratch3/workspace/bbentley_smith_edu-lynx/aligned"
OUTPUT_DIR="./outputs"
THREADS=16
MIN_MAPQ=20
MIN_BASEQ=20

# Sample lists
ALL_SAMPLES="all_samples_bam.txt"           # All samples for ANGSD
HIGH_COV_SAMPLES="high_coverage_samples.txt" # High coverage samples for GATK
POPULATION_FILE="populations.txt"       # Sample -> Population mapping
ENVIRONMENT_FILE="environment.txt"      # Environmental data

echo "=== GATK Analysis (High Coverage Samples) ==="

# 1. Create sample-specific GVCFs (assuming this is done; see haplotype_caller.sh for array)
# This step would typically be done per sample:
# gatk HaplotypeCaller -R ${REFERENCE} -I sample.bam -O sample.g.vcf -ERC GVCF

# 2. Joint genotyping (assuming GVCFs exist; run in an array)
# Do separately
#gatk GenomicsDBImport \
#    --genomicsdb-workspace-path ${OUTPUT_DIR}/gatk/genomicsdb_${CHR} \
#    --sample-name-map gvcf_samples.map \
#    --reader-threads ${THREADS} \
#    -L ${CHR} \
#    --batch-size 50 \
#    --tmp-dir ./tmp

#gatk GenotypeGVCFs \
#    -R ${REFERENCE} \
#    -V gendb://${OUTPUT_DIR}/gatk/genomicsdb \
#    -O ${OUTPUT_DIR}/gatk/joint_calls.vcf.gz

# Combine vcfs from each chromosome
#bcftools concat -f /scratch3/workspace/bbentley_smith_edu-lynx/outputs/gatk/chromosome.vcfs -Oz -o ${OUTPUT_DIR}/gatk/lynx_WGR_merged.vcf.gz
#tabix -p vcf ${OUTPUT_DIR}/gatk/lynx_WGR_merged.vcf.gz

#conda deactivate

module load gatk/4.5.0.0 # Issues with GATK in the Conda environment, load here for now.

# Keep only SNPs
gatk SelectVariants \
    -R $REFERENCE \
    -V ${OUTPUT_DIR}/gatk/lynx_WGR_merged.vcf.gz \
    -select-type SNP \
    -O ${OUTPUT_DIR}/gatk/lynx_WGR_SNPs.vcf.gz


# 3. Hard filtering
gatk VariantFiltration \
    -R ${REFERENCE} \
    -V ${OUTPUT_DIR}/gatk/lynx_WGR_SNPs.vcf.gz \
    -O ${OUTPUT_DIR}/gatk/lynx_WGR_filtered.vcf.gz \
    --filter-expression "QD < 2.0" --filter-name "QD2" \
    --filter-expression "QUAL < 30.0" --filter-name "QUAL30" \
    --filter-expression "SOR > 3.0" --filter-name "SOR3" \
    --filter-expression "FS > 60.0" --filter-name "FS60" \
    --filter-expression "MQ < 40.0" --filter-name "MQ40" \
    --filter-expression "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" \
    --filter-expression "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" \
    --filter-expression "DP < 15" --filter-name "DP15" \
    --filter-expression "DP > 1500" --filter-name "DP1500"

# 4. Select only PASS variants
gatk SelectVariants \
    -V ${OUTPUT_DIR}/gatk/lynx_WGR_filtered.vcf.gz \
    -O ${OUTPUT_DIR}/gatk/lynx_WGR_pass_only.vcf.gz \
    --exclude-filtered
    
module purge
module load conda/latest
conda activate bb_WGR

bcftools +setGT ${OUTPUT_DIR}/gatk/lynx_WGR_pass_only.vcf.gz \
  --output-type z \
  --output ${OUTPUT_DIR}/gatk/lynx_WGR_DP4.vcf.gz \
  -- -t q -n . -i 'FORMAT/DP<4'

# 5. Additional filtering for population genetics
vcftools --gzvcf ${OUTPUT_DIR}/gatk/lynx_WGR_DP4.vcf.gz \
    --maf 0.05 \
    --max-missing 0.7 \
    --minQ 30 \
    --recode --recode-INFO-all \
    --out ${OUTPUT_DIR}/gatk/lynx_WGR_popgen

bgzip ${OUTPUT_DIR}/gatk/lynx_WGR_popgen.recode.vcf
tabix -p vcf ${OUTPUT_DIR}/gatk/lynx_WGR_popgen.recode.vcf.gz

# No unique locus identifiers, causes issues for PLINK. Add w/ bcftools
bcftools annotate --set-id '%CHROM\_%POS' -Oz -o lynx_WGR_popgen.recode.vcf.gz lynx_WGR_popgen.annot.vcf.gz
tabix -p vcf ${OUTPUT_DIR}/gatk/lynx_WGR_popgen.annot.vcf.gz

# Prune the VCF for LD - REMOVED TOO MANY VARIANTS, using PLINK instead
#bcftools +prune ${OUTPUT_DIR}/gatk/lynx_WGR_popgen.recode.vcf.gz \
#  -Oz --output ${OUTPUT_DIR}/gatk/lynx_WGR_popgen.LDpruned.vcf.gz \
#  -w 50000 -m 0.2

echo "=== Population Structure Analysis ==="

# Convert VCF to PLINK format for additional analyses
plink --vcf ${OUTPUT_DIR}/gatk/lynx_WGR_popgen.annot.vcf.gz \
    --make-bed \
    --out ${OUTPUT_DIR}/popgen/lynx_WGR_snps \
    --allow-extra-chr

# LD pruning for structure analysis 
plink --bfile ${OUTPUT_DIR}/popgen/lynx_WGR_snps \
    --indep-pairwise 50 10 0.2 \
    --out ${OUTPUT_DIR}/popgen/lynx_WGR_ld_pruned \
    --allow-extra-chr

plink --bfile ${OUTPUT_DIR}/popgen/lynx_WGR_snps \
    --extract ${OUTPUT_DIR}/popgen/lynx_WGR_ld_pruned.prune.in \
    --make-bed \
    --out ${OUTPUT_DIR}/popgen/lynx_WGR_pruned_snps \
    --allow-extra-chr
    
# Convert back to a VCF to use in R for popgen
plink --bfile ${OUTPUT_DIR}/popgen/lynx_WGR_pruned_snps \
      --recode vcf \
      --keep-allele-order \
      --out ${OUTPUT_DIR}/popgen/lynx_WGR_pruned_snps \
       --allow-extra-chr
    
conda deactivate