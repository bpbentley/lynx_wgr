#!/bin/bash
#SBATCH -J lynx_WGR_ANGSD_Sep25
#SBATCH -o ./logs/%x_%a_%A.log
#SBATCH -e ./logs/%x_%a_%A.err
#SBATCH -p cpu
#SBATCH -t 48:00:00
#SBATCH -c 32  # Number of Cores per Task
#SBATCH --mem=240G  # Requested Memory

module load conda/latest
conda activate bb_WGR

# Population Genetics and GEA Analysis Pipeline

# Configuration
REFERENCE="./ref/GCF_007474595.2_mLynCan4.pri.v2_genomic.fna"
BAM_DIR="/scratch3/workspace/bbentley_smith_edu-lynx/aligned"
OUTPUT_DIR="./angsd"
THREADS=32
MIN_MAPQ=30
MIN_BASEQ=20

# Sample lists
ALL_SAMPLES="all_samples_bam.txt"           # All samples for ANGSD

#mkdir -p ${OUTPUT_DIR}/{angsd,gatk,popgen,gea}

echo "=== ANGSD Analysis (All Samples) ==="

# 1. Generate site allele frequency likelihoods
#angsd -bam ${ALL_SAMPLES} \
#    -ref ${REFERENCE} \
#    -out ${OUTPUT_DIR}/angsd/all_samples \
#    -nThreads ${THREADS} \
#    -GL 2 \
#    -doSaf 1 \
#    -anc ${REFERENCE} \
#    -minMapQ ${MIN_MAPQ} \
#    -minQ ${MIN_BASEQ} \
#    -doCounts 1 \
#    -setMinDepth 2 \
#    -setMaxDepth 1500

# 2. Optimize SFS and calculate theta
#realSFS ${OUTPUT_DIR}/angsd/all_samples.saf.idx \
#    -P ${THREADS} \
#    > ${OUTPUT_DIR}/angsd/all_samples.sfs

#realSFS saf2theta ${OUTPUT_DIR}/angsd/all_samples.saf.idx \
#    -sfs ${OUTPUT_DIR}/angsd/all_samples.sfs \
#    -outname ${OUTPUT_DIR}/angsd/all_samples

# 3. Calculate Tajima's D
#thetaStat do_stat ${OUTPUT_DIR}/angsd/all_samples.thetas.idx \
#    -win 50000 -step 10000 \
#    -outnames ${OUTPUT_DIR}/angsd/all_samples.thetasWindow

# 4. Generate genotype likelihoods for population structure
angsd -bam ${ALL_SAMPLES} \
    -ref ${REFERENCE} \
    -out ${OUTPUT_DIR}/all_samples_geno \
    -nThreads ${THREADS} \
    -GL 2 \
    -doGlf 2 \
    -doMajorMinor 1 \
    -SNP_pval 2e-6 \
    -doMaf 1 \
    -minMaf 0.05 \
    -minMapQ ${MIN_MAPQ} \
    -minQ ${MIN_BASEQ} \
    -minInd 30 \
    -P ${THREADS}

conda deactivate

# 5. PCA with PCAngsd

#conda activate pcangsd
#    pcangsd \
#    -b ${OUTPUT_DIR}/angsd/all_samples_geno.beagle.gz \
#    -o ${OUTPUT_DIR}/angsd/pca \
#    -t ${THREADS}
#conda deactivate