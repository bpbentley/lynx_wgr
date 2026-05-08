#!/bin/bash
#SBATCH -J lynx_WGR_GenomicsDBimport
#SBATCH -o ./logs/GenomicsDBimport/geno_all/%x_%a_%A.log
#SBATCH -e ./logs/GenomicsDBimport/geno_all/%x_%a_%A.err
#SBATCH -p cpu
#SBATCH -t 36:00:00
#SBATCH -c 16  # Number of Cores per Task
#SBATCH --mem=200G  # Requested Memory
#SBATCH --array=20-54

module load gatk/4.5.0.0

CHR=$(sed -n ${SLURM_ARRAY_TASK_ID}p ./ref/chromosome.list)
THREADS=16
OUTPUT_DIR="./outputs"
REFERENCE="./ref/GCF_007474595.2_mLynCan4.pri.v2_genomic.fna"

gatk GenomicsDBImport \
    --genomicsdb-workspace-path ${OUTPUT_DIR}/gatk/genomicsdb_${CHR} \
    --sample-name-map gvcf_samples2.map \
    --reader-threads ${THREADS} \
    -L ${CHR} \
    --batch-size 50 \
    --tmp-dir ./tmp
    
gatk GenotypeGVCFs \
    -R ${REFERENCE} \
    -V gendb://${OUTPUT_DIR}/gatk/genomicsdb_${CHR} \
    -O ${OUTPUT_DIR}/gatk/chr/joint_calls_${CHR}.vcf.gz
    
