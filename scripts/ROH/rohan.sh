#!/bin/bash
#SBATCH -J lynx_WGR_ROHan_1mb_1e8
#SBATCH -o ./logs/ROHan/100kb/%x_%a_%A.log
#SBATCH -e ./logs/ROHan/100kb/%x_%a_%A.err
#SBATCH -p cpu
#SBATCH -t 24:00:00
#SBATCH -c 24  # Number of Cores per Task
#SBATCH --mem=120G  # Requested Memory
#SBATCH --array=1-57

module load gsl/2.7.1
module load samtools/1.19.2

# Ts/Tv ratio - 2.424 from vcftools
#vcftools --gzvcf outputs/gatk/lynx_WGR_popgen.annot.vcf.gz --TsTv 1000
# Using 1mb windows and a rhomu of 2e-4 per Taylor et al. 2024

BAM=$(sed -n ${SLURM_ARRAY_TASK_ID}p ./all_samples_bam.txt)
ID=$(basename "$BAM" | sed 's/_dedup\.bam$//')

/home/bbentley_smith_edu/rohan/bin/rohan --tstv 2.4 --rohmu 2e-4 --size 1000000 --auto ./ref/autosome.list.txt -o ./rohan/${ID}_1mb_autosome_2e4 -t 24 /scratch3/workspace/bbentley_smith_edu-lynx/ref/GCF_007474595.2_mLynCan4.pri.v2_genomic.fna ${BAM}