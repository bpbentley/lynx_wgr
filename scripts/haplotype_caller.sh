#!/bin/bash
#SBATCH -J lynx_WGR_haplotype_caller
#SBATCH -o ./logs/hap_call/remain/%x_%a_%A.log
#SBATCH -e ./logs/hap_call/remain/%x_%a_%A.err
#SBATCH -p cpu
#SBATCH -t 12:00:00
#SBATCH -c 20  # Number of Cores per Task
#SBATCH --mem=200G  # Requested Memory
#SBATCH --array=31

module load gatk/4.5.0.0

BAM=$(sed -n ${SLURM_ARRAY_TASK_ID}p ./remaining_samples.txt)
SAMPLE=$(basename "$BAM" | sed 's/_dedup\.bam$//')
OUTDIR=/scratch3/workspace/bbentley_smith_edu-lynx/outputs/gatk
REFERENCE="./ref/GCF_007474595.2_mLynCan4.pri.v2_genomic.fna"

# Make sure you create the sequence dictionary first
# gatk CreateSequenceDictionary -R reference.fasta

gatk HaplotypeCaller -R ${REFERENCE} -I $BAM -O $OUTDIR/${SAMPLE}.g.vcf -ERC GVCF