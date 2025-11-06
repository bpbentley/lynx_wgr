#!/bin/bash
#SBATCH -J lynx_WGR_trim_align
#SBATCH -o ./logs/trim_align/%x_%a_%A.log
#SBATCH -e ./logs/trim_align/%x_%a_%A.err
#SBATCH -p cpu
#SBATCH -t 05:00:00
#SBATCH -c 20  # Number of Cores per Task
#SBATCH --mem=200G  # Requested Memory
#SBATCH --array=62

# Whole Genome Resequencing Pipeline
# For ~4X coverage diploid data
# Includes quality trimming and reference alignment
# With Claude AI

module load conda/latest

conda activate bb_WGR
module load openjdk/21

# Set variables
REFERENCE="ref/GCF_007474595.2_mLynCan4.pri.v2_genomic.fna"
SAMPLE=$(sed -n ${SLURM_ARRAY_TASK_ID}p ./sample_list.txt)
RAW_R1="data/raw/${SAMPLE}_R1.fastq.gz"
RAW_R2="data/raw/${SAMPLE}_R2.fastq.gz"
THREADS=20

# Create output directories
#mkdir -p trimmed aligned logs

#==============================================================================
# STEP 1: Quality trimming with BBDuk (BBTools)
#==============================================================================

echo "Starting quality trimming with BBDuk..."

~/bbmap/bbduk.sh \
    in1=${RAW_R1} \
    in2=${RAW_R2} \
    out1=trimmed/${SAMPLE}_R1_trimmed.fastq.gz \
    out2=trimmed/${SAMPLE}_R2_trimmed.fastq.gz \
    ref=adapters \
    ktrim=r \
    k=23 \
    mink=11 \
    hdist=1 \
    tpe \
    tbo \
    qtrim=rl \
    trimq=20 \
    minlen=50 \
    threads=${THREADS} \
    stats=logs/${SAMPLE}_bbduk_stats.txt \
    2>&1 | tee logs/${SAMPLE}_bbduk.log

# BBDuk parameter annotations:
# ref=adapters         : Use built-in adapter sequences for trimming
# ktrim=r             : Trim adapters from right end of reads
# k=23                : Kmer length for adapter detection (23 is good default)
# mink=11             : Minimum kmer length at read ends (allows shorter matches)
# hdist=1             : Maximum Hamming distance for kmer matches (1 mismatch allowed)
# tpe                 : Trim paired reads to same length (maintains pairing)
# tbo                 : Trim adapters based on pair overlap detection
# qtrim=rl            : Quality trim from right and left ends
# trimq=20            : Trim bases with quality below Q20
# minlen=50           : Discard reads shorter than 50bp after trimming
# threads=8           : Use 8 CPU threads for processing
# stats=file          : Output detailed trimming statistics

#==============================================================================
# Reference genome indexing (run once per reference)
#==============================================================================

#echo "Checking BWA index..."

# Check if BWA index exists, create if not
#if [ ! -f "${REFERENCE}.bwt" ]; then
#    echo "Creating BWA index for reference genome..."
#    bwa index ${REFERENCE}
#    # Creates .amb, .ann, .bwt, .pac, .sa files
#fi

# Also create samtools faidx index (needed for many downstream tools)
#if [ ! -f "${REFERENCE}.fai" ]; then
#    echo "Creating samtools faidx index..."
#    samtools faidx ${REFERENCE}
#fi

#==============================================================================
# Read alignment with BWA-MEM
#==============================================================================

#echo "Starting alignment with BWA-MEM..."

bwa mem \
    -t ${THREADS} \
    -M \
    -R "@RG\tID:${SAMPLE}\tSM:${SAMPLE}\tPL:ILLUMINA\tLB:${SAMPLE}_lib" \
    ${REFERENCE} \
    trimmed/${SAMPLE}_R1_trimmed.fastq.gz \
    trimmed/${SAMPLE}_R2_trimmed.fastq.gz \
    2>logs/${SAMPLE}_bwa.log | \
samtools view -@ ${THREADS} -bS - | \
samtools sort -@ ${THREADS} -o aligned/${SAMPLE}_sorted.bam -

# BWA-MEM parameter annotations:
# -t 8                : Use 8 threads for alignment
# -M                  : Mark shorter split hits as secondary (for Picard compatibility)
# -R "@RG\t..."      : Add read group information to BAM header
#   ID:sample_name    : Read group identifier
#   SM:sample_name    : Sample name
#   PL:ILLUMINA      : Platform (sequencing technology)
#   LB:library_name  : Library identifier

# Samtools piping annotations:
# samtools view -bS  : Convert SAM to BAM format (-b=binary, -S=input is SAM)
# samtools sort      : Sort BAM by coordinate for efficient access
# -@ 8               : Use 8 threads for sorting
# -o output.bam      : Specify output filename

#==============================================================================
# Remove PCR duplicates (recommended for GATK and ANGSD)
#==============================================================================

echo "Marking and removing PCR duplicates..."

# Option 1: Using Picard MarkDuplicates (recommended)
picard MarkDuplicates \
    INPUT=aligned/${SAMPLE}_sorted.bam \
    OUTPUT=aligned/${SAMPLE}_dedup.bam \
    METRICS_FILE=logs/${SAMPLE}_duplicate_metrics.txt \
    REMOVE_DUPLICATES=true \
    CREATE_INDEX=true \
    VALIDATION_STRINGENCY=SILENT \
    2>&1 | tee logs/${SAMPLE}_dedup.log

# Picard MarkDuplicates parameters:
# REMOVE_DUPLICATES=true    : Actually remove duplicates (not just mark them)
# CREATE_INDEX=true         : Create BAM index automatically
# METRICS_FILE             : Output duplicate statistics
# VALIDATION_STRINGENCY    : SILENT avoids warnings for minor format issues

#==============================================================================
# Index the final BAM file
#==============================================================================

echo "Indexing final BAM file..."

FINAL_BAM="aligned/${SAMPLE}_dedup.bam"
samtools index ${FINAL_BAM}
# Creates .bai index file for random access to BAM

#==============================================================================
# Generate alignment statistics
#==============================================================================

echo "Generating alignment statistics..."

# Basic alignment stats
samtools stats ${FINAL_BAM} > logs/${SAMPLE}_alignment_stats.txt

# Flagstat summary
samtools flagstat ${FINAL_BAM} > logs/${SAMPLE}_flagstat.txt

# Coverage depth (useful for 4X coverage data)
samtools depth ${FINAL_BAM} | \
    awk '{sum+=$3; count++} END {print "Average depth:", sum/count}' \
    > logs/${SAMPLE}_coverage.txt

#==============================================================================
# Quality control checks specific for low-coverage data
#==============================================================================

echo "Performing QC checks for low-coverage data..."

# Check for proper pairing (important for downstream variant calling)
samtools view -f 2 -c ${FINAL_BAM} > logs/${SAMPLE}_properly_paired_count.txt

# Calculate percentage of genome covered at different depths
# Useful for 4X data to see coverage distribution
samtools depth ${FINAL_BAM} | \
    awk '{
        depth[$3]++; 
        total++
    } 
    END {
        for(d=1; d<=10; d++) {
            covered = 0;
            for(i=d; i<=50; i++) covered += depth[i];
            print "Positions with depth >=" d ":", (covered/total)*100 "%"
        }
    }' > logs/${SAMPLE}_depth_coverage.txt

#==============================================================================
# Additional QC and preparation for ANGSD/GATK
#==============================================================================

echo "Performing additional QC for ANGSD/GATK preparation..."

# 1. Calculate insert size distribution (important for paired-end data)
samtools stats ${FINAL_BAM} | grep "insert size" > logs/${SAMPLE}_insert_size.txt

# 2. Generate more detailed coverage statistics
# Bedtools genomecov - useful for ANGSD depth filtering
if command -v bedtools &> /dev/null; then
    bedtools genomecov -ibam ${FINAL_BAM} -d > logs/${SAMPLE}_per_base_coverage.txt
    bedtools genomecov -ibam ${FINAL_BAM} > logs/${SAMPLE}_coverage_histogram.txt
fi

# 3. Check mapping quality distribution
samtools view ${FINAL_BAM} | \
    awk '{if($5 >= 0) mapq[$5]++} END {for(q in mapq) print q, mapq[q]}' | \
    sort -n > logs/${SAMPLE}_mapq_distribution.txt

# 4. Count reads in different categories for ANGSD
echo "Counting reads by category for ANGSD input assessment..."
echo "Total reads:" > logs/${SAMPLE}_read_categories.txt
samtools view -c ${FINAL_BAM} >> logs/${SAMPLE}_read_categories.txt
echo "Properly paired reads:" >> logs/${SAMPLE}_read_categories.txt
samtools view -c -f 2 ${FINAL_BAM} >> logs/${SAMPLE}_read_categories.txt
echo "High quality mapped reads (MAPQ>=20):" >> logs/${SAMPLE}_read_categories.txt
samtools view -c -q 20 ${FINAL_BAM} >> logs/${SAMPLE}_read_categories.txt
echo "Unique mapped reads (MAPQ>=1, not secondary):" >> logs/${SAMPLE}_read_categories.txt
samtools view -c -q 1 -F 256 ${FINAL_BAM} >> logs/${SAMPLE}_read_categories.txt

#==============================================================================
# Create ANGSD-ready file list
#==============================================================================

echo "Creating ANGSD input file list..."

# ANGSD typically works with lists of BAM files
echo "${PWD}/${FINAL_BAM}" > ${SAMPLE}_bamlist.txt

echo "ANGSD bamlist created: ${SAMPLE}_bamlist.txt"

#==============================================================================
# Validate final BAM for GATK compatibility
#==============================================================================

echo "Validating BAM file for GATK compatibility..."

# Check if GATK ValidateSamFile is available
if command -v gatk &> /dev/null; then
    gatk ValidateSamFile \
        --java-options "-Xmx4g" \
        -I ${FINAL_BAM} \
        -MODE SUMMARY \
        -O logs/${SAMPLE}_validation_report.txt \
        2>&1 | tee logs/${SAMPLE}_validation.log
    
    echo "GATK validation completed. Check logs/${SAMPLE}_validation_report.txt for any issues."
fi

echo "Pipeline completed successfully!"
echo "Output files:"
echo "  - Trimmed reads: trimmed/${SAMPLE}_R*_trimmed.fastq.gz"
echo "  - Final processed BAM: ${FINAL_BAM}"
echo "  - ANGSD input list: ${SAMPLE}_bamlist.txt"
echo "  - Logs and stats: logs/${SAMPLE}_*"
echo ""
echo "Ready for downstream analysis with:"
echo "  - ANGSD: Use ${SAMPLE}_bamlist.txt as input"
echo "  - GATK: Use ${FINAL_BAM} for variant calling"