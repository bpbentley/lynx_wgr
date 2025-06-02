#!/bin/bash
#SBATCH -J lynx_WGR_trim_align
#SBATCH -o ./logs/trim_align/%x_%a_%A.log
#SBATCH -e ./logs/trim_align/%x_%a_%A.err
#SBATCH -p cpu
#SBATCH -t 06:00:00
#SBATCH -c 20  # Number of Cores per Task
#SBATCH --mem=200G  # Requested Memory
#SBATCH --array=1-61

# Whole Genome Resequencing Pipeline
# For ~4X coverage diploid data
# Includes quality trimming and reference alignment
# With Claude AI

module load conda/latest
conda activate bat1k_cons

# Set variables
REFERENCE="ref/GCF_007474595.2_mLynCan4.pri.v2_genomic.fna"
SAMPLE=$(sed -n ${SLURM_ARRAY_TASK_ID}p ./sample_list.txt)
RAW_R1="data/raw/${SAMPLE}_R1.fastq.gz"
RAW_R2="data/raw/${SAMPLE}_R2.fastq.gz"
THREADS=20

# Create output directories
#mkdir -p trimmed aligned logs

#==============================================================================
# OPTION 1: Quality trimming with BBDuk (BBTools)
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
# OPTION 2: Quality trimming with Trimmomatic (alternative to BBDuk)
#==============================================================================

# Uncomment the following section if using Trimmomatic instead of BBDuk

# echo "Starting quality trimming with Trimmomatic..."
# 
# java -jar trimmomatic-0.39.jar PE \
#     -threads ${THREADS} \
#     -phred33 \
#     ${RAW_R1} ${RAW_R2} \
#     trimmed/${SAMPLE}_R1_trimmed.fastq.gz \
#     trimmed/${SAMPLE}_R1_unpaired.fastq.gz \
#     trimmed/${SAMPLE}_R2_trimmed.fastq.gz \
#     trimmed/${SAMPLE}_R2_unpaired.fastq.gz \
#     ILLUMINACLIP:TruSeq3-PE.fa:2:30:10:2:keepBothReads \
#     LEADING:3 \
#     TRAILING:3 \
#     SLIDINGWINDOW:4:15 \
#     MINLEN:50 \
#     2>&1 | tee logs/${SAMPLE}_trimmomatic.log

# Trimmomatic parameter annotations:
# PE                           : Paired-end mode
# -phred33                     : Quality scores are Phred+33 encoded
# ILLUMINACLIP                 : Remove Illumina adapters
#   TruSeq3-PE.fa:2:30:10:2   : adapter_file:seed_mismatches:palindrome_clip:simple_clip:min_adapter_length
#   keepBothReads             : Keep both reads even if only one passes filter
# LEADING:3                    : Remove low quality bases from beginning (Q<3)
# TRAILING:3                   : Remove low quality bases from end (Q<3)
# SLIDINGWINDOW:4:15          : Sliding window of 4bp, cut when average quality <15
# MINLEN:50                   : Drop reads shorter than 50bp

#==============================================================================
# Reference genome indexing (run once per reference)
#==============================================================================

echo "Checking BWA index..."

# Check if BWA index exists, create if not
if [ ! -f "${REFERENCE}.bwt" ]; then
    echo "Creating BWA index for reference genome..."
    bwa index ${REFERENCE}
    # Creates .amb, .ann, .bwt, .pac, .sa files
fi

# Also create samtools faidx index (needed for many downstream tools)
if [ ! -f "${REFERENCE}.fai" ]; then
    echo "Creating samtools faidx index..."
    samtools faidx ${REFERENCE}
fi

#==============================================================================
# Read alignment with BWA-MEM
#==============================================================================

echo "Starting alignment with BWA-MEM..."

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
java -Xmx8g -jar picard.jar MarkDuplicates \
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
# -Xmx8g                   : Allocate 8GB RAM (adjust based on available memory)

# Option 2: Using samtools markdup (alternative, faster but less detailed)
# samtools markdup -r -@ ${THREADS} aligned/${SAMPLE}_sorted.bam aligned/${SAMPLE}_dedup.bam
# samtools index aligned/${SAMPLE}_dedup.bam

#==============================================================================
# Base Quality Score Recalibration (BQSR) - Important for GATK
#==============================================================================

echo "Performing Base Quality Score Recalibration..."

# Step 1: Build recalibration model (requires known variant sites)
# You'll need a VCF of known variants for your species (e.g., dbSNP for humans)
# For non-model organisms, you may need to create this from high-confidence variants

if [ -f "known_variants.vcf.gz" ]; then
    echo "Running BQSR with known variants..."
    
    # Create recalibration table
    gatk BaseRecalibrator \
        --java-options "-Xmx8g" \
        -R ${REFERENCE} \
        -I aligned/${SAMPLE}_dedup.bam \
        --known-sites known_variants.vcf.gz \
        -O logs/${SAMPLE}_recal_data.table \
        2>&1 | tee logs/${SAMPLE}_bqsr1.log
    
    # Apply recalibration
    gatk ApplyBQSR \
        --java-options "-Xmx8g" \
        -R ${REFERENCE} \
        -I aligned/${SAMPLE}_dedup.bam \
        --bqsr-recal-file logs/${SAMPLE}_recal_data.table \
        -O aligned/${SAMPLE}_recalibrated.bam \
        2>&1 | tee logs/${SAMPLE}_bqsr2.log
    
    # Optional: Generate post-recalibration table for comparison
    gatk BaseRecalibrator \
        --java-options "-Xmx8g" \
        -R ${REFERENCE} \
        -I aligned/${SAMPLE}_recalibrated.bam \
        --known-sites known_variants.vcf.gz \
        -O logs/${SAMPLE}_post_recal_data.table
    
    # Generate recalibration plots
    gatk AnalyzeCovariates \
        --java-options "-Xmx8g" \
        -before logs/${SAMPLE}_recal_data.table \
        -after logs/${SAMPLE}_post_recal_data.table \
        -plots logs/${SAMPLE}_recalibration_plots.pdf
    
    FINAL_BAM="aligned/${SAMPLE}_recalibrated.bam"
    echo "BQSR completed. Using recalibrated BAM for downstream analysis."
    
else
    echo "No known variants file found. Skipping BQSR."
    echo "For non-model organisms, consider running initial variant calling"
    echo "to create a high-confidence variant set for BQSR."
    FINAL_BAM="aligned/${SAMPLE}_dedup.bam"
fi

#==============================================================================
# Index the final BAM file
#==============================================================================

echo "Indexing final BAM file..."

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

#==============================================================================
# Recommended next steps for ANGSD and GATK
#==============================================================================

echo ""
echo "=== RECOMMENDED NEXT STEPS ==="
echo ""
echo "For ANGSD analysis (population genetics):"
echo "  1. Create a comprehensive bamlist with all samples"
echo "  2. Filter by mapping quality: -minMapQ 20"
echo "  3. Filter by base quality: -minQ 20"
echo "  4. Set minimum depth: -minInd [number] (e.g., 50% of samples)"
echo "  5. For 4X data, consider -setMinDepth 1 -setMaxDepth [2*coverage]"
echo ""
echo "For GATK variant calling:"
echo "  1. Use HaplotypeCaller in GVCF mode for each sample"
echo "  2. Combine GVCFs with GenomicsDBImport or CombineGVCFs"
echo "  3. Joint genotype with GenotypeGVCFs"
echo "  4. Apply hard filters or VQSR for variant quality control"
echo "  5. For low coverage: consider --standard-min-confidence-threshold-for-calling 10"

#==============================================================================
# Notes for downstream analysis with 4X coverage diploid data:
#==============================================================================

# 1. ANGSD-specific considerations for low coverage:
#    - Use genotype likelihoods instead of called genotypes
#    - Consider -GL 2 (GATK model) for compatibility
#    - Use -doMaf 1 -doMajorMinor 1 for allele frequency estimation
#    - For PCA: -doIBS 2 -doCov 1 -makeMatrix 1
#    - For Fst: -doSaf 1 for site allele frequency spectra

# 2. GATK considerations:
#    - Use --emit-ref-confidence GVCF for better joint calling
#    - Consider --min-base-quality-score 20
#    - For low coverage: --standard-min-confidence-threshold-for-calling 10
#    - Apply appropriate filters: QD < 2.0, FS > 60.0, MQ < 40.0
#    - Consider using --max-alternate-alleles 3 for efficiency

# 3. General low-coverage recommendations:
#    - Imputation may be beneficial (BEAGLE, IMPUTE2)
#    - Consider population-based calling strategies
#    - Use strict quality filters but not overly restrictive depth filters
#    - Validate results with higher coverage subset if available
