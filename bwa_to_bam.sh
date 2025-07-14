#!/bin/bash

#SBATCH --job-name=bwa_to_bam
#SBATCH --partition=cascade
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=200G
#SBATCH --nodes=1
#SBATCH --time=12:00:00
#SBATCH --output=bwa_to_bam_%j.out
#SBATCH --error=bwa_to_bam_%j.err


# Load samtools module (if required)
module load SAMtools/1.16.1

module load GCCcore/11.3.0

module load BWA/0.7.17


# Path to the reference genome
REFERENCE_GENOME="/data/scratch/projects/punim1214/nicoleRNAseq/N2419466_30-1070643705_RNA_2024-09-23/refs/ncbi_dataset/data/GCF_000001735.4/GCF_000001735.4_TAIR10.1_genomic.fna"

# Directory containing FASTQ files
FASTQ_DIR="/data/scratch/projects/punim1214/nicoleRNAseq/N2419466_30-1070643705_RNA_2024-09-23/raw_data/trimmed"

# Output directory for BAM files
OUTPUT_DIR="/data/scratch/projects/punim1214/nicoleRNAseq/N2419466_30-1070643705_RNA_2024-09-23/raw_data/bams"

# Ensure output directory exists
mkdir -p "$OUTPUT_DIR"

# Loop over all FASTQ files and run bwa mem for paired-end reads
for fastq1 in "$FASTQ_DIR"/*_1_paired.fq.gz; do
    # Get the corresponding fastq2 file (replace _1 with _2)
    fastq2="${fastq1/_1_paired.fq.gz/_2_paired.fq.gz}"

    # Extract the sample name (e.g., from "sample1_1.fastq" get "sample1")
    sample_name=$(basename "$fastq1" "_1_paired.fq.gz")

    # Define the output SAM and BAM files
    sam_output="$OUTPUT_DIR/${sample_name}.sam"
    bam_output="$OUTPUT_DIR/${sample_name}.bam"

    # Run bwa mem for paired-end alignment
    bwa mem "$REFERENCE_GENOME" "$fastq1" "$fastq2" > "$sam_output"

    # Convert SAM to BAM and sort it
    samtools view -S -b "$sam_output" | samtools sort -o "$bam_output"

    # Index the BAM file
    samtools index "$bam_output"

    # Remove the intermediate SAM file
    rm "$sam_output"

    echo "Finished processing $sample_name"
done

echo "All samples processed!"
