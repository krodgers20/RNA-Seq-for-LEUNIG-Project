#!/bin/bash

#SBATCH --job-name=sort_and_index_bams
#SBATCH --partition=cascade
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --nodes=1
#SBATCH --time=48:00:00
#SBATCH --mem=200G
#SBATCH --output=sortbam_%j.out
#SBATCH --error=sortbam_%j.err

#load samtools module
module load SAMtools/1.16.1

# Directory containing BAM files
BAM_DIR="/data/scratch/projects/punim1214/nicoleRNAseq/N2419466_30-1070643705_RNA_2024-09-23/raw_data/bams"

#loop through all BAM files in the directory
for bam_file in "$BAM_DIR"/*.bam; do
    # Get the base filename (without extension)
    base_name=$(basename "$bam_file" .bam)

    # Sort the BAM file and output to a new sorted BAM file
    sorted_bam_file="${BAM_DIR}/${base_name}_sorted.bam"
    echo "Sorting $bam_file..."
    samtools sort -o "$sorted_bam_file" "$bam_file"

    # Index the sorted BAM file
    echo "Indexing $sorted_bam_file..."
    samtools index "$sorted_bam_file"

    echo "Done processing $bam_file"
done

echo "All BAM files have been sorted and indexed."
