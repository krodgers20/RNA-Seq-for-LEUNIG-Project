#!/bin/bash

#SBATCH --job-name=stringtie_bams
#SBATCH --partition=cascade
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --nodes=1
#SBATCH --time=48:00:00
#SBATCH --mem=200G
#SBATCH --output=stringtie_%j.out
#SBATCH --error=stringtie_%j.err

#load modules
module load GCC/11.3.0
module load StringTie/2.2.1

# Path to the annotation file (GTF format)
annotation_file="/data/scratch/projects/punim1214/nicoleRNAseq/N2419466_30-1070643705_RNA_2024-09-23/refs/Arabidopsis_thaliana.TAIR10.59_RENAMED2.gff3"

# Path to the directory containing sorted BAM files
bam_dir="/data/scratch/projects/punim1214/nicoleRNAseq/N2419466_30-1070643705_RNA_2024-09-23/raw_data/bams"

# Output directory to store StringTie results
output_dir="/data/scratch/projects/punim1214/nicoleRNAseq/N2419466_30-1070643705_RNA_2024-09-23/raw_data/stringtie_out"

# Create output directory if it doesn't exist
mkdir -p "$output_dir"

# Loop through all sorted BAM files in the bam_dir
for bam_file in "$bam_dir"/*_sorted.bam; do
  # Get the basename of the BAM file (without the extension)
  base_name=$(basename "$bam_file" _sorted.bam)

  # Define output file paths
  output_gtf="$output_dir/${base_name}.gtf"

  # Run StringTie on the current BAM file
  stringtie "$bam_file" -G "$annotation_file" -e -B -o "$output_gtf"

  # Print progress message
  echo "Processed $bam_file"
done
