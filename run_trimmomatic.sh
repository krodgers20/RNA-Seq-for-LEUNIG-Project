#!/bin/bash

#SBATCH --job-name=trimmomatic_processing
#SBATCH --partition=cascade
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=200G
#SBATCH --nodes=1
#SBATCH --time=12:00:00
#SBATCH --output=trimmomatic_%j.out
#SBATCH --error=trimmomatic_%j.err

#Load the Trimmomatic module
module load Trimmomatic/0.39-Java-11
java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.39.jar



#Define directories
FASTQ_DIR="/data/scratch/projects/punim1214/nicoleRNAseq/N2419466_30-1070643705_RNA_2024-09-23/raw_data"
TRIMMED_DIR="/data/scratch/projects/punim1214/nicoleRNAseq/N2419466_30-1070643705_RNA_2024-09-23/raw_data/trimmed"
ADAPTERS_FILE="/data/scratch/projects/punim1214/nicoleRNAseq/N2419466_30-1070643705_RNA_2024-09-23/raw_data/adapters.fa"

#create output directory if it doesn't exist
mkdir -p $TRIMMED_DIR

#Download the Illumina adapters file if it doesn't exist
#if [ ! -f $ADAPTERS_FILE ]; then
#    mkdir -p $(dirname $ADAPTERS_FILE)
#    wget -O $ADAPTERS_FILE "https://raw.githubusercontent.com/usadellab/Trimmomatic/master/adapters/TruSeq3-PE.fa"
#fi

#Loop through FASTQ files and run Trimmomatic
for fastq1 in $FASTQ_DIR/*_1.fq.gz; do
    #Get the base name of the file (without path and extension)
    base_name=$(basename "$fastq1" _1.fq.gz)

    #Define the corresponding paired file
    fastq2="${FASTQ_DIR}/${base_name}_2.fq.gz"

    #Define the output trimmed FASTQ files
    trimmed_paired1="${TRIMMED_DIR}/${base_name}_1_paired.fq.gz"
    trimmed_unpaired1="${TRIMMED_DIR}/${base_name}_1_unpaired.fq.gz"
    trimmed_paired2="${TRIMMED_DIR}/${base_name}_2_paired.fq.gz"
    trimmed_unpaired2="${TRIMMED_DIR}/${base_name}_2_unpaired.fq.gz"

    echo "Processing $base_name with Trimmomatic"

       # Run Trimmomatic
    java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.39.jar PE -phred33 \
        $fastq1 $fastq2 \
        $trimmed_paired1 $trimmed_unpaired1 \
        $trimmed_paired2 $trimmed_unpaired2 \
        ILLUMINACLIP:$ADAPTERS_FILE:2:30:10 \
        SLIDINGWINDOW:25:25 MINLEN:70

    echo "Trimming completed for $base_name"
done

echo "your shit just got trimmed."
