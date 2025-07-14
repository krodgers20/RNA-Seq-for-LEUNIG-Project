# Workflow
## 1. md5 sum checks on downloaded files
```
md5sum -c md5.md5
```
## 2. fastQC for quality checking 
```
#load module
module load FastQC/0.12.1-Java-11

#run for all fq.gz files in raw data directoy
fastqc *fq.gz
```
## 3. Download Arabidposis reference genome 
```
#make refs dir for reference genome: /data/scratch/projects/punim1214/nicoleRNAseq/refs
#download arabidopsis reference genome - used Col-0 version 10.1. manually downloaded from internet and scp file transferred from local machine. download link: https://www.ncbi.nlm.nih.gov/datasets/taxonomy/3702/ accessed 25/09/2024
```
## 4. Trim off adapters and low quality bases
```
#run "run_trimmomatic.sh" script. ensure to change directories and file names for this specific purpose.
```
## 5. Index reference genome
```
#index the reference genome using bwa
module load BWA/0.7.17
module load bwa-mem2/2.2.1
module load GCCcore/11.3.0

bwa index GCF_000001735.4_TAIR10.1_genomic.fna
```
## 6. Align reads to indexed reference genome with BWA
```
#run "bwa_to_bam.sh" script to align fastq read files to indexed reference genome. output in .bam format
```
## 7. Sort and index bam files
```
#run "sort_and_index_bam.sh" script to sort and index the bam files. This gives alignment and quantitative data in a file which is readily useable in downstream packages. (output sorted.bam)
```
## 8. Get quantitative data with stringtie 
```
#run "stringtie_batch.sh" script to run stringtie on all sorted bam files. stringtie directly quantifies aligned reads. including the -B option in the command prepares the output for the downstream R package Ballgown.
#stringtie requires a .gtf annotation file for the reference genome to run in this way
```
## 9. Move to R, load packages and set working directory
```
library(data.table)
library(dplyr)
library(ggplot2)
library(stringr)
library(tidyverse)

# Define the main directory containing sample-specific subdirectories
data_dir <- "/data/scratch/projects/punim1214/nicoleRNAseq/N2419466_30-1070643705_RNA_2024-09-23/raw_data/stringtie_out"
```
## 10. Extract quantitative data from stringtie output files
```
library(data.table)
library(dplyr)
library(ggplot2)
library(stringr)
library(tidyverse)

# Define the main directory containing sample-specific subdirectories
data_dir <- "/data/scratch/projects/punim1214/nicoleRNAseq/N2419466_30-1070643705_RNA_2024-09-23/raw_data/stringtie_out"

# Find all e_data.ctab and t_data.ctab files recursively in subdirectories
e_data_files <- list.files(data_dir, pattern = "e_data.ctab$", recursive = TRUE, full.names = TRUE)
t_data_files <- list.files(data_dir, pattern = "t_data.ctab$", recursive = TRUE, full.names = TRUE)

# Initialize a list to store combined results
all_samples_data <- list()  # Store all samples data

# Loop through each e_data file
for (e_file in e_data_files) {
  sample_name <- sub(".*/(.*?)-LGJ\\d+_L2/e_data.ctab$", "\\1", e_file)  # Get sample name
  cat("Processing sample:", sample_name, "\n")  # Debug print

  e_data <- fread(e_file)  # Read e_data

  # Load the corresponding e2t.ctab file from the same subdirectory
  e2t_file <- sub("e_data.ctab$", "e2t.ctab", e_file)  # Get corresponding e2t file
  e2t <- fread(e2t_file)  # Read e2t data

  # Merge e_data with e2t to add the t_id column
  merged_data <- merge(e_data, e2t, by = "e_id", all.x = TRUE)

  # Aggregate raw counts by t_id
  transcript_counts <- merged_data[, .(total_rcount = sum(rcount)), by = t_id]

  # Load the corresponding t_data file
  t_data_file <- sub("e_data.ctab$", "t_data.ctab", e_file)  # Get corresponding t_data file
  if (file.exists(t_data_file)) {
    t_data <- fread(t_data_file)  # Read t_data

    # Merge t_data with the transcript counts to get t_name
    t_data <- t_data[, .(t_id, t_name)]  # Select relevant columns from t_data
    final_data <- merge(transcript_counts, t_data, by = "t_id", all.x = TRUE)

    # Add sample information to final_data
    final_data[, sample := sample_name]  # Assign sample name to each row

    # Store the final aggregated data in the list
    all_samples_data[[sample_name]] <- final_data
  } else {
    cat("Corresponding t_data file not found for sample:", sample_name, "\n")
  }
}

# Combine all samples data into one data.table
all_samples_combined <- rbindlist(all_samples_data)

# Display the combined data
print(head(all_samples_combined))
```
## 11. Get data ready for use with DESeq2, normalize data
```
library(DESeq2)

# Convert the data.table into a wide format for DESeq2
count_matrix <- dcast(all_samples_combined, t_id ~ sample, value.var = "total_rcount")

# Remove the first column (t_id) for the counts matrix
rownames(count_matrix) <- count_matrix$t_id
count_matrix$t_id <- NULL  # Remove t_id column

# Create sample metadata

sample_info <- data.table(
  sample = colnames(count_matrix),
  treatment = rep(c("Col", "luh", "ABI"), each = 3)
)

# Convert 'treatment' to a factor
sample_info$treatment <- as.factor(sample_info$treatment)

# Create DESeqDataSet object
dds <- DESeqDataSetFromMatrix(countData = count_matrix, colData = sample_info, design = ~ treatment)

#estimate size factors
dds <- estimateSizeFactors(dds)

# Get normalized counts
normalized_counts <- counts(dds, normalized = TRUE)

# Ensure the row names are set to t_id
rownames(normalized_counts) <- dds$t_id


plot(as.matrix(count_matrix) %>% log10,as.matrix(normalized_counts) %>% log10,pch=20,col="#00000055",cex=.2)


# Convert to data.table
normalized_counts_dt <- as.data.table(normalized_counts)

# Create a new column 't_id' filled with the row numbers
normalized_counts_dt[, t_id := .I]

#add t_name column
# Merge to add t_name based on t_id
normalized_counts_dt <- merge(normalized_counts_dt,
                              all_samples_combined[, .(t_id, t_name)],
                              by = "t_id",
                              all.x = TRUE)
normalized_counts_dt <- unique(normalized_counts_dt)

saveRDS(normalized_counts_dt, "normalized_counts_dt.RDS")


# Set "Col" as the reference level for the treatment factor
dds$treatment <- relevel(dds$treatment, ref = "Col")

# Run the DESeq analysis
#dds <- DESeq(dds)

# Get normalized counts
normalized_counts <- counts(dds, normalized = TRUE)
```
## 12. ANOVA on normalized data 
```
# Set "Col" as the reference level for the treatment factor
dds$treatment <- relevel(dds$treatment, ref = "Col")

# Run the DESeq analysis
#dds <- DESeq(dds)

# Get normalized counts
normalized_counts <- counts(dds, normalized = TRUE)


# Define the treatment columns
treatment_columns <- c("ABI4", "ABI5", "ABI6", "Col5", "Col6", "Col7", "luh4", "luh5", "luh6")

# Reshape the data for ANOVA calculations
normalized_counts_dt_long <- melt(normalized_counts_dt, id.vars = "t_name",
                                  measure.vars = treatment_columns,
                                  variable.name = "replicate",
                                  value.name = "expression")

# Create a new column for treatment based on the replicate name
normalized_counts_dt_long[, treatment := fifelse(replicate %in% c("ABI4", "ABI5", "ABI6"), "ABI",
                                                 fifelse(replicate %in% c("Col5", "Col6", "Col7"), "Col",
                                                         "luh"))]

# Calculate p-values and differences using data.table syntax
normalized_counts_dt_unTFed <- normalized_counts_dt_long[, .(
  pColVsABI = getAnovaPval(expression, treatment),
  dColVsABI = mean(expression[treatment == "Col"]) - mean(expression[treatment == "ABI"]),

  pColVSluh = getAnovaPval(expression[treatment %in% c("Col", "luh")],
                           treatment[treatment %in% c("Col", "luh")]),
  dColVSluh = mean(expression[treatment == "Col"]) - mean(expression[treatment == "luh"]),

  pluhVSABI = getAnovaPval(expression[treatment %in% c("luh", "ABI")],
                           treatment[treatment %in% c("luh", "ABI")]),
  dluhVsABI = mean(expression[treatment == "luh"]) - mean(expression[treatment == "ABI"])
), by = .(t_name)]

# Save results to an RDS file
saveRDS(normalized_counts_dt_unTFed, "normalized_counts_dtunTFed.RDS")
```
## 13. Annotate genes of interest in data table
```
#add geneName and purpose
#remove "transcript:" substring from t_name values
normalized_counts_dt_unTFed[, t_name := gsub("transcript:", "", t_name)]


# Create a new data.table to hold the results
matched_data <- normalized_counts_dt_unTFed[, .(t_name, geneName = NA_character_, purpose = NA_character_)]

# Loop through gene_data and find substring matches
for (i in 1:nrow(gene_data)) {
  # Get the current geneID
  gene_id_substring <- gene_data$geneID[i]

  # Check for matches in pvals
  matches <- str_detect(normalized_counts_dt_unTFed$t_name, gene_id_substring)

  # Update matched_data with corresponding geneName and purpose
  matched_data[matches, geneName := gene_data$geneName[i]]
  matched_data[matches, purpose := gene_data$purpose[i]]
}

# Merge the matched data with pvals
normalized_counts_dt_unTFed <- merge(normalized_counts_dt_unTFed, matched_data[, .(t_name, geneName, purpose)],
                               by = "t_name", all.x = TRUE)
```
## 14. Get columns for each treatment and comparisons between treatments, save as .csv 
```
#add mean tscript tpm value column for each transcript in each sample into a new dt
normalized_counts_dt_unTFed2 <- normalized_counts_dt_long[,.(
  ColMeanCount = .SD[treatment=="Col",mean(expression)],
  luhMeanCount = .SD[treatment=="luh",mean(expression)],
  ABIMeanCount = .SD[treatment=="ABI",mean(expression)]
),by=.(t_name)]

#remove transcript: from the t_names for merging
normalized_counts_dt_unTFed2[, t_name := gsub("transcript:", "", t_name)]

#merge on t_name
normalized_counts_dt_unTFed <- merge(normalized_counts_dt_unTFed,
                       normalized_counts_dt_unTFed2[, .(t_name, ColMeanCount, luhMeanCount, ABIMeanCount)],
                       by = "t_name",
                       all.x = TRUE)

#save some RDSes
saveRDS(normalized_counts_dt_unTFed, "normalized_counts_dt_unTFed.RDS")


#col based on significance and regulation
normalized_counts_dt_unTFed[, regulation_CvA := ifelse(pColVsABI < 0.05, # Check if p-value is significant
                                                 ifelse(dColVsABI > 1, "Downregulated",   # Downregulated
                                                        ifelse(dColVsABI < -1, "Upregulated", NA)),  # Upregulated
                                                 "Not Significant")]  # Not significant

normalized_counts_dt_unTFed[, regulation_CvL := ifelse(pColVSluh < 0.05, # Check if p-value is significant
                                                       ifelse(dColVSluh > 1, "Downregulated",   # Downregulated
                                                              ifelse(dColVSluh < -1, "Upregulated", NA)),  # Upregulated
                                                       "Not Significant")]  # Not significant


normalized_counts_dt_unTFed[, regulation_LvA := ifelse(pluhVSABI < 0.05, # Check if p-value is significant
                                                       ifelse(dluhVsABI > 1, "Downregulated",   # Downregulated
                                                              ifelse(dluhVsABI < -1, "Upregulated", NA)),  # Upregulated
                                                       "Not Significant")]  # Not significant

#specify cols of interest
#col vs luh
ColVluh_columns <- c("t_name", "geneName", "purpose", "ColMeanCount", "luhMeanCount", "dColVSluh", "pColVSluh", "regulation_CvL")  # specify the columns by name

# Save the specified columns to a CSV file
fwrite(normalized_counts_dt_unTFed[, ..ColVluh_columns], file = "ColVsluh_DEG.csv")

#col vs abi
ColVABI_columns <- c("t_name", "geneName", "purpose", "ColMeanCount", "ABIMeanCount", "dColVsABI", "pColVsABI", "regulation_CvA")  # specify the columns by name

# Save the specified columns to a CSV file
fwrite(normalized_counts_dt_unTFed[, ..ColVABI_columns], file = "ColVsABI_DEG.csv")

#luh vs abi
luhVABI_columns <- c("t_name", "geneName", "purpose", "luhMeanCount", "ABIMeanCount", "dluhVsABI", "pluhVSABI", "regulation_LvA")  # specify the columns by name

# Save the specified columns to a CSV file
fwrite(normalized_counts_dt_unTFed[, ..luhVABI_columns], file = "luhVsABI_DEG.csv")
```
