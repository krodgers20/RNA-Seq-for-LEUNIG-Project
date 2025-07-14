#pulling out the gene list of ups and downs in abi only 
#pull the genes which are not sig in CvL but are up or down in A vs both 

library(dplyr)
ABIGuys <- normalized_counts_dt_forEnhancedVolcano %>% filter(regulation_CvL == "Not Significant" & regulation_CvA != "Not Significant" & regulation_LvA != "Not Significant")

upABIGuys <- ABIGuys %>% filter(regulation_CvA == "Upregulated" & regulation_LvA == "Upregulated")
downABIGuys <- ABIGuys %>% filter(regulation_CvA == "Downregulated" & regulation_LvA == "Downregulated")
nrow(downABIGuys)
nrow(upABIGuys)
nrow(ABIGuys)

library(Biostrings)
output_file <- "ABI_downprots.fasta"
fasta_sequences <- readAAStringSet("/data/scratch/projects/punim1214/nicoleRNAseq/TAIR10_proteins.fasta", format="fasta")
for (id in downABIGuys$t_name) {
  # Find the sequences where the header matches the gene_id
  matching_sequences <- fasta_sequences[grep(id, names(fasta_sequences))]
  
  # If matches are found, write them to the output FASTA file
  if (length(matching_sequences) > 0) {
    writeXStringSet(matching_sequences, filepath = output_file, append = TRUE)
  }
}

#define the ABI prot fatsas
ABI_upfasta <- readAAStringSet("/data/scratch/projects/punim1214/nicoleRNAseq/N2419466_30-1070643705_RNA_2024-09-23/raw_data/ABI_prots.fasta", format="fasta")

# Extract gene IDs up to the first '.' character
ABIup_gene_ids <- sub("\\..*$", "", names(ABI_upfasta))

# Get indices of the first occurrence of each unique gene ID
unique_indices <- !duplicated(ABIup_gene_ids)

# Filter sequences to retain only the first instance of each unique gene ID
unique_UPfasta_sequences <- ABI_upfasta[unique_indices]
writeXStringSet(unique_UPfasta_sequences, filepath = "/data/scratch/projects/punim1214/nicoleRNAseq/N2419466_30-1070643705_RNA_2024-09-23/raw_data/ABI_UPprots_nodupes.fasta", append = TRUE)



#downregs
#define the ABI prot fatsas
ABI_downfasta <- readAAStringSet("/data/scratch/projects/punim1214/nicoleRNAseq/N2419466_30-1070643705_RNA_2024-09-23/raw_data/ABI_downprots.fasta", format="fasta")

# Extract gene IDs up to the first '.' character
ABIdown_gene_ids <- sub("\\..*$", "", names(ABI_downfasta))

# Get indices of the first occurrence of each unique gene ID
unique_indices <- !duplicated(ABIdown_gene_ids)

# Filter sequences to retain only the first instance of each unique gene ID
unique_DOWNfasta_sequences <- ABI_downfasta[unique_indices]
writeXStringSet(unique_DOWNfasta_sequences, filepath = "/data/scratch/projects/punim1214/nicoleRNAseq/N2419466_30-1070643705_RNA_2024-09-23/raw_data/ABI_DOWNprots_nodupes.fasta", append = TRUE)



#sort the up and down guys on pval and foldchange for abi 
upABIGuys_sorted <- upABIGuys %>%
  arrange(log2FC_LvA, log2FC_CvA, pluhVSABI, pColVsABI)  
downABIGuys_sorted <- downABIGuys %>%
  arrange(log2FC_LvA, log2FC_CvA, pluhVSABI, pColVsABI) 

# Remove duplicates based on substring in col1 up to the first period
upABIGuys_sorted_nodupes <- upABIGuys_sorted %>%
  mutate(substring_col1 = sub("\\..*$", "", t_name)) %>%  # Extract substring before the first period
  distinct(substring_col1, .keep_all = TRUE) %>%  # Remove duplicates based on this substring
  select(-substring_col1)  # Optionally remove the temporary substring column

downABIGuys_sorted_nodupes <- downABIGuys_sorted %>%
  mutate(substring_col1 = sub("\\..*$", "", t_name)) %>%  # Extract substring before the first period
  distinct(substring_col1, .keep_all = TRUE) %>%  # Remove duplicates based on this substring
  select(-substring_col1)  # Optionally remove the temporary substring column

#take only the top 2000 based on the above rankings for input into stringdb 
top2k_downABI <- downABIGuys_sorted_nodupes[1:2000, ]
top2k_upABI <- upABIGuys_sorted_nodupes[1:2000, ]

fasta_sequences <- readAAStringSet("/data/scratch/projects/punim1214/nicoleRNAseq/TAIR10_proteins.fasta", format="fasta")
output_file <- "top2k_downABI.fasta"
for (id in top2k_downABI$t_name) {
  # Find the sequences where the header matches the gene_id
  matching_sequences <- fasta_sequences[grep(id, names(fasta_sequences))]
  
  # If matches are found, write them to the output FASTA file
  if (length(matching_sequences) > 0) {
    writeXStringSet(matching_sequences, filepath = output_file, append = TRUE)
  }
}

output_file <- "top2k_upABI.fasta"
for (id in top2k_upABI$t_name) {
  # Find the sequences where the header matches the gene_id
  matching_sequences <- fasta_sequences[grep(id, names(fasta_sequences))]
  
  # If matches are found, write them to the output FASTA file
  if (length(matching_sequences) > 0) {
    writeXStringSet(matching_sequences, filepath = output_file, append = TRUE)
  }
}

writeLines(as.character(normalized_counts_dt_forEnhancedVolcano$t_name), "nicoleRNAseq_allt_names.txt")

# Remove duplicates based on substring in col1 up to the first period
normalized_counts_dt_forEnhancedVolcano_nodupes <- normalized_counts_dt_forEnhancedVolcano %>%
  mutate(substring_col1 = sub("\\..*$", "", t_name)) %>%  # Extract substring before the first period
  distinct(substring_col1, .keep_all = TRUE) %>%  # Remove duplicates based on this substring
  select(-substring_col1)  # Optionally remove the temporary substring column

normalized_counts_dt_forEnhancedVolcano_nodupes$t_name <- gsub("\\.*", "", normalized_counts_dt_forEnhancedVolcano_nodupes$t_name)

writeLines(as.character(normalized_counts_dt_forEnhancedVolcano_nodupes$t_name), "nicoleRNAseq_allt_names_3.txt")

write.csv(x = normalized_counts_dt_forEnhancedVolcano_nodupes, file = "RNASeqDataFromKelly.csv", row.names = FALSE, col.names = TRUE)


#getting non duped lists of DEGs for each of the three comparisons - will actually need two per comparisn because up and dow erg should be seaprate
library(dplyr)
CvL_upreg <- normalized_counts_dt_forEnhancedVolcano_nodupes %>% filter(regulation_CvL == "Upregulated")
CvL_downreg <- normalized_counts_dt_forEnhancedVolcano_nodupes %>% filter(regulation_CvL == "Downregulated")

CvA_upreg <- normalized_counts_dt_forEnhancedVolcano_nodupes %>% filter(regulation_CvA == "Upregulated")
CvA_downreg <- normalized_counts_dt_forEnhancedVolcano_nodupes %>% filter(regulation_CvA == "Downregulated")

LvA_upreg <- normalized_counts_dt_forEnhancedVolcano_nodupes %>% filter(regulation_LvA == "Upregulated")
LvA_downreg <- normalized_counts_dt_forEnhancedVolcano_nodupes %>% filter(regulation_LvA == "Downregulated")
