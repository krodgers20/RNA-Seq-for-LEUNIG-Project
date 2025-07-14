### for making volcano plots using the enhanced volcano package in r

#volcano plotting with EnhancedVolcano
if (!requireNamespace('BiocManager', quietly = TRUE))
  install.packages('BiocManager')

BiocManager::install('EnhancedVolcano')
library(EnhancedVolcano)
library(magrittr)

#dont run this#################
# Function to calculate log2 fold change between two specified columns
calculate_log2_fold_change <- function(col1, col2) {
  # Calculate fold change, adding a small constant to avoid division by zero
  fold_change <- (col2 + 1) / (col1 + 1)  # Avoid division by zero
  
  # Calculate log2 fold change
  log2_fold_change <- log2(fold_change)
  
  return(log2_fold_change)
}

#Col vs luh
normalized_counts_dt_unTFed$log2FC_CvL <- calculate_log2_fold_change(
  normalized_counts_dt_unTFed$ColMeanCount, 
  normalized_counts_dt_unTFed$luhMeanCount
)


#basic v plot - requires log2fold change, adjusted or unadjusted pvals, point labels 
EnhancedVolcano(normalized_counts_dt_unTFed,
                lab = normalized_counts_dt_unTFed$geneName,
                x = 'log2FC_CvL',
                y = 'pColVSluh',
                FCcutoff = 1,
                pCutoff = 0.05,
                title = 'Col vs luh')



#Col vs ABI
normalized_counts_dt_unTFed$log2FC_CvA <- calculate_log2_fold_change(
  normalized_counts_dt_unTFed$ColMeanCount, 
  normalized_counts_dt_unTFed$ABIMeanCount
)


#basic v plot - requires log2fold change, adjusted or unadjusted pvals, point labels 
EnhancedVolcano(normalized_counts_dt_unTFed,
                lab = normalized_counts_dt_unTFed$geneName,
                x = 'log2FC_CvA',
                y = 'pColVsABI',
                FCcutoff = 1,
                pCutoff = 0.05,
                title = 'Col vs ABI')

#luh vs ABI
normalized_counts_dt_unTFed$log2FC_LvA <- calculate_log2_fold_change(
  normalized_counts_dt_unTFed$luhMeanCount, 
  normalized_counts_dt_unTFed$ABIMeanCount
)


#basic v plot - requires log2fold change, adjusted or unadjusted pvals, point labels 
EnhancedVolcano(normalized_counts_dt_unTFed,
                lab = normalized_counts_dt_unTFed$geneName,
                x = 'log2FC_LvA',
                y = 'pluhVSABI',
                FCcutoff = 1,
                pCutoff = 0.05,
                title = 'luh vs ABI')

####finished dont run this zone###############
saveRDS(normalized_counts_dt_unTFed, "normalized_counts_dt_forEnhancedVolcano.RDS")

#safe to run below##############
#doing it with the enhanced volcano object
EnhancedVolcano(normalized_counts_dt_forEnhancedVolcano,
                lab = normalized_counts_dt_forEnhancedVolcano$geneName,
                x = 'log2FC_CvL',
                y = 'pColVSluh',
                FCcutoff = 1,
                pCutoff = 0.05,
                title = 'Col vs luh')

EnhancedVolcano(normalized_counts_dt_forEnhancedVolcano,
                lab = normalized_counts_dt_forEnhancedVolcano$geneName,
                x = 'log2FC_CvA',
                y = 'pColVsABI',
                FCcutoff = 1,
                pCutoff = 0.05,
                title = 'Col vs ABI')

EnhancedVolcano(normalized_counts_dt_forEnhancedVolcano,
                lab = normalized_counts_dt_forEnhancedVolcano$geneName,
                x = 'log2FC_LvA',
                y = 'pluhVSABI',
                FCcutoff = 1,
                pCutoff = 0.05,
                title = 'luh vs ABI')

#change gene name value based on t_name (AT id)
#normalized_counts_dt_forEnhancedVolcano$geneName[normalized_counts_dt_forEnhancedVolcano$t_name == "AT4G39950.1"] <- "CYP79B2"
