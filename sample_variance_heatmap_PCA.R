#measure sample distances
sampleDists <- dist(t(normalized_counts))

#visualize as a heatmap

library(pheatmap)
library(RColorBrewer)
library(DESeq2)
library(ggplot2)

sampleDistMatrix <- as.matrix(sampleDists)
colours <- colorRampPalette(rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         col = colours)


#PCA
# Remove non-numeric columns (like gene names or other annotations)
PCA_matrix <- normalized_counts_dt[, .SD, .SDcols = patterns("ABI|Col|luh")]

# Transpose the matrix if genes are rows and samples are columns
# PCA expects samples as rows, so transpose the matrix if needed
t_PCA_matrix <- t(as.matrix(PCA_matrix))

# Remove columns (genes) with zero variance
t_PCA_matrix <- t_PCA_matrix[, apply(t_PCA_matrix, 2, var) > 0]

# Perform PCA
pca_result <- prcomp(t_PCA_matrix, scale. = TRUE)  # scale. = TRUE standardizes the data

# Extract PCA components for plotting
pca_data <- as.data.frame(pca_result$x)
# Create the 'treatment' column for a data frame
pca_data$treatment <- substr(rownames(pca_data), 1, 3)
pca_data$sample <- colnames(t_PCA_matrix)  # Add sample names if desired


# Calculate the percentage of variance explained by each PC
explained_variance <- pca_result$sdev^2
percent_variance <- explained_variance / sum(explained_variance) * 100

# Plot the PCA
ggplot(pca_data, aes(x = PC1, y = PC2, color = treatment)) +
  geom_point(size = 3) +
  xlab(paste0("PC1: ", round(percent_variance[1], 2), "% variance")) +
  ylab(paste0("PC2: ", round(percent_variance[2], 2), "% variance")) +
  ggtitle("PCA Plot") +
  theme_minimal()
