rm(list = ls())

library(DESeq2)

data <- readRDS("rds_objects/filtered_data.RDS")
classes <- readRDS("rds_objects/filtered_classes.RDS")

samples_info <- data.frame( 
  condition = factor(classes,levels = c("control","ipf"))
  )
# dataframe made to pass to DESeq2 to tell it this info
# that we have the label of each sample
# the order (control, ipf) is meant to be there for DESeq2
# to recognize healthy and diseased samples

rownames(samples_info) <- colnames(data)

dds <- DESeqDataSetFromMatrix(countData = data,
                              colData = samples_info,
                              design = ~ condition)

# DESeq2 data should have to be raw count data
# and not already normalized data
# the package normalizes the data on its own

dds <- DESeq(dds)
results <- results(dds)
head(results[order(results$padj),]) # order padj in descending order

significance_threshold <- 0.01 # for p-adjusted

significant_results <- results[which(results$padj < significance_threshold), ]

upregulated <- rownames(significant_results[significant_results$log2FoldChange > 0,])
downregulated <- rownames(significant_results[significant_results$log2FoldChange < 0,])

# Plotting the results ----

library(pheatmap)
N <- 100 # To select the 100 most significant differential expressed genes

top_genes <- head(rownames(significant_results[order(significant_results$padj),]), N)
# To get the order of adjusted P-values from lowest to highest
# i.e., from most significant to least significant
# then we get the rownames (gene names) of the most significant genes
# then head(rownames, 100) for top 100 genes in a dataframe

norm_counts <- counts(dds, normalized=TRUE)

heatmap_data <- norm_counts[top_genes,]

pdf("plots/heatmap_plot.pdf", width = 10, height = 8)
my_colors <- colorRampPalette(c("blue","white","red"))(50)

annotation_colors <- list(
  condition = c(control = "yellow", ipf = "orange")
)

pheatmap(heatmap_data, scale = "row", clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean", clustering_method = "complete",
         annotation_col = samples_info, annotation_colors = annotation_colors,
         color = my_colors, fontsize_row = 6, fontsize_col = 3)

dev.off()