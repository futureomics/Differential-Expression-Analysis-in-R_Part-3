# Step 3: Differential Expression Analysis (DEA)

if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install("DESeq2")
library(DESeq2)

# Prepare DESeq2 data
dds <- DESeqDataSetFromMatrix(countData = round(expression_data), 
                              colData = data.frame(condition = subtype), 
                              design = ~ condition)
dds <- DESeq(dds)
results <- results(dds, contrast = c("condition", "Her2", "TNBC"))

# Filter significant DEGs
significant_genes <- results[!is.na(results$padj) & results$padj < 0.05, ]
cat("Number of significant DEGs:", nrow(significant_genes), "\n")
write.csv(as.data.frame(significant_genes), file = file.path("significant_genes_her2_vs_tnbc.csv"))

# Volcano plot
library(ggplot2)
#Install the package
install.packages("ggrepel")

# Load the package into your current session
library(ggrepel)
install.packages("magrittr") # package installations are only needed the first time you use it
install.packages("dplyr")    # alternative installation of the %>%
library(magrittr) # needs to be run every time you start R and want to use %>%
library(dplyr)    # alternatively, this also loads %>%

volcano_df <- as.data.frame(results)
volcano_df$gene <- rownames(volcano_df)
volcano_df$significant <- ifelse(volcano_df$padj < 0.05, "Significant", "Not significant")

# Volcano plot with top 10 DEGs labeled
top_genes <- volcano_df %>%
  filter(padj < 0.005) %>%
  arrange(padj) %>%
  head(10)

volcano_with_labels_plot <- ggplot(volcano_df, aes(x = log2FoldChange, y = -log10(padj), color = significant)) +
  geom_point(alpha = 0.7) +
  scale_color_manual(values = c("gray", "red")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "blue") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "blue") +
  geom_label_repel(data = top_genes, aes(label = gene), box.padding = 0.5, max.overlaps = Inf) +
  labs(title = "Volcano Plot with Top 10 DEGs",
       x = "log2 Fold Change (Her2/TNBC)",
       y = "-log10(Adjusted p-value)") +
  theme_minimal()

ggsave(file.path("volcano_with_labels_plot.png"), plot = volcano_with_labels_plot)

# MA Plot (Log2 fold change vs mean expression)
ma_plot <- ggplot(volcano_df, aes(x = log2(baseMean), y = log2FoldChange)) +
  geom_point(alpha = 0.7, color = "gray") +
  geom_point(data = subset(volcano_df, padj < 0.05), aes(color = significant), size = 2) +
  scale_color_manual(values = c("red")) +
  labs(title = "MA Plot", x = "Log2 Mean Expression", y = "Log2 Fold Change") +
  theme_minimal()

ggsave(file.path("ma_plot.png"), plot = ma_plot)
