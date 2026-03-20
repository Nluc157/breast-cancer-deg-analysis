# ------------------------------------------------------------------------------
# Differential Gene Expression Analysis of TCGA BRCA Dataset 
#
# Author: Nicholas Lucido
# Description:
# This script performs differential gene expression analysis between
# tumor and control samples using limma. Downstream analyses include:
# - Volcano plot visualization
# - GO enrichment analysis
# - KEGG enrichment analysis
# - GSEA pathway analysis
# - Tumor stage comparison
#
# Input files:
#   data/filtered_counts_v2.csv
#   data/metadata_v2.csv
#
# Output files:
#   figures/
#   results/
# ------------------------------------------------------------------------------

library(limma)
library(dplyr)

filtered_counts <- read.csv("data/filtered_counts_v2.csv", row.names = 1)
metadata <- read.csv("data/metadata_v2.csv", row.names = 1)

head(metadata)
head(filtered_counts)

dim(filtered_counts)
dim(metadata)

# Check the first few column names (sample IDs) of filtered_counts (ignoring the first column)
head(colnames(filtered_counts)[2:ncol(filtered_counts)])

# Check the first few row names of metadata
head(rownames(metadata))

# Replace periods with dashes in filtered_counts column names (from the 2nd column onward)
colnames(filtered_counts)[2:ncol(filtered_counts)] <- gsub("\\.", "-", colnames(filtered_counts)[2:ncol(filtered_counts)])

# Replace periods with dashes in metadata row names
rownames(metadata) <- gsub("\\.", "-", rownames(metadata))

# Check if the sample IDs match
all(colnames(filtered_counts)[2:ncol(filtered_counts)] %in% rownames(metadata))




# -------------------------------------------------------------------------------




# Differential Gene Expression Analysis (Limma) of Groups 

# 1. Set 'Control' as the reference level in the Group column
metadata$Group <- factor(metadata$Group)
metadata$Group <- relevel(metadata$Group, ref = "Control")

# 2. Create design matrix
design <- model.matrix(~ Group, data = metadata)

# 3. Prepare expression matrix
expr_matrix <- as.matrix(filtered_counts[, -1])  # Assuming column 1 is gene name
rownames(expr_matrix) <- filtered_counts$Gene    # Or use rownames(filtered_counts) if preferred

# 4. Fit the linear model and compute statistics
fit <- lmFit(expr_matrix, design)
fit <- eBayes(fit)

# 5. Extract top differentially expressed genes comparing Cancer vs Control
results <- topTable(fit, coef = "GroupCancer", number = Inf)

# View the results
head(results)

# Number of significantly differentially expressed genes using FDR < 0.05
sum(results$adj.P.Val < 0.05)

# Filter for adjusted p-value < 0.05 and absolute logFC > 1
signif_genes <- results[results$adj.P.Val < 0.05 & abs(results$logFC) > 1, ]
nrow(signif_genes)

#Save to figures & results folder 
png("figures/logFC_distribution.png", width = 800, height = 600)

hist(results$logFC, breaks = 100, main = "Distribution of Log Fold Changes", xlab = "logFC")

dev.off()

write.csv(signif_genes, "results/filtered_DEGs.csv", row.names = TRUE)




# -------------------------------------------------------------------------------

  
  
  
# Volcano Plot Visual 
  
library(EnhancedVolcano)


# Sort the genes by adjusted p-value and logFC, and select the top 10 for labeling
top_genes_to_label <- head(signif_genes[order(signif_genes$adj.P.Val, abs(signif_genes$logFC)), ], 10)

# Save volcano plot to figures folder 
png("figures/volcano_plot.png", width = 1000, height = 800)

# Create volcano plot with all significant genes and only top 10 labeled, with bold labels
EnhancedVolcano(signif_genes,
                lab = ifelse(signif_genes$ID %in% top_genes_to_label$ID, signif_genes$ID, ""),  # Only label the top 10 genes
                x = 'logFC',
                y = 'adj.P.Val',
                title = 'Differentially Expressed Genes',
                pCutoff = 0.05,
                FCcutoff = 1,
                drawConnectors = TRUE,
                widthConnectors = 0.5,
                max.overlaps = 10,
                labFace = "bold"  # Makes labels bold
)

# Close the graphics device
dev.off()




#-------------------------------------------------------------------------------

  
  

# GO Analysis

library(clusterProfiler)
library(org.Hs.eg.db)     # Annotation database for humans
library(enrichplot)       # Nice enrichment visualizations
library(DOSE)             # Needed for enrichplot

# Get gene symbols
gene_symbols <- signif_genes$ID

# Convert gene symbols to Entrez IDs
entrez_ids <- bitr(gene_symbols,
                   fromType = "SYMBOL",
                   toType = "ENTREZID",
                   OrgDb = org.Hs.eg.db)

go_results <- enrichGO(gene         = entrez_ids$ENTREZID,
                       OrgDb        = org.Hs.eg.db,
                       keyType      = "ENTREZID",
                       ont          = "BP",            # Can be "BP", "MF", or "CC"
                       pAdjustMethod = "BH",
                       pvalueCutoff = 0.05,
                       qvalueCutoff = 0.2,
                       readable     = TRUE)            # Converts back to gene symbols


# Export GO bar plot 
png("figures/GO_barplot.png", width = 1000, height = 800)
barplot(go_results, showCategory = 10, title = "Top GO Biological Processes")
dev.off()

# Export GO dot plot
png("figures/GO_dotplot.png", width = 1000, height = 800)
dotplot(go_results, showCategory = 10, title = "GO Enrichment Dotplot")
dev.off()




# -------------------------------------------------------------------------------

  
  

# KEGG Enrichment Analysis

  kegg_results <- enrichKEGG(
    gene = entrez_ids$ENTREZID,
    organism = "hsa", 
    pvalueCutoff = 0.05
  )

# Export KEGG enrichment results
write.csv(as.data.frame(kegg_results),
          "results/KEGG_enrichment_results.csv")

# Export KEGG dot plot
png("figures/KEGG_dotplot.png", width = 1000, height = 800)
dotplot(kegg_results, showCategory = 10, title = "KEGG Pathway Enrichment")
dev.off()

# Export KEGG bar plot
png("figures/KEGG_barplot.png", width = 1000, height = 800)
barplot(kegg_results, showCategory = 10, title = "Top KEGG Pathways")
dev.off()




# -------------------------------------------------------------------------------

  


# GSEA: KEGG Pathway Analysis

library(clusterProfiler)
library(org.Hs.eg.db)
library(DOSE)
library(enrichplot)

# Prepare ranked gene list
gene_list <- signif_genes$logFC
names(gene_list) <- signif_genes$ID
gene_list <- sort(gene_list, decreasing = TRUE)

# Convert gene symbols to Entrez IDs
gene_df <- bitr(
  names(gene_list),
  fromType = "SYMBOL",
  toType   = "ENTREZID",
  OrgDb    = org.Hs.eg.db
)

# Merge Entrez IDs with logFC values
gene_df <- merge(
  gene_df,
  data.frame(SYMBOL = names(gene_list), logFC = gene_list),
  by = "SYMBOL"
)

# Keep highest absolute logFC per Entrez ID
gene_df <- gene_df[order(abs(gene_df$logFC), decreasing = TRUE), ]
gene_df <- gene_df[!duplicated(gene_df$ENTREZID), ]

# Create ranked vector for GSEA
gene_ranks <- gene_df$logFC
names(gene_ranks) <- gene_df$ENTREZID
gene_ranks <- sort(gene_ranks, decreasing = TRUE)

# Run KEGG GSEA
gsea_kegg <- gseKEGG(
  geneList     = gene_ranks,
  organism     = "hsa",
  minGSSize    = 10,
  maxGSSize    = 500,
  pvalueCutoff = 0.05,
  verbose      = FALSE
)

# Export results
write.csv(
  as.data.frame(gsea_kegg@result),
  "results/GSEA_KEGG_results.csv",
  row.names = FALSE
)

# Save KEGG GSEA dot plot
png("figures/GSEA_KEGG_dotplot.png", width = 1000, height = 800)
dotplot(gsea_kegg, showCategory = 10, title = "Top Enriched KEGG Pathways")
dev.off()

# Save enrichment plot for Cell Cycle pathway
png("figures/GSEA_CellCycle_Enrichment.png", width = 1000, height = 800)
gseaplot2(gsea_kegg, geneSetID = "hsa04110", title = "Cell Cycle Pathway")
dev.off()




# --------------------------------------------------------------------------------

  
  

# Expression Differences Across Tumor Stages


library(limma)

# Standardize sample IDs for matching
colnames(filtered_counts)[2:ncol(filtered_counts)] <- gsub(
  "\\.", "-",
  colnames(filtered_counts)[2:ncol(filtered_counts)]
)
rownames(metadata) <- gsub("\\.", "-", rownames(metadata))

# Create expression matrix from filtered counts
expr_matrix <- as.matrix(filtered_counts[, -1])
rownames(expr_matrix) <- filtered_counts$Gene
mode(expr_matrix) <- "numeric"

# Subset metadata to tumor samples only
tumor_metadata <- metadata[metadata$Group != "Control", ]

# Define tumor stage groups from Roman numeral stage labels
tumor_metadata$StageGroup <- ifelse(
  grepl("^III", tumor_metadata$Tumor.Stage), "Stage III",
  ifelse(
    grepl("^II", tumor_metadata$Tumor.Stage), "Stage II",
    ifelse(
      grepl("^IV", tumor_metadata$Tumor.Stage), "Stage IV",
      ifelse(
        grepl("^I", tumor_metadata$Tumor.Stage), "Stage I",
        NA
      )
    )
  )
)

# Remove samples with missing stage assignment
tumor_metadata <- tumor_metadata[!is.na(tumor_metadata$StageGroup), ]

# Exclude Stage IV due to insufficient sample size
tumor_metadata_subset <- tumor_metadata[tumor_metadata$StageGroup != "Stage IV", ]

# Match metadata and expression matrix sample IDs
tumor_samples_subset <- intersect(rownames(tumor_metadata_subset), colnames(expr_matrix))
tumor_metadata_subset <- tumor_metadata_subset[tumor_samples_subset, , drop = FALSE]
expr_tumor <- expr_matrix[, tumor_samples_subset, drop = FALSE]

# Set factor order for stage comparison
tumor_metadata_subset$StageGroup <- factor(
  tumor_metadata_subset$StageGroup,
  levels = c("Stage I", "Stage II", "Stage III")
)

# Check sample counts per stage
print(table(tumor_metadata_subset$StageGroup))

# Create design matrix
design_stage <- model.matrix(~ 0 + StageGroup, data = tumor_metadata_subset)
colnames(design_stage) <- levels(tumor_metadata_subset$StageGroup)

# Fit linear model and run empirical Bayes moderation
fit_stage <- lmFit(expr_tumor, design_stage)
fit_stage <- eBayes(fit_stage)

# Perform overall test for differences across Stage I, II, and III
stage_diff <- topTable(
  fit_stage,
  number = Inf,
  adjust.method = "BH",
  sort.by = "F"
)

# Export stage comparison results
write.csv(
  stage_diff,
  "results/tumor_stage_overall_F_test_results.csv",
  row.names = TRUE
)




# -------------------------------------------------------------------------------

  
  

# Export Table of Top Stage-Associated Genes

library(gridExtra)
library(grid)

# Select the top 15 stage-associated genes
top15 <- head(stage_diff, 15)

# Round numeric columns except p-values
for (col in names(top15)) {
  if (is.numeric(top15[[col]]) && !(col %in% c("P.Value", "adj.P.Val"))) {
    top15[[col]] <- round(top15[[col]], 2)
  }
}

# Format p-values in scientific notation
top15$P.Value <- formatC(top15$P.Value, format = "e", digits = 1)
top15$adj.P.Val <- formatC(top15$adj.P.Val, format = "e", digits = 1)

# Reformat table for clearer presentation
top15_display <- data.frame(
  Gene = top15$ID,
  EntrezID = rownames(top15),
  `Stage I` = top15$Stage.I,
  `Stage II` = top15$Stage.II,
  `Stage III` = top15$Stage.III,
  AvgExpr = top15$AveExpr,
  F = top15$F,
  P.Value = top15$P.Value,
  Adj.P.Val = top15$adj.P.Val,
  row.names = NULL
)

# Create graphical table
table_grob <- tableGrob(top15_display, rows = NULL)

# Save table as PNG with padding around the edges
png(
  filename = "figures/Top15_Stage_Associated_Genes.png",
  width = 2200,
  height = 900,
  res = 200
)

grid.newpage()

pushViewport(
  viewport(
    x = 0.5,
    y = 0.5,
    width = 0.95,
    height = 0.90
  )
)

grid.draw(table_grob)
popViewport()

dev.off()




# --------------------------------------------------------------------------------


  
  
# Boxplots of Top 6 Stage-Associated Genes Across Tumor Stages

library(ggplot2)
library(reshape2)
library(dplyr)

# Select the top 6 genes based on adjusted p-value
top6_genes <- head(stage_diff[order(stage_diff$adj.P.Val), ], 6)$ID

# Extract expression values for the top 6 genes
subset_expr <- expr_tumor[top6_genes, , drop = FALSE]

# Convert expression matrix to long format
expr_long <- melt(
  subset_expr,
  varnames = c("Gene", "Sample"),
  value.name = "Expression"
)

# Add stage information
expr_long <- expr_long %>%
  mutate(Sample = as.character(Sample)) %>%
  left_join(
    tumor_metadata_subset %>%
      mutate(Sample = rownames(tumor_metadata_subset)) %>%
      dplyr::select(Sample, StageGroup),
    by = "Sample"
  )

# Ensure correct stage order
expr_long$StageGroup <- factor(
  expr_long$StageGroup,
  levels = c("Stage I", "Stage II", "Stage III")
)

# Create faceted boxplots
p <- ggplot(expr_long, aes(x = StageGroup, y = Expression)) +
  
  geom_boxplot(
    fill = "#6BAED6",
    color = "black",
    outlier.shape = NA,
    width = 0.6
  ) +
  
  geom_jitter(
    width = 0.12,
    size = 1.2,
    alpha = 0.6,
    color = "black"
  ) +
  
  facet_wrap(~ Gene, scales = "free_y", ncol = 3) +
  
  labs(
    title = "Expression of Top 6 Stage-Associated Genes Across Tumor Stages",
    subtitle = "TCGA BRCA cohort",
    x = "Tumor Stage",
    y = "Normalized Expression"
  ) +
  
  theme_bw(base_size = 13) +
  
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    strip.background = element_rect(fill = "white", color = "black"),
    strip.text = element_text(face = "bold", size = 12),
    panel.grid.major = element_line(color = "gray90"),
    panel.grid.minor = element_blank()
  )

# Display plot
print(p)

# Export plot
ggsave(
  filename = "figures/Top6_Genes_Boxplot.png",
  plot = p,
  width = 10,
  height = 6,
  units = "in",
  dpi = 300
)




# -------------------------------------------------------------------------------

  
  

# Export Top 10 Differentially Expressed Genes to Word Document


library(officer)
library(flextable)

# Select top 10 differentially expressed genes from tumor vs control results
top10 <- head(results[, c("ID", "logFC", "AveExpr", "t", "P.Value", "adj.P.Val", "B")], 10)

# Rename columns for presentation
colnames(top10) <- c("Gene", "logFC", "AvgExpr", "t", "P.Value", "Adj.P.Val", "B")

# Round numeric columns except p-values
for (col in c("logFC", "AvgExpr", "t", "B")) {
  top10[[col]] <- round(top10[[col]], 2)
}

# Format p-values in scientific notation
top10$P.Value <- formatC(top10$P.Value, format = "e", digits = 1)
top10$Adj.P.Val <- formatC(top10$Adj.P.Val, format = "e", digits = 1)

# Create formatted flextable
ft <- flextable(top10)
ft <- autofit(ft)
ft <- set_caption(
  ft,
  caption = "Table 1. Top 10 Differentially Expressed Genes (Tumor vs Control) in the TCGA BRCA Cohort"
)

# Create Word document
doc <- read_docx()
doc <- body_add_par(doc, "Results", style = "heading 1")
doc <- body_add_par(doc, "Top 10 Differentially Expressed Genes", style = "heading 2")
doc <- body_add_flextable(doc, value = ft)

# Save Word document
print(doc, target = "results/Top10_Genes.docx")

  
  
  
  








