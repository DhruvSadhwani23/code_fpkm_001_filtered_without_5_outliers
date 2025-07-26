load("/mnt/808e843c-fb08-404a-9f45-efd817d81a77/srip2025_dhruv_fshd_wgcna/srip2025_dhruv_fshd_wgcna_115650_001_fpkm/srip2025_dhruv_fshd_wgcna_115650_001_downstream_analysis/115650_FPKM_001_DGE_analysis/volcano_log2_fold_change_final/WGCNA_115650_001_fpkm_log2_fold_change_workspace.RData")
#this is for loading this particular workspace

#loading the wgcna workspace here 
load("/mnt/808e843c-fb08-404a-9f45-efd817d81a77/srip2025_dhruv_fshd_wgcna/srip2025_dhruv_fshd_wgcna_115650_001_fpkm/srip2025_dhruv_fshd_wgcna_115650_001_WGCNA_analysis/WGCNA_115650_001_fpkm_workspace.RData")
samples.to.be.excluded <- c('GSM3186678', 'GSM3186679', 'GSM3186673', 'GSM3186676', 'GSM3186709')
log_dataOPsubset_INPUT <- log_dataOPsubset[, !(colnames(log_dataOPsubset) %in% samples.to.be.excluded)]
log_dataOPsubset_INPUT <- t(log_dataOPsubset_INPUT)
getwd()
############  PERFORMING DGE USING VIMMA   ############

# Load libraries
library(limma)
library(ggplot2)

# 1. Align sample order
log_dataOPsubset_INPUT <- log_dataOPsubset_INPUT[rownames(phenoDataOP), ]

# 2. Design matrix
design <- model.matrix(~ status_binary, data = phenoDataOP)

# 3. Transpose for limma
expr_for_limma <- t(log_dataOPsubset_INPUT)

# 4. Fit model
fit <- lmFit(expr_for_limma, design)
fit <- eBayes(fit)

# 5. Get all results
res <- topTable(fit, coef = "status_binary", number = Inf, sort.by = "none")
res$Gene <- rownames(res)

# 6. Mark significance
res$Significance <- ifelse(
  res$adj.P.Val < 0.05 & abs(res$logFC) > 1,
  ifelse(res$logFC > 0, "Upregulated", "Downregulated"),
  "Not Significant"
)

# 7. Volcano plot
ggplot(res, aes(x = logFC, y = -log10(adj.P.Val), color = Significance)) +
  geom_point(alpha = 0.7, size = 1.5) +
  scale_color_manual(values = c(
    "Upregulated" = "red",
    "Downregulated" = "blue",
    "Not Significant" = "grey"
  )) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  labs(
    title = "Volcano Plot: FSHD vs Control (limma)",
    x = "log2 Fold Change (FSHD vs Control)",
    y = "-log10 Adjusted p-value",
    color = "Significance"
  ) +
  theme_minimal()


#the results df still has some missing genes, no ensembl ids as well, checking for why such is happening
#saving the results' file as a csv, we need to check for missing genes and also look at the results properly
write.csv(res, "significant_genes_volcano_plot.csv", row.names = TRUE)


#HERE WHAT WE HAVE DONE IS FIRSTLY WE HAVE PLOTTED THE VOLCANO PLOT FOR THE LOG2 FOLD CHANGE USING GGPLOT, NOW WE'LL 
#CONVERT THE ENSEMBL IDS TO GENE SYMBOLS
#EARLIER MAYBE THERE WAS SOME PROBLEM DURING CONVERSION WHICH LED TO SOME MISSING GENE IDS HENCE WE HAVE CHANGED THE ORDER HERE 

# 1. Extract and clean Ensembl IDs from res rownames
clean_ids <- sub("\\..*", "", rownames(res))

# 2. Load your Ensembl-to-symbol mapping file
mapping <- read.csv("mart_export.txt", stringsAsFactors = FALSE)

# 3. Rename mapping columns to standard names
colnames(mapping) <- c("Ensembl.Gene.ID", "Gene.name")

# 4. Create a mapping table from cleaned Ensembl IDs and original rownames
id_table <- data.frame(
  Ensembl.Gene.ID = clean_ids,
  original_rownames = rownames(res),
  stringsAsFactors = FALSE
)

# 5. Merge with mapping
annotated_ids <- merge(
  id_table,
  mapping,
  by = "Ensembl.Gene.ID",
  all.x = TRUE
)

# 6. Restore the original order to match res
annotated_ids <- annotated_ids[match(clean_ids, annotated_ids$Ensembl.Gene.ID), ]

# 7. Create final gene name column: use gene symbol if available, otherwise cleaned Ensembl ID
annotated_ids$final_gene_name <- ifelse(
  is.na(annotated_ids$Gene.name),
  annotated_ids$Ensembl.Gene.ID,
  annotated_ids$Gene.name
)

# 8. Add BOTH columns to the res dataframe
res$Ensembl_ID <- annotated_ids$Ensembl.Gene.ID
res$Gene       <- annotated_ids$final_gene_name

# 9. (Optional) Save updated results
write.csv(res, "limma_results_with_gene_names.csv", row.names = FALSE)


#now in these results as well, we will separate the csvs with upregulated, downregulated and the non significant genes
#and for genes whose names have again vanished due to some reason, we could  refer the ensembl ids csv to get their ensembl ids  


# Split
upregulated_genes     <- subset(res, Significance == "Upregulated")
downregulated_genes   <- subset(res, Significance == "Downregulated")
non_significant_genes <- subset(res, Significance == "Not Significant")

# Sort
upregulated_genes     <- upregulated_genes[order(-upregulated_genes$logFC), ]
downregulated_genes   <- downregulated_genes[order(downregulated_genes$logFC), ]
non_significant_genes <- non_significant_genes[order(-abs(non_significant_genes$logFC)), ]

# Check in console
head(upregulated_genes)
head(downregulated_genes)
head(non_significant_genes)

# Optional RStudio viewer
View(upregulated_genes)
View(downregulated_genes)
View(non_significant_genes)

# Save
write.csv(upregulated_genes, "Upregulated_genes_sorted.csv", row.names = FALSE)
write.csv(downregulated_genes, "Downregulated_genes_sorted.csv", row.names = FALSE)
write.csv(non_significant_genes, "NonSignificant_genes_sorted.csv", row.names = FALSE)


save.image("WGCNA_115650_001_fpkm_log2_fold_change_workspace.RData")

