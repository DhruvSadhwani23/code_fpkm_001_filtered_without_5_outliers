load("/mnt/808e843c-fb08-404a-9f45-efd817d81a77/srip2025_dhruv_fshd_wgcna/srip2025_dhruv_fshd_wgcna_115650_001_fpkm/srip2025_dhruv_fshd_wgcna_115650_001_downstream_analysis/115650_FPKM_001_DGE_analysis/violin_plots/outputs/violin_plots_115650_fpkm_001_workspace.RData")

#now we'll be plotting the violin plots for the expression of top 10 hub genes for each of  red and brown( all top 25 for turquoise)
# 2 genes of interest of green as well
#also for clock genes (CLOCK, BMAL1 (ARNTL), PER1/2/3, CRY1/2, NR1D1/2 (REV-ERBα/β), RORA/B/C)

# and few other imp genes: MYOD1, TCAP, GRB10, PAX7, GSK3 beta, wnt pathway molecules
#WNT pathway molecules include: WNT3A (ligand), FZD (Frizzled receptor), LRP5/6 (co-receptor), DVL (Dishevelled), AXIN, APC, GSK3β, β-catenin, TCF/LEF (transcription factors), DKK1, SFRP1, WIF1 (extracellular inhibitors) 

#firstly loading the wgcna workspace and the log2 fold change workspace for getting list of genes with gene symbols

# 1️⃣ Strip version from expression matrix colnames
clean_ensembl_ids <- sub("\\..*", "", colnames(log_dataOPsubset_INPUT))

# 2️⃣ Build mapping dataframe
mapping <- data.frame(
  Ensembl_ID = sub("\\..*", "", rownames(res)),
  Gene_Symbol = res$Gene,
  stringsAsFactors = FALSE
)

# 3️⃣ Match colnames to Ensembl IDs
col_info <- data.frame(
  Original_Colname = colnames(log_dataOPsubset_INPUT),
  Clean_Ensembl_ID = clean_ensembl_ids,
  stringsAsFactors = FALSE
)

# 4️⃣ Merge colnames with mapping
annotated_cols <- merge(
  col_info,
  mapping,
  by.x = "Clean_Ensembl_ID",
  by.y = "Ensembl_ID",
  all.x = TRUE
)

# 5️⃣ Assign fallback names where mapping fails
annotated_cols$Final_Name <- ifelse(
  is.na(annotated_cols$Gene_Symbol) | annotated_cols$Gene_Symbol == "",
  annotated_cols$Clean_Ensembl_ID,  # fallback to cleaned Ensembl ID
  annotated_cols$Gene_Symbol
)

# 6️⃣ Restore original column order
annotated_cols <- annotated_cols[match(colnames(log_dataOPsubset_INPUT), annotated_cols$Original_Colname), ]

# 7️⃣ Ensure no duplicates
annotated_cols$Final_Name <- make.unique(annotated_cols$Final_Name)

# 8️⃣ Apply to matrix
log_dataOPsubset_INPUT_annotated <- log_dataOPsubset_INPUT
colnames(log_dataOPsubset_INPUT_annotated) <- annotated_cols$Final_Name


### CODE FOR PLOTTING VIOLIN PLOTS N CALC P VALS USING WILCOXON RANK SUM TEST FOR GENES OF INTEREST

library(ggplot2)
library(reshape2)
library(dplyr)

# 1️⃣ Match common samples
common_samples <- intersect(rownames(log_dataOPsubset_INPUT_annotated), rownames(phenoDataOP))
log_dataOPsubset_INPUT_annotated <- log_dataOPsubset_INPUT_annotated[common_samples, ]
phenoDataOP <- phenoDataOP[common_samples, ]

# 2️⃣ Genes of interest
my_genes_of_interest <- c(  "MMP14", "SFXN3","ARF3","CALHM2","STING1","C1R","GNB1","C1S","OS9","FKBP10",
                            "PPIB","GBA1","ARHGAP1","DDAH2","RUSC1","PLXNB2","MARVELD1","RAB5C","TMEM214",
                            "PDIA4","NECAP2","MRGPRF","PLOD3","CD44","COLGALT1"
                            
)

# 3️⃣ Filter for present genes
present_genes <- my_genes_of_interest[my_genes_of_interest %in% colnames(log_dataOPsubset_INPUT_annotated)]
if (length(present_genes) == 0) stop("None of the requested genes are in the expression matrix!")

# 4️⃣ Subset and reshape
expr_subset <- as.data.frame(log_dataOPsubset_INPUT_annotated[, present_genes, drop = FALSE])
expr_subset$SampleID <- rownames(expr_subset)
expr_long <- melt(expr_subset, id.vars = "SampleID", variable.name = "Gene", value.name = "Expression")

# 5️⃣ Merge phenotype
phenoDataOP$SampleID <- rownames(phenoDataOP)
expr_long <- merge(expr_long, phenoDataOP[, c("SampleID", "status_binary")], by = "SampleID")
expr_long$Status <- factor(expr_long$status_binary, levels = c(0, 1), labels = c("Control", "FSHD"))

# 6️⃣ Wilcoxon test + label
pvals <- expr_long %>%
  group_by(Gene) %>%
  summarise(
    p_value = wilcox.test(Expression ~ Status)$p.value
  ) %>%
  mutate(
    sig_label = ifelse(p_value < 0.05, 
                       paste0("Significant (p = ", signif(p_value, 3), ")"),
                       paste0("Not Significant (p = ", signif(p_value, 3), ")")),
    Gene_labeled = paste0(Gene, "\n", sig_label)
  )

# 7️⃣ Join back to use labeled gene names in facet
expr_long <- left_join(expr_long, pvals[, c("Gene", "Gene_labeled")], by = "Gene")

# 8️⃣ Plot with new facet labels
p <- ggplot(expr_long, aes(x = Status, y = Expression, fill = Status)) +
  geom_violin(trim = FALSE, alpha = 0.6) +
  geom_jitter(width = 0.1, alpha = 0.4, size = 0.6) +
  facet_wrap(~ Gene_labeled, scales = "fixed") +
  ylim(-1, 7.5) +
  labs(
    title = "Expression Comparison: Control vs FSHD",
    x = "",
    y = "Log-normalized Expression"
  ) +
  theme_minimal() +
  theme(legend.position = "none",
        strip.text = element_text(size = 10, face = "bold"))

# 9️⃣ Save
print(p)
ggsave("ViolinPlots_TurquoiseGenes.png", plot = p, width = 10, height = 6, dpi = 300)
save.image("violin_plots_115650_fpkm_001_workspace.RData")
