# --------------------------------------
# INSTALLING REQUIRED PACKAGES & LIBRARIES
# --------------------------------------

# BiocManager is required to install Bioconductor packages like clusterProfiler, org.Hs.eg.db, etc.
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# Installing 'clusterProfiler' for GO/KEGG enrichment and other functional analyses
# Initially this may fail due to unmet dependencies (e.g. fgsea, DOSE)
BiocManager::install("clusterProfiler")

# Installing 'AnnotationDbi' - needed for mapping and querying biological annotation databases
BiocManager::install("AnnotationDbi")

# Installing 'org.Hs.eg.db' - the human genome annotation database (Entrez ID, symbols, etc.)
BiocManager::install("org.Hs.eg.db")

# --------------------------------------
# WORKAROUND FOR INSTALLATION ERRORS
# --------------------------------------

# The initial attempt to install clusterProfiler and its dependencies (fgsea, DOSE, enrichplot) may fail due to:
# - Compilation issues (missing system libraries or permissions)
# - Inaccessible system-wide library paths (/usr/lib/R/library)

# To fix these, we switched to using GitHub versions via 'remotes', which downloads source builds
# into user-accessible paths and avoids system library permission issues.

# Make sure 'remotes' is installed to fetch packages directly from GitHub
if (!requireNamespace("remotes", quietly = TRUE))
  install.packages("remotes")

# fgsea (fast GSEA) initially failed with a "ScoreRuler.o" compilation error
# This was resolved by installing necessary system dependencies in the terminal:
#firstly ran the command ssh -X tanuja@10.0.68.137 then entered the password of the system(Compbio@iitgn), to get the terminal of the system setup in my laptop
# sudo apt-get update
# sudo apt-get install libcurl4-openssl-dev libssl-dev libxml2-dev build-essential
#from the pass n stuff, we got the access to sudo user then ran the above commands to fix the permission errors we were getting


# Then we successfully installed it from GitHub:
remotes::install_github("ctlab/fgsea")

# DOSE and enrichplot had failed earlier due to unmet fgsea and system-level dependencies.
# Installing them *after* fgsea was installed fixed the problem.
# They were then successfully installed using:
remotes::install_github("YuLab-SMU/DOSE")
remotes::install_github("YuLab-SMU/enrichplot")

# Just to be sure, re-install clusterProfiler from GitHub now that all dependencies are in place.
remotes::install_github("YuLab-SMU/clusterProfiler")

# --------------------------------------
# LOADING LIBRARIES
# --------------------------------------

# If no errors appear below, all libraries are successfully installed.
library(clusterProfiler)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(fgsea)        # Now works after GitHub install
library(enrichplot)   # Now loads correctly
library(DOSE)        # Now loads correctly

# --------------------------------------
# OPTIONAL: Drop Ensembl Version Numbers
# --------------------------------------

# Many Ensembl gene IDs come in the format "ENSG00000123456.1"
# Some databases/tools don't recognize IDs with version numbers.
# This line strips the version number after the period.

# Example input vector:
# gene_ids <- c("ENSG00000123456.1", "ENSG00000234567.2")

# Code to remove version numbers:

#for this, firstly loading the dataset and copying the pre-processing code from the WGCNA pipeline to get our list of genes (emsembl ids first which we'll then convert
# into gene symbols for getting the gene ontology enrichment analysis results in gene symbols which are easy to interpret innstead of ensembl ids)

### Loading the workspace:
load("/mnt/808e843c-fb08-404a-9f45-efd817d81a77/srip2025_dhruv_fshd_wgcna/srip2025_dhruv_fshd_wgcna_115650_001_fpkm/GO_analysis_115650_001_fpkm_workspace.RData")


################### LOADING DATASET AND NECESSARY LIBRARIES
setwd("/home/tanuja/Documents/srip2025_dhruv_fshd_wgcna/srip2025_dhruv_fshd_wgcna_115650_001_fpkm")
getwd()
dataOP <- read.csv("~/Documents/srip2025_dhruv_fshd_wgcna/srip2025_dhruv_fshd_wgcna_115650_001_fpkm/srip2025_dhruv_fshd_wgcna_115650_001_dataset.csv", header=T, stringsAsFactors = FALSE)
library(WGCNA)  #we did that by writing necessary linux code in terminal
library(DESeq2)
library(GEOquery)
library(tidyverse)
library(gridExtra)
library(CorLevelPlot)
#ADDING PHENODATA TO GET THE TRAITS AS WELL, WHICH WILL HELP US IDENTIFY MODULE TRAITS RELATIONSHIP AND ALSO GENE OF SIGNIIFICANCE
geo_id <- "GSE115650"
gse <- getGEO(geo_id, GSEMatrix = TRUE)   #stored the gse id in a variable and we got the phenodata using the getGEO function
phenoData <- pData(phenoData(gse[[1]]))   
phenoDataOP <- phenoData[,c(1,2,55,57, 58,59,60:63, 65:67)]    

#got the code from gpt for adding the col as a vector
# Create a named vector with sample IDs and severity group labels (from the image you provided)

severity_groups <- c(
  "01-0036" = 1, "01-0037" = 1, "01-0048" = 1, "01-0045" = 1, "01-0044" = 1,
  "01-0043" = 1, "01-0042" = 1, "01-0041" = 1, "01-0035" = 1, "01-0034" = 1,
  "01-0033" = 1, "01-0030" = 1, "01-0029" = 1, "01-0027" = 1, "01-0026" = 1,
  "01-0025" = 1, "01-0024" = 1, "01-0023" = 1, "01-0022-1" = 1,
  
  "01-0049" = 2, "01-0047" = 2, "01-0046" = 2, "01-0032" = 2, "01-0021" = 2,
  
  "32-0010" = 3, "32-0011" = 3, "32-0012" = 3, "32-0013" = 3, "32-0014" = 3,
  "32-0015" = 3, "32-0016" = 3, "32-0017" = 3, "32-0018" = 3, "32-0019" = 3,
  "01-0038" = 3,
  
  "32-0002" = 4, "32-0003" = 4, "32-0004" = 4, "32-0005" = 4, "32-0006" = 4,
  "32-0007" = 4, "32-0008" = 4, "32-0009" = 4
)
phenoDataOP <- phenoDataOP %>%
  mutate(RNA_Seq_severity_group = severity_groups[title])
head(phenoDataOP)
colnames(phenoDataOP)



library(dplyr)

dataOP <- dataOP %>%
  gather(key = "samples", value = "counts", -X ) %>%
  mutate(samples= gsub('\\.', '-', samples), samples = gsub('^X', '', samples)) %>%
  inner_join(phenoDataOP, by = c("samples" = "title")) %>%
  dplyr::select(1,3,4)  %>%
  spread(key = "geo_accession", value = "counts") %>%
  column_to_rownames(var = "X")

dataOPsubset <- dataOP[rowSums(dataOP >= 0.5) >= 22, ]
log_dataOPsubset <- log2(dataOPsubset + 1)
samples.to.be.excluded <- c('GSM3186678', 'GSM3186679', 'GSM3186673', 'GSM3186676', 'GSM3186709')
log_dataOPsubset_clean <- log_dataOPsubset[, !(colnames(log_dataOPsubset) %in% samples.to.be.excluded)]
log_dataOPsubset_cleanOP <- log_dataOPsubset_clean[!colnames(log_dataOPsubset_clean) %in% ("GSM3186698")]


genes_to_test <- rownames(log_dataOPsubset_cleanOP)
gene_ids_clean <- gsub("\\.\\d+$", "", genes_to_test)   #removed the versions from the ensembl ids

#converting the ensembl ids to gene symbols:
# Load required libraries
library(AnnotationDbi)
library(org.Hs.eg.db)

# Map Ensembl IDs to gene symbols using the mapIDs function
gene_symbols <- mapIds(org.Hs.eg.db, 
                       keys = gene_ids_clean, 
                       column = "SYMBOL", 
                       keytype = "ENSEMBL", 
                       multiVals = "first")

# gene_symbols is a named vector with Ensembl IDs as names and gene symbols as values
print(gene_symbols)
gene_symbols_list <- as.character(gene_symbols)  #we got the dictionery with ensembl id n gene symbol for each list from the mapIDs function
# we converted into a list containing only the gene symbols

print(gene_symbols_list)

# Now 'gene_ids_clean' will contain:
# [1] "ENSG00000123456" "ENSG00000234567"

# You can now safely use gene_ids_clean for mapping, enrichment, or conversion.

#now beginning with the code for GO enrichment analysis


GO_results <- enrichGO(gene = gene_symbols_list, OrgDb = "org.Hs.eg.db", keyType = "SYMBOL", ont = "BP")
as.data.frame(GO_results) #converting the results to a df to view, used enrichGO funtion to perform GO enrichment analysis
# as inputs, we gave it the gene symbols list from our cleaned n preprocessed dataset, the lib to get entire human species genome so that it can compare, key type is symbol because
#we provided the gene symbols, if it would have been ensembl ids, then it would have been ensembl
#ont = bp (biological process) because we are interested in that only here because that helps us the most in making interpretations.
library(enrichplot)



# Print the barplot (ggplot object)
print(barplot(GO_results, showCategory = 15))



#now the output of the enrichGO function, ie the GO_results, no of rows gives us the no of significantly enriched GO terms (it describes something abt the gene's fn)
# could be biological process, molecular fn or the cellular component it is involved in.
#it shows the statistically sig GO terms, which means they occur more freq than expected n their occurence is not random

head(GO_results)
print(GO_results)
# now we are checking the GO enrichment terms for our top 25 hub genes
#starting with mmp14
gene_of_interest <- "MMP14"

# Find rows where this gene appears in the GO term gene list
rows_with_gene <- GO_results[grep(gene_of_interest, GO_results$geneID), ]

# View GO terms associated with your gene
rows_with_gene[, c("Description")]

#here we get a list of biological processes
#so if i check for the description of a gene in Gene ontology results, n it gives multiple biological processes, 
#then it means that the group of genes which might be associated with each of those processes, the gene which 
#i am checking is present in the groups of those genes which are associated with the processes which are displayed in the output list above

save.image("GO_analysis_115650_001_fpkm_workspace.RData")




