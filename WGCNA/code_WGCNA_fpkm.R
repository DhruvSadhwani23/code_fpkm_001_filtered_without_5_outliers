################### LOADING DATASET AND NECESSARY LIBRARIES
setwd("/home/tanuja/Documents/srip2025_dhruv_fshd_wgcna/srip2025_dhruv_fshd_wgcna_115650_001_fpkm")
getwd()
dataOP <- read.csv("~/Documents/srip2025_dhruv_fshd_wgcna/srip2025_dhruv_fshd_wgcna_115650_001_fpkm/srip2025_dhruv_fshd_wgcna_115650_001_dataset.csv", header=T, stringsAsFactors = FALSE)
library(WGCNA)  #we did that by writing necessary linux code in terminal
library(DESeq2)
library(GEOquery)
library(tidyverse)
library(gridExtra)
library(CorLevelPlot) # we loaded the necessary libraries
# did same stuff, set the working directory, loaded the csv file and the libraries 

load("/mnt/808e843c-fb08-404a-9f45-efd817d81a77/srip2025_dhruv_fshd_wgcna/srip2025_dhruv_fshd_wgcna_115650_001_fpkm/srip2025_dhruv_fshd_wgcna_115650_001_WGCNA_analysis/WGCNA_115650_001_fpkm_workspace.RData")


dataOP[1:10,1:5]  #just seeing the dataset once
colnames(dataOP)
# we'll see if there is any problem with the col names as they are non uniform, n if so,we will modify them further..


# WE WILL BE FILTERING LOW EXPRESSED GENES USING 2 METHODS: (PART OF DATA CLEANING)
#First is using the goodSampleGenes function (used for like a final cleanup of data) in wgcna used for filtering out missing or problematic genes in our expression matrix
# like missing values(NAs) in genes or samples, const genes/samples and low expressed genes
# another method we'll be using is manually filtering out the less expressed genes (helps to remove noise)
# and in wgcna, an expression matrix is simply a numeric matrix that contains rows as genes or transcripts and samples as columns
#with the entries as expression values, which could be raw counts, normalised counts (like the vst we applied to the prev dataset), or the FPKM/TPM/RPKM values(which are again normalisation methods)
# raw counts mean the reads of that particular gene were present in our RNA seq data, which represent how high was a gene expressed in that particular sample
#so yes, we have the expression matrix here with normalised fpkm values and we will use both the methods to filter out less expressed n other problematic genes 
#then we'll be using hierarchical clustering to remove outlier genes/samples from our data 

#this is our entire process for data cleaning and preprocessing
#after doing all of this, we will begin with the WGCNA steps, 
#like how we did in the tutorial + taking a ref of what they have done in the paper + help from chatgpt regarding the pipeline's structure
# we will also log transform this data before performing wgcna analysis.

#########################################DATA PREPROCESSING:

# from this we see that this is the counts for the first sample and the first col includes the gene ids, so we would want to set it as rownames, so they 
# would only be as row identifiers and not part of the actual expression matrix
# after this,R keeps the row names internally (they don’t appear as a column), but you can view them(code below)
#hence we are converting the first column to rownames, so that R knows this is not data and just identifiers for diff genes
#else they could also interfere with the operations we apply to columns which has the values.

dataOP[1,1] # see, the first col contains ensembl ids
head(dataOP)
rownames(dataOP)
colnames(dataOP)

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

# also now we will add the relevant columns to our data from the pheodata as per our understanding of the tutorial and
# taking reference with what they have done in the paper
# added columns title(to merge it with dataOP as it contains sample ids,geo accession ids to get proper geo ids of samples)
#rem columns are related to traits like active, subject status, fat fraction , pathological level etc..
# the whole list: [1] "title"              "geo_accession"      "active:ch1"         "fat.fraction:ch1"   "inflam:ch1"        
#[6] "path.score:ch1"     "rnaseq.rank:ch1"    "rnaseq.score:ch1"   "stir:ch1"           "subject status:ch1"
#[11] "t1:ch1"            
#they have divided into groups of 4 on the basisof sum of rna expression, as 1,2,3,4 with increasing levels, 
# so Ill also be manually adding that row into my traits columns for further wgcna analysis.

#merging phenodata and our dataset to get a clean matrix in wgcna format ready for analysis: genes as rows n proper sample names as columns (from GEO, as GSM ids)

library(dplyr)

dataOP <- dataOP %>%
  gather(key = "samples", value = "counts", -X ) %>%
  mutate(samples= gsub('\\.', '-', samples), samples = gsub('^X', '', samples)) %>%
  inner_join(phenoDataOP, by = c("samples" = "title")) %>%
  dplyr::select(1,3,4)  %>%
  spread(key = "geo_accession", value = "counts") %>%
  column_to_rownames(var = "X")
head(dataOP)
colnames(dataOP)
rownames(dataOP)  

#keeping all the entries numeric for the operations we need to perform ahead( lets say filtering for example),Sometimes when you read the data from a file, it might import some or all of the columns as factors or characters, especially if there are mixed types or formatting issues.
#so need to account for those by converting all of such into a matrix of numeric values.
#now here we will convert it into numeric, we were doing it earlier, hence it was erasing the ensembl ids, so 
#first we converted them into rownames n then they are not present as a column, just R knows that they are as identifiers for diff genes


dataOP[] <- sapply(dataOP, as.numeric)  # although this wasnt needed here as we saw from str(dataOP) r was already reading all of the values as numeric but here its just for our understanding.


# DATA CLEANING + LOG TRANSFORMATION
#removing low expressed genes manually first
# we are keeping only those genes that have fpkm>1 in atleast half of the samples , hence filtering out the less expressed genes
# then goodSampleGenes function n  clustering as discussed above to remove ouliers

dataOPsubset <- dataOP[rowSums(dataOP >= 0.5) >= 22, ]   #total number of samples are 43, keeping those genes with >0.5 value in
#more than half of the samples , a common threshold
nrow(dataOPsubset)  # left with 10310 genes
head(dataOPsubset)


#log tranformation, stabilizes variance n keeps data normalised from the highly skewed fpkm values(they are in general)
#normally distributed data is an assumption of wgcna and also for calculating peasron coeff, which is again used for clustering.
# also without log tranformation, a single highly expressed gene can dominate distance calculations, after the tranform, highly expressed n lowly expressed genes contribute more equally to similarity measures. 
#also for samplegenes function , it is extremely sensitivie to extreme values, hence it is also applied after log transform to ensure sample n gene variance is more stable.
#hence the order is 1) Log tranform, 2) goodsamplegenes fn 3) remove outliers using clustering

##LOG TRANSFORM:
log_dataOPsubset <- log2(dataOPsubset + 1)

#goodSampleGenes function:

gsg <- goodSamplesGenes(t(log_dataOPsubset))
summary(gsg)  # no bad genes

#now using PCA and hclust to remove outlier samples
#HCLUST

sampleTree <- hclust(dist(t(log_dataOPsubset)), method = "average")
plot(sampleTree, main = "Sample clustering", xlab = "", sub = "")
abline(h = 50, col = "red")

# we did check with a couple of values, we'll check with pca also once 
library(ggplot2)
#PCA
# Perform PCA on transposed expression matrix
pca <- prcomp(t(log_dataOPsubset))

# Extract PCA coordinates and calculate percentage variance explained
pca.dat <- as.data.frame(pca$x)
pca.var <- pca$sdev^2
pca.var.percent <- (pca.var / sum(pca.var)) * 100

# Add sample names as a column
pca.dat$Sample <- rownames(pca.dat)

# Simplify status labels from phenodata
# Assumes phenodata has rownames as sample names and a column "subject status:ch1"
pca.dat$Status <- ifelse(
  phenoDataOP[rownames(pca.dat), "subject status:ch1"] == "Control",
  "Control",
  "FSHD"
)

# PCA plot with colored groups
ggplot(pca.dat, aes(x = PC1, y = PC2, color = Status)) +
  geom_point(size = 3) +
  geom_text(aes(label = Sample), vjust = -1.2, size = 3) +
  labs(
    x = paste0('PC1: ', round(pca.var.percent[1], 2), '%'),
    y = paste0('PC2: ', round(pca.var.percent[2], 2), '%'),
    color = 'Condition'
  ) +
  theme_minimal()

# apart from checking visually, we also checked using the z score method, which gave no outliers
pc1_scores <- pca$x[,1]           # PC1 values for each sample
pc1_z <- scale(pc1_scores)        # Convert to z-scores (mean 0, sd 1)

outliers_pca <- names(pc1_z[abs(pc1_z) > 2.5])  # Threshold: z > ±2.5
outliers_pca

#NOW EXCLUDING OUTLIER SAMPLES
#THOSE ARE VISUALLY SEEN AND COMMON IN BOTH
samples.to.be.excluded <- c('GSM3186678', 'GSM3186679', 'GSM3186673', 'GSM3186676', 'GSM3186709')
log_dataOPsubset_clean <- log_dataOPsubset[, !(colnames(log_dataOPsubset) %in% samples.to.be.excluded)]


#from hclust:  673,709,676,677,678,679 (outliers seen)
#from PCA: 678,676,709, 673,679  (")

head(log_dataOPsubset_clean)

##########    WGCNA CORE CODE  #######


#network construction, firstly selecting a set of soft threshold powers
power <- c(c(1:10), seq(from = 12, to = 50, by = 2 ))
sft <- pickSoftThreshold(t(log_dataOPsubset_clean), powerVector = power,networkType = "signed", verbose = 5)
sft.data <- sft$fitIndices

#visualisation to pick correct power

a1 <- ggplot(sft.data, aes(Power,SFT.R.sq, label= power)) + geom_point() +geom_text(nudge_y = 0.1) + geom_hline(yintercept = 0.8, color="red") + labs(x = 'Power', y = 'scaled free topology model fit, signedR^2') + theme_classic()

a2 <- ggplot(sft.data, aes(Power, mean.k., label= Power)) + geom_point() + geom_text(nudge_y = 0.1) + labs(x='Power', y='Mean connectivity') + theme_classic()
grid.arrange(a1,a2,nrow=2)
print(a1)
print(a2)


#now we'll convert the matrix to numeric
log_dataOPsubset_clean <- t(log_dataOPsubset_clean)
# here log_dataOPsubset initially had genes as rows and samples as columns (genes x samples form), we have converted it into its transposed matrix
#as for the core wgcna code and its functions like picksoftThreshold, blockwisemodules requires genes as columns and samples as rows(samples x genes)
log_dataOPsubset_clean[] <- sapply(log_dataOPsubset_clean, as.numeric)
soft_power <- 20
temp_cor <- cor
cor <- WGCNA::cor
#writing the fn to identify gene modules:

bwnet <- blockwiseModules(log_dataOPsubset_clean, maxBlockSize = 14000, TOMType = "signed",power = soft_power, mergeCutHeight = 0.25, numericLabels = FALSE, randomSeed = 1234, verbose = 3 )   # most imp wgcna code

cor <- temp_cor



############ FINDING MODULE EIGENGENES:

module_eigengenes <- bwnet$MEs
head(module_eigengenes)
table(bwnet$colors)
plotDendroAndColors(bwnet$dendrograms[[1]], cbind(bwnet$unmergedColors, bwnet$colors),
                    c("unmerged", "merged"),
                    dendroLabels = FALSE,
                    addGuide = TRUE,
                    hang= 0.03,
                    guideHang = 0.05)

#after this, we will calculate the module trait relationships and visualize them as a heatmap
colnames(phenoDataOP)
phenoDataOP <- phenoDataOP %>%
  mutate(
    gender_binary = ifelse(`gender:ch1` == "male", 1, 0),
    status_binary = ifelse(`subject status:ch1` == "Control", 0, 1)
  )

colnames(phenoDataOP)


traits <- phenoDataOP[, !colnames(phenoDataOP) %in% c("gender:ch1", "title", "geo_accession", "subject status:ch1", "batch:ch1")]  #just keeping the numeric columns in our traits matrix for finding correlation between module n traits
nrow(traits)
nrow(log_dataOPsubset_clean)
#here we have forgot to drop the outlier samples from the phenodataOP and also the traits so we are getting the error that as no of rows, we are unable to merge the 2 datasets and hence find the correlation
#dropping them from both traits n phenodataOP


phenoDataOP <- phenoDataOP[!rownames(phenoDataOP) %in% (samples.to.be.excluded), ]
traits <- traits[!rownames(traits) %in% (samples.to.be.excluded), ]
nrow(traits)

#there are some NA values in traits (also the phenodataOP), due to which we are unable to calculate the correlation
#those are mainly due to thw missing values of fshd characteristics in mri samples , n few missing values as well, we'll deal with that.
#dropping those in control(making them 0 ofc) and remooving the one sample (fshd) with missing values
#also we dont know the gender for control so removing it from our traits for correlation analysis

traits <- traits[!rownames(traits) %in% ("GSM3186698"), ]
phenoDataOP <- phenoDataOP[!rownames(traits) %in% ("GSM3186698"), ]
traits[is.na(traits)] <- 0
#also dropping that subset from log_dataOPsubset_clean and also module_eigengenes
module_eigengenes <- module_eigengenes[!rownames(module_eigengenes) %in% ("GSM3186698"), ]
log_dataOPsubset_clean <- log_dataOPsubset_clean[!rownames(log_dataOPsubset_clean) %in% ("GSM3186698"), ]

#calculating the module trait relationships

nSamples <- nrow(log_dataOPsubset_clean)
nGenes <- ncol(log_dataOPsubset_clean)
str(traits)
traits[] <- lapply(traits, function(x) as.numeric(as.character(x)))  #converting the character values to numeric, they are numbers but internally they were stored as characters 
#which gave us an error while calculating correlation, we fixed that.
str(traits)   #checking now
module.trait.corr <- cor(module_eigengenes, traits, use = 'p')
module.trait.corr.pvals <- corPvalueStudent(module.trait.corr, nSamples)
str(traits)
traits[] <- lapply(traits, function(x) as.numeric(as.character(x)))  #converting the character values to numeric, they are numbers but internally they were stored as characters 
#which gave us an error while calculating correlation, we fixed that.
str(traits)   #checking now


#now lets visualize the module trait relationship as a heatmap

heatmap.data <- merge(module_eigengenes, traits, by = 'row.names')
head(heatmap.data)
colnames(heatmap.data)
heatmap.data <- heatmap.data[, !colnames(heatmap.data) %in% c("gender_binary")]  
# just removed the gender col as we didnt have the gender values for the control samples
heatmap.data <- heatmap.data %>% 
  column_to_rownames(var = 'Row.names')
colnames(heatmap.data)[colnames(heatmap.data) %in% c("active:ch1", "fat.fraction:ch1", "inflam:ch1", "path.score:ch1", "rnaseq.rank:ch1","rnaseq.score:ch1","stir:ch1","t1:ch1","RNA_Seq_severity_group","status_binary")] <- c("Active", "Fat Fraction", "Inflammation", "Pathology Score", "RNA Seq rank", "RNA Seq Score", "STIR Signal", "T1 signal", "RNA Seq Severity Group", "Disease State")



##png(file="/home/tanuja/Documents/srip2025_dhruv_fshd_wgcna/srip2025_dhruv_fshd_wgcna_115650_001_fpkm.png", width=14,height=8)
CorLevelPlot(heatmap.data,
             y = names(heatmap.data)[1:10],
             x = names(heatmap.data)[11:20],
             titleX = "Trait Data",
             titleY = "Modules",
             rotTitleY = 90,
            
             rotLabX = 45,
             signifSymbols = c("***", "**", "*", ""),
             signifCutpoints = c(0, 0.001, 0.01, 0.05, 1),
             col = c("blue1", "skyblue", "white", "pink", "red"))

##dev.off()

module.gene.mapping <- as.data.frame(bwnet$colors)   #extracting the genes for the module turqoise, as we can see that it is the most relevant module
module.gene.mapping %>% 
  filter(`bwnet$colors` == 'turquoise') %>% 
  rownames()

#just removing that gender binary column from traits, as its incorrect because the NAs columns were also converted into zero, which is for female
#we were plotting the heatmap using heatmap.data n in that we removed the row before plotting n for the earlier part where we 
#calculated the module trait correaltions and their p values using module eigengenes and traits, that part of the code we ran again after this 
#so basically theres no problem changing it here later.
traits <- traits[, !colnames(traits) %in% c("gender_binary")]  


# 6B. Intramodular analysis: Identifying driver genes
# Calculate the module membership and the associated p-values
# The module membership/intramodular connectivity is calculated as the correlation of the module eigengenes and the gene expression profile. 
# This quantifies the similarity of all genes on the array to every module.

module.membership.measure <- cor(module_eigengenes, log_dataOPsubset_clean, use = 'p')
module.membership.measure.pvals <- corPvalueStudent(module.membership.measure, nSamples)
head(module.membership.measure.pvals)

# Calculate the gene significance and associated p-values
#our top most trait of interest being ofc whether the sample is affected or unaffected 

gene.signf.corr <- cor(log_dataOPsubset_clean, traits$status_binary, use = 'p')
gene.signf.corr.pvals <- corPvalueStudent(gene.signf.corr, nSamples)
gene.signf.corr.pvals %>% 
  as.data.frame() %>% 
  arrange(V1) %>% 
  head(25)
# Using the gene significance you can identify genes that have a high significance for trait of interest 
# Using the module membership measures you can identify genes with high module membership in interesting modules.
colnames(heatmap.data)

#so what we'll do is that for our topmost trait of interest, ie status_binary (affected/ unaffected), from the heatmap
#we can see that the module turquoise has highest correlation for that trait
#so now we'll plot a scatterplot of the gene significance v/s the module membership for the turquoise module, where the points would be representing the genes in that module
#we'll see if theres any correlation between (the significant genes for our trait of interest and the genes having
#high module membership within that particular module(which are basically the genes in a particular module that have high correlation with the eigengene of that module))

#we'll also plot a barplot for each module with its no of genes so that its easy to visualise them (at the end after the scatterplot)

library(ggplot2)
library(dplyr)

# Prepare data
genes_per_module <- as.data.frame(table(bwnet$colors)) %>%
  rename(Module = Var1, Count = Freq) %>%
  arrange(desc(Count))

# Get WGCNA-style module colors (ensures color names match)
module_colors <- sort(unique(bwnet$colors))

# Ensure correct names for the scale
fill_colors <- setNames(module_colors, module_colors)

# Bar plot
module_plot <- ggplot(genes_per_module, aes(x = reorder(Module, -Count), y = Count, fill = Module)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = fill_colors) +  # This keeps WGCNA module colors
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(hjust = 0.5)
  ) +
  labs(x = "Modules", y = "Number of Genes") +
  geom_text(aes(label = Count), vjust = -0.5, size = 5) +
  guides(fill = "none")

# Show plot
print(module_plot)

#firstly, we will extract the genes corresponding to the module turquoise, which is most correlated with our gene of interest
#we did this earlier as well, from bwnet object extracting colors and filtering out turquoise genes.
rownames(module.membership.measure)

# Step 0: Transpose if needed (only run this once if matrices were in wrong orientation)
module.membership.measure <- t(module.membership.measure)
# gene.signf.corr is already fine, so don't transpose it

# Step 1: Extract genes belonging to the each module of interest
red_genes <- module.gene.mapping %>%
  filter(`bwnet$colors` == 'red') %>%
  rownames()

# Step 2: Plot Module Membership vs Gene Significance 

verboseScatterplot(
  abs(module.membership.measure[red_genes, "MEred"]),
  abs(gene.signf.corr[red_genes, "MEred" ]),
  xlab = "Module Membership in red module",
  ylab = "Gene significance for disease state",
  main = "Module membership vs. gene significance\n",
  cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2,
  col = "red"
)

#syntax explanation:  


#hence the plot clearly shows that the genes of significance and genes with high module membership are highly correlated which shows that 
#the genes high associated with a trait are often also the most important elements of modules associated with that trait.



#NOW WE CAN SAY THAT THE ACTUAL NETWORK ANALYSIS BEGINS
#FROM THIS, WE'LL GET THE DATA FOR THE MODULES/(GENE NETWORKS N SUBNETWORKS) WHICH WE WILL BE EXTRACTING HERE
#THEN WE'LL PERFORM DOWNSTREAM ANALYSIS ON THE DATA TO EXPLORE OTHER GENES WHICH COULD BE USED TO FIND OTHER GENES 
#WHICH ARE REGULATED BY DUX4 OR INVOLVED IN SOME COMMON PATHWAY/ CO-EXPRESSED WITH DUX4

#EXTRACTING MODULES/ NETWORKS(OR SUBNETWORKS OF GENES) OF INTEREST: which we have got taking reference from the paper n the our traits of interest
#then have extracted those modules which have high correlation with those traits:
table(bwnet$colors)


#our traits of interest: status_binary
#those referred in the paper:  (dux4 expression levels)rnaseq.score:ch1 (turquoise), RNA_Seq_severity_group(yellow), (stir1---turquoise, t1----turquoise, pathology----turquoise, fat----turquoise)---these are the mri resuts found to be a reliable predictor of dux4 target expression as well as muscle pathology
print("Calculating TOM similarity matrix...")
TOM = TOMsimilarityFromExpr((log_dataOPsubset_clean), power = 20)
print("TOM calculation completed.\n")

# Output base directory
baseDir <-  "."
# Get all module colors
allModules = unique(bwnet$colors)
print(paste("Modules detected:", paste(allModules, collapse = ", ")))

# Step 2: Calculate module eigengenes
moduleColors <- bwnet$colors
print("Calculating module eigengenes...")
MEs = moduleEigengenes((log_dataOPsubset_clean),moduleColors)$eigengenes
print("Module eigengenes calculated.\n")

# Step 3: Compute kME values (correlation of gene expression with module eigengenes)
print("Calculating kME (module membership) for all genes...")
kMEall = as.data.frame(cor((log_dataOPsubset_clean), MEs, use = "p"))
print("kME calculation completed.\n")


# Step 4: Loop over modules with high correlation with trait of interest
#we are doing for all the modules(extracting the networks), we'll analyse later

for (module in allModules) {
  
  print(paste0("Processing module: ", module, "..."))
  
  # Genes in the module (original IDs with possible versions)
  inModule <- moduleColors == module
  modProbes <- rownames(t(log_dataOPsubset_clean))[inModule]  # original IDs with version suffixes
  print(modProbes)
  # Strip version suffix from gene IDs (e.g. ENSG00000123456.12 -> ENSG00000123456)
  clean_probes <- sub("\\..*$", "", modProbes)
  
  # Subset TOM matrix (genes in module)
  modTOM <- TOM[inModule, inModule]
  
  # Create output directory for module
  moduleDir <- file.path(baseDir, module)
  if (!dir.exists(moduleDir)) {
    dir.create(moduleDir, recursive = TRUE)
    print(paste("Created directory:", moduleDir))
  }
  
  # Extract kME values (module membership) for current module genes
  # kMEall rownames should be original modProbes (with versions)
  kMEcol <- paste0("ME", module)
  if (!(kMEcol %in% colnames(kMEall))) {
    warning(paste("Warning: kME column", kMEcol, "not found. Skipping module."))
    next
  }
  
  # Extract kME using original modProbes as rownames
  modkME <- kMEall[modProbes, kMEcol, drop = FALSE]  # ensure result is a vector
  
  # Export network to Cytoscape
  print(paste("Exporting Cytoscape files for module:", module))
  exportNetworkToCytoscape(modTOM,
                           edgeFile = file.path(moduleDir, paste0("CytoscapeInput-edges-", module, ".txt")),
                           nodeFile = file.path(moduleDir, paste0("CytoscapeInput-nodes-", module, ".txt")),
                           weighted = TRUE,
                           threshold = 0.1,
                           nodeNames = modProbes,
                           
                           nodeAttr = data.frame(module = module, kME = modkME)
  )
  
  print(paste("Export completed for module:", module, "\n"))
}

print("Wohoooooooo! All modules processed successfully.")



save.image("WGCNA_115650_001_fpkm_workspace.RData")

