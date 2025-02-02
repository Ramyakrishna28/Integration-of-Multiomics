BiocManager::install(version = '3.20')

# Install the 'multiGSEA' package
BiocManager::install("multiGSEA")
BiocManager::install("multiGSEA", force = TRUE)

library(multiGSEA)
# load human organism library

BiocManager::install("org.Hs.eg.db")
library("org.Hs.eg.db")

# To run the analysis of this vignette, load the installed version of multiGSEA.
library(multiGSEA)
library(magrittr)

# load transcriptomic features
data(transcriptome)

# load proteomic features
data(proteome)

# load metabolomic features
data(metabolome)

head(transcriptome)  # Display the first few rows of transcriptome data
head(proteome)       # Display the first few rows of proteome data
head(metabolome)     # Display the first few rows of metabolome data

2.Load Data into multiGSEA
# create data structure
omics_data <- initOmicsDataStructure(layer = c(
  "transcriptome",
  "proteome",
  "metabolome"
))

## add transcriptome layer
omics_data$transcriptome <- rankFeatures(
  transcriptome$logFC,
  transcriptome$pValue
)
names(omics_data$transcriptome) <- transcriptome$Symbol

## add proteome layer
omics_data$proteome <- rankFeatures(proteome$logFC, proteome$pValue)
names(omics_data$proteome) <- proteome$Symbol

## add metabolome layer
## HMDB features have to be updated to the new HMDB format
omics_data$metabolome <- rankFeatures(metabolome$logFC, metabolome$pValue)
names(omics_data$metabolome) <- metabolome$HMDB
names(omics_data$metabolome) <- gsub(
  "HMDB", "HMDB00",
  names(omics_data$metabolome)
)

# prediction of the ls values for each omics layer 
omics_short <- lapply(names(omics_data), function(name) {
  head(omics_data[[name]])
})
names(omics_short) <- names(omics_data)
Omics_short


MultiGSEA Data Structure:

The input data for multiGSEA must be in the form of a nested list, where each sublist corresponds to a specific omics layer.

omics_data <- list( transcriptome = data.frame(Symbol = c("Gene1", "Gene2", "Gene3"), logFC = c(1.2, -0.5, 0.8), pValue = c(0.001, 0.05, 0.02)), 
proteome = data.frame(Symbol = c("Protein1", "Protein2", "Protein3"), logFC = c(0.9, -0.3, 1.5), pValue = c(0.01, 0.04, 0.001)), 
metabolome = data.frame(Symbol = c("Metabolite1", "Metabolite2", "Metabolite3"), logFC = c(1.1, -0.2, 0.5), pValue = c(0.05, 0.1, 0.03)) )

3. Gene Set Enrichment Analysis (GSEA):
# Download and customize pathway definitions
databases <- c("kegg", "reactome")
layers <- names(omics_data)

pathways <- getMultiOmicsFeatures(
  dbs = databases, layer = layers,
  returnTranscriptome = "SYMBOL",
  returnProteome = "SYMBOL",
  returnMetabolome = "HMDB",
  useLocal = FALSE
)

pathways_short <- lapply(names(pathways), function(name) {
  head(pathways[[name]], 2)
})
names(pathways_short) <- names(pathways)
pathways_short

4.Combined multiomics enrichment : 
measure pathway responses more comprehensively by aggregating p-values from different omics layers 

# use the multiGSEA function to calculate the enrichment scores
# for all omics layer at once.
enrichment_scores <- multiGSEA(pathways, omics_data)

df <- extractPvalues(
  enrichmentScores = enrichment_scores,
  pathwayNames = names(pathways[[1]])
)

df$combined_pval <- combinePvalues(df)
df$combined_padj <- p.adjust(df$combined_pval, method = "BH")

df <- cbind(data.frame(pathway = names(pathways[[1]])), df)

head(df)
