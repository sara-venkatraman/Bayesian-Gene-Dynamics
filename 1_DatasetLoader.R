# --- Script 1: Read data and generate a subset of genes ---

# "Processed Data for Clustering" contains the following subdirectories. The files
# in each subdirectory have the same names.
# - "/Differentially-Expressed/" contains gene expression data for differentially-expressed 
#   genes identified in Schlamp et al., 2020
# - "/Neighbors (STRING > 950)/" contains gene expression data for the "neighbors" of the 
#   differentially-expressed genes, i.e. genes that have a composite STRING score of at least
#   950 with at least one differentially-expressed gene
# - "/Combined Genes/" contains the union of the genes in the previous two directories.

# Set the desired subdirectory:
subdirectory <- "Differentially-Expressed"

# Read gene names and Flybase IDs and format them as character vectors
geneNames <- read.csv(paste("../Processed Data/", subdirectory, "/GeneNames.csv", sep=""), header=T)
geneIDs <- read.csv(paste("../Processed Data/", subdirectory, "/GeneIDs.csv", sep=""), header=T)  
geneNames <- as.character(geneNames[,1])
geneIDs <- as.character(geneIDs[,1])

# Dataframe of expression measurements (as log_2 fold change) 
geneData <- read.csv(paste("../Processed Data/", subdirectory, "/LogChange.csv", sep=""), row.names=1)

# Define the hours corresponding to each time point in the dataset
hours <- c(0, 1, 2, 4, 5, 6, 8, 10, 12, 14, 16, 20, 24, 30, 36, 42, 48)

# (If desired) dataframes of normalized counts from two replicates
# normCountsRep1 <- read.csv(paste("../Processed Data/", subdirectory, "/NormCountsRep1.csv", sep=""), row.names=1)
# normCountsRep2 <- read.csv(paste("../Processed Data/", subdirectory, "/NormCountsRep2.csv", sep=""), row.names=1)

# Read prior matrix
priorMatrix <- read.csv(paste("../Processed Data/", subdirectory, "/PriorMatrix.csv", sep=""), row.names=1)
rownames(priorMatrix) <- geneNames;  colnames(priorMatrix) <- geneNames

# --- Subset selection ---

# Set subset size
subsetSize <- 150

# If using the combined set of genes (differentially-expressed and neighbors), set
# the desired proportion of each set to use
DEsubsetProp <- 0.6
neighborSubsetProp <- 1 - DEsubsetProp

if(subdirectory == "Combined Genes") {
  DEgeneNames <- read.csv("../Processed Data/Differentially-Expressed/GeneNames.csv")[,1]
  neighborGeneNames <- read.csv("../Processed Data/Neighbors (STRING > 950)/GeneNames.csv")[,1]
  DEsubset <- as.character(sample(DEgeneNames, size=DEsubsetProp*subsetSize))
  neighborSubset <- as.character(sample(neighborGeneNames, size=neighborSubsetProp*subsetSize))
  geneSubset <- unique(c(DEsubset, neighborSubset)) # There should be no duplicates anyway
} else {
  geneSubset <- as.character(sample(geneNames, size=subsetSize))
}
