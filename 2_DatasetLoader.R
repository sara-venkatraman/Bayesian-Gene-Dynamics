# --- Script 2: Read data ---

# "Processed Data for Clustering" contains the following subdirectories. The files
# in each subdirectory have the same names.
# - "/Differentially-Expressed/" only contains gene data of differentially-expressed genes
# - "/Neighbors (STRING > 900)/" only contains gene data for the "neighbors" (identified 
#   by STRING-DB) of the differentially-expressed genes, with composite STRING scores >= 900
# - "/Combined Genes/" contains the union of the genes in the previous two directories.

# First, set the desired subdirectory:
subdirectory <- "Combined Genes"

# List of names and Flybase IDs of differentially-expressed genes
geneNames <- read.csv(paste("../Processed Data for Clustering/", subdirectory, "/GeneNames.csv", sep=""), header=T)
geneIDs <- read.csv(paste("../Processed Data for Clustering/", subdirectory, "/GeneIDs.csv", sep=""), header=T)  

# Dataframe of expression measurements (in normalized counts) from two replicates  
normCountsRep1 <- read.csv(paste("../Processed Data for Clustering/", subdirectory, "/NormCountsRep1.csv", sep=""), row.names=1)
normCountsRep2 <- read.csv(paste("../Processed Data for Clustering/", subdirectory, "/NormCountsRep2.csv", sep=""), row.names=1)
geneData <- read.csv(paste("../Processed Data for Clustering/", subdirectory, "/LogChange.csv", sep=""), row.names=1)

# Rename/reformat dataframes and vectors to use with functions in EmpiricalBayesFunctions.R
geneNames <- as.character(geneNames$geneName)
geneIDs <- as.character(geneIDs$geneID)

# --- Subset selection: adjust as needed depending on choice of subdirectory ---
# Define a subset of genes on which to run the R^2 analysis: 60% of the subset is 
# differentially-expressed, the other 40% is neighbors
DEgeneNames <- read.csv("../Processed Data for Clustering/Differentially-Expressed/GeneNames.csv")
neighborGeneNames <- read.csv("../Processed Data for Clustering/Neighbors (STRING > 900)/GeneNames.csv")

subsetSize <- 200
DEsubset <- as.character(sample(DEgeneNames$geneName, size=0.6*subsetSize))
neighborSubset <- as.character(sample(neighborGeneNames$geneName, size=0.4*subsetSize))
genesSubset <- unique(c(DEsubset, neighborSubset)) # There should be no duplicates anyway
