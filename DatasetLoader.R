# --- Script 2: Read data ---

# Note: "Processed Data for Clustering" contains the following subdirectories. The files
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

# Define a subset of genes on which to run the R^2 analysis
subsetSize <- 250
genesSubset <- sample(geneNames, size=subsetSize)
