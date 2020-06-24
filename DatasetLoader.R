# --- Script 2: Read data ---

# List of names and Flybase IDs of differentially-expressed genes
geneNames <- read.csv("../Processed Data for Clustering/DEgeneNames.csv", header=T)
geneIDs <- read.csv("../Processed Data for Clustering/DEgeneIDs.csv", header=T)  

# Dataframe of expression measurements (in normalized counts) from two replicates  
normCountsRep1 <- read.csv("../Processed Data for Clustering/DEnormCountsRep1.csv", row.names=1)
normCountsRep2 <- read.csv("../Processed Data for Clustering/DEnormCountsRep2.csv", row.names=1)
geneData <- read.csv("../Processed Data for Clustering/DElogChange.csv", row.names=1)

# Rename/reformat dataframes and vectors to use with functions in EmpiricalBayesFunctions.R
geneNames <- as.character(geneNames$geneName)
geneIDs <- as.character(geneIDs$geneID)

# Define a subset of genes on which to run the R^2 analysis
subsetSize <- 250
genesSubset <- sample(geneNames, size=subsetSize)
