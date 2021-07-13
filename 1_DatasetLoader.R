# --- Script 1: Read time-course gene expression data from CSV ---

# Read dataframe of expression measurements (as log_2 fold change) 
geneData <- read.csv("Data/geneData.csv", row.names=1)

# Save the list of gene names
geneNames <- rownames(geneData)

# Define the hours corresponding to each time point in the dataset
hours <- c(0, 1, 2, 4, 5, 6, 8, 10, 12, 14, 16, 20, 24, 30, 36, 42, 48)

# Read prior matrix
priorMatrix <- read.csv("Data/priorMatrix.csv", row.names=1)

# Make sure row and column names of the prior matrix are the same
colnames(priorMatrix) <- rownames(priorMatrix)
