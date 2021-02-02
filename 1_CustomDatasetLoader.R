# Dataframe of expression measurements (as log_2 fold change) 
# geneData <- read.table("../Processed Data/Post-Mating Response/WvsV_2methods2fold_Wlog2fc.txt", row.names=1, header=T)
geneData <- read.table("../Processed Data/Post-Mating Response/WvsV_3methods_Wlog2fc.txt", row.names=1, header=T)

# Define the hours corresponding to each time point in the dataset
hours <- c(0.5, 1, 2, 3, 4, 5, 6, 8, 12, 24)

# Set gene names to be rownames of the expression dataframe
geneNames <- rownames(geneData)
geneSubset <- geneNames
