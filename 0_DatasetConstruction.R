# --- Script 1: Dataset construction ---

# Create CSV dataset of differentially-expressed genes with gene names and Flybase IDs
# (no need to run if DEGenes.csv already exists)
load("../Drosophila Data (Processed)/full_results.Rdata")
DEGenes <- read.csv("../Drosophila Data (Processed)/nonFlatGeneData.csv")
DEGeneIDs <- as.character(full_results[DEGenes$X, "gene_ID"])
DEGenes <- cbind(DEGeneIDs, DEGenes)
colnames(DEGenes)[1:2] <- c("geneID", "geneName") 
write.csv(DEGenes, file="DEGenes.csv", row.names=F)

# ------------------------------------------------------------------

# Create three CSV datasets of differentially-expressed genes: one with 
# log-fold change data, one with normalized count data from replicate 1,
# and one with normalized count data from replicate 2
load("../Drosophila Data (Processed)/full_results.Rdata")
DEgeneNames <- read.csv("../Drosophila Data (Processed)/DEGeneNames.csv", stringsAsFactors=F, header=F)
DEgeneIDs <- as.data.frame(as.character(full_results[DEgeneNames$V1, "gene_ID"]))
colnames(DEgeneNames) <- c("geneName");  colnames(DEgeneIDs) <- c("geneID")

DEnormCountsRep1 <- full_results[DEgeneNames$geneName, 75:95]
DEnormCountsRep2 <- full_results[DEgeneNames$geneName, 96:115] # Missing time point 4
DElogChange <- full_results[DEgeneNames$geneName, 3:19]        # Short dataset

write.csv(DEgeneNames, "../Processed Data for Clustering/DEgeneNames.csv", row.names=F)
write.csv(DEgeneIDs, "../Processed Data for Clustering/DEgeneIDs.csv", row.names=F)
write.csv(DEnormCountsRep1, "../Processed Data for Clustering/DEnormCountsRep1.csv")
write.csv(DEnormCountsRep2, "../Processed Data for Clustering/DEnormCountsRep2.csv")
write.csv(DElogChange, "../Processed Data for Clustering/DElogChange.csv")

# ------------------------------------------------------------------

# Create CSV dataset of STRING-DB scores in the form of a similarity matrix
# (no need to run if Dme_STRING_matrix.csv already exists)

# Read STRING-DB data and remove the protein IDs in the first two columns
stringDB <- read.table("../Drosophila Data (Processed)/Dme_STRING_full.txt", header=T)
stringDB <- stringDB[, -c(1,2)]

# Reshape the STRING-DB into a matrix with one row per Flybase ID (ID_A)
library(tidyr)
stringDBMatrix <- spread(stringDB, FBgn_item_id_b, combined_score)
rownames(stringDBMatrix) <- stringDBMatrix$FBgn_item_id_a
stringDBMatrix <- stringDBMatrix[,-1]
write.csv(stringDBMatrix, file="Dme_STRING_matrix.csv")

# ------------------------------------------------------------------

# Create CSV dataset of *experimental* STRING-DB scores in the form of a 
# similarity matrix (no need to run if Dme_STRING_exp_matrix.csv already exists)

# Read STRING-DB data and remove all columns but the gene names and experimental scores
stringDBexp <- read.table("../Drosophila Data (Processed)/Dme_STRING_individual_scores.txt", header=T)
stringDBexp <- stringDBexp[,c("gene1","gene2","experimental")]
stringDBexp <- stringDBexp[complete.cases(stringDBexp),]

# Reshape the STRING-DB into a matrix with one row per Flybase ID (gene1)
library(tidyr)
stringDBMatrix <- spread(stringDBexp, gene2, experimental)
rownames(stringDBMatrix) <- stringDBMatrix$gene1
stringDBMatrix <- stringDBMatrix[,-1]
write.csv(stringDBMatrix, file="Dme_STRING_exp_matrix.csv")

# ------------------------------------------------------------------

# Create CSV dataset of differentially-expressed genes and their neighbors as
# identified by composite STRING-DB scores. To run this code, first run the functions 
# Gene.Name.To.Flybase.ID and Flybase.ID.To.Gene.Name in EmpiricalBayesFunctions.R.

# Returns the Flybase ID of a given gene name
Gene.Name.To.Flybase.ID <- function(geneName) {
  index <- match(geneName, geneNames)
  flybaseID <- geneIDs[index]
  return(flybaseID)
}

# Returns the gene name of a given Flybase ID
Flybase.ID.To.Gene.Name <- function(geneID) {
  index <- match(geneID, geneIDs)
  geneName <- geneNames[index]
  return(geneName)
}

# Read gene names and Flybase IDs (converted to character) and read composite STRING scores
geneNames <- read.csv("../Processed Data/Differentially-Expressed/GeneNames.csv", header=T)
geneIDs <- read.csv("../Processed Data/Differentially-Expressed/GeneIDs.csv", header=T)  
geneNames <- as.character(geneNames$geneName)
geneIDs <- as.character(geneIDs$geneID)
stringDBMatrix <- read.csv("../Processed Data/Prior Matrices/Dme_STRING_matrix.csv", header=T, row.names=1)

# Initialize vector for storing names of neighbors not in the DE gene set and the number
# of new additions (i.e. length of DEgenesNeighbors)
DEgenesNeighbors <- c()
numNewNeighbors <- 0

# For each DE gene, find its row index in the STRING matrix, then search that row for other
# genes with which it is associated (that are not already DE genes)
for(i in 1:length(geneNames)) {
  # Find Flybase ID of current gene in STRING-DB matrix 
  index <- match(Gene.Name.To.Flybase.ID(geneNames[i]), rownames(stringDBMatrix))
  
  # If the index is not NA, i.e. if the current gene is represented in the matrix, then...
  if(!is.na(index)) {
    # Get the IDs of the associated genes
    STRINGscores <- stringDBMatrix[index,]
    
    # Find the column indices of genes scores above some threshold
    associatedIndices <- which(STRINGscores >= 950)
    associatedIDs <- colnames(stringDBMatrix)[associatedIndices]
    
    # Find which of the associated genes are already in the set of DE genes
    notDEgenes <- !(associatedIDs %in% geneIDs) # indicates whether or not each gene is DE
    indicesOfNeighbors <- which(notDEgenes)
    neighborIDs <- associatedIDs[indicesOfNeighbors]
    numAdditions <- length(indicesOfNeighbors)
    if(numAdditions > 0) {
      DEgenesNeighbors[(numNewNeighbors+1):(numNewNeighbors+numAdditions)] <- neighborIDs
    }
    numNewNeighbors <- length(DEgenesNeighbors)
  }
}

DEgenesNeighbors <- unique(DEgenesNeighbors)

# Sanity check: these numbers should all be high
# stringDBMatrix[index, associatedIDs[1:5]]

# How many new neighbors are added?
# - With STRING score threshold of >= 200, we add 9328 new genes to the set of 951.
# - With STRING score threshold of >= 500, we add 5607 new genes.
# - With STRING score threshold of >= 700, we add 3539 new genes.
# - With STRING score threshold of >= 900, we add 2147 new genes (2013 after excluding NAs).

# Construct a dataset containing the log-fold change expression values from both the
# differentially-expressed genes and the neighbors. Also save the normalized counts from
# both replicates

# Load the full dataset from the experiment and set row names of full dataset to Flybase ID for easier access to certain rows
load("../Drosophila Data (Processed)/full_results.Rdata")

rownames(full_results) <- full_results$gene_ID

# Extract the expression data for the neighbor genes identified above and exclude rows with missing data
geneDataNewNeighbors <- full_results[DEgenesNeighbors,]
geneDataNewNeighbors <- geneDataNewNeighbors[complete.cases(geneDataNewNeighbors),]

# Set row names back to gene names and save only the log-fold change data from early timeframe
rownames(geneDataNewNeighbors) <- geneDataNewNeighbors$gene_name

# Extract log-fold changes and data from two replicates
neighborsLogChange <- geneDataNewNeighbors[, 3:19]
neighborsNormCountsRep1 <- geneDataNewNeighbors[, 75:95]
neighborsNormCountsRep2 <- geneDataNewNeighbors[, 96:115] # Missing time point 4

# Write CSV files for log-fold changes, two replicates, gene names, and gene IDs
write.csv(neighborsLogChange, "../Processed Data/DEneighbors950LogChange.csv")
write.csv(neighborsNormCountsRep1, "../Processed Data/DEneighbors950normCountsRep1.csv")
write.csv(neighborsNormCountsRep2, "../Processed Data/DEneighbors950normCountsRep2.csv")
write.csv(geneDataNewNeighbors$gene_name, "../Processed Data/DEneighbors950geneNames.csv", row.names=F)
write.csv(geneDataNewNeighbors$gene_ID, "../Processed Data/DEneighbors950geneIDs.csv", row.names=F)

# ------------------------------------------------------------------

# Create CSV dataset (including log-fold change, replicate data, gene IDs, and gene names) 
# which combines data from differentially-expressed genes *and* their neighbors as
# identified by STRING-DB

# Read data for DE genes 
DElogChange <- read.csv("../Processed Data/Differentially-Expressed/LogChange.csv")
DEnormCountsRep1 <- read.csv("../Processed Data/Differentially-Expressed/NormCountsRep1.csv")
DEnormCountsRep2 <- read.csv("../Processed Data/Differentially-Expressed/NormCountsRep2.csv")
DEnames <- read.csv("../Processed Data/Differentially-Expressed/GeneNames.csv")
DEids <- read.csv("../Processed Data/Differentially-Expressed/GeneIDs.csv")

# Read data for DE gene neighbors
neighborsLogChange <- read.csv("../Processed Data/Neighbors (STRING > 950)/LogChange.csv")
neighborsNormCountsRep1 <- read.csv("../Processed Data/Neighbors (STRING > 950)/NormCountsRep1.csv")
neighborsNormCountsRep2 <- read.csv("../Processed Data/Neighbors (STRING > 950)/NormCountsRep2.csv")
neighborsNames <- read.csv("../Processed Data/Neighbors (STRING > 950)/GeneNames.csv")
neighborsIDs <- read.csv("../Processed Data/Neighbors (STRING > 950)/GeneIDs.csv")

# First combine gene names from both sets and save alphabetical ordering
colnames(DEnames) <- "geneName";  colnames(neighborsNames) <- "geneName"
combinedNames <- rbind(DEnames, neighborsNames)
ordering <- order(combinedNames)

# Combine remaining datasets
colnames(DEids) <- "geneID";  colnames(neighborsIDs) <- "geneID"
colnames(DElogChange)[1] <- "geneName";  colnames(neighborsLogChange)[1] <- "geneName" 
combinedIDs <- rbind(DEids, neighborsIDs)
combinedLogChange <- rbind(DElogChange, neighborsLogChange)
combinedNormCountsRep1 <- rbind(DEnormCountsRep1, neighborsNormCountsRep1)
combinedNormCountsRep2 <- rbind(DEnormCountsRep2, neighborsNormCountsRep2)

# Sort all combined datasets
combinedNames$geneName <- combinedNames$geneName[ordering]
combinedIDs$geneID <- combinedIDs$geneID[ordering]
combinedLogChange <- combinedLogChange[ordering,]
combinedNormCountsRep1 <- combinedNormCountsRep1[ordering,]
combinedNormCountsRep2 <- combinedNormCountsRep2[ordering,]

# Adjust row names
rownames(combinedLogChange) <- combinedNames$geneName
rownames(combinedNormCountsRep1) <- combinedNames$geneName
rownames(combinedNormCountsRep2) <- combinedNames$geneName

# Remove redundant gene name column
combinedLogChange <- combinedLogChange[,-1]
combinedNormCountsRep1 <- combinedNormCountsRep1[,-1]
combinedNormCountsRep2 <- combinedNormCountsRep2[,-1]

# Write CSVs
write.csv(combinedLogChange, "../Processed Data/CombinedLogChange.csv")
write.csv(combinedNormCountsRep1, "../Processed Data/CombinedNormCountsRep1.csv")
write.csv(combinedNormCountsRep2, "../Processed Data/CombinedNormCountsRep2.csv")
write.csv(combinedNames$geneName, "../Processed Data/CombinedGeneNames.csv", row.names=F)
write.csv(combinedIDs$geneID, "../Processed Data/CombinedGeneIDs.csv", row.names=F)

# --- Create a prior information matrix just for the genes in the combined dataset ---

stringDBMatrix <- read.csv("../Processed Data/Prior Matrices/Dme_STRING_matrix.csv", header=T, row.names=1)
geneData <- read.csv("../Processed Data/Combined Genes/LogChange.csv", row.names=1)
geneIDs <- read.csv("../Processed Data/Combined Genes/GeneIDs.csv")[,1]
geneNames <- read.csv("../Processed Data/Combined Genes/GeneNames.csv")[,1]
riorMatrix <- matrix(0L, nrow(geneData), nrow(geneData))
rownames(priorMatrix) <- geneIDs;  colnames(priorMatrix) <- geneIDs
availableStringPrior <- as.matrix(stringDBMatrix[geneIDs[geneIDs %in% rownames(stringDBMatrix)], geneIDs[geneIDs %in% rownames(stringDBMatrix)]])
availableStringPrior[is.na(availableStringPrior)] <- 0
priorMatrix[geneIDs[geneIDs %in% rownames(stringDBMatrix)], geneIDs[geneIDs %in% rownames(stringDBMatrix)]] <- availableStringPrior
priorMatrix[priorMatrix <= 500] <- 0
rownames(priorMatrix) <- geneNames;  colnames(priorMatrix) <- geneNames
write.csv(priorMatrix, "../Processed Data/Prior Matrices/CombinedGenesPriorMatrix.csv")

# --- Create a prior information matrix just for the genes in the DE dataset ---

stringDBMatrix <- read.csv("../Processed Data/Prior Matrices/Dme_STRING_matrix.csv", header=T, row.names=1)
geneData <- read.csv("../Processed Data/Differentially-Expressed/LogChange.csv", row.names=1)
geneIDs <- read.csv("../Processed Data/Differentially-Expressed/GeneIDs.csv")[,1]
geneNames <- read.csv("../Processed Data/Differentially-Expressed/GeneNames.csv")[,1]
priorMatrix <- matrix(0L, nrow(geneData), nrow(geneData))
rownames(priorMatrix) <- geneIDs;  colnames(priorMatrix) <- geneIDs
availableStringPrior <- as.matrix(stringDBMatrix[geneIDs[geneIDs %in% rownames(stringDBMatrix)], geneIDs[geneIDs %in% rownames(stringDBMatrix)]])
availableStringPrior[is.na(availableStringPrior)] <- 0
priorMatrix[geneIDs[geneIDs %in% rownames(stringDBMatrix)], geneIDs[geneIDs %in% rownames(stringDBMatrix)]] <- availableStringPrior
priorMatrix[priorMatrix <= 500] <- 0
rownames(priorMatrix) <- geneNames;  colnames(priorMatrix) <- geneNames
write.csv(priorMatrix, "../Processed Data/Prior Matrices/DEPriorMatrix.csv")

# --- Create a prior information matrix for DE genes with entries 0,1,NA ---

# 0 = neither STRING nor replicate data agrees that there's an association
# 1 = one of the two sources indicates an association

subdirectory <- "Combined Genes"

stringDBMatrix <- read.csv("../Processed Data/Prior Matrices/Dme_STRING_matrix.csv", header=T, row.names=1)
geneData <- read.csv(paste("../Processed Data/", subdirectory, "/LogChange.csv", sep=""), row.names=1)
geneIDs <- read.csv(paste("../Processed Data/", subdirectory, "/GeneIDs.csv", sep=""))[,1]
geneNames <- read.csv(paste("../Processed Data/", subdirectory, "/GeneNames.csv", sep=""))[,1]
normCountsRep1 <- read.csv(paste("../Processed Data/", subdirectory, "/NormCountsRep1.csv", sep=""), row.names=1)
normCountsRep2 <- read.csv(paste("../Processed Data/", subdirectory, "/NormCountsRep2.csv", sep=""), row.names=1)

# First get the prior from STRING
priorMatrixString <- matrix(NA, nrow(geneData), nrow(geneData))
rownames(priorMatrixString) <- geneIDs;  colnames(priorMatrixString) <- geneIDs
availableStringPrior <- as.matrix(stringDBMatrix[geneIDs[geneIDs %in% rownames(stringDBMatrix)], geneIDs[geneIDs %in% rownames(stringDBMatrix)]])
availableStringPrior[is.na(availableStringPrior)] <- 0
priorMatrixString[geneIDs[geneIDs %in% rownames(stringDBMatrix)], geneIDs[geneIDs %in% rownames(stringDBMatrix)]] <- availableStringPrior
priorMatrixString[priorMatrixString <= 500] <- 0

# Next get the prior from the replicate data (average of correlations from 
# the two replicates of normalized counts)
priorMatrixRep <- (cor(t(normCountsRep1)) + cor(t(normCountsRep2)))/2
priorMatrixRep[abs(priorMatrixRep) < 0.8] <- 0

# Combine the priors (1 = associated, NA = unknown, 0 = not associated)
priorMatrix <- (priorMatrixRep != 0) + 0      # Store 1's from replicate data
priorMatrix[priorMatrixString > 0] <- 1       # Store 1's from STRING scores
priorMatrix[is.na(priorMatrixString)] <- NA   # Replace STRING unknowns with NA
rownames(priorMatrix) <- geneNames;  colnames(priorMatrix) <- geneNames

# Combine the priors (equal weight)
# priorMatrix <- ((priorMatrixString != 0) + 0) + ((priorMatrixRep != 0) + 0)
# rownames(priorMatrix) <- geneNames;  colnames(priorMatrix) <- geneNames

# Write to CSV
write.csv(priorMatrix, paste("../Processed Data/", subdirectory, "/PriorMatrix.csv", sep=""))

# --- Create a prior information matrix for the post-mating response data ---

stringDBMatrix <- read.csv("../Processed Data/Prior Matrices/Dme_STRING_matrix.csv", header=T, row.names=1)
load("../Processed Data/Post-Mating Response/C_expr_matrix.RData")
geneData <- Cdata;  remove(Cdata)
geneNames <- rownames(geneData)

priorMatrixString <- matrix(NA, nrow(geneData), nrow(geneData))
rownames(priorMatrixString) <- geneNames;  colnames(priorMatrixString) <- geneNames
availableStringPrior <- as.matrix(stringDBMatrix[geneNames[geneNames %in% rownames(stringDBMatrix)], geneNames[geneNames %in% rownames(stringDBMatrix)]])
priorMatrixString[geneNames[geneNames %in% rownames(stringDBMatrix)], geneNames[geneNames %in% rownames(stringDBMatrix)]] <- availableStringPrior
priorMatrixString[priorMatrixString <= 500] <- 0
priorMatrixString[priorMatrixString > 0] <- 1

write.csv(priorMatrixString, "../Processed Data/Post-Mating Response/PriorMatrix.csv")



