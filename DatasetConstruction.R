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

# Create a CSV dataset of differentially-expressed genes and their neighbors as
# identified by composite STRING-DB scores. To run this code, first run the functions 
# Gene.Name.To.Flybase.ID and Flybase.ID.To.Gene.Name in EmpiricalBayesFunctions.R.

# Read gene names and Flybase IDs (converted to character) and read composite STRING scores
geneNames <- read.csv("../Processed Data for Clustering/DEgeneNames.csv", header=T)
geneIDs <- read.csv("../Processed Data for Clustering/DEgeneIDs.csv", header=T)  
geneNames <- as.character(geneNames$geneName)
geneIDs <- as.character(geneIDs$geneID)
stringDBMatrix <- read.csv("../Drosophila Data (Processed)/Dme_STRING_matrix.csv", header=T, row.names=1)

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
    associatedIndices <- which(STRINGscores >= 900)
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

# Write CSV files
write.csv(neighborsLogChange, "../Processed Data for Clustering/DEneighbors900LogChange.csv")
write.csv(neighborsNormCountsRep1, "../Processed Data for Clustering/DEneighbors900normCountsRep1.csv")
write.csv(neighborsNormCountsRep2, "../Processed Data for Clustering/DEneighbors900normCountsRep2.csv")
