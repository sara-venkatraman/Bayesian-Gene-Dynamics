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
