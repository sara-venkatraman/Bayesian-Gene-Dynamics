# --- Script 4: Prior matrix construction ---

# Load datasets needed for constructing priors
source("DatasetLoader.R")

# Load functions in EmpiricalBayesFunctions.R
source("EmpiricalBayesFunctions.R")

# --- 1. Construct prior matrix from full STRING-DB scores ---

# Read STRING-DB similarity matrix
stringDBMatrix <- read.csv("../Processed Data for Clustering/Prior Matrices/Dme_STRING_matrix.csv", header=T, row.names=1)

# Construct a prior similarity matrix (binary) for a subset of DE genes using STRING-DB.
# For now, treat NAs in the STRING-DB matrix as either unknown or non-similar.
# Entry (i,j) of the prior matrix is 1 if genes i, j have a similarity score above 500 and
# is 0 otherwise. Entry (i,j) = 0 also if gene i and/or j was not included in STRING-DB.
# Matrix is formatted according to requirements in EmpiricalBayesFunctions.R.
stringMatrixSubset <- matrix(rep(0, subsetSize^2), nrow=subsetSize, ncol=subsetSize)
rownames(stringMatrixSubset) <- genesSubset;  colnames(stringMatrixSubset) <- genesSubset
for(i in 1:subsetSize) {
  for(j in i:subsetSize) {
    STRINGscore <- stringDBMatrix[Gene.Name.To.Flybase.ID(genesSubset[i]), Gene.Name.To.Flybase.ID(genesSubset[j])]
    if(!is.na(STRINGscore) && !is.null(STRINGscore)) { 
      stringMatrixSubset[i,j] <- STRINGscore 
      stringMatrixSubset[j,i] <- stringMatrixSubset[i,j]
    }
  }
}

priorMatrix <- (stringMatrixSubset >= 500) + 0

# Print the number of 0's and 1's in the binary prior matrix
table(priorMatrix)

# --- 2. Construct prior matrix from experimental STRING-DB scores ---

# Read STRING-DB similarity matrix
stringDBMatrix <- read.csv("../Processed Data for Clustering/Prior Matrices/Dme_STRING_exp_matrix.csv", header=T, row.names=1)

# Construct a prior similarity matrix (binary) for a subset of DE genes using STRING-DB.
# For now, treat NAs in the STRING-DB matrix as either unknown or non-similar.
# Entry (i,j) of the prior matrix is 1 if genes i, j have a similarity score above 500 and
# is 0 otherwise. Entry (i,j) = 0 also if gene i and/or j was not included in STRING-DB.
# Matrix is formatted according to requirements in EmpiricalBayesFunctions.R.
stringMatrixSubset <- matrix(rep(0, subsetSize^2), nrow=subsetSize, ncol=subsetSize)
rownames(stringMatrixSubset) <- genesSubset;  colnames(stringMatrixSubset) <- genesSubset
for(i in 1:subsetSize) {
  for(j in i:subsetSize) {
    STRINGscore <- stringDBMatrix[Gene.Name.To.Flybase.ID(genesSubset[i]), Gene.Name.To.Flybase.ID(genesSubset[j])]
    if(!is.na(STRINGscore) && !is.null(STRINGscore)) { 
      stringMatrixSubset[i,j] <- STRINGscore 
      stringMatrixSubset[j,i] <- stringMatrixSubset[i,j]
    }
  }
}

priorMatrix <- (stringMatrixSubset >= 100) + 0

# Print the number of 0's and 1's in the binary prior matrix
table(priorMatrix)

# --- 3. Construct prior matrix from second replicate of normalized counts ---

# Method 1: Set prior to 1 if the correlation between two genes is > 0.85
coefficientMatrix <- matrix(rep(0, subsetSize^2), nrow=subsetSize, ncol=subsetSize)
rownames(coefficientMatrix) <- genesSubset;  colnames(coefficientMatrix) <- genesSubset
for(i in 1:(subsetSize-1)) {
  for(j in (i+1):subsetSize) {
    countsProfile1 <- as.numeric(normCountsRep2[genesSubset[i],])
    countsProfile2 <- as.numeric(normCountsRep2[genesSubset[j],])
    geneCorrelation <- cor(countsProfile1, countsProfile2)
    coefficientMatrix[i,j] <- geneCorrelation
    coefficientMatrix[j,i] <- coefficientMatrix[i,j]
  }
}
priorMatrix <- (abs(coefficientMatrix) >= 0.85) + 0

# Print the number of 0's and 1's in the binary prior matrix
table(priorMatrix)

# Method 2: Set prior to 1 if the R^2 from regressing one gene on the other (in 
# either direction) is > 0.85
coefficientMatrix <- matrix(rep(0, subsetSize^2), nrow=subsetSize, ncol=subsetSize)
rownames(coefficientMatrix) <- genesSubset;  colnames(coefficientMatrix) <- genesSubset
for(i in 1:(subsetSize-1)) {
  for(j in (i+1):subsetSize) {
    countsProfile1 <- as.numeric(normCountsRep2[genesSubset[i],])
    countsProfile2 <- as.numeric(normCountsRep2[genesSubset[j],])
    model1 <- lm(countsProfile1 ~ countsProfile2)
    model2 <- lm(countsProfile2 ~ countsProfile1)
    R2 <- max(summary(model1)$r.squared, summary(model2)$r.squared)
    coefficientMatrix[i,j] <- R2
    coefficientMatrix[j,i] <- coefficientMatrix[i,j]
  }
}
priorMatrix <- (abs(coefficientMatrix) >= 0.85) + 0

# Print the number of 0's and 1's in the binary prior matrix
table(priorMatrix)

# --- Construct prior matrix using both experimental STRING scores *and* replicates ---
#     To construct this, run both section 2 and (one method from) section 3.

# Step 1: Construct binary prior matrix from STRING experimental scores as in section 2
stringDBMatrix <- read.csv("../Processed Data for Clustering/Prior Matrices/Dme_STRING_exp_matrix.csv", header=T, row.names=1)
stringMatrixSubset <- matrix(rep(0, subsetSize^2), nrow=subsetSize, ncol=subsetSize)
rownames(stringMatrixSubset) <- genesSubset;  colnames(stringMatrixSubset) <- genesSubset
for(i in 1:subsetSize) {
  for(j in i:subsetSize) {
    STRINGscore <- stringDBMatrix[Gene.Name.To.Flybase.ID(genesSubset[i]), Gene.Name.To.Flybase.ID(genesSubset[j])]
    if(!is.na(STRINGscore) && !is.null(STRINGscore)) { 
      stringMatrixSubset[i,j] <- STRINGscore 
      stringMatrixSubset[j,i] <- stringMatrixSubset[i,j]
    }
  }
}
priorMatrixString <- (stringMatrixSubset >= 200) + 0
table(priorMatrixString)

# Step 2: Construct binary prior matrix from replicate data as in section 3, method 1
coefficientMatrix <- matrix(rep(0, subsetSize^2), nrow=subsetSize, ncol=subsetSize)
rownames(coefficientMatrix) <- genesSubset;  colnames(coefficientMatrix) <- genesSubset
for(i in 1:(subsetSize-1)) {
  for(j in (i+1):subsetSize) {
    countsProfile1 <- as.numeric(normCountsRep2[genesSubset[i],])
    countsProfile2 <- as.numeric(normCountsRep2[genesSubset[j],])
    geneCorrelation <- cor(countsProfile1, countsProfile2)
    coefficientMatrix[i,j] <- geneCorrelation
    coefficientMatrix[j,i] <- coefficientMatrix[i,j]
  }
}
priorMatrixReplicates <- (abs(coefficientMatrix) >= 0.85) + 0
table(priorMatrixReplicates)

# Step 3: Combine prior matrices
priorMatrix <- priorMatrixString + priorMatrixReplicates

# Check counts of 0, 1, 2 entries
table(priorMatrix)

# Check which entries were 2
which(priorMatrix == 2, arr.ind = TRUE)
