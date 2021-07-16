# --- Non-Bayesian lead-lag R^2 computations ---

# This script implements the non-Bayesian lead-lag R^2 calculations 
# described in Section 2.1 of our paper. The code here is structured 
# similarly to the script "2_LLR2Calculations.R" in the parent directory. 
# Approximate timing: 3.20 minutes for 1735 genes on a 2017 3.1 GHz Intel 
# Core i5 MacBook Pro.

# Start the timer for this script
startTime <- Sys.time()

# Load parallel computing package
library(parallel)

# Read the gene expression data
geneData <- read.csv("../Data/geneData.csv", row.names=1)
geneNames <- rownames(geneData)
hours <- c(0, 1, 2, 4, 5, 6, 8, 10, 12, 14, 16, 20, 24, 30, 36, 42, 48)

# Reshape gene expression dataframe into a list
geneDataList <- as.list(as.data.frame(t(geneData)))

# Define function for numerically integrating between two points
Integrator <- function(f, lower, upper) integrate(f, lower, upper)$value

# Integrate spline interpolants constructed from each gene expression 
# time profile. Each interpolant is integrated up to each time point.
intLowerBound <- head(hours, -1)
intUpperBound <- hours[-1]
integrals <- mclapply(geneDataList, function(y) 
  c(0, mapply(Integrator, list(splinefun(x=hours, y=y, method="natural")),
              intLowerBound, intUpperBound)))

# Define matrix specifying the genes comprising each unique pair 
# (first row = gene 1, second row = gene 2)
genePairs <- combn(seq_along(geneDataList), 2)

# Get a list of response vectors Y for all regressions (time profile of gene 1)
Y <- lapply(genePairs[1, ], function(i) matrix(geneDataList[[i]]))

# Reshape the genePairs matrix into a list
genePairsAsList <- unname(as.list(as.data.frame(genePairs)))
genePairsAsList <- lapply(genePairsAsList, as.list) 

# Get design matrices for all regressions
Get.Design.Matrix <- function(idx1, idx2) {
  cbind(geneDataList[[idx2]], integrals[[idx2]], integrals[[idx1]], hours[], 1)
}
X <- mclapply(genePairsAsList, function(pair) do.call(Get.Design.Matrix, pair))

# Define function to compute the three non-Bayesian lead-lag R^2 values 
# given the design matrix X and response vector Y for one gene pair
Get.R2 <- function(x, y) {
  # Fit the three models (equations 4, 7, and 8 in our paper)
  LLR2model <- lm.fit(x, y)
  LLR2model.other <- lm.fit(x[,c(1,2,5)], y)
  LLR2model.own <- lm.fit(x[,3:5], y)
  
  # Compute main LLR2 value
  MSS <- sum(LLR2model$fitted.values^2)
  RSS <- sum(LLR2model$residuals^2)
  LLR2 <- MSS/(MSS + RSS)
  
  # Compute LLR2.other
  MSS <- sum(LLR2model.other$fitted.values^2)
  RSS <- sum(LLR2model.other$residuals^2)
  LLR2.other <- MSS/(MSS + RSS) 
  
  # Compute LLR2.own
  MSS <- sum(LLR2model.own$fitted.values^2)
  RSS <- sum(LLR2model.own$residuals^2)
  LLR2.own <- MSS/(MSS + RSS) 
  
  # Combine results into a list
  list(LLR2, LLR2.other, LLR2.own)
}

# Get the data needed for each regression (design matrix, response, prior) as a list
XY_list <- Map(list, X, Y)

# Compute the three R^2 values for each pair of genes
R2NonBayes <- mclapply(XY_list, function(x) do.call(Get.R2, x))

# Define function to get R^2 similarity matrix from R^2 list
Get.R2.Matrix.From.List <- function(pairs, R2list) {
  # Initialize matrix
  n <- length(geneDataList)
  R2Matrix <- matrix(0L, n, n)
  
  # Get indices of R^2 list that will fill both triangular halves of the matrix
  all_pairs <- cbind(pairs, pairs[2:1, ])
  rows <- as.list(all_pairs[1, ])
  cols <- as.list(all_pairs[2, ])
  
  # Fill in matrix
  Fill.Matrix.Entry <- function(row.i, col.i, R2.i) R2Matrix[row.i, col.i] <<- R2.i
  invisible(Map(Fill.Matrix.Entry, rows, cols, R2list))
  R2Matrix
}

# Get matrix for each of the three LLR2 values
nonBayesLLR2Mat <- Get.R2.Matrix.From.List(genePairs, lapply(R2NonBayes, `[[`, 1))
nonBayesLLR2Mat.other <- Get.R2.Matrix.From.List(genePairs, lapply(R2NonBayes, `[[`, 2))
nonBayesLLR2Mat.own <- Get.R2.Matrix.From.List(genePairs, lapply(R2NonBayes, `[[`, 3))

# Add in the diagonal entries
for (i in seq_along(geneDataList)) {
  nonBayesLLR2Mat[i,i] <- 1
  nonBayesLLR2Mat.other[i,i] <- 1
  nonBayesLLR2Mat.own[i,i] <- 1
}

# Use gene names as the row and column names of each similarity matrix
colnames(nonBayesLLR2Mat) <- rownames(nonBayesLLR2Mat) <- names(geneDataList)
colnames(nonBayesLLR2Mat.other) <- rownames(nonBayesLLR2Mat.other) <- names(geneDataList)
colnames(nonBayesLLR2Mat.own) <- rownames(nonBayesLLR2Mat.own) <- names(geneDataList)

# Write R^2 similarity matrices to CSV
write.csv(nonBayesLLR2Mat, "NonBayesLLR2.csv")
write.csv(nonBayesLLR2Mat.other, "NonBayesLLR2_other.csv")
write.csv(nonBayesLLR2Mat.own, "NonBayesLLR2_own.csv")

# Stop timer and print runtime
Sys.time() - startTime
