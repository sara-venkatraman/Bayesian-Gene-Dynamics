# --- Script 2: Non-Bayesian lead-lag R^2 computations ---

# Requirements: Run the script "1_DatasetLoader.R" first.
# Timing: 3.56 seconds for 156 genes, 54.49 seconds for 951 genes, 3.20 minutes for 1735 genes

library(parallel)

# Start the timer for this script
startTime <- Sys.time()

# Format gene expression dataframe into a list (produces one list element
# per row of the dataframe)
geneDataList <- as.list(as.data.frame(t(geneData)))

# Define function for spline integrations
intLowerBound <- head(hours, -1)
intUpperBound <- hours[-1]
Integrator <- function(f, lower, upper) integrate(f, lower, upper)$value

# Integrate the splines
integrals <- mclapply(geneDataList, function(y) {
  c(0, mapply(Integrator, list(splinefun(x=hours, y=y, method = "natural")), intLowerBound, intUpperBound))
})

# Matrix specifying the genes comprising each unique pair 
# (first row = gene 1, second row = gene 2)
genePairs <- combn(seq_along(geneDataList), 2)

# Get response vectors for all regressions (time profile of gene 1)
Y <- lapply(genePairs[1, ], function(i) matrix(geneDataList[[i]]))

# Get design matrices for all regressions
Get.Design.Matrix <- function(idx1, idx2) {
  cbind(geneDataList[[idx2]], integrals[[idx2]], integrals[[idx1]], hours[], 1)
}
genePairsAsList <- unname(as.list(as.data.frame(genePairs)))
genePairsAsList <- lapply(genePairsAsList, as.list) 
X <- mclapply(genePairsAsList, function(pair) do.call(Get.Design.Matrix, pair))

# Define function to compute R^2
Get.R2 <- function(x, y) {
  # Fit the three models
  LLR2model <- lm.fit(x, y)
  LLR2model.other <- lm.fit(x[,c(1,2,5)], y)
  LLR2model.own <- lm.fit(x[,3:5], y)
  
  # Compute LLR2
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

# Compute all R^2 values
XY_list <- Map(list, X, Y)
R2NonBayes <- mclapply(XY_list, function(x) do.call(Get.R2, x))

# Define function to get R^2 similarity matrix from R^2 list
Get.R2.Matrix.From.List <- function(pairs, R2list) {
  # Initialize matrix
  n <- length(geneDataList)
  R2Matrix <- matrix(0L, n, n)
  
  # Get x, y value pairs (array indices) to fill both (redundant) halves
  all_pairs <- cbind(pairs, pairs[2:1, ])
  rows <- as.list(all_pairs[1, ])
  cols <- as.list(all_pairs[2, ])
  
  # Fill in matrix
  Fill.Matrix.Entry <- function(row.i, col.i, R2.i) R2Matrix[row.i, col.i] <<- R2.i
  invisible(Map(Fill.Matrix.Entry, rows, cols, R2list))
  R2Matrix
}

# Get matrix for each of the three R^2 values
nonBayesLLR2Mat <- Get.R2.Matrix.From.List(genePairs, lapply(R2NonBayes, `[[`, 1))
nonBayesLLR2Mat.other <- Get.R2.Matrix.From.List(genePairs, lapply(R2NonBayes, `[[`, 2))
nonBayesLLR2Mat.own <- Get.R2.Matrix.From.List(genePairs, lapply(R2NonBayes, `[[`, 3))

# Add in the diagonals
for (i in seq_along(geneDataList)) {
  nonBayesLLR2Mat[i,i] <- 1
  nonBayesLLR2Mat.other[i,i] <- 1
  nonBayesLLR2Mat.own[i,i] <- 1
}

# Add gene names as row and column names
colnames(nonBayesLLR2Mat) <- rownames(nonBayesLLR2Mat) <- names(geneDataList)
colnames(nonBayesLLR2Mat.other) <- rownames(nonBayesLLR2Mat.other) <- names(geneDataList)
colnames(nonBayesLLR2Mat.own) <- rownames(nonBayesLLR2Mat.own) <- names(geneDataList)

# Stop timer and print runtime
Sys.time() - startTime

# Write R^2 similarity matrices to CSV, if desired
write.csv(nonBayesLLR2Mat, "LLR2.csv")
write.csv(nonBayesLLR2Mat.other, "LLR2_other.csv")
write.csv(nonBayesLLR2Mat.own, "LLR2_own.csv")


