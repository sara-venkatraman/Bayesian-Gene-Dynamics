# --- Script 2: Bayesian lead-lag R^2 computations ---

# This script implements an optimized version of Algorithm 1 (Appendix A.2)
# in our paper. Approximate timing: 22 minutes for 1735 genes on a 2017 3.1
# GHz Intel Core i5 MacBook Pro.

# Run the script "1_DatasetLoader.R" first to load gene expression data
source("1_DatasetLoader.R")

# Load parallel computing package
library(parallel)

# Start the timer for this script
startTime <- Sys.time()

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
              intLowerBound, intUpperBound))
)

# Define matrix specifying the genes comprising each unique pair 
# (first row = gene 1, second row = gene 2)
genePairs <- combn(seq_along(geneDataList), 2)

# Get a list of response vectors Y for all regressions (time profile of gene 1)
Y <- lapply(genePairs[1, ], function(i) matrix(geneDataList[[i]]))

# Reshape the genePairs matrix into a list
genePairsAsList <- unname(as.list(as.data.frame(genePairs)))
genePairsAsList <- lapply(genePairsAsList, as.list) 

# Get a list of design matrices X for all regressions
Get.Design.Matrix <- function(idx1, idx2) {
  cbind(geneDataList[[idx2]], integrals[[idx2]], integrals[[idx1]], hours[], 1)
}
X <- mclapply(genePairsAsList, function(pair) do.call(Get.Design.Matrix, pair))

# Get a list of the prior matrix entries for all regressions
Get.Prior.Indicator <- function(idx1, idx2) {
  priorMatrix[idx1, idx2]
}
priorList <- mclapply(genePairsAsList, function(pair) do.call(Get.Prior.Indicator, pair))

# Define function to compute the three Bayesian lead-lag R^2 values given the
# design matrix X, response vector Y, and prior (0, 1, or NA) for one gene pair
Get.R2.Bayes <- function(x, y, prior) {
  # Get dimensions of design matrix x
  n <- nrow(x)
  p <- ncol(x)
  
  # Fit the three models (equations 4, 7, and 8 in our paper)
  LLR2model <- lm.fit(x, y)
  LLR2model.other <- lm.fit(x[,c(1,2,5)], y)
  LLR2model.own <- lm.fit(x[,3:5], y)
  
  # Set prior mean of regression coefficients
  if(is.na(prior))
    priorMean <- matrix(0, nrow=5, ncol=1)
  else
    priorMean <- matrix(c(prior > 0, prior > 0, 0, 0, 0), ncol=1)

  # Compute main LLR2 value
  LScoefs <- matrix(LLR2model$coefficients, ncol=1)
  LSfit <- matrix(LLR2model$fitted.values, ncol=1)
  if(is.na(prior) || prior > 0) {
    sigmaSq <- norm(y - LSfit, "2")^2 / (n-p)
    g <- ((norm(LSfit - x %*% priorMean,"2")^2 - p*sigmaSq) / (p*sigmaSq))[1]
  } else
    g <- 1
  posteriorMean <- (1/(1+g))*priorMean + (g/(1+g))*LScoefs
  posteriorFit <- x %*% posteriorMean
  LLR2 <- var(posteriorFit)/(var(posteriorFit) + var(y-posteriorFit))
  
  # Compute LLR2.other
  LScoefs <- matrix(LLR2model.other$coefficients, ncol=1)
  LSfit <- matrix(LLR2model.other$fitted.values, ncol=1)
  if(is.na(prior) || prior > 0) {
    sigmaSq <- norm(y - LSfit, "2")^2 / (n-p)
    g <- ((norm(LSfit - x[,c(1,2,5)] %*% priorMean[c(1,2,5),],"2")^2 - p*sigmaSq) / (p*sigmaSq))[1]
  } else
    g <- 1
  posteriorMean <- (1/(1+g))*priorMean[c(1,2,5),] + (g/(1+g))*LScoefs
  posteriorFit <- x[,c(1,2,5)] %*% posteriorMean
  LLR2.other <- var(posteriorFit)/(var(posteriorFit) + var(y-posteriorFit))

  # Compute LLR2.own
  LScoefs <- matrix(LLR2model.own$coefficients, ncol=1)
  LSfit <- matrix(LLR2model.own$fitted.values, ncol=1)
  if(is.na(prior) || prior > 0) {
    sigmaSq <- norm(y - LSfit, "2")^2 / (n-p)
    g <- ((norm(LSfit,"2")^2 - p*sigmaSq) / (p*sigmaSq))[1]
  } else
    g <- 1
  posteriorMean <- (g/(1+g))*LScoefs
  posteriorFit <- x[,3:5] %*% posteriorMean
  LLR2.own <- var(posteriorFit)/(var(posteriorFit) + var(y-posteriorFit))

  # Combine results into a list
  list(LLR2, LLR2.other, LLR2.own) 
}

# Get the data needed for each regression (design matrix, response, prior) as a list
XYprior_list <- Map(list, X, Y, priorList)

# Compute the three R^2 values for each pair of genes
R2Bayes <- mclapply(XYprior_list, function(x) do.call(Get.R2.Bayes, x))

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
bayesLLR2Mat <- Get.R2.Matrix.From.List(genePairs, lapply(R2Bayes, `[[`, 1))
bayesLLR2Mat.other <- Get.R2.Matrix.From.List(genePairs, lapply(R2Bayes, `[[`, 2))
bayesLLR2Mat.own <- Get.R2.Matrix.From.List(genePairs, lapply(R2Bayes, `[[`, 3))

# Add in the diagonal entries
for (i in seq_along(geneDataList)) {
  bayesLLR2Mat[i,i] <- 1
  bayesLLR2Mat.other[i,i] <- 1
  bayesLLR2Mat.own[i,i] <- 1
}

# Use gene names as the row and column names of each similarity matrix
colnames(bayesLLR2Mat) <- rownames(bayesLLR2Mat) <- names(geneDataList)
colnames(bayesLLR2Mat.other) <- rownames(bayesLLR2Mat.other) <- names(geneDataList)
colnames(bayesLLR2Mat.own) <- rownames(bayesLLR2Mat.own) <- names(geneDataList)

# Write R^2 similarity matrices to CSV
write.csv(bayesLLR2Mat, "BayesLLR2.csv")
write.csv(bayesLLR2Mat.other, "BayesLLR2_Other.csv")
write.csv(bayesLLR2Mat.own, "BayesLLR2_Own.csv")

# Stop timer and print runtime
Sys.time() - startTime

