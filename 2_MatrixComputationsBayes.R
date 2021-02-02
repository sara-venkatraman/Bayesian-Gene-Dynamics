# --- Script 2: Non-Bayesian lead-lag R^2 computations ---

# Requirements: Run the script "1_DatasetLoader.R" first.
# Timing: 2.58 minutes for 951 genes, 21.8 minutes for 1735 genes

library(parallel)
# options(mc.cores = 3L)

# Start the timer for this script
startTime <- Sys.time()

# Format gene data into a list
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

Get.Prior.Indicator <- function(idx1, idx2) {
  (priorMatrix[idx1, idx2] > 0) + 0
}

priorList <- mclapply(genePairsAsList, function(pair) do.call(Get.Prior.Indicator, pair))

# Define function to compute R^2
Get.R2.Bayes <- function(x, y, prior) {
  # Get dimensions of design matrix x
  n <- nrow(x)
  p <- ncol(x)

  # Fit the three models
  LLR2model <- lm.fit(x, y)
  LLR2model.other <- lm.fit(x[,c(1,2,5)], y)
  LLR2model.own <- lm.fit(x[,3:5], y)

  # Set prior mean of regression coefficients
  priorMean <- matrix(c(prior, prior, 0, 0, 0), ncol=1)
  
  # Compute LLR2
  LScoefs <- matrix(LLR2model$coefficients, ncol=1)
  LSfit <- matrix(LLR2model$fitted.values, ncol=1)
  sigmaSq <- norm(y - LSfit, "2")^2 / (n-p)
  g <- ((norm(LSfit - x %*% priorMean,"2")^2 - p*sigmaSq) / (p*sigmaSq))[1]
  posteriorMean <- (1/(1+g))*priorMean + (g/(1+g))*LScoefs
  posteriorFit <- x %*% posteriorMean
  LLR2 <- var(posteriorFit)/(var(posteriorFit) + var(y-posteriorFit))

  # Compute LLR2.other
  LScoefs <- matrix(LLR2model.other$coefficients, ncol=1)
  LSfit <- matrix(LLR2model.other$fitted.values, ncol=1)
  sigmaSq <- norm(y - LSfit, "2")^2 / (n-p)
  g <- ((norm(LSfit - x[,c(1,2,5)] %*% priorMean[c(1,2,5),],"2")^2 - p*sigmaSq) / (p*sigmaSq))[1]
  posteriorMean <- (1/(1+g))*priorMean[c(1,2,5),] + (g/(1+g))*LScoefs
  posteriorFit <- x[,c(1,2,5)] %*% posteriorMean
  LLR2.other <- var(posteriorFit)/(var(posteriorFit) + var(y-posteriorFit))

  # Compute LLR2.own
  LScoefs <- matrix(LLR2model.own$coefficients, ncol=1)
  LSfit <- matrix(LLR2model.own$fitted.values, ncol=1)
  sigmaSq <- norm(y - LSfit, "2")^2 / (n-p)
  g <- ((norm(LSfit - x[,3:5] %*% priorMean[3:5,],"2")^2 - p*sigmaSq) / (p*sigmaSq))[1]
  posteriorMean <- (1/(1+g))*priorMean[3:5,] + (g/(1+g))*LScoefs
  posteriorFit <- x[,3:5] %*% posteriorMean
  LLR2.own <- var(posteriorFit)/(var(posteriorFit) + var(y-posteriorFit))

  # Combine results into a list
  list(LLR2, LLR2.other, LLR2.own) 
}

# Compute all R^2 values
XYprior_list <- Map(list, X, Y, priorList)
R2Bayes <- mclapply(XYprior_list, function(x) do.call(Get.R2.Bayes, x))

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
bayesLLR2Mat <- Get.R2.Matrix.From.List(genePairs, lapply(R2Bayes, `[[`, 1))
bayesLLR2Mat.other <- Get.R2.Matrix.From.List(genePairs, lapply(R2Bayes, `[[`, 2))
bayesLLR2Mat.own <- Get.R2.Matrix.From.List(genePairs, lapply(R2Bayes, `[[`, 3))

# Add in the diagonals
for (i in seq_along(geneDataList)) {
  bayesLLR2Mat[i,i] <- 1
  bayesLLR2Mat.other[i,i] <- 1
  bayesLLR2Mat.own[i,i] <- 1
}

# Add gene names as row and column names
colnames(bayesLLR2Mat) <- rownames(bayesLLR2Mat) <- names(geneDataList)
colnames(bayesLLR2Mat.other) <- rownames(bayesLLR2Mat.other) <- names(geneDataList)
colnames(bayesLLR2Mat.own) <- rownames(bayesLLR2Mat.own) <- names(geneDataList)

# Stop timer and print runtime
Sys.time() - startTime