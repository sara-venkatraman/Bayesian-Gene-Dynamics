# ======= Script 1: Lead-lag R^2 computations =======

# Load parallel computing package
library(parallel)

#' @description
#' Compute a matrix of Bayesian or ordinary least-squares lead-lag \eqn{R^2} 
#' (LLR2) values for a given gene expression dataset.
#' 
#' Approximate timing: 22 minutes for 1735 genes with the Bayesian method, 
#' or 5.56 minutes with the non-Bayesian (least-squares) method, on a 2017
#' 3.1 GHz Intel Core i5 MacBook Pro.
#'
#' @param geneData N-by-n numeric matrix of N genes' expressions measured 
#' at n time points. The row names should be the names of the genes.
#' @param hours Length-n numeric vector containing the hours corresponding
#' to each time point.
#' @param bayes Logical indicating whether or not the Bayesian method
#' should be used to calculate the LLR2 values. Default TRUE.
#' @param priorMatrix (Bayesian method only) N-by-N prior adjacency 
#' matrix whose (i,j) entry is 1 if the two corresponding genes are likely 
#' to be associated, 0 if unlikely, and NA for unknown.
#' @param writeToCSV Logical indicating whether the LLR2 matrices should be 
#' written to CSV files in the current working directory. Default TRUE.
#'
#' @return Length-3 list containing symmetric N-by-N matrices \code{LLR2Mat},
#' \code{LLR2Mat.other}, and \code{LLR2Mat.own}. For details of how the 
#' three LLR2 values are calculated, see Sections 2.1 and 3 of Venkatraman 
#' et al. 2021 or source code for the functions \code{Compute.LLR2.Bayes} and
#' \code{Compute.LLR2.OLS}.
#' @export
LLR2 <- function(geneData, hours, bayes=TRUE, priorMatrix=NULL, writeToCSV=TRUE) {
  # For Bayesian LLR2 calculations, do not proceed if priorMatrix is missing
  if(bayes == TRUE && is.null(priorMatrix)) {
    stop("Missing prior adjacency matrix needed for Bayesian LLR2 calculations")
  }
  
  # Start the timer
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
                intLowerBound, intUpperBound)))
  
  # Define matrix specifying the genes comprising each unique pair 
  # (first row = gene 1, second row = gene 2)
  genePairs <- combn(seq_along(geneDataList), 2)
  
  # Get a list of response vectors Y for all regressions
  Y <- lapply(genePairs[1, ], function(i) matrix(geneDataList[[i]]))
  
  # Reshape the genePairs matrix into a list
  genePairsAsList <- unname(as.list(as.data.frame(genePairs)))
  genePairsAsList <- lapply(genePairsAsList, as.list) 
  
  # Get a list of design matrices X for all regressions
  Get.Design.Matrix <- function(idx1, idx2) {
    cbind(geneDataList[[idx2]], integrals[[idx2]], integrals[[idx1]], hours[], 1)
  }
  X <- mclapply(genePairsAsList, function(pair) do.call(Get.Design.Matrix, pair))
  
  # Compute the three R^2 values for each pair of genes
  if(bayes == TRUE) {
    Get.Prior.Indicator <- function(idx1, idx2) priorMatrix[idx1, idx2]
    priorList <- mclapply(genePairsAsList, function(pair) do.call(Get.Prior.Indicator, pair))
    regressionDataList <- Map(list, X, Y, priorList)
    R2 <- mclapply(regressionDataList, function(x) do.call(Compute.LLR2.Bayes, x))
  } else {
    regressionDataList <- Map(list, X, Y)
    R2 <- mclapply(regressionDataList, function(x) do.call(Compute.LLR2.OLS, x))
  }
    
  # Get a symmetric matrix for each of the three LLR2 values
  LLR2Mat <- Get.R2.Matrix.From.List(genePairs, lapply(R2, `[[`, 1))
  LLR2Mat.other <- Get.R2.Matrix.From.List(genePairs, lapply(R2, `[[`, 2))
  LLR2Mat.own <- Get.R2.Matrix.From.List(genePairs, lapply(R2, `[[`, 3))
  
  # Add in the diagonal entries
  for (i in seq_along(geneDataList))
    LLR2Mat[i,i] <- LLR2Mat.other[i,i] <- LLR2Mat.own[i,i] <- 1
  
  # Use gene names as the row and column names of each matrix
  colnames(LLR2Mat) <- rownames(LLR2Mat) <- colnames(LLR2Mat.other) <- 
    rownames(LLR2Mat.other) <- colnames(LLR2Mat.own) <- rownames(LLR2Mat.own) <-
    names(geneDataList)
  
  # Stop timer and print runtime
  print(Sys.time() - startTime)
  
  # Write R^2 similarity matrices to CSV if desired
  if(writeToCSV == TRUE) {
    write.csv(LLR2Mat, ifelse(bayes, "BayesLLR2.csv", "LLR2.csv"))
    write.csv(LLR2Mat.other, ifelse(bayes, "BayesLLR2_Other.csv", "LLR2_Other.csv"))
    write.csv(LLR2Mat.own, ifelse(bayes, "BayesLLR2_Own.csv", "LLR2_Own.csv"))
  }
  
  # Combine results into a list and return
  list(LLR2Mat=LLR2Mat, LLR2Mat.other=LLR2Mat.other, LLR2Mat.own=LLR2Mat.own)
}

# ========================================================================

#' @description 
#' Computes the Bayesian LLR2, LLR2_own, and LLR2_other for a single gene
#' pair via Bayesian regression, given the design matrix X, response
#' vector Y, and prior indicator for the corresponding regression of gene A
#' on gene B. The returned LLR2 value is the maximum of LLR2(gene A, gene B)
#' and LLR2(gene B, gene A); the latter can be computed even with the X and 
#' Y constructed for the former.  
#'
#' @param X n-by-5 design matrix, where n is the number of time points. 
#' This design matrix is defined according to equation 5 in Venkatraman et
#' al. 2021.
#' @param Y n-by-1 response vector, defined according to equation 5. 
#' @param prior 0, 1, or NA depending on whether the two genes have an 
#' unlikely, likely, or unknown relationship according to prior biological
#' information.
#'
#' @return Length-3 list containing the Bayesian LLR2, LLR2_own, and 
#' LLR2_other values. For details of how these three quantities are 
#' calculated, see Section 3 of Venkatraman et al. 2021. The returned LLR2
#' value is the maximum of LLR2(gene A, gene B) and LLR2(gene B, gene A). 
#' Whichever one of these two orderings was used determines the order for 
#' the LLR2_own and LLR2_other calculations.
#' @export
Compute.LLR2.Bayes <- function(X, Y, prior) {
  # Get dimensions of design matrix X
  n <- nrow(X)
  p <- ncol(X)
  
  # Define design matrix and response vector for the other direction
  X2 <- cbind(Y, X[,3], X[,2], X[,4], X[,5])
  Y2 <- X[,1]
  
  # Fit the LLR2 model (equation 4 in our paper) in both directions
  LLR2model.dir1 <- lm.fit(X, Y)
  LLR2model.dir2 <- lm.fit(X2, Y2)
  
  # Set prior mean of regression coefficients
  if(is.na(prior))
    priorMean <- matrix(0, nrow=5, ncol=1)
  else
    priorMean <- matrix(c(prior > 0, prior > 0, 0, 0, 0), ncol=1)
  
  # Compute main LLR2 value, first direction
  LScoefs <- matrix(LLR2model.dir1$coefficients, ncol=1)
  LSfit <- matrix(LLR2model.dir1$fitted.values, ncol=1)
  if(is.na(prior) || prior > 0) {
    sigmaSq <- norm(Y - LSfit, "2")^2 / (n-p)
    g <- ((norm(LSfit - X %*% priorMean,"2")^2 - p*sigmaSq) / (p*sigmaSq))[1]
  } else
    g <- 1
  posteriorMean <- (1/(1+g))*priorMean + (g/(1+g))*LScoefs
  posteriorFit <- X %*% posteriorMean
  LLR2.dir1 <- var(posteriorFit)/(var(posteriorFit) + var(Y - posteriorFit))
  
  # Compute main LLR2 value, second direction
  LScoefs <- matrix(LLR2model.dir2$coefficients, ncol=1)
  LSfit <- matrix(LLR2model.dir2$fitted.values, ncol=1)
  if(is.na(prior) || prior > 0) {
    sigmaSq <- norm(Y2 - LSfit, "2")^2 / (n-p)
    g <- ((norm(LSfit - X2 %*% priorMean,"2")^2 - p*sigmaSq) / (p*sigmaSq))[1]
  } else
    g <- 1
  posteriorMean <- (1/(1+g))*priorMean + (g/(1+g))*LScoefs
  posteriorFit <- X2 %*% posteriorMean
  LLR2.dir2 <- var(posteriorFit)/(var(posteriorFit) + var(Y2 - posteriorFit))
  
  # Fit the two sub-models (equations 18, 19 in our paper) with direction 
  # indicated by the main LLR2 value
  if(LLR2.dir2 > LLR2.dir1) {
    X <- X2
    Y <- Y2
  }
  LLR2model.other <- lm.fit(X[,c(1,2,5)], Y)
  LLR2model.own <- lm.fit(X[,3:5], Y)
  
  # Compute LLR2.other
  LScoefs <- matrix(LLR2model.other$coefficients, ncol=1)
  LSfit <- matrix(LLR2model.other$fitted.values, ncol=1)
  if(is.na(prior) || prior > 0) {
    sigmaSq <- norm(Y - LSfit, "2")^2 / (n-p)
    g <- ((norm(LSfit - X[,c(1,2,5)] %*% priorMean[c(1,2,5),],"2")^2 - p*sigmaSq) / (p*sigmaSq))[1]
  } else
    g <- 1
  posteriorMean <- (1/(1+g))*priorMean[c(1,2,5),] + (g/(1+g))*LScoefs
  posteriorFit <- X[,c(1,2,5)] %*% posteriorMean
  LLR2.other <- var(posteriorFit)/(var(posteriorFit) + var(Y - posteriorFit))
  
  # Compute LLR2.own
  LScoefs <- matrix(LLR2model.own$coefficients, ncol=1)
  LSfit <- matrix(LLR2model.own$fitted.values, ncol=1)
  if(is.na(prior) || prior > 0) {
    sigmaSq <- norm(Y - LSfit, "2")^2 / (n-p)
    g <- ((norm(LSfit,"2")^2 - p*sigmaSq) / (p*sigmaSq))[1]
  } else
    g <- 1
  posteriorMean <- (g/(1+g))*LScoefs
  posteriorFit <- X[,3:5] %*% posteriorMean
  LLR2.own <- var(posteriorFit)/(var(posteriorFit) + var(Y - posteriorFit))
  
  # Combine results into a list and return
  LLR2 <- max(LLR2.dir1, LLR2.dir2)
  list(LLR2, LLR2.other, LLR2.own) 
}

# ========================================================================

#' @description
#' Computes the non-Bayesian LLR2, LLR2_own, and LLR2_other for a single gene
#' pair via ordinary least-squares regression, given the design matrix X 
#' and response vector Y for the corresponding regression of gene A on gene
#' B. The returned LLR2 value is the maximum of LLR2(gene A, gene B) and 
#' LLR2(gene B, gene A); the latter can be computed even with the X and Y
#' constructed for the former.  

#' @param X n-by-5 design matrix, where n is the number of time points. 
#' This design matrix is defined according to equation 5 in Venkatraman et
#' al. 2021.
#' @param Y n-by-1 response vector, defined according to equation 5. 
#'
#' @return Length-3 list containing the non-Bayesian LLR2, LLR2_own, and 
#' LLR2_other values. For details of how these three quantities are 
#' calculated, see Sections 2.1 and 3.5 of Venkatraman et al. 2021. The 
#' returned LLR2 value is the maximum of LLR2(gene A, gene B) and 
#' LLR2(gene B, gene A). Whichever one of these two orderings was used 
#' determines the order for the LLR2_own and LLR2_other calculations.
#' @export
Compute.LLR2.OLS <- function(X, Y) {
  # Define design matrix and response vector for the other direction
  X2 <- cbind(Y, X[,3], X[,2], X[,4], X[,5])
  Y2 <- X[,1]
  
  # Fit the LLR2 model (equation 4 in our paper) in both directions
  LLR2model.dir1 <- lm.fit(X, Y)
  LLR2model.dir2 <- lm.fit(X2, Y2)
  
  # Compute main LLR2 value, first direction
  MSS <- sum(LLR2model.dir1$fitted.values^2)
  RSS <- sum(LLR2model.dir1$residuals^2)
  LLR2.dir1 <- MSS/(MSS + RSS)
  
  # Compute main LLR2 value, second direction
  MSS <- sum(LLR2model.dir2$fitted.values^2)
  RSS <- sum(LLR2model.dir2$residuals^2)
  LLR2.dir2 <- MSS/(MSS + RSS)
  
  # Fit the two sub-models (equations 18, 19 in our paper) with direction 
  # indicated by the main LLR2 value
  if(LLR2.dir2 > LLR2.dir1) {
    X <- X2
    Y <- Y2
  }
  LLR2model.other <- lm.fit(X[,c(1,2,5)], Y)
  LLR2model.own <- lm.fit(X[,3:5], Y)
  
  # Compute LLR2.other
  MSS <- sum(LLR2model.other$fitted.values^2)
  RSS <- sum(LLR2model.other$residuals^2)
  LLR2.other <- MSS/(MSS + RSS) 
  
  # Compute LLR2.own
  MSS <- sum(LLR2model.own$fitted.values^2)
  RSS <- sum(LLR2model.own$residuals^2)
  LLR2.own <- MSS/(MSS + RSS) 
  
  # Combine results into a list and return
  LLR2 <- max(LLR2.dir1, LLR2.dir2)
  list(LLR2, LLR2.other, LLR2.own)
}

# ========================================================================

#' @description
#' Converts a list of pairwise lead-lag R^2 (LLR2) values into a symmetric 
#' matrix. Called within \code{LLR2}; not meant for standalone use.
#'
#' @param genePairs 2-by-(N(N-1)/2) matrix specifying the indices of genes  
#' comprising  each unique pair (first row is gene 1, second row is gene 2).
#' (see \code{LLR2} source for definition).
#' @param R2list Length-(N(N-1)/2) list of the LLR2 values that are to be
#' turned into a symmetric matrix. The ordering should match the ordering 
#' of the columns. 
#'
#' @return N-by-N symmetric matrix of LLR2 values.
#' @export
Get.R2.Matrix.From.List <- function(genePairs, R2list) {
  # Get total number of gene pairs M
  M <- ncol(genePairs)
  
  # Use M to compute total number of genes N, initialize size-N matrix
  N <- (1 + sqrt(1 + 8*M))/2
  R2Matrix <- matrix(0L, N, N)
  
  # Get indices of R^2 list that will fill both triangular halves of the matrix
  allPairs <- cbind(genePairs, genePairs[2:1, ])
  rows <- as.list(allPairs[1, ])
  cols <- as.list(allPairs[2, ])
  
  # Fill in matrix
  Fill.Matrix.Entry <- function(row.i, col.i, R2.i) R2Matrix[row.i, col.i] <<- R2.i
  invisible(Map(Fill.Matrix.Entry, rows, cols, R2list))
  R2Matrix
}
