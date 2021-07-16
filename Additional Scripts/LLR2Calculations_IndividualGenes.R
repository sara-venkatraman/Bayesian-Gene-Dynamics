# --- Compute the non-Bayesian and Bayesian LLR2 values for a specific gene pair ---

# This is in contrast to the script "2_LLR2Calculations.R" in the parent
# directory and the script "LLR2Calculations_NonBayes" in this directory,
# which respectively compute the Bayesian and non-Bayesian (ordinary least
# squares) lead-lag R^2 values for all possible gene pairs. This script is
# used to obtain the entries of Tables 2 and 3 in Section 4.2 of our paper.

# Read the gene expression data
geneData <- read.csv("../Data/geneData.csv", row.names=1)
geneNames <- rownames(geneData)
hours <- c(0, 1, 2, 4, 5, 6, 8, 10, 12, 14, 16, 20, 24, 30, 36, 42, 48)
priorMatrix <- read.csv("../Data/priorMatrix.csv", row.names=1)
colnames(priorMatrix) <- rownames(priorMatrix)

# Input: geneName, the name of a gene as a string
# Output: an natural cubic spline interpolant of the gene's expression data;
# the interpolant is a function that can be evaluated at any real number.
Expression.Profile.Interpolant <- function(geneName) {
  profile <- geneData[geneName,]
  numPoints <- length(profile)
  splinefun(x=hours, y=profile, method="natural")
}  

# Input: gene1 and gene2, the names of two genes as strings.
# Output: design matrix X for the lead-lag model shown in equation 5 of 
# our paper.
Get.Design.Matrix <- function(gene1, gene2) {
  # Initialize design matrix
  X <- matrix(0, nrow=length(hours), ncol=5)
  
  # Get the temporal expression data for the two genes
  profile1 <- geneData[gene1,]
  profile2 <- geneData[gene2,]
  
  # Compute the interpolants for each gene
  interp1 <- Expression.Profile.Interpolant(gene1)
  interp2 <- Expression.Profile.Interpolant(gene2)
  
  # Fill in first row of design matrix manually
  X[1,] <- c(profile2[[1]], 0, 0, hours[1], 1)
  
  # Fill in remaining entries using numerical integration
  for(i in 2:length(hours)) {
    X[i,1] <- profile2[[i]]
    X[i,2] <- integrate(interp2, lower=hours[i-1], upper=hours[i])$value
    X[i,3] <- integrate(interp1, lower=hours[i-1], upper=hours[i])$value
    X[i,4] <- hours[i]
    X[i,5] <- 1
  }
  X
}

# Computes the three non-Bayesian lead-lag R^2 values via ordinary least
# squares regression, given the names of two genes.
# Input: gene1 and gene2, the names of two genes as strings. gene1 is used
# as the "response" gene ("gene A") and is regressed on gene2 ("gene B").
# Output: list containing the LLR2, LLR2.other, and LLR2.own values.
Compute.Gene.Pair.R2 <- function(gene1, gene2) {
  # Get the response vector Y and the design matrix X
  y <- geneData[gene1,]
  y <- t(unname(y))
  x <- Get.Design.Matrix(gene1, gene2)
  
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

# Computes the three Bayesian lead-lag R^2 values via Bayesian regression,
# as described in Section 3 of our paper, given the names of two genes.
# Input: gene1 and gene2, the names of two genes as strings. gene1 is used
# as the "response" gene ("gene A") and is regressed on gene2 ("gene B").
# Output: list containing the Bayesian LLR2, LLR2.other, and LLR2.own values.
Compute.Gene.Pair.R2.Bayes <- function(gene1, gene2, prior) {
  # Get the response vector Y and the design matrix X
  y <- geneData[gene1,]
  y <- t(unname(y))
  x <- Get.Design.Matrix(gene1, gene2)
  
  # Get dimensions of X
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
  
  # Compute main LLR2
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

# --- Tables 2-3: Bayesian LLR2 calculations for immune response/metabolism case study ---

# Define the set of six immune response and metabolic genes
immMetGenes <- c("IM1", "IM2", "FASN1", "UGP", "mino", "fbp")

# Initialize LLR2 and LLR2 - LLR2.own matrices
LLR2Mat <- matrix(1, nrow=length(immMetGenes), ncol=length(immMetGenes))
LLR2DiffMat <- matrix(1, nrow=length(immMetGenes), ncol=length(immMetGenes))
rownames(LLR2Mat) <- colnames(LLR2Mat) <- immMetGenes
rownames(LLR2DiffMat) <- colnames(LLR2DiffMat) <- immMetGenes

# Fill in entries
for(i in 1:length(immMetGenes)) {
  for(j in 1:length(immMetGenes)) {
    if(i != j) {
      gene1 <- immMetGenes[i]
      gene2 <- immMetGenes[j]
      LLR2 <- Compute.Gene.Pair.R2.Bayes(gene1, gene2, priorMatrix[gene1, gene2])
      LLR2Mat[i,j] <- LLR2[[1]]
      LLR2DiffMat[i,j] <- LLR2[[1]] - LLR2[[3]]
    }
  }
}

# Print table 2
round(LLR2Mat, 2)

# Print table 3
round(LLR2DiffMat, 2)
