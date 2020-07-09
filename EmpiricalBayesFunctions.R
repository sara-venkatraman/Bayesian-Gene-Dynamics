# --- Script 3: Function definitions ---

# To use the functions in this script, the following objects need to have already been constructed.
# Examples of how to construct these objects from raw data are in DatasetConstruction.R, and examples
# of how to read these objects from CSV files are in EmpiricalBayesDataAnalysis.R.
# - geneData: A G-by-n dataframe containing measurements (e.g. log-fold change 
#   or normalized counts) of gene expression for G genes over n time points. 
#   Row names of this dataframe are the gene names.
# - geneNames: A length-G vector containing names of each gene in geneData
# - geneIDs: A length-G vector containing Flybase IDs of each gene in geneData
# - genesSubset: A length-N vector of gene names from geneNames, where N <= G, for which the
#   R^2 analysis will be conducted
# - priorMatrix; A N-by-N prior matrix of binary similarities amongst every pair of N genes
#   in genesSubset. Row names of this matrix are the gene names in genesSubset.

hours <- c(0, 1, 2, 4, 5, 6, 8, 10, 12, 14, 16, 20, 24, 30, 36, 42, 48)

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

# Construct the spline interpolant for a gene given its name
Expression.Profile.Interpolant <- function(geneName) {
  profile <- geneData[geneName,]
  numPoints <- length(profile)
  interpolant <- splinefun(x=hours, y=profile, method="natural")
  return(interpolant)
}

# Plot the expression profiles (observed points and spline interpolants)
# of two genes given their names
Plot.Gene.Pair <- function(gene1, gene2) {
  interp1 <- Expression.Profile.Interpolant(gene1)
  interp2 <- Expression.Profile.Interpolant(gene2)
  profile1 <- geneData[gene1,]
  profile2 <- geneData[gene2,]
  ylim1 <- min(min(profile1), min(profile2))
  ylim2 <- max(max(profile1), max(profile2))
  curve(interp1, from=0, to=45, col="red", main=paste(gene1, ", ", gene2, sep=""), xlab="Time", ylab="Expression (log-fold)", ylim=c(ylim1, ylim2), lwd=1.5)
  curve(interp2, from=0, to=45, col="blue", lwd=1.5, add=T)
  points(hours, profile1, pch=20, col="red")
  points(hours, profile2, pch=20, col="blue")
  legend("bottomright", legend=c(gene1, gene2), col=c("red", "blue"), fill=c("red", "blue"), horiz=TRUE, cex=0.6) 
}

# Plot the expression profiles (observed points and spline interpolants)
# of up to 16 genes
library(RColorBrewer)
Plot.Gene.Group <- function(genesToPlot) {
  plotColors <- c(brewer.pal(n=8, name="Dark2"), brewer.pal(n=8, name="Set2"))[1:length(genesToPlot)]
  plotTitle <- strwrap(paste(genesToPlot, collapse=", "))
  plotExtrema <- Time.Profile.Extrema(genesToPlot)
  par(xpd=TRUE, mar = par()$mar + c(0,0,0,1))
  interp <- Expression.Profile.Interpolant(genesToPlot[1])
  profile <- geneData[genesToPlot[1],]
  curve(interp, from=0, to=48, col=plotColors[1], main=plotTitle, xlab="Time", ylab="Expression (log-fold)", ylim=c(plotExtrema$min, plotExtrema$max), lwd=1.5)
  points(hours, profile, pch=20, col=plotColors[1])
  for(i in 2:length(genesToPlot)) {
    interp <- Expression.Profile.Interpolant(genesToPlot[i])
    profile <- geneData[genesToPlot[i],]
    curve(interp, from=0, to=48, col=plotColors[i], add=T, xlab="Time", ylab="Expression (log-fold)", lwd=1.5)
    points(hours, profile, pch=20, col=plotColors[i])
  }
  legend("topright", legend=genesToPlot, col=plotColors, fill=plotColors, cex=0.6, xpd=TRUE, inset=c(-.18,0)) 
  par(xpd=FALSE, mar = par()$mar-c(0,0,0,1))
}

# For plotting purposes only: returns the largest and smallest expression
# levels observed over a set of genes
Time.Profile.Extrema <- function(genesToPlot) {
  minima <- c();  maxima <- c()
  for(i in 1:length(genesToPlot)) {
    profile <- geneData[genesToPlot[i],]
    minima[i] <- min(profile);  maxima[i] <- max(profile)
  }
  return(list(min=min(minima), max=max(maxima)))
}

# Construct the design matrix for the lead-lag R^2 model given two gene names
Construct.Design.Matrix <- function(gene1, gene2) {
  # Initialize design matrix
  X <- matrix(0, nrow=length(hours), ncol=5)
  
  # Get vectors of log-fold change expression measurements for the two genes
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
  return(X)
}

# Given a gene name, return the time profile of that gene as a column vector
Get.Gene.Data.As.Vector <- function(geneName) {
  Y <- geneData[geneName,]
  Y <- t(unname(Y))
  return(Y)
}

# Given two gene names, compute the following three (Bayesian) R^2 values:
# 1. The R^2 from regressing gene1 only on gene2 (simultaneous R^2 from simple
#    linear regression)
# 2. The lead-lag R^2 from regressing gene1 on gene2 and itself using an empirical
#    Bayes approach (see writeup for more details -- regression coefficients are 
#    taken to be the posterior mean from a normal-inverse gamma model with the g-prior)
# 3. The R^2 from regressing gene1 only on time and its own time integrals
# If g is not provided (default), compute the optimal g
Compute.Gene.Pair.R2.Bayes <- function(gene1, gene2, g, priorMean) {
  # Get the response vector and the design matrix for lead-lag regression (model 2)
  response <- Get.Gene.Data.As.Vector(gene1)
  X <- Construct.Design.Matrix(gene1, gene2)
  
  # If prior information is available for this gene pair, choose non-zero prior mean.
  # Otherwise, choose prior mean of zero.
  if(missing(priorMean)) {
    if(gene1 %in% rownames(priorMatrix) && gene2 %in% rownames(priorMatrix)) {
      if(priorMatrix[gene1, gene2] >= 1) { priorMean <- c(1,1,0,0,0) }
      else { priorMean <- c(0,0,0,0,0) }
    }
    else { priorMean <- c(0,0,0,0,0) }
  }
  
  # If g is not specified, set g to its optimal value
  if(missing(g)) {
    g <- Optimize.g(X, response, priorMean)
  }
  
  # For the first R^2 value: set up design matrix (first and last column of full
  # design matrix), prior mean (first and last entry of full prior mean), posterior mean
  # of Model 1 coefficients, least-squares estimates for Model 1, and fitted values
  XModel1 <- X[,c(1,5)]
  priorMeanModel1 <- priorMean[c(1,5)]
  leastSquaresModel1 <- solve(t(XModel1) %*% XModel1) %*% t(XModel1) %*% response
  posteriorMeanModel1 <- (1/(1+g))*priorMeanModel1 + (g/(1+g))*leastSquaresModel1
  fitModel1 <- XModel1 %*% posteriorMeanModel1
  R2Model1 <- var(fitModel1)/(var(fitModel1) + var(response-fitModel1))
  
  # For the second R^2 value: compute least-squares estimates for Model 2, posterior mean
  # of Model 2 coefficients, and fitted values
  leastSquaresModel2 <- solve(t(X) %*% X) %*% t(X) %*% response
  posteriorMeanModel2 <- (1/(1+g))*priorMean + (g/(1+g))*leastSquaresModel2
  fitModel2 <- X %*% posteriorMeanModel2
  R2Model2 <- var(fitModel2)/(var(fitModel2) + var(response-fitModel2))
  
  # For the third R^2 value: set up design matrix (last three columns of full
  # design matrix), prior mean (last three entries of full prior mean), posterior mean
  # of Model 3 coefficients, least-squares estimates for Model 3, and fitted values
  XModel3 <- X[,3:5]
  priorMeanModel3 <- priorMean[3:5]
  leastSquaresModel3 <- solve(t(XModel3) %*% XModel3) %*% t(XModel3) %*% response
  posteriorMeanModel3 <- (1/(1+g))*priorMeanModel3 + (g/(1+g))*leastSquaresModel3
  fitModel3 <- XModel3 %*% posteriorMeanModel3
  R2Model3 <- var(fitModel3)/(var(fitModel3) + var(response-fitModel3))
  
  # Return the R^2 values
  return(list(R2Model1=R2Model1, R2Model2=R2Model2, R2Model3=R2Model3, g=g))
}

# Given two gene names, compute the following three R^2 values:
# 1. The R^2 from regressing gene1 only on gene2 (simultaneous R^2 from simple
#    linear regression)
# 2. The lead-lag R^2 from regressing gene1 on gene2 and itself (see 
#    writeup for more details)
# 3. The R^2 from regressing gene1 only on time and its own time integrals
Compute.Gene.Pair.R2 <- function(gene1, gene2) {
  # Get the response vector and the design matrix for lead-lag regression (model 2)
  response <- geneData[gene1,]
  response <- t(unname(response))
  X <- Construct.Design.Matrix(gene1, gene2)
  
  # Compute the R^2 from model 1
  model1 <- lm(response ~ X[,c(1,5)] + 0)
  R2Model1 <- summary(model1)$r.squared
  
  # Compute the R^2 from model 2
  model2 <- lm(response ~ X + 0)
  R2Model2 <- summary(model2)$r.squared
  
  # Compute the R^2 from model 3
  model3 <- lm(response ~ X[,3:5] + 0)
  R2Model3 <- summary(model3)$r.squared
  
  # Return the R^2 values
  return(list(R2Model1=R2Model1, R2Model2=R2Model2, R2Model3=R2Model3, g=NA))
}

# Given a vector of gene names, construct the following three matrices:
# (see previous two functions for definitions of the three models)
# Matrix 1: (i,j) is the model 1 R^2 between genes i and j
# Matrix 2: (i,j) is the model 2 R^2 (lead-lag) between genes i and j
# Matrix 3: (i,j) is the model 3 R^2 for gene i
# These R^2 values can be computed with or without the Bayesian approach.
# While the three R^2 metrics are not symmetric, the matrices are made to be
# symmetric depending on the value of the model 2 R^2 (see implementation).
# Additionally return a matrix whose (i,j) value is the value of g used in 
# the g-prior for the regression between genes i and j if bayes=TRUE. Otherwise
# all entries of the matrix are 1.
Compute.R2.Matrices <- function(genesSubset, bayes=TRUE) {
  # Initialize the similarity matrices and set row/column names
  subsetSize <- length(genesSubset)
  matrix1 <- matrix(1, nrow=subsetSize, ncol=subsetSize)
  matrix2 <- matrix(1, nrow=subsetSize, ncol=subsetSize)
  matrix3 <- matrix(1, nrow=subsetSize, ncol=subsetSize)
  gMatrix <- matrix(1, nrow=subsetSize, ncol=subsetSize)
  
  # Set row/column names
  rownames(matrix1) <- genesSubset; colnames(matrix1) <- genesSubset
  rownames(matrix2) <- genesSubset; colnames(matrix2) <- genesSubset
  rownames(matrix3) <- genesSubset; colnames(matrix3) <- genesSubset
  rownames(gMatrix) <- genesSubset; colnames(gMatrix) <- genesSubset

  # Fill in the off-diagonal entries (i,j) with R^2 value between genes i and j
  for(i in 1:subsetSize) {
    for(j in i:subsetSize) {
      if(i != j) {
        # Compute similarity between two genes (both directions of dependence)
        gene1 <- genesSubset[i];  gene2 <- genesSubset[j]
        
        if(bayes == TRUE) {
          # Get design matrices and response vectors for both orderings of the gene pair
          X1 <- Construct.Design.Matrix(gene1, gene2)
          X2 <- Construct.Design.Matrix(gene2, gene1)
          Y1 <- Get.Gene.Data.As.Vector(gene1)
          Y2 <- Get.Gene.Data.As.Vector(gene2)
          
          # Find optimal value of g for both orderings of the gene pair
          if(gene1 %in% rownames(priorMatrix) && gene2 %in% rownames(priorMatrix)) {
            if(priorMatrix[gene1, gene2] == 1) { priorMean <- c(1,1,0,0,0) }
            else { priorMean <- c(0,0,0,0,0) }
          }
          else { priorMean <- c(0,0,0,0,0) }
          
          # Compute R^2 values for both orderings of the gene pair
          R2.1 <- Compute.Gene.Pair.R2.Bayes(gene1, gene2)
          R2.2 <- Compute.Gene.Pair.R2.Bayes(gene2, gene1)
        }
        else {
          R2.1 <- Compute.Gene.Pair.R2(gene1, gene2)
          R2.2 <- Compute.Gene.Pair.R2(gene2, gene1)
          # g1 <- 1;  g2 <- 1;
        }
        
        if(R2.1$R2Model2-R2.1$R2Model3 > R2.2$R2Model2-R2.2$R2Model3) {
          R2Model1 <- R2.1$R2Model1
          R2Model2 <- R2.1$R2Model2
          R2Model3 <- R2.1$R2Model3
          g <- R2.1$g
          # g <- g1
        } else {
          R2Model1 <- R2.2$R2Model1
          R2Model2 <- R2.2$R2Model2
          R2Model3 <- R2.2$R2Model3
          g <- R2.2$g
          # g <- g2
        }
        
        # Fill in upper-triangular entries of similarity matrices
        matrix1[i,j] <- R2Model1
        matrix2[i,j] <- R2Model2
        matrix3[i,j] <- R2Model3
        gMatrix[i,j] <- g

        # Fill in lower-triangular entries of similarity matrices
        matrix1[j,i] <- R2Model1
        matrix2[j,i] <- R2Model2
        matrix3[j,i] <- R2Model3
        gMatrix[j,i] <- g
      }
    }
  }
  return(list(matrix1=matrix1, matrix2=matrix2, matrix3=matrix3, gMatrix=gMatrix))
}

# Given a vector of gene names, and the output of the Compute.R2.Matrices function,
# turn the upper-triangular halves of each R^2 matrix into a vector (this is 
# primarily for use with ggplot). Row names are the gene pairs. Optional parameter 
# for prior separation: if TRUE, the returned "priorVector" will contain "0" if the 
# corresponding entry of the prior matrix = 0, "1A" if the STRING entry = 1, "1B" if
# the replicate entry = 1, and "2" if both the STRING and replicate entries = 1.
Vectorize.R2.Matrices <- function(genesSubset, R2Matrices, separatePriors=FALSE) {
  matrix1 <- R2Matrices$matrix1
  matrix2 <- R2Matrices$matrix2
  matrix3 <- R2Matrices$matrix3
  gMatrix <- R2Matrices$gMatrix
  
  subsetSize <- length(genesSubset)
  vector1 <- matrix(0, nrow=sum(1:(subsetSize-1)))
  vector2 <- matrix(0, nrow=sum(1:(subsetSize-1)))
  vector3 <- matrix(0, nrow=sum(1:(subsetSize-1)))
  priorVector <- matrix(0, nrow=sum(1:(subsetSize-1)))
  gVector <- matrix(0, nrow=sum(1:(subsetSize-1)))
  pairLabels <- c();  k <- 1
  for(i in 1:(subsetSize-1)) {
    for(j in (i+1):subsetSize) {
      vector1[k] <- matrix1[i,j]
      vector2[k] <- matrix2[i,j]
      vector3[k] <- matrix3[i,j]
      gVector[k] <- gMatrix[i,j]
      if(separatePriors == TRUE) {
        if(priorMatrix[i,j] == 0) { priorVector[k] <- 0 }
        else if(priorMatrix[i,j] == 2) { priorVector[k] <- 2 }
        else if(priorMatrixString[i,j] == 1) { priorVector[k] <- 0.67 }
        else if(priorMatrixReplicates[i,j] == 1) { priorVector[k] <- 1.33 }
      } 
      else { priorVector[k] <- priorMatrix[i,j] }
      pairLabels[k] <-  paste(rownames(matrix1)[i], ", ", colnames(matrix1)[j], sep="")
      k <- k+1
    } 
  }
  rownames(vector1) <- pairLabels
  rownames(vector2) <- pairLabels
  rownames(vector3) <- pairLabels
  rownames(priorVector) <- pairLabels
  rownames(gVector) <- pairLabels
  return(list(vector1=vector1, vector2=vector2, vector3=vector3, priorVector=priorVector, gVector=gVector))
}

# Given a vector of gene names and the output of the Compute.R2.Matrices function,
# draw an interactive scatterplot of the model 1 R^2 against the difference of the 
# R^2 values from models 2 and 3. Color points differently if prior information
# is available for the corresponding gene pairs. Optional argument to specify whether
# points should be colored according to whether zero, one, or two prior sources 
# indicated an association.
library(plotly);  library(ggplot2)
Draw.Metric.Scatterplot.For.Binary.Prior <- function(R2Matrices, bayes, colorPriors=FALSE) {
  R2Vectors <- Vectorize.R2.Matrices(genesSubset, R2Matrices, separatePriors=colorPriors)
  vec1 <- R2Vectors$vector1
  vec2 <- R2Vectors$vector2
  vec3 <- R2Vectors$vector3
  gVec <- R2Vectors$gVector
  
  if(colorPriors == TRUE) { priorVec <- R2Vectors$priorVector } 
  else { priorVec <- R2Vectors$priorVector > 0 }
  
  plotData <- as.data.frame(cbind(round(vec1, 4), round(abs(vec2-vec3), 4), priorVec))
  colnames(plotData) <- c("x.axis", "y.axis", "prior")
  plotData$prior <- as.factor(plotData$prior)
  
  if(bayes == TRUE) {
    hoverText <- paste(row.names(plotData), paste(rep("; g = ", nrow(vec1)), round(gVec, 3), sep=""), sep="")
  } else {
    hoverText <- row.names(plotData)
  }
  
  p <- ggplot(plotData, aes(x=x.axis, y=y.axis, color=prior, text=hoverText)) + 
    geom_point(shape=1, size=0.9) + xlab('Model 1 R^2') + ylab('Difference between model 2 and model 3 R^2') + 
    ggtitle("Comparison of Empirical Bayes R^2 Values") + theme_light() + theme(legend.position="none")
  if(colorPriors == TRUE) {
    p <- p + scale_color_manual(values=c("navy","orangered3","goldenrod1","forestgreen"))
  }
  else {
    p <- p + scale_color_manual(values=c("navy","orangered3"))
  }
  if(bayes == TRUE) { p <- p + ggtitle("Comparison of Empirical Bayes R^2 Values") }
  else { p <- p + ggtitle("Comparison of R^2 Values (Non-Bayesian)")}
  ggplotly(p)
}

# Given a design matrix X and a response vector Y, find the value of g which minimizes
# the empirical risk of the fitted values obtained from Bayesian regression of Y on X 
# with a g-prior placed on the regression coefficients with prior mean mu0
Optimize.g <- function(X, Y, mu0) {
  H <- X %*% solve(t(X) %*% X) %*% t(X)
  Yh <- H %*% Y
  objective <- function(g) {
    term1 <- 1/((1+g)^2) * (2*t(Y) %*% X %*% mu0 + 2*t(mu0) %*% t(X) %*% Yh)
    term2 <- -1/(1+g) * t(Y) %*% Yh
    term3 <- g/((1+g)^2) * (t(Y) %*% Yh + 2*t(Yh) %*% Yh)
    term4 <- -2/((1+g)^3) * t(mu0) %*% t(X) %*% X %*% mu0
    term5 <- -2*g/((1+g)^3) * (2*t(mu0) %*% t(X) %*% Yh)
    term6 <- -2*(g^2)/((1+g)^3) %*% t(Yh) %*% Yh
    result <- (1/length(Y)) * (term1 + term2 + term3 + term4 + term5 + term6)
    return(result)
  }
  g <- tryCatch({ uniroot(objective, interval=c(0.01, 10))$root }, error=function(cond) { 1 })
  return(g)
}



# # Given two gene names, compute the following three (Bayesian) R^2 values:
# # 1. The R^2 from regressing gene1 only on gene2 (simultaneous R^2 from simple
# #    linear regression)
# # 2. The lead-lag R^2 from regressing gene1 on gene2 and itself using an empirical
# #    Bayes approach (see writeup for more details -- regression coefficients are 
# #    taken to be the posterior mean from a normal-inverse gamma model)
# # 3. The R^2 from regressing gene1 only on time and its own time integrals
# Compute.Gene.Pair.R2.Bayes <- function(gene1, gene2, inversePriorCovMatrix) {
#   # Get the response vector and the design matrix for lead-lag regression (model 2)
#   response <- Get.Gene.Data.As.Vector(gene1)
#   X <- Construct.Design.Matrix(gene1, gene2)
#   
#   # If inverse prior covariance matrix is not provided, set it to X'X (this is equivalent
#   # to using the Zellner prior with g=1)
#   if(missing(inversePriorCovMatrix)) { inversePriorCovMatrix <- t(X) %*% X }
#   
#   # If prior information is available for this gene pair, choose non-zero prior mean.
#   # Otherwise, choose prior mean of zero.
#   if(gene1 %in% rownames(priorMatrix) && gene2 %in% rownames(priorMatrix)) {
#     if(priorMatrix[gene1, gene2] == 1) { priorMean <- c(1,1,0,0,0) }
#     else { priorMean <- c(0,0,0,0,0) }
#   }
#   else { priorMean <- c(0,0,0,0,0) }
#   
#   # For the first R^2 value: set up design matrix (first and last column of full
#   # design matrix), inverse of the prior covariance matrix (the four corner entries
#   # of the full inverse prior covariance), prior mean (first and last entry of full
#   # prior mean), posterior mean, fitted values
#   XModel1 <- X[,c(1,5)]
#   inversePriorCovModel1 <- inversePriorCovMatrix[c(1,5),c(1,5)]
#   priorMeanModel1 <- priorMean[c(1,5)]
#   posteriorCovModel1 <- solve(t(XModel1) %*% XModel1 + inversePriorCovModel1)
#   posteriorMeanModel1 <- posteriorCovModel1 %*% (inversePriorCovModel1 %*% priorMeanModel1 + t(XModel1) %*% response)
#   fitModel1 <- XModel1 %*% posteriorMeanModel1
#   R2Model1 <- var(fitModel1)/(var(fitModel1) + var(response-fitModel1))
#   
#   # For the second R^2 value: compute posterior mean of the vector of five  
#   # regression coefficients, fitted values, and lead-lag R^2
#   posteriorCovModel2 <- solve(t(X) %*% X + inversePriorCovMatrix)
#   posteriorMeanModel2 <- posteriorCovModel2 %*% (inversePriorCovMatrix %*% priorMean + t(X) %*% response)
#   fitModel2 <- X %*% posteriorMeanModel2
#   R2Model2 <- var(fitModel2)/(var(fitModel2) + var(response-fitModel2))
#   
#   # For the third R^2 value: set up design matrix (last three columns of full
#   # design matrix), inverse of the prior covariance matrix (lower-right 3x3 block 
#   # of full inverse prior covariance), prior mean (last three entries of full
#   # prior mean), posterior mean, fitted values
#   XModel3 <- X[,3:5]
#   inversePriorCovModel3 <- inversePriorCovMatrix[3:5,3:5]
#   priorMeanModel3 <- priorMean[3:5]
#   posteriorCovModel3 <- solve(t(XModel3) %*% XModel3 + inversePriorCovModel3)
#   posteriorMeanModel3 <- posteriorCovModel3 %*% (inversePriorCovModel3 %*% priorMeanModel3 + t(XModel3) %*% response)
#   fitModel3 <- XModel3 %*% posteriorMeanModel3
#   R2Model3 <- var(fitModel3)/(var(fitModel3) + var(response-fitModel3))
#   
#   # Return the R^2 values
#   return(list(R2Model1=R2Model1, R2Model2=R2Model2, R2Model3=R2Model3))
# }
