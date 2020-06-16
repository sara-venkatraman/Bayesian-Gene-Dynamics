# --- Script 5: Exploratory data analysis ---

# Load datasets needed for constructing priors
source("DatasetLoader.R")

# Load functions in EmpiricalBayesFunctions.R
source("EmpiricalBayesFunctions.R")

# Before running code below, also need to run the relevant section of PriorConstruction.R.

# --- Initial plots --- 

# There seems to be a lead-lag pattern here
Plot.Gene.Pair("Galk", "CG12158") # prior = 0
Compute.Gene.Pair.R2.Bayes("Galk", "CG12158") # coefficients shrink towards 0
Compute.Gene.Pair.R2("Galk", "CG12158")

# Try reversing the response gene and predictor gene from previous test
Compute.Gene.Pair.R2.Bayes("CG12158", "Galk")
Compute.Gene.Pair.R2("CG12158", "Galk")

# Prior = 1 for this pair
Plot.Gene.Pair("cDIP", "CG12158")
Compute.Gene.Pair.R2.Bayes("cDIP", "CG12158")
Compute.Gene.Pair.R2("cDIP", "CG12158")

# Two immunity genes, no prior available
Plot.Gene.Pair("IM14", "IM2")
Compute.Gene.Pair.R2.Bayes("IM14", "IM2")
Compute.Gene.Pair.R2("IM14", "IM2")

# Another pair of immunity genes, no prior available
Plot.Gene.Pair("IM4", "IM2")
Compute.Gene.Pair.R2.Bayes("IM4", "IM2")
Compute.Gene.Pair.R2("IM4", "IM2")

# The pattern here is inverted; not sure how to set g here to get Bayes LLR2 
# larger. The R^2 from model1 and model3 are very close.
Plot.Gene.Pair("IM4", "GstD7") # prior = 1
X <- Construct.Design.Matrix("IM4", "GstD7")
Compute.Gene.Pair.R2.Bayes("IM4", "GstD7", 3 * solve(t(X) %*% X))
Compute.Gene.Pair.R2("IM4", "GstD7")

# Try reversing the response gene and predictor gene from previous test
X <- Construct.Design.Matrix("GstD7", "IM4")
Compute.Gene.Pair.R2.Bayes("GstD7", "IM4", 3 * solve(t(X) %*% X))
Compute.Gene.Pair.R2("GstD7", "IM4")

# These genes have very similar profiles; use small g to increase LLR2
Plot.Gene.Pair("IM4", "IMPPP") # prior = 1
X <- Construct.Design.Matrix("IM4", "IMPPP")
Compute.Gene.Pair.R2.Bayes("IM4", "IMPPP", 3 * solve(t(X) %*% X))
Compute.Gene.Pair.R2("IM4", "IMPPP")

# Try reversing the response gene and predictor gene from previous test
X <- Construct.Design.Matrix("IMPPP", "IM4")
Compute.Gene.Pair.R2.Bayes("IMPPP", "IM4", 3 * solve(t(X) %*% X))
Compute.Gene.Pair.R2("IMPPP", "IM4")

# --- Construct scatterplots, using STRING-DB prior information ---

# To do: use Bayes estimates for reference model parameters as well. Ideally,
# genes that are very similar should show large difference in LLR2 and reference
# R2, indicating that adding a co-regulated gene to the model improves similarity
# metric significantly (i.e. want to push model2 and model3 R^2 further apart)

bayesR2Matrices <- Compute.R2.Matrices(Flybase.ID.To.Gene.Name(genesSubset), bayes=T)
bayesMatrix1 <- bayesR2Matrices$matrix1
bayesMatrix2 <- bayesR2Matrices$matrix2
bayesMatrix3 <- bayesR2Matrices$matrix3

# Vectorize upper-triangular halves of the R^2 similarity matrices - rownames
# should be the corresponding pair of genes
bayesMatrix1Vec <- matrix(0, nrow=sum(1:(N-1)))
bayesMatrix2Vec <- matrix(0, nrow=sum(1:(N-1)))
bayesMatrix3Vec <- matrix(0, nrow=sum(1:(N-1)))
priorMatrixVec <- matrix(0, nrow=sum(1:(N-1)))
pairLabels <- c();  k <- 1
for(i in 1:(N-1)) {
  for(j in (i+1):N) {
    bayesMatrix1Vec[k] <- bayesMatrix1[i,j]
    bayesMatrix2Vec[k] <- bayesMatrix2[i,j]
    bayesMatrix3Vec[k] <- bayesMatrix3[i,j]
    priorMatrixVec[k] <- priorMatrix[i,j]
    pairLabels[k] <-  paste(rownames(bayesMatrix1)[i], ", ", colnames(bayesMatrix1)[j], sep="")
    k <- k+1
  } 
}
rownames(bayesMatrix1Vec) <- pairLabels; bayesMatrix1Vec <- round(bayesMatrix1Vec, 5)
rownames(bayesMatrix2Vec) <- pairLabels; bayesMatrix2Vec <- round(bayesMatrix2Vec, 5)
rownames(bayesMatrix3Vec) <- pairLabels; bayesMatrix3Vec <- round(bayesMatrix3Vec, 5)
rownames(priorMatrixVec) <- pairLabels

# Scatterplot of metric1 vs. metric2-metric3, with points colored by whether or
# not prior information was available
library(plotly);  library(ggplot2)
plot1data <- as.data.frame(cbind(bayesMatrix1Vec, abs(bayesMatrix2Vec-bayesMatrix3Vec), priorMatrixVec))
colnames(plot1data) <- c("Model1", "Model2Model3Diff", "prior1")
plot1data$prior1 <- as.character(plot1data$prior1)
p1 <- ggplot(plot1data, aes(Model1, Model2Model3Diff, color=prior1, text=row.names(plot1data))) + 
  geom_point(size=0.9) + xlab('Model 1 R^2') + 
  ylab('Difference between model 2 and model 3 R^2') + 
  ggtitle("Comparison of Empirical Bayes R^2 Values") +
  theme_light() + theme(legend.position="none") +
  scale_color_manual(values=c("navy","orangered3"))
ggplotly(p1)

# Re-do the above scatterplot without Bayesian computations

nonBayesR2Matrices <- Compute.R2.Matrices(Flybase.ID.To.Gene.Name(genesSubset), bayes=F)
nonBayesMatrix1 <- nonBayesR2Matrices$matrix1
nonBayesMatrix2 <- nonBayesR2Matrices$matrix2
nonBayesMatrix3 <- nonBayesR2Matrices$matrix3

# Vectorize upper-triangular halves of the R^2 similarity matrices - rownames
# should be the corresponding pair of genes
nonBayesMatrix1Vec <- matrix(0, nrow=sum(1:(N-1)))
nonBayesMatrix2Vec <- matrix(0, nrow=sum(1:(N-1)))
nonBayesMatrix3Vec <- matrix(0, nrow=sum(1:(N-1)))
pairLabels <- c();  k <- 1
for(i in 1:(N-1)) {
  for(j in (i+1):N) {
    nonBayesMatrix1Vec[k] <- nonBayesMatrix1[i,j]
    nonBayesMatrix2Vec[k] <- nonBayesMatrix2[i,j]
    nonBayesMatrix3Vec[k] <- nonBayesMatrix3[i,j]
    pairLabels[k] <-  paste(rownames(nonBayesMatrix1)[i], ", ", colnames(nonBayesMatrix1)[j], sep="")
    k <- k+1
  } 
}
rownames(nonBayesMatrix1Vec) <- pairLabels; nonBayesMatrix1Vec <- round(nonBayesMatrix1Vec, 5)
rownames(nonBayesMatrix2Vec) <- pairLabels; nonBayesMatrix2Vec <- round(nonBayesMatrix2Vec, 5)
rownames(nonBayesMatrix3Vec) <- pairLabels; nonBayesMatrix3Vec <- round(nonBayesMatrix3Vec, 5)

plot2data <- as.data.frame(cbind(nonBayesMatrix1Vec, nonBayesMatrix2Vec-nonBayesMatrix3Vec, priorMatrixVec))
colnames(plot2data) <- c("Model1", "Model2Model3Diff", "prior1")
plot2data$prior1 <- as.character(plot2data$prior1)
p2 <- ggplot(plot2data, aes(Model1, Model2Model3Diff, color=prior1, text=row.names(plot2data))) + 
  geom_point(size=0.9) + xlab('Model 1 R^2') + 
  ylab('Difference between model 2 and model 3 R^2') + 
  ggtitle("Comparison of (non-Bayesian) R^2 Values") + 
  theme_light() + theme(legend.position="none") + 
  scale_color_manual(values=c("navy","orangered3"))
ggplotly(p2)

# First scatterplot shows that prior information may have been weighted too heavily.
# Gene pairs with prior = 1 do not really look similar.

Plot.Gene.Group(c("CG32751", "per", "Cbs", "Cpr62Ba", "salt"))

par(mfrow=c(1,2))
hist(bayesMatrix1Vec, col="lightblue", main="Model 1 R^2 (Bayes)", xlab="R^2")
hist(nonBayesMatrix1Vec, col="lightblue", main="Model 1 R^2 (OLS)", xlab="R^2")

par(mfrow=c(1,2))
hist(bayesMatrix2Vec, col="forestgreen", main="Model 2 R^2 (Bayes)", xlab="R^2")
hist(nonBayesMatrix2Vec, col="forestgreen", main="Model 2 R^2 (OLS)", xlab="R^2")

par(mfrow=c(1,2))
hist(bayesMatrix3Vec, col="goldenrod", main="Model 2 R^2 (Bayes)", xlab="R^2")
hist(nonBayesMatrix3Vec, col="goldenrod", main="Model 2 R^2 (OLS)", xlab="R^2")

# --- Construct scatterplots, using data-driven prior information ---

# Some interesting plots from this prior
Plot.Gene.Group(c("IM14","IM1","IM2","CG5791", "Drs")) # This group was identified using the priors
Plot.Gene.Pair("Slob", "per")
nonBayesR2Matrices <- Compute.R2.Matrices(genesSubset, bayes=FALSE)
bayesR2Matrices <- Compute.R2.Matrices(genesSubset, bayes=TRUE)
Draw.Metric.Scatterplot.For.Binary.Prior(nonBayesR2Matrices, bayes=FALSE)
Draw.Metric.Scatterplot.For.Binary.Prior(bayesR2Matrices, bayes=TRUE)

# --- Hierarchical clustering using empirical Bayes R^2 ---

hierClust <- hclust(as.dist(1-nonBayesR2Matrices$matrix2), method="ward.D")
clusters <- cutree(hierClust, k=15)
table(clusters)
plot(hierClust, cex=0.3);  rect.hclust(hierClust, k=15, border="red")
Plot.Gene.Group((genesSubset[clusters == 9]))

# --- Analysis on simulated data: matrices/vectors formatted as required by EmpiricalBayesFunctions.R
#     How to use this code: first run one of the following lines, then run all the code in the subsequent part ---

# Simulation: two genes that are visually similar (one is translation of the other), prior=1
gene1Profile <- sin(0.8*hours);  gene2Profile <- gene1Profile-0.5 + rnorm(length(hours), mean=0, sd=0.1);  prior <- 1
# Simulation: same time profiles as above, but prior=0
gene1Profile <- sin(0.8*hours);  gene2Profile <- gene1Profile-0.5 + rnorm(length(hours), mean=0, sd=0.1);  prior <- 0
# Simulation: two genes that are visually dissimilar, prior=0
gene1Profile <- sin(0.8*hours);  gene2Profile <- exp(-((hours-10)^2)/8);  prior <- 1
# Simulation: two genes that are visually dissimilar, prior=1
gene1Profile <- sin(0.8*hours);  gene2Profile <- exp(-((hours-10)^2)/8);  prior <- 0

# Draw plot of simulated profiles, construct prior matrix, compute non-Bayesian and Bayesian R^2 values 
geneData <- as.data.frame(rbind(gene1Profile, gene2Profile)); rownames(geneData) <- c("gene1","gene2");
priorMatrix <- matrix(c(1,prior,prior,1), nrow=2, ncol=2); rownames(priorMatrix) <- c("gene1","gene2"); colnames(priorMatrix) <- c("gene1","gene2")
Plot.Gene.Pair("gene1","gene2"); Compute.Gene.Pair.R2("gene1","gene2"); Compute.Gene.Pair.R2.Bayes("gene1","gene2")

# --- Signal retention ---

# How many gene pairs are in each quadrant of the plot, and what do the pair plots look like for those in the upper-right quadrant?
bayesR2MatricesVec <- Vectorize.R2.Matrices(genesSubset, bayesR2Matrices)
bayesVec1 <- bayesR2MatricesVec$vector1;  bayesVec2 <- bayesR2MatricesVec$vector2;  bayesVec3 <- bayesR2MatricesVec$vector3
upperRight <- rownames(bayesVec1)[bayesVec1 > 0.7 & bayesVec2-bayesVec3 > 0.4]
lowerRight <- rownames(bayesVec1)[bayesVec1 > 0.7 & bayesVec2-bayesVec3 <= 0.4] 
lowerLeft <- rownames(bayesVec1)[bayesVec1 <= 0.7 & bayesVec2-bayesVec3 <= 0.4] 
upperLeft <- rownames(bayesVec1)[bayesVec1 <= 0.7 & bayesVec2-bayesVec3 > 0.4]  # upper-left quadrant - interesting patterns here, but not always close
lapply(list(upperRight, lowerRight, lowerLeft, upperLeft), function(x) { return(length(x)/length(bayesVec1)) }) # Print the above in proportions

# With 951 differentially-expressed genes, there are choose(951,2)=451725 unique pairs of genes. So roughly 3475 (= 0.007 * 451725) 
# pairs are "retained" in the upper-right quadrant of the plot.

# Plot all the gene pairs that are in the "high X, high Y" category
pdf("../Graphs/GenePairsUpperRightQuadrant.pdf", height=9, width=11)
for(i in seq(1, length(upperRight), by=9)) {
  par(mfrow=c(3,3), mai=c(0.6,0.6,0.6,0.6))
  for(j in 1:9) {
    if(i+(j-1) <= length(upperRight)) {
      genePair <- strsplit(upperRight[i+(j-1)], ", ")[[1]]
      Plot.Gene.Pair(genePair[1], genePair[2]) 
    }
  }
}
dev.off(); par(mfrow=c(1,1))

# Try the above again, but reduce the search space further by requiring higher values on X and Y axes. This is now 0.019% of all gene pairs.
highXhighYStrong <- rownames(bayesVec1)[bayesVec1 > 0.7 & bayesVec2-bayesVec3 > 0.4] # upper-right quadrant
pdf("../Graphs/GenePairsUpperRightQuadrantOptimizedStrong.pdf", height=9, width=11)
for(i in seq(1, length(highXhighYStrong), by=9)) {
  par(mfrow=c(3,3), mai=c(0.6,0.6,0.6,0.6))
  for(j in 1:9) {
    genePair <- strsplit(highXhighYStrong[i+(j-1)], ", ")[[1]]
    Plot.Gene.Pair(genePair[1], genePair[2])
  }
}
dev.off();  par(mfrow=c(1,1))

# Some interesting plots that come out of the above PDF:
Plot.Gene.Group(c("Sodh-1","UGP","glob1","Gip","Pgi", "GstE9", "CG6084"))
Plot.Gene.Group(c("vri","Slob","per","CG33511","Pdp1"))
Plot.Gene.Pair("CG33511","vri")  # In STRING-DB, no annotation available for CG33511. vri (clock-controlled) affects hair and cell growth

# Try hierarchical clustering on the genes that appeared in highXhighYStrong
highXhighYStrongGeneNames <- unique(strsplit(paste(highXhighYStrong, collapse=", "), ", ")[[1]])
bayesMatrix2SubsetStrong <- bayesR2Matrices$matrix2[highXhighYStrongGeneNames, highXhighYStrongGeneNames]

hierClustStrong <- hclust(as.dist(1-bayesMatrix2SubsetStrong), method="ward.D")
clustersStrong <- cutree(hierClustStrong, k=10);  table(clustersStrong)
plot(hierClustStrong, cex=0.3);  rect.hclust(hierClustStrong, k=10, border="red")
Plot.Gene.Group((highXhighYStrongGeneNames[clustersStrong == 3])) # this cluster looks good ("Drs", "CG11842", "IM14", "spirit", "CG18563", "SPH93")
Plot.Gene.Group((highXhighYStrongGeneNames[clustersStrong == 8])) # this cluster looks good ("CG30031", "deltaTry", "CG30025")
Plot.Gene.Group((highXhighYStrongGeneNames[clustersStrong == 9])) # this cluster looks ok
Plot.Gene.Group((highXhighYStrongGeneNames[clustersStrong == 10])) # this cluster looks ok
Plot.Gene.Group((highXhighYStrongGeneNames[clustersStrong == 6])) # this cluster looks ok
Plot.Gene.Group((highXhighYStrongGeneNames[clustersStrong == 5])) # this cluster looks ok

# Were these genes assigned to the same cluster? Almost!
Plot.Gene.Group(c("CG43202","IM3","Hayan", "IM14","IM23", "CG33470"))
clustersStrong[match(c("CG43202","IM3","Hayan", "IM14","IM23", "CG33470"), names(clustersStrong))]

# Experimenting with other agglomerations; ward.D or ward.D2 seem to give the most balancec clusters
hierClustStrong <- hclust(as.dist(1-bayesMatrix2SubsetStrong), method="ward.D2")
clustersStrong <- cutree(hierClustStrong, k=10);  table(clustersStrong)
plot(hierClustStrong, cex=0.3);  rect.hclust(hierClustStrong, k=10, border="red")
Plot.Gene.Group((highXhighYStrongGeneNames[clustersStrong == 10])) 

# Interesting inverted patterns: prior = 1 (strong relationship just based on simple linear model), but very low
# R^2 from model 1 and model2-model 3
par(mfrow=c(2,3), mai=c(0.4,0.4,0.4,0.4))
Plot.Gene.Pair("CG9692","CG4594")
Plot.Gene.Pair("CG31189","CG32751")
Plot.Gene.Pair("CG42329","SpdS")
Plot.Gene.Pair("EAChm","CG31102")
Plot.Gene.Pair("Smvt","vri")
par(mfrow=c(1,1))

# --- Test optimization of g ---

gene1 <- "CG33511";  gene2 <- "vri";  Plot.Gene.Pair(gene1, gene2)
X <- Construct.Design.Matrix(gene1, gene2)
Y <- Get.Gene.Data.As.Vector(gene1)
mu0 <- matrix(c(1,1,0,0,0), nrow=5)
gOpt <- Optimize.g(X, Y, mu0);  gOpt
Compute.Gene.Pair.R2.Bayes(gene1, gene2, gOpt)
Compute.Gene.Pair.R2.Bayes(gene1, gene2, 3)
Compute.Gene.Pair.R2(gene1, gene2)
# Comments: optimal g was quite small - not very good when the prior is 0, but helps separate models 2 and 3 much more when
# prior is 1

gene1 <- "Smvt";  gene2 <- "vri";  Plot.Gene.Pair(gene1, gene2)
X <- Construct.Design.Matrix(gene1, gene2)
Y <- Get.Gene.Data.As.Vector(gene1)
mu0 <- matrix(c(1,1,0,0,0), nrow=5)
gOpt <- Optimize.g(X, Y, mu0);  gOpt
Compute.Gene.Pair.R2.Bayes(gene1, gene2, gOpt)
Compute.Gene.Pair.R2.Bayes(gene1, gene2)
Compute.Gene.Pair.R2(gene1, gene2)
# Comments: same as above

gene1 <- "CG18577";  gene2 <- "Mal-A1";  Plot.Gene.Pair(gene1, gene2)
X <- Construct.Design.Matrix(gene1, gene2)
Y <- Get.Gene.Data.As.Vector(gene1)
mu0 <- matrix(c(1,1,0,0,0), nrow=5)
gOpt <- Optimize.g(X, Y, mu0);  gOpt
Compute.Gene.Pair.R2.Bayes(gene1, gene2, gOpt)
Compute.Gene.Pair.R2.Bayes(gene1, gene2, 3)
Compute.Gene.Pair.R2(gene1, gene2)
# Comments: same as above

gene1 <- "CG18577";  gene2 <- "scb";  Plot.Gene.Pair(gene1, gene2)
X <- Construct.Design.Matrix(gene1, gene2)
Y <- Get.Gene.Data.As.Vector(gene1)
mu0 <- matrix(c(0,0,0,0,0), nrow=5)
gOpt <- Optimize.g(X, Y, mu0);  gOpt
Compute.Gene.Pair.R2.Bayes(gene1, gene2, gOpt)
Compute.Gene.Pair.R2.Bayes(gene1, gene2, 3)
Compute.Gene.Pair.R2(gene1, gene2)

# Plot the objective function
gVec <- c(0.01, 1, 2, 3, 4, 5, 6, 50)
objVec <- rep(0, length(gVec))
H <- X %*% solve(t(X) %*% X) %*% t(X);  Yh <- H %*% Y
for(i in 1:length(gVec)) { objVec[i] <- objective(gVec[i]) }
plot(gVec, objVec, type="o")
objVec

# Low STRING score but similar time profiles
par(mfrow=c(2,3), mai=c(.57,.55,.57,.55))
Plot.Gene.Pair("Pgm","CG34331") # include
stringDBMatrix[Gene.Name.To.Flybase.ID("Pgm"), Gene.Name.To.Flybase.ID("CG34331")] # score: NA

Plot.Gene.Pair("CG6543","Nipsnap") #include
stringDBMatrix[Gene.Name.To.Flybase.ID("CG6543"), Gene.Name.To.Flybase.ID("Nipsnap")] # score: 151

Plot.Gene.Pair("Nipsnap","CG1907") # include
stringMatrixSubset["Nipsnap","CG1907"] # score: 154

Plot.Gene.Pair("CG6543","CG1907") # include
stringDBMatrix[Gene.Name.To.Flybase.ID("CG6543"), Gene.Name.To.Flybase.ID("CG1907")] # score: 354

Plot.Gene.Pair("CG1907", "lic") # include
stringDBMatrix[Gene.Name.To.Flybase.ID("lic"), Gene.Name.To.Flybase.ID("CG1907")] # score: NA

Plot.Gene.Pair("ORMDL", "Mdh1")  # include
stringDBMatrix[Gene.Name.To.Flybase.ID("ORMDL"), Gene.Name.To.Flybase.ID("Mdh1")] # score: NA
par(mfrow=c(1,1), mar=c(1,1,1,1))

Plot.Gene.Pair("CG6543","Ctr1B")
stringDBMatrix[Gene.Name.To.Flybase.ID("CG6543"), Gene.Name.To.Flybase.ID("Ctr1B")] # score: NA

Plot.Gene.Pair("CG6543","CG4594")
stringDBMatrix[Gene.Name.To.Flybase.ID("CG6543"), Gene.Name.To.Flybase.ID("CG4594")] # score: 429

Plot.Gene.Pair("Nipsnap", "lic")
stringDBMatrix[Gene.Name.To.Flybase.ID("Nipsnap"), Gene.Name.To.Flybase.ID("lic")] # score: NA

# Group some of the above genes
Plot.Gene.Group(c("CG6543","CG4594","Ctr1B"))
Plot.Gene.Group(c("Pgm","CG34331","CG6543","Nipsnap","CG1907","ORMDL","Mdh1"))
Plot.Gene.Group(c("Pgm","CG6543","CG4594","Ctr1B", "CG34331", "CG6503", "Nipsnap", "Nplp2"))
stringMatrixSubset[c("Pgm","CG6543","CG4594","Ctr1B", "CG34331", "CG6503", "Nipsnap", "Nplp2"), c("Pgm","CG6543","CG4594","Ctr1B", "CG34331", "CG6503", "Nipsnap", "Nplp2")]
stringDBMatrix[Gene.Name.To.Flybase.ID(c("Pgm","CG34331","CG6543","Nipsnap","CG1907","ORMDL","Mdh1")),  Gene.Name.To.Flybase.ID(c("Pgm","CG34331","CG6543","Nipsnap","CG1907","ORMDL","Mdh1"))]

# High STRING score but different time profiles
par(mfrow=c(3,3), mai=c(.6,.55,.4,.5))
Plot.Gene.Pair("CG11459","Gal")
stringDBMatrix[Gene.Name.To.Flybase.ID("CG11459"), Gene.Name.To.Flybase.ID("Gal")] # score: 901

Plot.Gene.Pair("Tpi","Eip71CD")
stringDBMatrix[Gene.Name.To.Flybase.ID("Tpi"), Gene.Name.To.Flybase.ID("Eip71CD")] # score: 550

Plot.Gene.Pair("Jon25Biii","Jon65Aii")
stringDBMatrix[Gene.Name.To.Flybase.ID("Jon25Biii"), Gene.Name.To.Flybase.ID("Jon65Aii")] # score: 914

Plot.Gene.Pair("CG8837","ninaD")
stringDBMatrix[Gene.Name.To.Flybase.ID("ninaD"), Gene.Name.To.Flybase.ID("CG8837")] # score: 900

Plot.Gene.Pair("Got1","ImpL3") # Fold-change is small
stringDBMatrix[Gene.Name.To.Flybase.ID("Got1"), Gene.Name.To.Flybase.ID("ImpL3")] # score: 975

Plot.Gene.Pair("PyK","ImpL3") # Fold-change is small
stringDBMatrix[Gene.Name.To.Flybase.ID("PyK"), Gene.Name.To.Flybase.ID("ImpL3")] # score: 943

Plot.Gene.Pair("NimB3","NimC2")
stringDBMatrix[Gene.Name.To.Flybase.ID("NimB3"), Gene.Name.To.Flybase.ID("NimC2")] # score: 559

Plot.Gene.Pair("NimB3","NimB4")
stringDBMatrix[Gene.Name.To.Flybase.ID("NimB3"), Gene.Name.To.Flybase.ID("NimB4")] # score: 680

Plot.Gene.Pair("Pgm","PGRP-SD")
stringDBMatrix[Gene.Name.To.Flybase.ID("Pgm"), Gene.Name.To.Flybase.ID("PGRP-SD")] # score: 900
par(mfrow=c(1,1))


Plot.Gene.Pair("Tpi","Hsp70Ba") # Fold-change is small for Tpi
stringDBMatrix[Gene.Name.To.Flybase.ID("Tpi"), Gene.Name.To.Flybase.ID("Hsp70Ba")] # score: 513

Plot.Gene.Pair("CG11459", "PGRP-SD")
stringDBMatrix[Gene.Name.To.Flybase.ID("CG11459"), Gene.Name.To.Flybase.ID("PGRP-SD")] # score: 513

Plot.Gene.Pair("CG11459", "CtsB1")
stringDBMatrix[Gene.Name.To.Flybase.ID("CG11459"), Gene.Name.To.Flybase.ID("CtsB1")] # score: 513

Plot.Gene.Group(c("CG11459","PGRP-SD","CtsB1","Rel", "p38c")) # Aggregate some previous examples
Plot.Gene.Group(c("NimB3","NimB4","NimC2")) # Aggregate some previous examples
stringMatrixSubset[c("NimB3","NimB4","NimC2"), c("NimB3","NimB4","NimC2")]
Plot.Gene.Group(c("Got1","Mdh1","r","CG1640","ImpL3")) # Aggregate some previous examples
Plot.Gene.Group(c("CG11459","p38c","PGRP-SD","CtsB1","Rel")) # Aggregate some previous examples
Plot.Gene.Group(c("NimB4","CG31205","NimB3","CG1092","CG34331"))
Plot.Gene.Group(c("NimC2","NimB3","CG10659","CG12057","CG11893"))

# Cyclical patterns
Plot.Gene.Pair("CG5999","CG6910")
stringDBMatrix[Gene.Name.To.Flybase.ID("CG5999"), Gene.Name.To.Flybase.ID("CG6910")] # score: 653

Plot.Gene.Pair("CG6296","CG7025")
stringMatrixSubset["CG6296","CG7025"] # score: 583

Plot.Gene.Pair("CG6543","CG10131")
stringMatrixSubset["CG6543","CG10131"] # score: 965

Plot.Gene.Pair("mino","CG1946")
stringMatrixSubset["mino","CG1946"] # score: 714

