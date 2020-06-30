# --- Script 5: Exploratory data analysis ---

# Load datasets needed for constructing priors
source("DatasetLoader.R")

# Load functions in EmpiricalBayesFunctions.R
source("EmpiricalBayesFunctions.R")

# Before running code below, also need to run the relevant section of PriorConstruction.R.

# --- Construct scatterplots (with and without Bayesian regression) ---

nonBayesR2Matrices <- Compute.R2.Matrices(genesSubset, bayes=FALSE)
bayesR2Matrices <- Compute.R2.Matrices(genesSubset, bayes=TRUE)

Draw.Metric.Scatterplot.For.Binary.Prior(nonBayesR2Matrices, bayes=FALSE)
Draw.Metric.Scatterplot.For.Binary.Prior(bayesR2Matrices, bayes=TRUE)

# --- Hierarchical clustering using empirical Bayes R^2 ---

avgSimMatrix <- (bayesR2Matrices$matrix1 + (bayesR2Matrices$matrix2-bayesR2Matrices$matrix3))/2
hierClust <- hclust(as.dist(1-avgSimMatrix), method="ward.D")
clusters <- cutree(hierClust, k=20)
table(clusters)
plot(hierClust, cex=0.3);  rect.hclust(hierClust, k=20, border="red")
Plot.Gene.Group((genesSubset[clusters == 20]))

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
upperRight <- rownames(bayesVec1)[bayesVec1 > 0.5 & bayesVec2-bayesVec3 > 0.2]
lowerRight <- rownames(bayesVec1)[bayesVec1 > 0.5 & bayesVec2-bayesVec3 <= 0.2] 
lowerLeft <- rownames(bayesVec1)[bayesVec1 <= 0.5 & bayesVec2-bayesVec3 <= 0.2] 
upperLeft <- rownames(bayesVec1)[bayesVec1 <= 0.5 & bayesVec2-bayesVec3 > 0.2]  # upper-left quadrant - interesting patterns here, but not always close
lapply(list(upperRight, lowerRight, lowerLeft, upperLeft), function(x) { return(length(x)/length(bayesVec1)) }) # Print the above in proportions

# Plot all the gene pairs that are in the upper-right quadrant
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

# Some interesting plots that come out of the above PDF ("connected components")
Plot.Gene.Group(c("Sodh-1","UGP","glob1","Gip","Pgi", "GstE9", "CG6084"))
Plot.Gene.Group(c("vri","Slob","per","CG33511","Pdp1"))
Plot.Gene.Pair("CG33511","vri") # In STRING-DB, no annotation available for CG33511. vri (clock-controlled) affects hair and cell growth

# Try hierarchical clustering on the genes that appeared in the upper right quadrant
upperRightGeneNames <- unique(strsplit(paste(upperRight, collapse=", "), ", ")[[1]])
bayesMatrix2Subset <- bayesR2Matrices$matrix2[upperRightGeneNames, upperRightGeneNames]

hierClustUpperRight <- hclust(as.dist(1-bayesMatrix2Subset), method="ward.D")
clustersUpperRight <- cutree(hierClustUpperRight, k=10);  table(clustersUpperRight)
plot(hierClustUpperRight, cex=0.3);  rect.hclust(hierClustUpperRight, k=10, border="red")

# Plot genes in one cluster
Plot.Gene.Group((upperRightGeneNames[clustersUpperRight == 3]))

# Two clusters that looked good were: 
# ("Drs", "CG11842", "IM14", "spirit", "CG18563", "SPH93")
# ("CG30031", "deltaTry", "CG30025")

# Interesting inverted patterns: prior = 1 (strong relationship just based on simple linear model), 
# but very low R^2 from model 1 and model2-model 3
par(mfrow=c(2,3), mai=c(0.4,0.4,0.4,0.4))
Plot.Gene.Pair("CG9692","CG4594")
Plot.Gene.Pair("CG31189","CG32751")
Plot.Gene.Pair("CG42329","SpdS")
Plot.Gene.Pair("EAChm","CG31102")
Plot.Gene.Pair("Smvt","vri")
par(mfrow=c(1,1))

# --- Test optimization of g for the normal-inverse gamma model ---

gene1 <- "CG33511";  gene2 <- "vri";  Plot.Gene.Pair(gene1, gene2)
X <- Construct.Design.Matrix(gene1, gene2)
Y <- Get.Gene.Data.As.Vector(gene1)
mu0 <- matrix(c(1,1,0,0,0), nrow=5)
gOpt <- Optimize.g(X, Y, mu0);  gOpt
Compute.Gene.Pair.R2.Bayes(gene1, gene2, gOpt)
Compute.Gene.Pair.R2.Bayes(gene1, gene2, 3)
Compute.Gene.Pair.R2(gene1, gene2)
# Comments: optimal g was quite small - not very good when the prior is 0, but
# helps separate models 2 and 3 much more when prior is 1

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

# Plot the objective function
gVec <- c(0.01, 1, 2, 3, 4, 5, 6, 50)
objVec <- rep(0, length(gVec))
H <- X %*% solve(t(X) %*% X) %*% t(X);  Yh <- H %*% Y
for(i in 1:length(gVec)) { objVec[i] <- objective(gVec[i]) }
plot(gVec, objVec, type="o")
objVec

# --- Exploring inconsistency between associations suggested by priors and
#     associations suggested by data ---

# Low STRING score but similar time profiles
par(mfrow=c(2,3), mai=c(.57,.55,.57,.55))
Plot.Gene.Pair("Pgm","CG34331")
stringDBMatrix[Gene.Name.To.Flybase.ID("Pgm"), Gene.Name.To.Flybase.ID("CG34331")] # score: NA

Plot.Gene.Pair("CG6543","Nipsnap")
stringDBMatrix[Gene.Name.To.Flybase.ID("CG6543"), Gene.Name.To.Flybase.ID("Nipsnap")] # score: 151

Plot.Gene.Pair("Nipsnap","CG1907")
stringMatrixSubset["Nipsnap","CG1907"] # score: 154

Plot.Gene.Pair("CG6543","CG1907")
stringDBMatrix[Gene.Name.To.Flybase.ID("CG6543"), Gene.Name.To.Flybase.ID("CG1907")] # score: 354

Plot.Gene.Pair("CG1907", "lic")
stringDBMatrix[Gene.Name.To.Flybase.ID("lic"), Gene.Name.To.Flybase.ID("CG1907")] # score: NA

Plot.Gene.Pair("ORMDL", "Mdh1")
stringDBMatrix[Gene.Name.To.Flybase.ID("ORMDL"), Gene.Name.To.Flybase.ID("Mdh1")] # score: NA
par(mfrow=c(1,1), mar=c(1,1,1,1))

# Group some of the above genes
Plot.Gene.Group(c("CG6543","CG4594","Ctr1B"))
Plot.Gene.Group(c("Pgm","CG6543","CG4594","Ctr1B", "CG34331", "CG6503", "Nipsnap", "Nplp2"))
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

# Group some of the above genes
Plot.Gene.Group(c("CG11459","PGRP-SD","CtsB1","Rel", "p38c"))
Plot.Gene.Group(c("NimB3","NimB4","NimC2"))
Plot.Gene.Group(c("Got1","Mdh1","r","CG1640","ImpL3"))
Plot.Gene.Group(c("NimB4","CG31205","NimB3","CG1092","CG34331"))
Plot.Gene.Group(c("NimC2","NimB3","CG10659","CG12057","CG11893"))
stringMatrixSubset[c("NimB3","NimB4","NimC2"), c("NimB3","NimB4","NimC2")]

# Cyclical patterns
Plot.Gene.Pair("CG5999","CG6910")
Plot.Gene.Group(c("IM14","IM1","IM2","CG5791", "Drs"))
stringDBMatrix[Gene.Name.To.Flybase.ID("CG5999"), Gene.Name.To.Flybase.ID("CG6910")] # score: 653

# --- Soft clustering with MFuzz package ---

library(Mfuzz)
geneDataMFuzz <- table2eset("../Processed Data for Clustering/DElogChange_MFuzz.txt")
mClustering <- mestimate(geneDataMFuzz)
mFuzzClust <- mfuzz(geneDataMFuzz, c=16, m=mClustering)
mfuzz.plot(geneDataMFuzz, cl=mFuzzClust, mfrow=c(4,4), time.labels=hours, new.window=F)

# --- Soft clustering with lead-lag R^2 ---

library(fclust)

# Save components of output from R^2 computations
bayesR2MatricesVec <- Vectorize.R2.Matrices(genesSubset, bayesR2Matrices)
bayesMat1 <- bayesR2Matrices$matrix1; bayesMat2 <- bayesR2Matrices$matrix2; bayesMat3 <- bayesR2Matrices$matrix3
bayesVec1 <- bayesR2MatricesVec$vector1;  bayesVec2 <- bayesR2MatricesVec$vector2;  bayesVec3 <- bayesR2MatricesVec$vector3

# Extract names of gene pairs falling in upper-right area of scatterplot
upperRight <- rownames(bayesVec1)[bayesVec1 > 0.4 & bayesVec2-bayesVec3 > 0.2]
upperRightGeneNames <- unique(strsplit(paste(upperRight, collapse=", "), ", ")[[1]])

# Extract subsets of similarity matrices corresponding to genes in 'upperRight'
bayesMatrix1Subset <- bayesMat1[upperRightGeneNames, upperRightGeneNames]
bayesMatrix2Subset <- bayesMat2[upperRightGeneNames, upperRightGeneNames]
bayesMatrix3Subset <- bayesMat3[upperRightGeneNames, upperRightGeneNames]

# Run soft clustering method and extract/view the weight matrix of cluster membership scores
dissimilarityMatrix <- 1 - bayesMatrix2Subset
numClusters <- 15
softClustering <- Fclust(X=as.dist(dissimilarityMatrix), k=numClusters, type="standard", ent=FALSE, noise=FALSE, stand=0, distance=TRUE)
weightMatrix <- softClustering$U
# View(weightMatrix)

# Plot individual clusters
clusterIndex <- 1
table(softClustering$U[,clusterIndex] > (1/numClusters + 1e-7) + 0)
genesInCluster <- rownames(weightMatrix)[weightMatrix[,clusterIndex] > (1/numClusters + 1e-7)]
Plot.Gene.Group(genesInCluster)

# Plot clusters in a PDF
pdf("GeneClusters.pdf", height=8, width=12)
par(mfrow=c(2,2), mai=c(0.8,1.1,0.8,1.1))
for(i in seq(1, numClusters)) {
  tol <- 0.00000018
  if(sum(softClustering$U[,i] > (1/numClusters + tol) + 0) > 2) {
    genesInCluster <- rownames(weightMatrix)[weightMatrix[,i] > (1/numClusters + tol)]
    Plot.Gene.Group(genesInCluster)
  }
}
dev.off(); par(mfrow=c(1,1))

# --- Interesting plots that we obtain after adding neighbors of DE genes to dataset ---

# Prior matrix = 2 for all these genes
Plot.Gene.Group(c("RpL39","RpS6","RpS21","RpL31","RpL10Ab","RpL26"))
c("RpL39","RpS6","RpS21","RpL31","RpL10Ab","RpL26") %in% DEgeneNames # None were DE!

# Prior matrix = 2 for all these genes
Plot.Gene.Group(c("Vha36-1","Vha14-1","VhaM8.9"))

# Add some more genes which did not have entry 2 in prior matrix
Plot.Gene.Group(c("Vha36-1","Vha14-1","VhaM8.9","Vha100-3","RpL12","VhaM9.7-b","Vha16-5"))

# Prior matrix = 2 for all these genes
Plot.Gene.Group(c("Cul3","Cand1","mib2"))

# --- Soft clustering implementation (to be continued) ---

# # Initialize weight matrix: a gene falls uniformly in any one of the clusters
# numClusters <- 16
# weightMatrixOld <- matrix(rep(1/numClusters, subsetSize*numClusters), nrow=subsetSize, ncol=numClusters)
# rownames(weightMatrixOld) <- genesSubset
# 
# # Initialize matrix of centroids
# centroids <- matrix(0, nrow=numClusters, ncol=length(hours))
# 
# # Optimization
# iterations <- 0;  continue <- TRUE
# while(continue) {
#   # Update centroids
#   for(j in 1:numClusters) {
#     
#   }
#   # Update weights
#   
#   # Check convergence
#   iterations <- iterations + 1
#   if(iterations > 100 && norm(weightMatrixNew - weightMatrixOld, "F") > 0.3) {
#     continue <- FALSE
#   }
# }

