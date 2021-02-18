# --- Read in R^2 matrices ---

# If reading these files, run DatasetLoader.R and PlottingFunctions.R first.

nonBayesLLR2Mat <- read.csv(paste("../Processed Data/R-Squared (", subdirectory, ")/NonBayesLLR2.csv", sep=""), row.names=1)
nonBayesLLR2Mat.other <- read.csv(paste("../Processed Data/R-Squared (", subdirectory, ")/NonBayesLLR2_Other.csv", sep=""), row.names=1)
nonBayesLLR2Mat.own <- read.csv(paste("../Processed Data/R-Squared (", subdirectory, ")/NonBayesLLR2_Own.csv", sep=""), row.names=1)
rownames(nonBayesLLR2Mat) <- geneNames;  colnames(nonBayesLLR2Mat) <- geneNames
rownames(nonBayesLLR2Mat.other) <- geneNames;  colnames(nonBayesLLR2Mat.other) <- geneNames
rownames(nonBayesLLR2Mat.own) <- geneNames;  colnames(nonBayesLLR2Mat.own) <- geneNames

bayesLLR2Mat <- read.csv(paste("../Processed Data/R-Squared (", subdirectory, ")/BayesLLR2.csv", sep=""), row.names=1)
bayesLLR2Mat.other <- read.csv(paste("../Processed Data/R-Squared (", subdirectory, ")/BayesLLR2_Other.csv", sep=""), row.names=1)
bayesLLR2Mat.own <- read.csv(paste("../Processed Data/R-Squared (", subdirectory, ")/BayesLLR2_Own.csv", sep=""), row.names=1)
rownames(bayesLLR2Mat) <- geneNames;  colnames(bayesLLR2Mat) <- geneNames
rownames(bayesLLR2Mat.other) <- geneNames;  colnames(bayesLLR2Mat.other) <- geneNames
rownames(bayesLLR2Mat.own) <- geneNames;  colnames(bayesLLR2Mat.own) <- geneNames

# --- Network diagrams ---

# Non-Bayesian analysis
adjNonBayes <- (nonBayesLLR2Mat > 0.95) + 0
networkNonBayes <- graph_from_adjacency_matrix(adjNonBayes , mode='undirected', diag=F)
plot(networkNonBayes, layout=layout.fruchterman.reingold, vertex.label=NA, vertex.size=2, vertex.color=cluster_label_prop(networkNonBayes)$membership)

# Plot largest connected component from non-Bayesian analysis
# Contains 1731/1735 genes (from full gene set)
# Contains 871/951 genes (from DE gene set)
networkNonBayesComponents <- clusters(networkNonBayes, mode="strong")
subnetworkNonBayes <- induced_subgraph(networkNonBayes, V(networkNonBayes)[networkNonBayesComponents$membership == which.max(networkNonBayesComponents$csize)])
plot(subnetworkNonBayes, layout=layout.fruchterman.reingold, vertex.label=NA, vertex.size=4, vertex.frame.color="darkslategray",edge.width=0.5)#vertex.color=cluster_label_prop(subnetworkNonBayes)$membership, vertex.frame.color="black")

# Bayesian analysis
adjBayes <- (bayesLLR2Mat > 0.85) + 0
networkBayes <- graph_from_adjacency_matrix(adjBayes , mode='undirected', diag=F)
plot(networkBayes, layout=layout_nicely, vertex.label=NA, vertex.size=3, vertex.color=cluster_label_prop(networkBayes)$membership)

# Plot largest connected component from Bayesian analysis
# Contains 426/1735 genes (from full gene set), using 0.95 R^2 cutoff for edges
# Contains 290/951 genes (from DE gene set), using 0.95 R^2 cutoff for edges
networkBayesComponents <- clusters(networkBayes, mode="strong")
subnetworkBayes <- induced_subgraph(networkBayes, V(networkBayes)[networkBayesComponents$membership == which.max(networkBayesComponents$csize)])
plot(subnetworkBayes, layout=layout_with_kk, vertex.label=NA, vertex.size=4.5, edge.width=0.7, vertex.color=cluster_label_prop(subnetworkBayes)$membership, vertex.frame.color="darkslategray", edge.color="dimgray")

# --- Plotting clusters from Bayesian analysis ---

# Hierarchical clustering: compute distance matrix and cluster via Ward's method
distMatrix <- 1 - (bayesLLR2Mat)
hierClust <- hclust(as.dist(distMatrix), method="ward.D")

# Extract k=40 clusters and plot clusters with <= 37 genes in a PDF
k <- 40
clusters <- cutree(hierClust, k=k);  table(clusters)
pdf("Output/GeneClusters.pdf", height=9, width=11)
par(mfrow=c(3,3))
for(i in 1:k) {
  if(i <= 37) { Plot.Gene.Group(geneNames[clusters == i]) }
}
dev.off(); par(mfrow=c(1,1))

# Cut the dendrogram at ~ 9 or 10 (yields 15 clusters)
plot(hierClust)
subGroups <- cutree(hierClust, h=9)
table(subGroups)
Plot.Gene.Group(geneNames[subGroups == 16], T, F)

# For each of the 15 clusters, plot the network formed from the cluster
pdf("Output/15GeneClusterNetworks.pdf", height=11, width=16)
par(mfrow=c(3,5))
for(i in 1:length(table(subGroups))) {
  # Get the gene names in the i-th cluster
  subGroupNames <- names(subGroups)[subGroups == i]
  
  # Form an adjacency matrix and the corresponding graph from those genes.
  # Define an edge between genes if the Bayesian LLR2 > 0.9
  adjBayesSubGroup <- (bayesLLR2Mat[subGroupNames, subGroupNames] > 0.9) + 0
  subnetHierClust <- graph_from_adjacency_matrix(adjBayesSubGroup , mode='undirected', diag=F)
  
  # Form a dataframe out of the list of edges in the network, where each row
  # defines an edge between the genes in columns 1 and 2. Add a third column
  # containing the prior adjacency matrix value for each of those edges.
  subnetEdges <- data.frame(as_edgelist(subnetHierClust))
  colnames(subnetEdges) <- c("Gene1", "Gene2")
  subnetEdges$Prior <- 0
  for(j in 1:nrow(subnetEdges))
    subnetEdges$Prior[j] <- priorMatrix[subnetEdges[j,"Gene1"], subnetEdges[j,"Gene2"]]
  
  # Write the dataframe to a CSV for this cluster
  write.csv(subnetEdges, paste("Output/Cluster", i, "Edges.csv", sep=""), row.names=F)
  
  # Edges between genes with unknown associations will be blue, and edges 
  # between genes with known associations will be red.
  E(subnetHierClust)$color[is.na(subnetEdges$Prior)] <- 'blue'
  E(subnetHierClust)$color[!is.na(subnetEdges$Prior)] <- alpha('red', 0.7)
  
  # Plot the network for this cluster on the PDF
  plot(subnetHierClust, layout=layout_nicely, vertex.label=NA, vertex.size=5, 
       edge.width=1, vertex.color="gray", vertex.frame.color="darkslategray", 
       main=paste("Network from Cluster", i), cex.main=1.5)
}
dev.off(); par(mfrow=c(1,1))

# For each of the 15 clusters, plot the time profiles of all genes in the cluster
pdf("Output/15GeneClusters.pdf", height=11, width=23)
par(mfrow=c(3,5))
plotColors <- c(brewer.pal(n=8, name="Dark2"), brewer.pal(n=8, name="Dark2"))
for(i in 1:length(table(subGroups))) {
  Plot.Gene.Group(geneNames[subGroups == i], monochrome=plotColors[i], points=F, 
                  plotTitle=paste("Cluster ", i, " (", length(geneNames[subGroups == i]), " genes)", sep=""),
                  titleSize=1.5)
}
dev.off(); par(mfrow=c(1,1))

# --- R^2 scatterplots ---

Draw.R2.Scatterplot(nonBayesLLR2Mat.other, nonBayesLLR2Mat - nonBayesLLR2Mat.own, 
                    priorMatrix, geneSubset, interactive=F)

Draw.R2.Scatterplot(bayesLLR2Mat.other, bayesLLR2Mat - bayesLLR2Mat.own, 
                    priorMatrix, geneSubset, interactive=F)

# --- Inferring R^2 cutoffs from data ---

# Vectorize the LLR2 matrix and the prior matrix (saved into CSV)
bayesLLR2Vec <- data.frame(Vectorize.Labeled.Square.Matrix(bayesLLR2Mat),
                           Vectorize.Labeled.Square.Matrix(priorMatrix))
bayesLLR2Vec <- bayesLLR2Vec[order(bayesLLR2Vec[,1], decreasing=T),]
colnames(bayesLLR2Vec) <- c("LLR2", "Prior")

# How many edges exist are there in the prior matrix? Compute the proportion
# of "1"s in one of the triangular halves of priorMatrix
numPositiveEdges <- sum(priorMatrix[upper.tri(priorMatrix)] == 1, na.rm=T)
numPossibleEdges <- sum(upper.tri(priorMatrix))
priorEdge1Prob <- numPositiveEdges / numPossibleEdges 

# It's about 3.3%. Let prior parameters for Beta distribution on p be
# alpha = 1 and beta = 29.7. Then the posterior mean is:
alpha <- 1;  beta <- 29.7
postEdge1Prob <- (numPositiveEdges + alpha) / (numPossibleEdges + alpha + beta)

# Take the top postEdgeProb proportion of R^2 values
topConnections <- bayesLLR2Vec[1:round(postEdge1Prob*nrow(bayesLLR2Vec)),]
topGenePairs <- strsplit(rownames(topConnections), ", ")
Plot.Gene.Group(topGenePairs[[50]])

# Compute P(W=1 | R^2 = r). Rather than using the built-in 'density' function,
# use the kdensity package, which returns a density that can be evaluated at a point
library(kdensity)
LLR2WithPrior0 <- bayesLLR2Vec[is.na(bayesLLR2Vec$Prior) | bayesLLR2Vec$Prior == 0,]
LLR2WithPrior1 <- bayesLLR2Vec[!is.na(bayesLLR2Vec$Prior) & bayesLLR2Vec$Prior == 1,]
LLR2Density0 <- kdensity(LLR2WithPrior0$LLR2)
LLR2Density1 <- kdensity(LLR2WithPrior1$LLR2)
plot(LLR2Density0, main="Density of LLR2 (W=0)")
plot(LLR2Density1, main="Density of LLR2 (W=1)")
Compute.Posterior.Edge.Prob <- function(R2Density0, R2Density1, r, w, priorEdge1Prob) {
  priorEdge0Prob <- 1 - priorEdge1Prob
  denominator <- R2Density0(r)*priorEdge0Prob + R2Density1(r)*priorEdge1Prob
  if(w == 0)
    numerator <- R2Density0(r)*priorEdge0Prob
  else
    numerator <- R2Density1(r)*priorEdge1Prob
  numerator/denominator
}
Compute.Posterior.Edge.Prob(LLR2Density0, LLR2Density1, 0.5, 1, priorEdge1Prob)

postEdgeProbUncond <- 0
for(i in 1:nrow(LLR2WithPrior1)) {
  postEdgeProbUncond <- postEdgeProbUncond + LLR2Density1(LLR2WithPrior1$LLR2[i])*priorEdge1Prob
}
postEdgeProbUncond


