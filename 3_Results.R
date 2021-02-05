# --- Network diagrams ---

# Non-Bayesian analysis
adjNonBayes <- (nonBayesLLR2Mat > 0.95) + 0
networkNonBayes <- graph_from_adjacency_matrix(adjNonBayes , mode='undirected', diag=F)
plot(networkNonBayes, layout=layout.fruchterman.reingold, vertex.label=NA, vertex.size=2, vertex.color=cluster_label_prop(networkNonBayes)$membership)

# Plot largest connected component from non-Bayesian analysis
# Contains 1731/1735 genes (from full gene set)
networkNonBayesComponents <- clusters(networkNonBayes, mode="strong")
subnetworkNonBayes <- induced_subgraph(networkNonBayes, V(networkNonBayes)[networkNonBayesComponents$membership == which.max(networkNonBayesComponents$csize)])
plot(subnetworkNonBayes, layout=layout.fruchterman.reingold, vertex.label=NA, vertex.size=3, vertex.color=cluster_label_prop(subnetworkNonBayes)$membership)

# Bayesian analysis
adjBayes <- (bayesLLR2Mat > 0.95) + 0
networkBayes <- graph_from_adjacency_matrix(adjBayes , mode='undirected', diag=F)
plot(networkBayes, layout=layout_nicely, vertex.label=NA, vertex.size=3, vertex.color=cluster_label_prop(networkBayes)$membership)

# Plot largest connected component from Bayesian analysis
# Contains 426/1735 genes (from full gene set)
networkBayesComponents <- clusters(networkBayes, mode="strong")
subnetworkBayes <- induced_subgraph(networkBayes, V(networkBayes)[networkBayesComponents$membership == which.max(networkBayesComponents$csize)])
plot(subnetworkBayes, layout=layout_with_kk, vertex.label=NA, vertex.size=4, vertex.color=cluster_label_prop(subnetworkBayes)$membership)

# Plot gene communities
pdf("GeneCommunities.pdf", height=9, width=11)
par(mfrow=c(3,3))
for(i in 1:max(subnetworkBayesCommunities$membership)) {
  genesInCommunity <- V(subnetworkBayes)[subnetworkBayesCommunities$membership == i]
  if(i <= 37) {
    Plot.Gene.Group(genesInCommunity)
  }
}
dev.off(); par(mfrow=c(1,1))

# --- Plotting clusters from Bayesian analysis ---

# Hierarchical clustering

distMatrix <- 1 - (bayesLLR2Mat)
k <- 40
hierClust <- hclust(as.dist(distMatrix), method="ward.D")
clusters <- cutree(hierClust, k=k)
table(clusters)

pdf("GeneClusters.pdf", height=9, width=11)
par(mfrow=c(3,3))
for(i in 1:k) {
  genesInCluster <- geneNames[clusters == i]
  if(i <= 37) {
    Plot.Gene.Group(genesInCluster)
  }
}
dev.off(); par(mfrow=c(1,1))

# --- R^2 scatterplots ---

Draw.R2.Scatterplot(nonBayesLLR2Mat.other, nonBayesLLR2Mat - nonBayesLLR2Mat.own, 
                    priorMatrix, geneSubset)

Draw.R2.Scatterplot(bayesLLR2Mat.other, bayesLLR2Mat - bayesLLR2Mat.own, 
                    priorMatrix, geneSubset)

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
priorEdgeProb <- numPositiveEdges / numPossibleEdges 

# It's about 3.3%. Let prior parameters for Beta distribution on p be
# alpha = 1 and beta = 29.7. Then the posterior mean is:
alpha <- 1;  beta <- 29.7
postEdgeProb <- (numPositiveEdges + alpha) / (numPossibleEdges + alpha + beta)

# Take the top postEdgeProb proportion of R^2 values
topConnections <- bayesLLR2Vec[1:round(postEdgeProb*length(bayesLLR2Vec))]
topGenePairs <- strsplit(names(topConnections), ", ")
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
Compute.Posterior.Edge.Prob(LLR2Density0, LLR2Density1, 0.5, 1, priorEdgeProb)
