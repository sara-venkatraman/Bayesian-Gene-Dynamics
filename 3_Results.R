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


