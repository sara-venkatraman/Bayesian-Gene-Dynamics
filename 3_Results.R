# --- Network diagrams ---

# Non-Bayesian analysis
r2MatAdjNonBayes <- (r2MatNonBayes > 0.95) + 0
diag(r2MatAdjNonBayes) <- 0
networkNonBayes <- graph_from_adjacency_matrix(r2MatAdjNonBayes , mode='undirected', diag=F)
plot(networkNonBayes, layout=layout_nicely, vertex.label=NA, vertex.size=2)

# Bayesian analysis
r2MatAdjBayes <- (r2MatBayes > 0.95) + 0
diag(r2MatAdjBayes) <- 0
network <- graph_from_adjacency_matrix(r2MatAdjBayes , mode='undirected', diag=F)
plot(network, layout=layout_nicely, vertex.label=NA, vertex.size=3)

# --- Plotting clusters from Bayesian analysis ---

# Hierarchical clustering

distMatrix <- 1-r2MatBayes
k <- 30
hierClust <- hclust(as.dist(distMatrix), method="ward.D")
clusters <- cutree(hierClust, k=k)
table(clusters)

pdf("GeneClusters.pdf", height=9, width=11)
par(mfrow=c(3,3))
for(i in 1:k) {
  genesInCluster <- geneNames[clusters == i]
  if(i <= 28) {
    Plot.Gene.Group(genesInCluster)
  }
}
dev.off(); par(mfrow=c(1,1))

# --- R^2 scatterplots ---

Draw.R2.Scatterplot(nonBayesLLR2Mat.other, nonBayesLLR2Mat - nonBayesLLR2Mat.own, 
                    priorMatrix, geneSubset)

Draw.R2.Scatterplot(bayesLLR2Mat.other, bayesLLR2Mat - bayesLLR2Mat.own, 
                    priorMatrix, geneSubset)


