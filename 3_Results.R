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
adjBayes <- (bayesLLR2Mat > 0.95) + 0
networkBayes <- graph_from_adjacency_matrix(adjBayes , mode='undirected', diag=F)
plot(networkBayes, layout=layout_nicely, vertex.label=NA, vertex.size=3, vertex.color=cluster_label_prop(networkBayes)$membership)

# Plot largest connected component from Bayesian analysis
# Contains 426/1735 genes (from full gene set), using 0.95 R^2 cutoff for edges
# Contains 290/951 genes (from DE gene set), using 0.95 R^2 cutoff for edges
networkBayesComponents <- clusters(networkBayes, mode="strong")
subnetworkBayes <- induced_subgraph(networkBayes, V(networkBayes)[networkBayesComponents$membership == which.max(networkBayesComponents$csize)])
plot(subnetworkBayes, layout=layout_with_kk, vertex.label=NA, vertex.size=4.5, edge.width=0.7, vertex.color=cluster_label_prop(subnetworkBayes)$membership, vertex.frame.color="darkslategray", edge.color="dimgray")

# --- Hierarchical clustering from Bayesian analysis ---

# Compute distance matrix and cluster via Ward's method
distMatrix <- 1 - (bayesLLR2Mat)
hierClust <- hclust(as.dist(distMatrix), method="ward.D")
plot(hierClust)

# Cut the dendrogram at ~ 9 or 10 (yields 15 clusters)
subGroups <- cutree(hierClust, h=9)
table(subGroups)
Plot.Gene.Group(geneNames[subGroups == 14], T, F)

# For each of the 15 clusters, plot the network formed from the cluster
pdf("Output/GeneClusterNetworks.pdf", height=11, width=16)
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
  for(j in 1:nrow(subnetEdges)) {
    subnetEdges$Prior[j] <- priorMatrix[subnetEdges[j,"Gene1"], subnetEdges[j,"Gene2"]]
  }
  
  # Write the dataframe to a CSV for this cluster
  write.csv(subnetEdges, paste("Output/Cluster", i, "Edges.csv", sep=""), row.names=F)
  
  # Edges between genes with unknown associations will be blue, and edges 
  # between genes with known associations will be red.
  E(subnetHierClust)$color[is.na(subnetEdges$Prior)] <- 'blue'
  E(subnetHierClust)$color[!is.na(subnetEdges$Prior)] <- alpha('red', 0.7)
  
  # Plot the network for this cluster on the PDF
  numUnknownEdges <- sum(is.na(subnetEdges$Prior))
  numKnownEdges <- sum(!is.na(subnetEdges$Prior))
  plotTitle <- paste("Network from Cluster ", i, "\n(", length(subGroupNames), " genes; ", numKnownEdges, " known edges, ", numUnknownEdges, " unknown)", sep="")
  plot(subnetHierClust, layout=layout_nicely, vertex.label=NA, vertex.size=5, 
       edge.width=1, vertex.color="gray", vertex.frame.color="darkslategray", 
       main=plotTitle, cex.main=1.5)
}
dev.off(); par(mfrow=c(1,1))

# For each of the 15 clusters, plot the time profiles of all genes in the cluster
pdf("Output/GeneClusters.pdf", height=11, width=23)
par(mfrow=c(3,5))
plotColors <- c(brewer.pal(n=8, name="Dark2"), brewer.pal(n=8, name="Dark2"))
for(i in 1:length(table(subGroups))) {
  plotTitle=paste("Cluster ", i, " (", length(geneNames[subGroups == i]), " genes)", sep="")
  Plot.Gene.Group(geneNames[subGroups == i], monochrome=plotColors[i], points=F, 
                  plotTitle=plotTitle, titleSize=1.5)
}
dev.off(); par(mfrow=c(1,1))

# Create a list in which each entry is a vector of gene names in one cluster
clusterList <- list()
for(i in 1:length(table(subGroups))) { 
  clusterList[[i]] <- names(subGroups)[subGroups == i]
}
save(clusterList, file="Output/GeneClusters.RData")

# Alternatively, extract k=40 clusters and plot clusters with <= 37 genes
k <- 40
clusters <- cutree(hierClust, k=k);  table(clusters)
pdf("Output/GeneClustersWithColors.pdf", height=9, width=11)
par(mfrow=c(3,3))
for(i in 1:k) {
  if(i <= 37) { Plot.Gene.Group(geneNames[clusters == i]) }
}
dev.off(); par(mfrow=c(1,1))

# --- R^2 scatterplots ---

Draw.R2.Scatterplot(nonBayesLLR2Mat.other, nonBayesLLR2Mat - nonBayesLLR2Mat.own, 
                    priorMatrix, geneSubset, interactive=F)

Draw.R2.Scatterplot(bayesLLR2Mat.other, bayesLLR2Mat - bayesLLR2Mat.own, 
                    priorMatrix, geneSubset, interactive=F)

# --- Large heatmap of the Bayesian R^2 similarity matrix ---

Draw.Heatmap.Diagonal.Block <- function(k) {
  simMatrixSection <- melt(as.matrix(bayesLLR2Mat[k:(k+346), k:(k+346)]))
  plotTitle <- paste("Diagonal block of LLR2 similarity matrix: Genes", k, "to", k+346)
  ggplot(data=simMatrixSection, aes(x=Var2, y=Var1, fill=value)) + geom_tile() +
    theme(axis.title.x = element_blank(), axis.title.y = element_blank()) +
    theme(axis.text.x=element_text(angle=90, hjust=0.95,vjust=0.2)) +
    scale_fill_distiller(palette="Blues", direction=1) + 
    ggtitle(label=plotTitle) + theme(plot.title=element_text(size=40, face="bold")) +
    theme(legend.position="none")
}

pdf("Output/HeatmapFullGeneSet.pdf", height=38, width=34)
Draw.Heatmap.Diagonal.Block(1)
Draw.Heatmap.Diagonal.Block(348)
Draw.Heatmap.Diagonal.Block(695)
Draw.Heatmap.Diagonal.Block(1042)
Draw.Heatmap.Diagonal.Block(1389)
dev.off()

