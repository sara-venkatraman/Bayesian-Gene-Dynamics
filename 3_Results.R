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
distMatrix <- 1 - bayesLLR2Mat
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

# For each of the 15 clusters: take the set of genes with at least one blue edge,
# Add all their neighbors (genes with at least one edge), plot the network
pdf("Output/GeneClusterSubnetworks.pdf", height=12, width=12)
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
  colnames(subnetEdges) <- c("Gene1", "Gene2"); subnetEdges$Prior <- 0
  for(j in 1:nrow(subnetEdges)) {
    prior <- priorMatrix[subnetEdges[j,"Gene1"], subnetEdges[j,"Gene2"]]
    subnetEdges$Prior[j] <- prior
  }
  
  # Get the nodes (genes) with at least one unknown edge
  unknownEdges <- subnetEdges[is.na(subnetEdges$Prior),]
  nodesWithUnknownEdges <- unique(c(unknownEdges$Gene1, unknownEdges$Gene2))
  
  # Get the neighbors of the nodes with at least one unknown edge, store in allNodes
  neighborNodes <- adjacent_vertices(subnetHierClust, nodesWithUnknownEdges)
  allNodes <- c()
  for(j in 1:length(neighborNodes)) {
    allNodes <- c(allNodes, names(neighborNodes[j]), names(neighborNodes[[j]]))
  }
  allNodes <- unique(allNodes)
  
  # Form a new subnetwork out of allNodes
  newAdjBayesSubGroup <- (bayesLLR2Mat[allNodes, allNodes] > 0.9) + 0
  newSubnet <- graph_from_adjacency_matrix(newAdjBayesSubGroup , mode='undirected', diag=F)
  
  # Get a new subnetEdges dataframe corresponding to allNodes
  newSubnetEdges <- data.frame(as_edgelist(newSubnet))
  colnames(newSubnetEdges) <- c("Gene1", "Gene2"); newSubnetEdges$Prior <- 0
  for(j in 1:nrow(newSubnetEdges)) {
    prior <- priorMatrix[newSubnetEdges[j,"Gene1"], newSubnetEdges[j,"Gene2"]]
    newSubnetEdges$Prior[j] <- prior
  }
  
  # Edges between genes with unknown associations will be blue, and edges 
  # between genes with known associations will be red.
  E(newSubnet)$color[is.na(newSubnetEdges$Prior)] <- 'blue'
  E(newSubnet)$color[!is.na(newSubnetEdges$Prior)] <- alpha('red', 0.7)
  
  # Plot the new subnetwork for this cluster on the PDF
  numUnknownEdges <- sum(is.na(newSubnetEdges$Prior))
  numKnownEdges <- sum(!is.na(newSubnetEdges$Prior))
  plotTitle <- paste("Subnetwork from Cluster ", i, "\n(", length(allNodes), " genes; ", numKnownEdges, " known edges, ", numUnknownEdges, " unknown)", sep="")
  if(i %in% c(2,7,9,12,14))
    netLayout <- layout_with_lgl
  else if(i %in% c(4,5,10,11))
    netLayout <- layout_with_kk
  else if(i %in% c(1,3,6,8,13))
    netLayout <- layout_nicely
  plot(newSubnet, layout=netLayout, vertex.size=8, vertex.label.family="Helvetica",
       vertex.label.cex=0.5, edge.width=1.2, vertex.color=alpha("papayawhip", 0.85), vertex.frame.color="peachpuff3", 
       main=plotTitle, cex.main=1.5)
}
dev.off()

# For each of the 15 clusters, plot the time profiles of all genes in the cluster,
# but color only the genes which had at least one unknown edge
pdf("Output/GeneClustersWithUnknownAssociations.pdf", height=11, width=23)
par(mfrow=c(3,5))
plotColors <- c(brewer.pal(n=8, name="Dark2"), brewer.pal(n=8, name="Dark2"))
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
  colnames(subnetEdges) <- c("Gene1", "Gene2"); subnetEdges$Prior <- 0
  for(j in 1:nrow(subnetEdges)) {
    prior <- priorMatrix[subnetEdges[j,"Gene1"], subnetEdges[j,"Gene2"]]
    subnetEdges$Prior[j] <- prior
  }
  
  # Get the nodes (genes) with at least one edge
  nodesWithEdges <- unique(c(subnetEdges$Gene1, subnetEdges$Gene2))
  
  # Get the nodes (genes) with at least one unknown edge
  unknownEdges <- subnetEdges[is.na(subnetEdges$Prior),]
  nodesWithUnknownEdges <- unique(c(unknownEdges$Gene1, unknownEdges$Gene2))

  # Plot the time profiles
  plotTitle=paste("Cluster ", i, ":  ", length(nodesWithEdges), " genes with at least one edge\n(", length(nodesWithUnknownEdges), " with at least one unknown edge)", sep="")
  Plot.Gene.Group(nodesWithEdges[! nodesWithEdges %in% nodesWithUnknownEdges], monochrome="darkgray", points=F, 
                  plotTitle=plotTitle, titleSize=1.5, genesForExtrema=nodesWithEdges)
  if(length(nodesWithUnknownEdges) > 0) {
    Plot.Gene.Group(nodesWithUnknownEdges, monochrome="dodgerblue3", points=F, 
                    plotTitle=plotTitle, titleSize=1.5, add=T) 
  }
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

# --- Plots for manuscript ---

# --- Inflated LLR^2 (Figure 1) ---

# Using ggplot
plotColors <- c("dodgerblue2","orangered2")
grid.arrange(Plot.Gene.Group(c("Act87E", "bmm"), plotTitle="Gene A: Act87E;  Gene B: bmm", plotLegend=T, legendPos="bottom", gg=T, plotGrid=T, titleSize=12, plotColors=plotColors),
             Plot.Gene.Group(c("Act79B", "Mal-A7"), plotTitle="Gene A: Act79B;  Gene B: Mal-A7", plotLegend=T, legendPos="bottom", gg=T, plotGrid=T, titleSize=12, plotColors=plotColors),
             Plot.Gene.Group(c("IM3","tok"), plotTitle="Gene A: IM3;  Gene B: tok", plotLegend=T, legendPos="bottom", gg=T, plotGrid=T, titleSize=12, plotColors=plotColors), ncol=3)

# Not using ggplot
par(mfrow=c(1,3))
Plot.Gene.Group(c("Act87E", "bmm"), plotTitle="Gene A: Act87E;  Gene B: bmm", titleSize=1.4, plotLegend=T, legendPos="topleft", gg=F, grid=T, plotColors=plotColors)
Plot.Gene.Group(c("Mal-A7","Act79B"), plotTitle="Gene A: Mal-A7;  Gene B: Act79B", titleSize=1.4, plotLegend=T, legendPos="topleft", gg=F, grid=T, plotColors=plotColors)
Plot.Gene.Group(c("IM3","tok"), plotTitle="Gene A: IM3;  Gene B: tok", titleSize=1.4, plotLegend=T, legendPos="topright", gg=F, grid=T, plotColors=plotColors)
par(mfrow=c(1,1))

# Print non-Bayes LLR2 values
nonBayesLLR2Mat["bmm","Act87E"]
nonBayesLLR2Mat["Mal-A7","Act79B"]
nonBayesLLR2Mat["tok","IM3"]

# --- Histograms of R^2 ---

# LLR2 distribution without priors
histData <- as.data.frame(nonBayesLLR2Mat[upper.tri(nonBayesLLR2Mat)])
colnames(histData) <- "LLR2"
h1 <- ggplot(histData, aes(x=LLR2)) + 
  geom_histogram(mapping=aes(y=..count../sum(..count..)), binwidth=0.02, color="skyblue3", fill="lightsteelblue1") +
  theme_bw() + xlab(TeX("Lead-lag $R^2$")) + ylab("Density") +
  ggtitle(TeX("Distribution of lead-lag $R^2$ values (without priors)")) +
  theme(plot.title = element_text(hjust = 0.5))
h1

# 95th percentile of this distribution = 0.9479654 (0.95)
quantile(histData$LLR2, 0.95)

# LLR2 distribution with priors
histData <- as.data.frame(bayesLLR2Mat[upper.tri(bayesLLR2Mat)])
colnames(histData) <- "LLR2"
h2 <- ggplot(histData, aes(x=LLR2)) + 
  geom_histogram(mapping=aes(y=..count../sum(..count..)), binwidth=0.02, color="skyblue3", fill="lightsteelblue1") +
  theme_bw() + xlab(TeX("Bayesian lead-lag $R^2$")) + ylab("Density") +
  ggtitle(TeX("Distribution of Bayesian lead-lag $R^2$ values")) +
  theme(plot.title = element_text(hjust = 0.5))
h2

# 95th percentile of this distribution = 0.7239886 (0.72)
quantile(histData$LLR2, 0.95)

# Side-by-side histograms
grid.arrange(h1, h2, ncol=2)

# LLR2 - LLR2.own distribution, with priors
histData <- as.data.frame(bayesLLR2Mat[upper.tri(bayesLLR2Mat)] - bayesLLR2Mat.own[upper.tri(bayesLLR2Mat.own)])
colnames(histData) <- "LLR2_diff"
h3 <- ggplot(histData, aes(x=LLR2_diff)) + 
  geom_histogram(mapping=aes(y=..count../sum(..count..)), binwidth=0.02, color="skyblue3", fill="lightsteelblue1") +
  theme_bw() + xlab(TeX("Bayesian LL$R^2$ - LL$R^2_{own}$")) + ylab("Density") +
  ggtitle(TeX("Distribution of Bayesian LL$R^2$ - LL$R^2_{own}$ values")) +
  theme(plot.title = element_text(hjust = 0.5))
h3

# 95th quantile of this distribution = 0.2167662 (0.22)
quantile(histData$LLR2_diff, 0.95)

# --- Known and uncharacterized genes (Figure 2) ---

grid.arrange(Plot.Gene.Group(c("per", "tim", "to", "vri", "CG11854", "CG18609", "Pdp1", "CG33511"), gg=T, plotGrid=T,plotLegend=T, plotTitle="Known and uncharacterized genes\n with circadian rhythm patterns"),
             Plot.Gene.Group(c("AttC", "DptA", "DptB", "Dro", "edin", "Mtk", "CG43920", "CR44404", "CR45045"), gg=T, plotGrid=T,plotLegend=T, plotTitle="Known and uncharacterized genes\n with immune response functions"), ncol=2)

# Also involved in immune response: "IM1", "IM14", "IM2", "IM23", "IM3", "IM33", "IM4", "IMPPP"

# --- Time profiles in each cluster, with ggplot ---

clusterColors <- c("darkorange3", "dodgerblue3", "forestgreen", "darkmagenta", "indianred2", "orange4", "navy", "red2", "blueviolet", "turquoise4", "olivedrab", "darkslategray", "antiquewhite4", "coral4", "goldenrod2")
monochrome <- T;  points <- F;  plotGrid <- T;  gg <- T;  titleSize <- 12
plotList <- list()
for(i in 1:length(table(subGroups))) {
  plotList[[i]] <- Plot.Gene.Group(geneNames[subGroups == i], plotTitle=paste("Cluster ", i, " (", table(subGroups)[i], " genes)", sep=""), plotColors=clusterColors[i], monochrome=monochrome, points=points, gg=gg, plotGrid=plotGrid, titleSize=titleSize)
}
ggsave(file="Clusters.pdf", arrangeGrob(grobs=plotList, ncol=4), width=12, height=9, units="in")

# --- Temporal profile plot of Imd and Toll regulated genes in cluster 7 ---

imdGenes <- c("AttA", "AttB", "AttC", "AttD", "Dro", "CecA2", "DptA", "DptB", "PGRP-SC2", "PGRP-SB1")
tollGenes <- c("PGRP-SA", "IM33", "IMPPP", "IM23", "IM1", "IM2", "IM4", "IM14", "IM3")
newGenes <- c("CR44404", "CG43236", "CG43202", "CG43920")
negativeGenes <- c("Acp1", "CG7214")
C7colors <- c(rep(alpha("orangered2", 0.65), 10), rep(alpha("dodgerblue3", 0.65), 9), rep(alpha("black", 0.65), 4), rep(alpha("forestgreen", 0.65), 2))
Plot.Gene.Group(c(imdGenes, tollGenes, newGenes, negativeGenes), plotColors=C7colors, 
                plotGrid=T, gg=T, points=T, lineLabels=F, plotTitle="Selected genes from cluster 7")

# --- Network of unknown genes in cluster 7 ---

# Get the entire network
adjBayes <- (bayesLLR2Mat > 0.9) + 0
networkBayes <- graph_from_adjacency_matrix(adjBayes , mode='undirected', diag=F)

# Get all the neighbors of the unknown genes in cluster 7
unknownC7neighbors <- adjacent_vertices(networkBayes, newGenes)

# Collapse this list of neighbors into one vector of nodes
allNodes <- c()
for(j in 1:length(unknownC7neighbors))
  allNodes <- c(allNodes, names(unknownC7neighbors[j]), names(unknownC7neighbors[[j]]))
allNodes <- unique(allNodes)

# Form a new subnetwork out of allNodes
C7adj <- (bayesLLR2Mat[allNodes, allNodes] > 0.9) + 0
C7net <- graph_from_adjacency_matrix(C7adj , mode='undirected', diag=F)

# Get a new dataframe of edges corresponding to allNodes
C7edges <- data.frame(as_edgelist(C7net))
colnames(C7edges) <- c("Gene1", "Gene2"); C7edges$Prior <- 0
for(j in 1:nrow(C7edges)) {
  prior <- priorMatrix[C7edges[j,"Gene1"], C7edges[j,"Gene2"]]
  C7edges$Prior[j] <- prior
}

# Edges between genes with unknown associations will be blue, and edges 
# for known associations will be red
E(C7net)$color[is.na(C7edges$Prior)] <- alpha('blue', 0.7)
E(C7net)$color[!is.na(C7edges$Prior)] <- alpha('red', 0.7)

# The four novel genes will be colored differently
V(C7net)$color[as_ids(V(C7net)) %in% newGenes] <- alpha("navajowhite", 0.95)
V(C7net)$color[! as_ids(V(C7net)) %in% newGenes] <- alpha("linen", 0.85)

# Node frame colors
V(C7net)$frame.color[as_ids(V(C7net)) %in% newGenes] <- "peachpuff3"
V(C7net)$frame.color[! as_ids(V(C7net)) %in% newGenes] <- "gray66"

# Plot the new subnetwork for cluster 7 on a PDF
numUnknownEdges <- sum(is.na(C7edges$Prior))
numKnownEdges <- sum(!is.na(C7edges$Prior))
pdf("Output/Cluster7Subnetwork.pdf", height=6, width=6)
plotTitle <- paste("New relationships detected in cluster 7\n(", length(allNodes), " genes; ", numKnownEdges, " previously-known edges, ", numUnknownEdges, " newly-identified edges)", sep="")
plot(C7net, layout=layout_with_lgl, vertex.size=21, vertex.label.family="Helvetica",
     vertex.label.cex=0.5, edge.width=1.2)
title(plotTitle, cex.main=0.8, font.main=1)
dev.off()

# --- Temporal profile plot selected genes in cluster 12 ---

C12genes <- c("fbp", "to", "AGBE", "Galk", "Gba1b", "CG11594", "CG10469", "CG13315")
C12colors <- c("black", "dodgerblue2", rep("orangered2", 3), rep("springgreen4", 3))
Plot.Gene.Group(C12genes, plotColors=C12colors, plotGrid=T, lineLabels=T,
                plotTitle="Selected genes from cluster 12", pointSize=1,
                plotLegend=T)

# --- Network of neighbors of gene "fbp" in cluster 12 ---

# Get all the neighbors of fbp
fbpNeighbors <- adjacent_vertices(networkBayes, "fbp")
fbpNeighbors <- c("fbp", names(fbpNeighbors[[1]]))

# Form a new subnetwork out of fbp and neighborrs
fbpAdj <- (bayesLLR2Mat[fbpNeighbors, fbpNeighbors] > 0.9) + 0
fbpNet <- graph_from_adjacency_matrix(fbpAdj , mode='undirected', diag=F)

# Get a new dataframe of edges corresponding to fbpNeighbors
fbpEdges <- data.frame(as_edgelist(fbpNet))
colnames(fbpEdges) <- c("Gene1", "Gene2"); fbpEdges$Prior <- 0
for(j in 1:nrow(fbpEdges)) {
  prior <- priorMatrix[fbpEdges[j,"Gene1"], fbpEdges[j,"Gene2"]]
  fbpEdges$Prior[j] <- prior
}

# Edges between genes with unknown associations will be blue, and edges 
# for known associations will be red
E(fbpNet)$color[is.na(fbpEdges$Prior)] <- alpha('blue', 0.1)
E(fbpNet)$color[!is.na(fbpEdges$Prior)] <- alpha('red', 0.1)

# The fbp gene will be colored differently
V(fbpNet)$color[as_ids(V(fbpNet)) == "fbp"] <- "navajowhite"
V(fbpNet)$color[! as_ids(V(fbpNet)) == "fbp"] <- alpha("linen", 0.85)

# Node frame colors
V(fbpNet)$frame.color[as_ids(V(fbpNet)) == "fbp"] <- "peachpuff3"
V(fbpNet)$frame.color[! as_ids(V(fbpNet)) == "fbp"] <- "gray66"

# Node sizes (a few highlighted nodes will be larger)
fbpNetHighlightNodes <- c("fbp", "Galk", "AGBE", "Gba1b", "CG11594", "CG10469", "CG13315") 
V(fbpNet)$size[as_ids(V(fbpNet)) %in% fbpNetHighlightNodes] <- 17
V(fbpNet)$size[! as_ids(V(fbpNet)) %in% fbpNetHighlightNodes] <- 7

# Node names (text only for highlighted nodes)
V(fbpNet)$label[! as_ids(V(fbpNet)) %in% fbpNetHighlightNodes] <- NA
V(fbpNet)$label[as_ids(V(fbpNet)) %in% fbpNetHighlightNodes] <- fbpNetHighlightNodes

# Darken edges for the highlighted nodes
E(fbpNet)$color[which(fbpEdges$Gene1 == "fbp" & fbpEdges$Gene2 %in% fbpNetHighlightNodes)] <- "blue"

# Plot the new subnetwork for fbp on a PDF
numUnknownEdges <- sum(is.na(fbpEdges$Prior))
numKnownEdges <- sum(!is.na(fbpEdges$Prior))
pdf("Output/fbpSubnetwork.pdf", height=6, width=6)
plotTitle <- paste("Neighbors of gene \"fbp\" in cluster 12\n(", length(fbpNeighbors), " genes; ", numKnownEdges, " known edges, ", numUnknownEdges, " newly-identified edges)", sep="")
plot(fbpNet, layout=layout_with_kk, vertex.label.family="Helvetica",
     vertex.label.cex=0.5, edge.width=1.2)
title(plotTitle, cex.main=0.8, font.main=1)
dev.off()

# --- Small-scale case study (immune response/metabolism) ---

immMetGenes <- c("IM1", "IM2", "FASN1", "UGP", "mino", "fbp")
Plot.Gene.Group(immMetGenes, plotGrid=T, plotLegend=T, plotTitle="Subset of genes involved in immune response\n (IM1, IM2) and metabolism (FASN1, UGP, mino, fbp)", titleSize=13, legendPos="right")
priorMatrix[immMetGenes, immMetGenes]
round(bayesLLR2Mat[immMetGenes, immMetGenes], 2)

# --- R^2 scatterplots ---

Draw.R2.Scatterplot(bayesLLR2Mat.other, bayesLLR2Mat-bayesLLR2Mat.own, priorMatrix, geneSubset, T, F)
Draw.R2.Scatterplot(nonBayesLLR2Mat.other, nonBayesLLR2Mat-nonBayesLLR2Mat.own, priorMatrix, geneSubset, F, F)
  


