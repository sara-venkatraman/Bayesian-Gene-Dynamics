# ======= Script 3: Results (clustering and network analysis) =======

# This script produces all tables and figures in "An empirical Bayes 
# approach to estimating dynamic models of co-regulated gene expression" 
# (Venkatraman et al., 2021). Figure/table numbers refer to our preprint: 
# https://arxiv.org/abs/2112.15326 or 
# https://www.biorxiv.org/content/10.1101/2021.07.08.451684v2

# ======= Setup: load gene expression data and R scripts ======= 

# Load LLR2 and plotting scripts
source("1_LLR2.R")
source("2_PlottingFunctions.R")

# Read gene expression data and prior adjacency matrix. Save gene names
# and make sure row and column names of the adjacency matrix are the same.
geneData <- read.csv("Data/GeneData.csv", row.names=1)
geneNames <- rownames(geneData)
priorMatrix <- read.csv("Data/PriorMatrix.csv", row.names=1)
colnames(priorMatrix) <- rownames(priorMatrix)

# Define hours corresponding to each time point
hours <- c(0, 1, 2, 4, 5, 6, 8, 10, 12, 14, 16, 20, 24, 30, 36, 42, 48)

# ========================================================================
# Run the next two lines to produce lead-lag R^2 similarity matrices *only*
# if they have not already been generated. Otherwise, skip to the next block
# of code to read the LLR2 matrices from CSV files.
bayesLLR2 <- LLR2(geneData, hours, bayes=TRUE, priorMatrix=priorMatrix, writeToCSV=TRUE)
nonBayesLLR2 <- LLR2(geneData, hours, bayes=FALSE, writeToCSV=TRUE)
# ========================================================================

# Read the non-Bayesian lead-lag R^2 matrices
nonBayesLLR2Mat <- read.csv("LLR2 Matrices/NonBayesLLR2.csv", row.names=1)
nonBayesLLR2Mat.other <- read.csv("LLR2 Matrices/NonBayesLLR2_Other.csv", row.names=1)
nonBayesLLR2Mat.own <- read.csv("LLR2 Matrices/NonBayesLLR2_Own.csv", row.names=1)
colnames(nonBayesLLR2Mat) <- rownames(nonBayesLLR2Mat)
colnames(nonBayesLLR2Mat.other) <- rownames(nonBayesLLR2Mat.other)
colnames(nonBayesLLR2Mat.own) <- rownames(nonBayesLLR2Mat.own)

# Read the Bayesian lead-lag R^2 matrices
bayesLLR2Mat <- read.csv("LLR2 Matrices/BayesLLR2.csv", row.names=1)
bayesLLR2Mat.other <- read.csv("LLR2 Matrices/BayesLLR2_Other.csv", row.names=1)
bayesLLR2Mat.own <- read.csv("LLR2 Matrices/BayesLLR2_Own.csv", row.names=1)
colnames(bayesLLR2Mat) <- rownames(bayesLLR2Mat)
colnames(bayesLLR2Mat.other) <- rownames(bayesLLR2Mat.other)
colnames(bayesLLR2Mat.own) <- rownames(bayesLLR2Mat.own)

# ======= Figure 1: Non-Bayesian network statistics and plot =======

# Get the entire network of genes (edge defined for LLR2 > 0.9)
adjacencyNonBayes <- (nonBayesLLR2Mat > 0.9) + 0
networkNonBayes <- graph_from_adjacency_matrix(adjacencyNonBayes, mode='undirected', diag=F)
edgesNonBayes <- data.frame(as_edgelist(networkNonBayes))
colnames(edgesNonBayes) <- c("Gene1", "Gene2")

# Get total number of edges
nrow(edgesNonBayes)

# How many genes are in the largest connected component? 
max(clusters(networkNonBayes, mode="strong")$csize)

# Non-Bayesian network plot
pdf("NonBayesNetwork.pdf", height=4, width=4)
plot(networkNonBayes, vertex.label=NA, vertex.size=3, vertex.frame.color=alpha("darkblue",0.45), edge.width=0.2, edge.color=alpha('darkgray', 0.7))
title("Network of genes formed\n without prior information", cex.main=0.7, font.main=1)
dev.off()

# ======= Figure 1: Bayesian network statistics and plot =======

# Get the entire network of genes (edge defined for LLR2 > 0.9). Classify
# edges as 0, 1, or NA according to prior.
adjacencyBayes <- (bayesLLR2Mat > 0.9) + 0
networkBayes <- graph_from_adjacency_matrix(adjacencyBayes, mode='undirected', diag=F)
edgesBayes <- data.frame(as_edgelist(networkBayes))
colnames(edgesBayes) <- c("Gene1", "Gene2")
edgesBayes$Prior <- apply(edgesBayes, 1, function(x) priorMatrix[x[1], x[2]])

# Get total number of edges
nrow(edgesBayes)

# Get number of known, unknown, and unlikely edges, respectively
sum(edgesBayes$Prior == 1, na.rm=TRUE)
sum(is.na(edgesBayes$Prior))
sum(edgesBayes$Prior == 0, na.rm=TRUE)

# How many genes are in the largest connected component? 
max(clusters(networkBayes, mode="strong")$csize)

# Bayesian network plot
E(networkBayes)$color[is.na(edgesBayes$Prior)] <- 'blue'
E(networkBayes)$color[!is.na(edgesBayes$Prior)] <- 'firebrick3'
E(networkBayes)$width[is.na(edgesBayes$Prior)] <- 0.6
E(networkBayes)$width[!is.na(edgesBayes$Prior)] <- 0.4
pdf("BayesNetwork.pdf", height=4, width=4)
plot(networkBayes, layout=layout_with_graphopt, vertex.label=NA, vertex.size=2.5, vertex.frame.color=alpha("darkblue",0.35))
title("Network of genes formed\n with prior information", cex.main=0.6, font.main=1)
dev.off()

# ======= Figure 2: Inflated LLR2 values =======

p1 <- Plot.Genes(geneData[c("Act87E", "bmm"),], hours,
                 plotColors=c("dodgerblue2","orangered2"),
                 plotTitle="<b>Genes <i>Act87E, bmm</i></b>",
                 plotSubtitle=paste("Lead-lag R<sup>2</sup> = ", round(nonBayesLLR2Mat["bmm", "Act87E"], 3), "<br>Bayesian lead-lag R<sup>2</sup> = ", round(bayesLLR2Mat["bmm", "Act87E"], 3), sep=""),
                 axisLabels=list(x="Time (hours after infection)", y="Expression (log<sub>2</sub>-fold change)"))

p2 <- Plot.Genes(geneData[c("Act79B", "Mal-A7"),], hours,
                 plotColors=c("dodgerblue2","orangered2"),
                 plotTitle="<b>Genes <i>Act79B, Mal-A7</i></b>",
                 plotSubtitle=paste("Lead-lag R<sup>2</sup> = ", round(nonBayesLLR2Mat["Act79B", "Mal-A7"], 3), "<br>Bayesian lead-lag R<sup>2</sup> = ", round(bayesLLR2Mat["Act79B", "Mal-A7"], 3), "0", sep=""),
                 axisLabels=list(x="Time (hours after infection)", y="Expression (log<sub>2</sub>-fold change)"))

p3 <- Plot.Genes(geneData[c("BomS3", "tok"),], hours,
                 plotColors=c("dodgerblue2","orangered2"),
                 plotTitle="<b>Genes <i>BomS3, tok</i></b>",
                 plotSubtitle=paste("Lead-lag R<sup>2</sup> = ", round(nonBayesLLR2Mat["BomS3", "tok"], 3), "<br>Bayesian lead-lag R<sup>2</sup> = ", round(bayesLLR2Mat["BomS3", "tok"], 3), sep=""),
                 axisLabels=list(x="Time (hours after infection)", y="Expression (log<sub>2</sub>-fold change)"))

ggsave(file="InflatedLLR2.pdf", arrangeGrob(grobs=list(p1,p2,p3), ncol=3), width=12, height=4.25, units="in")

# ======= Figure 3: Circadian rhythm, immune response pathways =======

p1 <- Plot.Genes(geneData[c("per", "tim", "to", "vri", "CG11854", "CG18609", "Pdp1", "CG33511"),],
                 hours, plotTitle="Known and uncharacterized genes<br> with circadian rhythm patterns",
                 axisLabels=list(x="Time (hours after infection)", y="Expression (log<sub>2</sub>-fold change)"))

p2 <- Plot.Genes(geneData[c("AttC", "DptA", "DptB", "Dro", "edin", "Mtk", "CG43920", "CG44404", "CG45045"),],
                 hours, plotTitle="Known and uncharacterized genes<br> with immune response functions",
                 axisLabels=list(x="Time (hours after infection)", y="Expression (log<sub>2</sub>-fold change)"))

ggsave(file="CircadianAndImmunePatterns.pdf", arrangeGrob(grobs=list(p1,p2), ncol=2), width=10, height=4.5, units="in")

# ======= Figure 4: Immunity/metabolism case study =======

immMetGenes <- c("IM1", "IM2", "FASN1", "UGP", "mino", "fbp")
p <- Plot.Genes(geneData[immMetGenes,], hours, axisLabels=list(x="Time (hours after infection)", y="Expression (log<sub>2</sub>-fold change)"),
                plotTitle="Selected genes involved in immune response<br> (<i>IM1, IM2</i>) and metabolism (<i>FASN1, UGP, mino, fbp</i>)") + 
  guides(color=guide_legend(nrow=1))
ggsave(file="ImmuneMetabolic.pdf", p, width=5.8, height=4.48, units="in")

# ======= Tables 1-5 (4 and 5 in Appendix): Immunity/metabolism case study =======

immMetGenes <- c("IM1", "IM2", "FASN1", "UGP", "mino", "fbp")

# Table 1: print section of prior adjacency matrix
priorMatrix[immMetGenes, immMetGenes]

# Table 2: print section of Bayesian LLR2 matrix
round(bayesLLR2Mat[immMetGenes, immMetGenes], 2)

# Table 3: print section of non-Bayesian LLR2 matrix
round(nonBayesLLR2Mat[immMetGenes, immMetGenes], 2)

# Get 95th percentile of empirical distributions of the LLR2 values
quantile(bayesLLR2Mat[upper.tri(bayesLLR2Mat)], 0.95)
quantile(nonBayesLLR2Mat[upper.tri(nonBayesLLR2Mat)], 0.95)

# Table 4: print section of Bayesian LLR2 - LLR2.own matrix
round(bayesLLR2Mat[immMetGenes, immMetGenes] - bayesLLR2Mat.own[immMetGenes, immMetGenes], 2)

# Table 5: print section of non-Bayesian LLR2 - LLR2.own matrix
round(nonBayesLLR2Mat[immMetGenes, immMetGenes] - nonBayesLLR2Mat.own[immMetGenes, immMetGenes], 2)

# Get 95th percentile of empirical distributions of LLR2 - LLR2.own values
quantile((bayesLLR2Mat - bayesLLR2Mat.own)[upper.tri(bayesLLR2Mat)], 0.95)
quantile((nonBayesLLR2Mat - nonBayesLLR2Mat.own)[upper.tri(nonBayesLLR2Mat)], 0.95)

# ======= Figure 5: Hierarchical clustering with the Bayesian LLR2 =======

# Cluster via Ward's method and check dendrogram 
hierClust <- hclust(as.dist(1 - bayesLLR2Mat), method="ward.D")
plot(hierClust, cex=0.2, col="lightsteelblue4")

# Cutting the dendrogram at height 10 yields 12 clusters
subGroups <- cutree(hierClust, h=10)
table(subGroups)

# Find min, max, median, and mean cluster sizes
min(table(subGroups))
max(table(subGroups))
median(table(subGroups))
mean(table(subGroups))

# Plot the expression trajectories of genes in each cluster
plotList <- list()
plotColors <- c("darkorange3", "dodgerblue3", "forestgreen", "darkmagenta",
                "turquoise4", "orange4", "navy", "darkgoldenrod3", 
                "blueviolet", "red3", "olivedrab", "darkslategray", 
                "antiquewhite4", "coral4", "goldenrod2")
for(i in 1:length(table(subGroups))) {
  numGenes <- table(subGroups)[i]
  plotList[[i]] <- Plot.Genes(geneData[subGroups == i,], hours, points=FALSE,
                              plotColors=rep(plotColors[i], numGenes),
                              plotTitle=paste("Cluster ", i, " (", numGenes, " genes)", sep=""),
                              plotLegend=FALSE, lineOpacity=0.2, 
                              axisLabels=list(x="Time (hours)", y="Expression"))
}
ggsave(file="ClustersLLR2.pdf", arrangeGrob(grobs=plotList, ncol=4), width=12, height=6.75, units="in")

# ======= Figure 6: Cluster 10 expression trajectories =======

imdGenes <- c("AttA", "AttB", "AttC", "Dro", "CecA2", "DptA", "DptB", "PGRP-SB1")
tollGenes <- c("PGRP-SA", "IMPPP", "IM1", "IM2", "IM4", "IM14", "IM23", "BomS3")
newGenes <- c("CG44404", "CG43236", "CG43202", "CG43920")
C10colors <- c(rep(alpha("orangered2", 0.65), 8), rep(alpha("dodgerblue3", 0.65), 8), rep(alpha("black", 0.75), 4))
p <- Plot.Genes(geneData[c(imdGenes, tollGenes, newGenes),1:8], hours[1:8],
                plotColors=C10colors, plotLegend=FALSE,
                plotTitle="<b>Selected genes from cluster 10</b>",
                plotSubtitle="(<span style='color:orangered2;'>Imd-regulated</span><span style='color:white;'> l</span><span style='color:orangered2;'>genes</span>; <span style='color:dodgerblue3;'>Toll-regulated genes</span>;<br><span style='color:black;'>genes potentially associated with Imd signaling</span>)",
                axisLabels=list(x="Time (hours after infection)", y="Expression (log<sub>2</sub>-fold change)"))
ggsave(file="Cluster10.pdf", p, width=6, height=4.5, units="in")

# ======= Figure 7: Cluster 10 network =======
# Note: to run this, first run the hierarchical clustering code under Figure 5.

newGenes <- c("CG44404", "CG43236", "CG43202", "CG43920")

# Get the entire network of genes (edge defined for LLR2 > 0.9)
adjacencyBayes <- (bayesLLR2Mat > 0.9) + 0
networkBayes <- graph_from_adjacency_matrix(adjacencyBayes , mode='undirected', diag=FALSE)

# Get the names of all neighbors of the unknown genes in cluster 10
C10neighbors <- adjacent_vertices(networkBayes, newGenes)
C10neighbors <- unique(geneNames[unlist(C10neighbors, use.names=FALSE)])

# Form a new subnetwork out of C10neighbors
C10adjacency <- (bayesLLR2Mat[C10neighbors, C10neighbors] > 0.9) + 0
C10network <- graph_from_adjacency_matrix(C10adjacency, mode='undirected', diag=FALSE)

# Get the priors associated with each edge in this subnetwork
C10edges <- data.frame(as_edgelist(C10network))
colnames(C10edges) <- c("Gene1", "Gene2")
C10edges$Prior <- unlist(apply(C10edges, 1, function(x) priorMatrix[x[1], x[2]]))

# Blue edges between genes with unknown associations,
# red edges between genes with known associations
E(C10network)$color[is.na(C10edges$Prior)] <- alpha('blue', 0.2)
E(C10network)$color[!is.na(C10edges$Prior)] <- alpha('red', 0.2)

# The four genes of interest will be colored differently
V(C10network)$color[as_ids(V(C10network)) %in% newGenes] <- alpha("navajowhite", 0.95)
V(C10network)$color[! as_ids(V(C10network)) %in% newGenes] <- alpha("linen", 0.85)

# Set the node outline colors
V(C10network)$frame.color[as_ids(V(C10network)) %in% newGenes] <- "peachpuff3"
V(C10network)$frame.color[! as_ids(V(C10network)) %in% newGenes] <- "gray66"

# Darken the edges that connect genes within cluster 10
E(C10network)$color[C10edges$Gene1 %in% geneNames[subGroups == 10] & C10edges$Gene2 %in% geneNames[subGroups == 10] & is.na(C10edges$Prior)] <- alpha('blue', 0.8)
E(C10network)$color[C10edges$Gene1 %in% geneNames[subGroups == 10] & C10edges$Gene2 %in% geneNames[subGroups == 10] & !is.na(C10edges$Prior)] <- alpha('red', 0.8)

# Plot the subnetwork on a PDF (note: network layout is random)
numUnknownEdges <- sum(is.na(C10edges$Prior))
numKnownEdges <- sum(!is.na(C10edges$Prior))
pdf("Cluster10Network.pdf", height=6, width=6)
plotTitle <- paste("New relationships detected in cluster 10\n(", length(C10neighbors), " genes; ", numKnownEdges, " edges previously known, ", numUnknownEdges, " edges newly identified)", sep="")
plot(C10network, layout=layout_with_dh, vertex.size=21, vertex.label.family="Helvetica",
     vertex.label.cex=0.5, edge.width=1.2, vertex.label.font=2)
title(plotTitle, cex.main=0.8, font.main=1)
dev.off()

# ======= Figure 8: LLR2 scatterplots =======

# First run this if using LLR2 matrices that were read from CSVs:
bayesLLR2 <- list(LLR2Mat=bayesLLR2Mat, LLR2Mat.other=bayesLLR2Mat.other, LLR2Mat.own=bayesLLR2Mat.own)
nonBayesLLR2 <- list(LLR2Mat=nonBayesLLR2Mat, LLR2Mat.other=nonBayesLLR2Mat.other, LLR2Mat.own=nonBayesLLR2Mat.own)

# Select random sample of 150 genes and draw scatterplots
geneSubset <- sample(geneNames, size=150)
p1 <- Plot.LLR2.Scatterplot(bayesLLR2, priorMatrix, geneSubset, interactive=F,
                            plotTitle="Bayesian lead-lag R<sup>2</sup> values")
p2 <- Plot.LLR2.Scatterplot(nonBayesLLR2, priorMatrix, geneSubset, interactive=F,
                            plotTitle="Non-Bayesian lead-lag R<sup>2</sup> values")
ggsave(file="R2Scatterplots.pdf", arrangeGrob(grobs=list(p1,p2), ncol=2), width=10, height=4.3, units="in")

# ======= Figure 9 (Appendix): Unknown genes in middle/upper-right of Bayesian LLR2 scatterplot =======

plotColors <- c("dodgerblue2", "orangered2", "goldenrod", "forestGreen")

p1 <- Plot.Genes(geneData[c("CR42868", "AttD", "CG9616"),], hours, 
                 plotTitle="Genes: <i>CR42868, AttD, CG9616</i>", plotColors=plotColors,
                 axisLabels=list(x="Time (hours after infection)", y="Expression (log<sub>2</sub>-fold change)"))

p2 <- Plot.Genes(geneData[c("Spn28Dc", "CR43364", "CR42715", "scb"),], hours,
                 plotTitle="Genes: <i>Spn28Dc, CR43364, CR42715, scb</i>", plotColors=plotColors,
                 axisLabels=list(x="Time (hours after infection)", y="Expression (log<sub>2</sub>-fold change)"))

p3 <- Plot.Genes(geneData[c("ACC", "Idh", "GstE9"),], hours, 
                 plotTitle="Genes: <i>ACC, Idh, GstE9</i>", plotColors=plotColors,
                 axisLabels=list(x="Time (hours after infection)", y="Expression (log<sub>2</sub>-fold change)"))

ggsave(file="R2ScatterplotsMiddle.pdf", arrangeGrob(grobs=list(p1,p2,p3), ncol=3), width=13, height=3.9, units="in")

# ======= Figure 10 (Appendix): Genes in upper-right region of Bayesian LLR2 scatterplot ---

plotColors <- c("dodgerblue2", "orangered2", "goldenrod", "forestGreen")

p1 <- Plot.Genes(geneData[c("alphaTry", "gammaTry", "CG30025"),], hours, 
                 plotTitle="Genes: <i>alphaTry, gammaTry, CG30025</i>", plotColors=plotColors,
                 axisLabels=list(x="Time (hours after infection)", y="Expression (log<sub>2</sub>-fold change)"))

p2 <- Plot.Genes(geneData[c("RpS26", "RpS6", "RpL13", "RpL7"),], hours, 
                 plotTitle="Genes: <i>RpS26, RpS6, RpL13, RpL7</i>", plotColors=plotColors,
                 axisLabels=list(x="Time (hours after infection)", y="Expression (log<sub>2</sub>-fold change)"))

p3 <- Plot.Genes(geneData[c("CG13096", "l(1)G0020", "Nop56"),], hours, 
                 plotTitle="Genes: <i>CG13096, l(1)G0020, Nop56</i>", plotColors=plotColors,
                 axisLabels=list(x="Time (hours after infection)", y="Expression (log<sub>2</sub>-fold change)"))

ggsave(file="R2ScatterplotsRight.pdf", arrangeGrob(grobs=list(p1,p2,p3), ncol=3), width=13, height=3.9, units="in")

# ======= Figure 11 (Appendix): Hierarchical clustering with the Non-Bayesian LLR2 =======

# Cluster via Ward's method and check dendrogram 
hierClust <- hclust(as.dist(1 - nonBayesLLR2Mat), method="ward.D")
plot(hierClust, cex=0.2, col="lightsteelblue4")

# Divide into 12 clusters
subGroups <- cutree(hierClust, k=12)
table(subGroups)

# Plot the expression trajectories of genes in each cluster
plotList <- list()
plotColors <- c("darkorange3", "dodgerblue3", "forestgreen", "darkmagenta",
                "turquoise4", "orange4", "navy", "darkgoldenrod3", 
                "blueviolet", "red3", "olivedrab", "darkslategray", 
                "antiquewhite4", "coral4", "goldenrod2")
for(i in 1:length(table(subGroups))) {
  numGenes <- table(subGroups)[i]
  plotList[[i]] <- Plot.Genes(geneData[subGroups == i,], hours, points=FALSE,
                              plotColors=rep(plotColors[i], numGenes),
                              plotTitle=paste("Cluster ", i, " (", numGenes, " genes)", sep=""),
                              plotLegend=FALSE, lineOpacity=0.2,
                              axisLabels=list(x="Time (hours)", y="Expression"))
}
ggsave(file="ClustersNonBayesLLR2.pdf", arrangeGrob(grobs=plotList, ncol=4), width=12, height=6.75, units="in")

# ======= Figure 12 (Appendix): Hierarchical clustering with correlation =======

# Cluster via Ward's method and check dendrogram 
hierClust <- hclust(as.dist(1-cor(t(geneData))), method="ward.D")
plot(hierClust, cex=0.2, col="lightsteelblue4")

# Divide into 12 clusters
subGroups <- cutree(hierClust, k=12)
table(subGroups)

# Plot the expression trajectories of genes in each cluster
plotList <- list()
plotColors <- c("darkorange3", "dodgerblue3", "forestgreen", "darkmagenta",
                "turquoise4", "orange4", "navy", "darkgoldenrod3", 
                "blueviolet", "red3", "olivedrab", "darkslategray", 
                "antiquewhite4", "coral4", "goldenrod2")
for(i in 1:length(table(subGroups))) {
  numGenes <- table(subGroups)[i]
  plotList[[i]] <- Plot.Genes(geneData[subGroups == i,], hours, points=FALSE,
                              plotColors=rep(plotColors[i], numGenes),
                              plotTitle=paste("Cluster ", i, " (", numGenes, " genes)", sep=""),
                              plotLegend=FALSE, lineOpacity=0.2,
                              axisLabels=list(x="Time (hours)", y="Expression"))
}
ggsave(file="ClustersCorrelation.pdf", arrangeGrob(grobs=plotList, ncol=4), width=12, height=6.75, units="in")

# ======= Figure 13 (Appendix): Cluster 2 expression trajectories =======

circadian <- c("per", "Clk", "vri", "Pdp1")
visual <- c("Rh5", "Rh6", "Pdh")
sodium <- c("Smvt", "salt")
C2colors <- c(rep("black", 2), "darkorange", "red", rep("pink2", 2), "chartreuse4", rep("purple", 2))
p <- Plot.Genes(geneData[c(circadian, visual, sodium),], hours, plotColors=C2colors,
                plotTitle="<b>Selected genes from cluster 2<b>", plotLegend=FALSE,
                plotSubtitle="(Regulators of circadian clock: <i>per, Clk, <span style='color:darkorange;'>vri</span>, <span style='color:red;'>Pdp1</i></span>;<br><span style='color:lightpink3;'>rhodopsins: <i>Rh5, Rh6</i></span>; <span style='color:chartreuse4;'>Pdh</span>; <span style='color:purple;'>sodium transporters: <i>Smvt, salt</i></span>)",
                axisLabels=list(x="Time (hours after infection)", y="Expression (log<sub>2</sub>-fold change)"))
ggsave(file="Cluster2.pdf", p, width=6, height=4.5, units="in")

# ======= Figure 14 (Appendix): Cluster 4 expression trajectories =======
# Note: to run this, first run the hierarchical clustering code under Figure 5.

# Get genes in fbp network: cluster 4 genes X such that LLR2(X, fbp) > 0.9
C4highlights <- c("Gale", "AGBE", "Gba1b", "fit", "fbp")
C4fbpNeighbors <- geneNames[subGroups == 4 & bayesLLR2Mat[,"fbp"] > 0.9]
C4fbpNeighbors <- C4fbpNeighbors[! C4fbpNeighbors %in% C4highlights]

# Plot trajectories
C4colors <- c(rep("snow3", length(C4fbpNeighbors)), rep("springgreen4", 3), "blue3", "red3")
p <- Plot.Genes(geneData[c(C4fbpNeighbors, C4highlights),], hours, plotColors=C4colors,
                plotLegend=FALSE, plotTitle="<b>Selected genes from cluster 4</b>",
                plotSubtitle="(<span style='color:snow4;'>genes with which <i>fbp</i> has LLR<sup>2</sup> > 0.9</span>;<br><span style='color:springgreen4;'>genes involved in carbohydrate metabolism</span>; <span style='color:red3;'><i>fbp</i></span>; <span style='color:blue3;'><i>fit</i></span>)",
                axisLabels=list(x="Time (hours after infection)", y="Expression (log<sub>2</sub>-fold change)"))
ggsave(file="Cluster4.pdf", p, width=6, height=4.5, units="in")

# ======= Figure 15 (Appendix): Cluster 4 subnetwork - neighbors of gene "fbp" =======
# Note: to run this, first run the hierarchical clustering code under Figure 5.

# Get the entire network of genes (edge defined for LLR2 > 0.9)
adjacencyBayes <- (bayesLLR2Mat > 0.9) + 0
networkBayes <- graph_from_adjacency_matrix(adjacencyBayes , mode='undirected', diag=FALSE)

# Get all neighbors of fbp
fbpNeighbors <- adjacent_vertices(networkBayes, "fbp")
fbpNeighbors <- c("fbp", names(fbpNeighbors[[1]]))

# Form a new subnetwork out of fbp and its neighbors
fbpAdjacency <- (bayesLLR2Mat[fbpNeighbors, fbpNeighbors] > 0.9) + 0
fbpNetwork <- graph_from_adjacency_matrix(fbpAdjacency , mode='undirected', diag=F)

# Get the priors associated with each edge in this subnetwork
fbpEdges <- data.frame(as_edgelist(fbpNetwork))
colnames(fbpEdges) <- c("Gene1", "Gene2"); fbpEdges$Prior <- 0
fbpEdges$Prior <- unlist(apply(fbpEdges, 1, function(x) priorMatrix[x[1], x[2]]))

# Blue edges between genes with unknown associations,
# red edges between genes with known associations
E(fbpNetwork)$color[is.na(fbpEdges$Prior)] <- alpha('blue', 0.08)
E(fbpNetwork)$color[!is.na(fbpEdges$Prior)] <- alpha('red', 0.08)

# Darken the edges that connect fbp to cluster 4 genes
E(fbpNetwork)$color[fbpEdges$Gene1 %in% geneNames[subGroups == 4] & fbpEdges$Gene2 %in% geneNames[subGroups == 4] & is.na(fbpEdges$Prior)] <- alpha('blue', 0.6)
E(fbpNetwork)$color[fbpEdges$Gene1 %in% geneNames[subGroups == 4] & fbpEdges$Gene2 %in% geneNames[subGroups == 4] & !is.na(fbpEdges$Prior)] <- alpha('red', 0.3)

# The fbp gene will be colored differently
V(fbpNetwork)$color[as_ids(V(fbpNetwork)) == "fbp"] <- "navajowhite"
V(fbpNetwork)$color[! as_ids(V(fbpNetwork)) == "fbp"] <- alpha("linen", 0.85)

# Set the node outline colors
V(fbpNetwork)$frame.color[as_ids(V(fbpNetwork)) == "fbp"] <- "peachpuff3"
V(fbpNetwork)$frame.color[! as_ids(V(fbpNetwork)) == "fbp"] <- "gray66"

# Set the node sizes (a few highlighted nodes will be larger)
fbpNetHighlightNodes <- c("fbp", "Gale", "AGBE", "Gba1b") 
V(fbpNetwork)$size[as_ids(V(fbpNetwork)) %in% fbpNetHighlightNodes] <- 17
V(fbpNetwork)$size[! as_ids(V(fbpNetwork)) %in% fbpNetHighlightNodes] <- 7

# Set the node names (displayed only for highlighted nodes)
V(fbpNetwork)$label[! as_ids(V(fbpNetwork)) %in% fbpNetHighlightNodes] <- NA
V(fbpNetwork)$label[as_ids(V(fbpNetwork)) %in% fbpNetHighlightNodes] <- fbpNetHighlightNodes

# Plot the subnetwork on a PDF
numUnknownEdges <- sum(is.na(fbpEdges$Prior))
numKnownEdges <- sum(!is.na(fbpEdges$Prior))
pdf("Cluster4Network.pdf", height=6, width=6)
plotTitle <- paste("Neighbors of gene \"fbp\"", "\n(", length(fbpNeighbors), " genes; ", numKnownEdges, " edges previously known, ", numUnknownEdges, " edges newly identified)", sep="")
plot(fbpNetwork, layout=layout_with_kk, vertex.label.family="Helvetica",
     vertex.label.cex=0.5, edge.width=1.2, vertex.label.font=2)
title(plotTitle, cex.main=0.8, font.main=1)
dev.off()

# ======= Figure 16 (Appendix): Cluster 6 expression trajectories =======

mannosidases <- c("LManI", "LManIII", "LManIV", "LManV", "LManVI")
maltases <- c("Mal-A2", "Mal-A3", "Mal-A4")
hemocytes <- c("NimB4", "NimC1", "NimC2", "eater", "Hml", "Gs1")
C6colors <- c(rep("firebrick3", 5), rep("goldenrod3", 3), rep("blue3", 4), "purple", "seagreen3")
p <- Plot.Genes(geneData[c(mannosidases, maltases, hemocytes), 1:11], hours[1:11], plotColors=C6colors,
                plotLegend=FALSE, plotTitle="<b>Selected genes from cluster 6<b>",
                plotSubtitle="(genes involved in carbohydrate metabolism: <span style='color:firebrick3;'>mannosidases</span>; <span style='color:goldenrod3';>maltases</span>;<br>genes involved in phagocytosis: <span style='color:blue3;'>Nimrod family</span>; <span style='color:purple;'><i>Hml</i></span>; <span style='color:seagreen3;'><i>Gs1</i></span>)",
                axisLabels=list(x="Time (hours after infection)", y="Expression (log<sub>2</sub>-fold change)"))
ggsave(file="Cluster6.pdf", p, width=6, height=4.5, units="in")
