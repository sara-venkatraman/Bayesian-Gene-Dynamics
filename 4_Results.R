# --- Script 4: Results ---

# This script produces all figures in our paper.

# Run the data loader script to read gene expression data and prior adjacency matrix
source("1_DatasetLoader.R")

# Run the LLR2 calculations script to produce the three lead-lag R^2
# similarity matrices *only* if they have not already been generated.
# source("2_LLR2Calculations.R")

# Run the plotting functions script
source("3_PlottingFunctions.R")

# Load library for arranging multiple plots on a grid
library(gridExtra)

# Read in the non-Bayesian lead-lag R^2 matrices (these matrices are only
# used to generate Figure 9)
nonBayesLLR2Mat <- read.csv("../Processed Data/R-Squared (Combined Genes)/NonBayesLLR2.csv", row.names=1)
nonBayesLLR2Mat.other <- read.csv("../Processed Data/R-Squared (Combined Genes)/NonBayesLLR2_Other.csv", row.names=1)
nonBayesLLR2Mat.own <- read.csv("../Processed Data/R-Squared (Combined Genes)/NonBayesLLR2_Own.csv", row.names=1)
colnames(nonBayesLLR2Mat) <- rownames(nonBayesLLR2Mat)
colnames(nonBayesLLR2Mat.other) <- rownames(nonBayesLLR2Mat.other)
colnames(nonBayesLLR2Mat.own) <- rownames(nonBayesLLR2Mat.own)

# Read in the Bayesian lead-lag R^2 matrices
bayesLLR2Mat <- read.csv("../Processed Data/R-Squared (Combined Genes)/BayesLLR2.csv", row.names=1)
bayesLLR2Mat.other <- read.csv("../Processed Data/R-Squared (Combined Genes)/BayesLLR2_Other.csv", row.names=1)
bayesLLR2Mat.own <- read.csv("../Processed Data/R-Squared (Combined Genes)/BayesLLR2_Own.csv", row.names=1)
colnames(bayesLLR2Mat) <- rownames(bayesLLR2Mat)
colnames(bayesLLR2Mat.other) <- rownames(bayesLLR2Mat.other)
colnames(bayesLLR2Mat.own) <- rownames(bayesLLR2Mat.own)

# --- Figure 1: known and unknown circadian rhythm, immune response genes ---

p1 <- Plot.Gene.Group(c("per", "tim", "to", "vri", "CG11854", "CG18609", "Pdp1", "CG33511"), plotTitle="Known and uncharacterized genes<br> with circadian rhythm patterns")
p2 <- Plot.Gene.Group(c("AttC", "DptA", "DptB", "Dro", "edin", "Mtk", "CG43920", "CG44404", "CG45045"), plotTitle="Known and uncharacterized genes<br> with immune response functions")
ggsave(file="CircadianAndImmunePatterns.pdf", arrangeGrob(grobs=list(p1,p2), ncol=2), width=10, height=4.5, units="in")

# --- Figure 2: inflated lead-lag R^2 ---

p1 <- Plot.Gene.Group(c("Act87E", "bmm"), plotTitle="Gene A: <i>Act87E</i>, Gene B: <i>bmm</i>", plotColors=c("dodgerblue2","orangered2"))
p2 <- Plot.Gene.Group(c("Act79B", "Mal-A7"), plotTitle="Gene A: <i>Act79B</i>, Gene B: <i>Mal-A7</i>", plotColors=c("dodgerblue2","orangered2"))
p3 <- Plot.Gene.Group(c("BomS3", "tok"), plotTitle="Gene A: <i>BomS3</i>, Gene B: <i>tok</i>", plotColors=c("dodgerblue2","orangered2"))
ggsave(file="InflatedLLR2.pdf", arrangeGrob(grobs=list(p1,p2,p3), ncol=3), width=12.5, height=4, units="in")

# Print non-Bayes lead-lag R^2 values
nonBayesLLR2Mat["bmm", "Act87E"]
nonBayesLLR2Mat["Mal-A7", "Act79B"]
nonBayesLLR2Mat["tok", "IM3"]

# --- Figure 3 and Tables 1-3: immunity/metabolism case study ---

immuneMetabolicGenes <- c("IM1", "IM2", "FASN1", "UGP", "mino", "fbp")
p <- Plot.Gene.Group(immuneMetabolicGenes, plotTitle="Selected genes involved in immune response<br> (<i>IM1, IM2</i>) and metabolism (<i>FASN1, UGP, mino, fbp</i>)", titleSize=12) + guides(color=guide_legend(nrow=1))
ggsave(file="ImmuneMetabolic.pdf", p, width=5.8, height=4.48, units="in")

# Print corresponding section of prior adjacency matrix
priorMatrix[immuneMetabolicGenes, immuneMetabolicGenes]

# Print corresponding sections of Bayesian LLR2 and LLR2 - LLR2_own matrices
# TODO: include code for asymmetric calculations
round(bayesLLR2Mat[immuneMetabolicGenes, immuneMetabolicGenes], 2)
round(bayesLLR2Mat[immuneMetabolicGenes, immuneMetabolicGenes] - bayesLLR2Mat.own[immuneMetabolicGenes, immuneMetabolicGenes], 2)

# --- Figure 4: Hierarchical clustering ---

# Compute distance matrix and cluster via Ward's method
hierClust <- hclust(as.dist(1 - bayesLLR2Mat), method="ward.D")

# Cut the dendrogram at height (yields 15 clusters; use "plot(hierClust)" for dendogram
subGroups <- cutree(hierClust, h=9)
table(subGroups)

# For each of the 15 clusters, plot the time profiles of all genes in the cluster
plotList <- list()
plotColors <- c("darkorange3", "dodgerblue3", "forestgreen", "darkmagenta", "indianred2", "orange4", "navy", "red2", "blueviolet", "turquoise4", "olivedrab", "darkslategray", "antiquewhite4", "coral4", "goldenrod2")
for(i in 1:length(table(subGroups))) {
  numGenesInCluster <- length(geneNames[subGroups == i])
  plotList[[i]] <- Plot.Gene.Group(geneNames[subGroups == i], plotColors=rep(plotColors[i], numGenesInCluster), points=FALSE, plotTitle=paste("Cluster ", i, " (", numGenesInCluster, " genes)", sep=""), titleSize=12, plotLegend=FALSE, lineOpacity=0.2)
}
ggsave(file="Clusters.pdf", arrangeGrob(grobs=plotList, ncol=4), width=12, height=9, units="in")

# --- Figure 5: Cluster 7 time profiles ---

imdGenes <- c("AttA", "AttB", "AttC", "Dro", "CecA2", "DptA", "DptB", "PGRP-SB1")
tollGenes <- c("PGRP-SA", "IMPPP", "IM1", "IM2", "IM4", "IM14", "IM23", "BomS3")
newGenes <- c("CG44404", "CG43236", "CG43202", "CG43920")
negativeGenes <- c("Acp1", "CG7214")
C7colors <- c(rep(alpha("orangered2", 0.65), 8), rep(alpha("dodgerblue3", 0.65), 8), rep(alpha("black", 0.75), 4), rep(alpha("forestgreen", 0.65), 2))
p <- Plot.Gene.Group(c(imdGenes, tollGenes, newGenes, negativeGenes), plotColors=C7colors, plotTimes=hours[1:8],
                plotLegend=FALSE, plotTitle="<b>Selected genes from cluster 7</b>",
                plotSubtitle="(<span style='color:orangered2;'>Imd-regulated</span><span style='color:white;'> l</span><span style='color:orangered2;'>genes</span>; <span style='color:dodgerblue3;'>Toll-regulated genes</span>; <span style='color:forestgreen;'>genes that encode cuticle proteins</span>;<br><span style='color:black;'>genes potentially associated with Imd signaling</span>)")
ggsave(file="Cluster7Genes.pdf", p, width=6, height=4.5, units="in")

# --- Figure 6: Cluster 7 subnetwork ---

newGenes <- c("CG44404", "CG43236", "CG43202", "CG43920")

# Get the entire network of genes (edge defined for LLR2 > 0.9)
adjBayes <- (bayesLLR2Mat > 0.9) + 0
networkBayes <- graph_from_adjacency_matrix(adjBayes , mode='undirected', diag=FALSE)

# Get all neighbors of the unknown genes in cluster 7 as a named list
unknownC7neighbors <- adjacent_vertices(networkBayes, newGenes)

# Collapse this list of neighbors into one vector of nodes
allNodes <- c()
allNodes <- lapply(unknownC7neighbors, function(x) c(allNodes, names(x)))
allNodes <- unique(unlist(allNodes, use.names=FALSE))

# Form a new subnetwork out of allNodes
C7adj <- (bayesLLR2Mat[allNodes, allNodes] > 0.9) + 0
C7net <- graph_from_adjacency_matrix(C7adj , mode='undirected', diag=FALSE)

# Get the priors associated with each edge in this subnetwork
C7edges <- data.frame(as_edgelist(C7net))
colnames(C7edges) <- c("Gene1", "Gene2")
C7edges$Prior <- unlist(apply(C7edges, 1, function(x) priorMatrix[x[1], x[2]]))

# Blue edges between genes with unknown associations,
# red edges between genes with known associations
E(C7net)$color[is.na(C7edges$Prior)] <- alpha('blue', 0.7)
E(C7net)$color[!is.na(C7edges$Prior)] <- alpha('red', 0.7)

# The four genes of interest will be colored differently
V(C7net)$color[as_ids(V(C7net)) %in% newGenes] <- alpha("navajowhite", 0.95)
V(C7net)$color[! as_ids(V(C7net)) %in% newGenes] <- alpha("linen", 0.85)

# Set the node outline colors
V(C7net)$frame.color[as_ids(V(C7net)) %in% newGenes] <- "peachpuff3"
V(C7net)$frame.color[! as_ids(V(C7net)) %in% newGenes] <- "gray66"

# Plot the subnetwork on a PDF (note: network layout is random)
numUnknownEdges <- sum(is.na(C7edges$Prior))
numKnownEdges <- sum(!is.na(C7edges$Prior))
pdf("Cluster7Subnetwork.pdf", height=6, width=6)
plotTitle <- paste("New relationships detected in cluster 7\n(", length(allNodes), " genes; ", numKnownEdges, " previously-known edges, ", numUnknownEdges, " newly-identified edges)", sep="")
plot(C7net, layout=layout_with_lgl, vertex.size=21, vertex.label.family="Helvetica",
     vertex.label.cex=0.5, edge.width=1.2)
title(plotTitle, cex.main=0.8, font.main=1)
dev.off()

# --- Figure 7: Cluster 12 time profiles ---

C12genes <- c("fbp", "to", "AGBE", "Galk", "Gba1b", "CG11594", "CG10469", "CG13315")
C12colors <- c("black", "dodgerblue2", rep("orangered2", 3), rep("springgreen4", 3))
p <- Plot.Gene.Group(C12genes, plotColors=C12colors, plotLegend=FALSE,
                plotTitle="<b>Selected genes from cluster 12</b>",
                plotSubtitle="(<span style='color:black;'>gene <i>fbp</i></span>; <span style='color:orangered2;'>genes involved in carbohydrate metabolism</span>;<br> <span style='color:springgreen4;'>uncharacterized genes</span>; <span style='color:dodgerblue2;'>gene <i>takeout</i></span>)")
ggsave(file="Cluster12Genes.pdf", p, width=6, height=4.5, units="in")

# --- Figure 8: Cluster 12 subnetwork - neighbors of gene 'fbp' ---

# Get the entire network of genes (edge defined for LLR2 > 0.9)
adjBayes <- (bayesLLR2Mat > 0.9) + 0
networkBayes <- graph_from_adjacency_matrix(adjBayes , mode='undirected', diag=FALSE)

# Get all neighbors of fbp
fbpNeighbors <- adjacent_vertices(networkBayes, "fbp")
fbpNeighbors <- c("fbp", names(fbpNeighbors[[1]]))

# Form a new subnetwork out of fbp and its neighbors
fbpAdj <- (bayesLLR2Mat[fbpNeighbors, fbpNeighbors] > 0.9) + 0
fbpNet <- graph_from_adjacency_matrix(fbpAdj , mode='undirected', diag=F)

# Get the priors associated with each edge in this subnetwork
fbpEdges <- data.frame(as_edgelist(fbpNet))
colnames(fbpEdges) <- c("Gene1", "Gene2"); fbpEdges$Prior <- 0
fbpEdges$Prior <- unlist(apply(fbpEdges, 1, function(x) priorMatrix[x[1], x[2]]))

# Blue edges between genes with unknown associations,
# red edges between genes with known associations
E(fbpNet)$color[is.na(fbpEdges$Prior)] <- alpha('blue', 0.1)
E(fbpNet)$color[!is.na(fbpEdges$Prior)] <- alpha('red', 0.1)

# The fbp gene will be colored differently
V(fbpNet)$color[as_ids(V(fbpNet)) == "fbp"] <- "navajowhite"
V(fbpNet)$color[! as_ids(V(fbpNet)) == "fbp"] <- alpha("linen", 0.85)

# Set the node outline colors
V(fbpNet)$frame.color[as_ids(V(fbpNet)) == "fbp"] <- "peachpuff3"
V(fbpNet)$frame.color[! as_ids(V(fbpNet)) == "fbp"] <- "gray66"

# Set the node sizes (a few highlighted nodes will be larger)
fbpNetHighlightNodes <- c("fbp", "Galk", "AGBE", "Gba1b", "CG11594", "CG10469", "CG13315") 
V(fbpNet)$size[as_ids(V(fbpNet)) %in% fbpNetHighlightNodes] <- 17
V(fbpNet)$size[! as_ids(V(fbpNet)) %in% fbpNetHighlightNodes] <- 7

# Set the node names (displayed only for highlighted nodes)
V(fbpNet)$label[! as_ids(V(fbpNet)) %in% fbpNetHighlightNodes] <- NA
V(fbpNet)$label[as_ids(V(fbpNet)) %in% fbpNetHighlightNodes] <- fbpNetHighlightNodes

# Darken edges for the highlighted nodes
E(fbpNet)$color[which(fbpEdges$Gene1 == "fbp" & fbpEdges$Gene2 %in% fbpNetHighlightNodes)] <- "blue"

# Plot the subnetwork on a PDF
numUnknownEdges <- sum(is.na(fbpEdges$Prior))
numKnownEdges <- sum(!is.na(fbpEdges$Prior))
pdf("fbpSubnetwork.pdf", height=6, width=6)
plotTitle <- paste("Neighbors of gene \"fbp\" in cluster 12\n(", length(fbpNeighbors), " genes; ", numKnownEdges, " known edges, ", numUnknownEdges, " newly-identified edges)", sep="")
plot(fbpNet, layout=layout_with_kk, vertex.label.family="Helvetica",
     vertex.label.cex=0.5, edge.width=1.2)
title(plotTitle, cex.main=0.8, font.main=1)
dev.off()

# --- Figure 9: R^2 scatterplots ---

geneSubset <- sample(geneNames, size=150)
p1 <- Draw.R2.Scatterplot(nonBayesLLR2Mat, nonBayesLLR2Mat.other, nonBayesLLR2Mat.own, priorMatrix, geneSubset, interactive=F, plotTitle=latex2exp::TeX("Non-Bayesian lead-lag $R^2$ values"))
p2 <- Draw.R2.Scatterplot(bayesLLR2Mat, bayesLLR2Mat.other, bayesLLR2Mat.own, priorMatrix, geneSubset, interactive=F)
ggsave(file="R2Scatterplots.pdf", arrangeGrob(grobs=list(p1,p2), ncol=2), width=10, height=4, units="in")

# --- Figure 10: Genes in middle/upper-right region of Bayesian LLR2 scatterplot ---

plotColors <- c("dodgerblue2", "orangered2", "goldenrod", "forestGreen", "magenta", "navy", "chocolate1")
p1 <- Plot.Gene.Group(c("CR42868", "AttD", "CG9616"), plotTitle="Genes: <i>CR42868, AttD, CG9616</i>", plotColors=plotColors)
p2 <- Plot.Gene.Group(c("Spn28Dc", "CR43364", "CR42715", "scb"), plotTitle="Genes: <i>Spn28Dc, CR43364, CR42715, scb</i>", plotColors=plotColors)
p3 <- Plot.Gene.Group(c("ACC", "Idh", "GstE9"), plotTitle="Genes: <i>ACC, Idh, GstE9</i>", plotColors=plotColors)
ggsave(file="R2ScatterplotsMiddle.pdf", arrangeGrob(grobs=list(p1,p2,p3), ncol=3), width=13, height=3.9, units="in")

# --- Figure 11: Genes in upper-right region of Bayesian LLR2 scatterplot ---

plotColors <- c("dodgerblue2", "orangered2", "goldenrod", "forestGreen", "magenta", "navy", "chocolate1")
p1 <- Plot.Gene.Group(c("alphaTry", "gammaTry", "CG30025"), plotTitle="Genes: <i>alphaTry, gammaTry, CG30025</i>", plotColors=plotColors)
p2 <- Plot.Gene.Group(c("RpS26", "RpS6", "RpL13", "RpL7"), plotTitle="Genes: <i>RpS26, RpS6, RpL13, RpL7</i>", plotColors=plotColors)
p3 <- Plot.Gene.Group(c("CG13096", "l(1)G0020", "Nop56"), plotTitle="Genes: <i>CG13096, l(1)G0020, Nop56</i>", plotColors=plotColors)
ggsave(file="R2ScatterplotsRight.pdf", arrangeGrob(grobs=list(p1,p2,p3), ncol=3), width=13, height=3.9, units="in")

# --- Figure 12: Cluster 2 time profiles ---

circadian <- c("tim", "per", "Clk", "vri", "Pdp1")
cuticle <- c("Cpr49Ab", "Cpr49Ae", "Cpr62Ba", "Cpr72Ec")
dopSynth <- c("e", "ple")
C2colors <- c(rep("black", 3), "orange", "red", rep("pink2", 4), rep("purple", 2))
p <- Plot.Gene.Group(c(circadian, cuticle, dopSynth), plotColors=C2colors,
                plotTitle="<b>Selected genes from cluster 2<b>", plotLegend=FALSE,
                plotSubtitle="(Regulators of circadian clock: <i>tim, per, Clk, <span style='color:darkorange;'>vri</span>, <span style='color:red;'>Pdp1</i></span>;<br> <span style='color:lightpink3;'>genes that encode cuticle proteins</span>; <span style='color:purple;'>genes involved in dopamine synthesis</span>)")
ggsave(file="Cluster2Genes.pdf", p, width=6, height=4.5, units="in")

# --- Figure 13: Cluster 5 time profiles ---

ribosome <- c("nop5", "Fib", "Nop60B", "CG12301", "CG32409", "U3-55K")
fattyAcid <- c("FASN1", "ACC", "AcCoAS", "ATPCL", "mino")
fasn1New <- c("CG3756", "CG3940", "CG8036", "eRF1", "aralar1", "CG32409", "CG31904", "CG15120", "CG1640", "CG16926") # CG32904 is CG32409?
lipidCat <- c("dob", "Lip2", "Lsd-1")
C5colors <- c(rep(alpha("orangered3", 0.7), 6), rep(alpha("blue2", 0.7), 5), rep(alpha("chartreuse3", 0.75), 10), rep(alpha("black", 0.7), 3))
p <- Plot.Gene.Group(c(ribosome, fattyAcid, fasn1New, lipidCat), plotColors=C5colors, plotLegend=FALSE,
                plotTitle="<b>Selected genes from cluster 5</b>", plotTimes=hours[1:12],
                plotSubtitle="(<span style='color:orangered3;'>ribosome biogenesis</span>; <span style='color:black;'>lipid catabolism</span>; <span style='color:blue2;'>fatty acid biosynthesis</span>;<br><span style='color:chartreuse3;'>genes with uncharacterized relationships to <i>FASN1</i></span>)")
ggsave(file="Cluster5Genes.pdf", p, width=6, height=4.5, units="in")

# --- Figure 14: Cluster 9 time profiles ---

maltaseUp <- c("Mal-A1", "Mal-A6", "Mal-A7", "Mal-A8")
hemoDown <- c("NimC1", "NimB4", "eater", "Hml")
other <- c("Galphaf", "tobi")
C9colors <- c(rep(alpha("orangered3", 0.7), 4), rep(alpha("black", 0.8), 4), "green3", "dodgerblue2")
p <- Plot.Gene.Group(c(maltaseUp, hemoDown, other), plotColors=C9colors, plotLegend=FALSE,
                plotTitle="<b>Selected genes from cluster 9<b>", plotTimes=hours[1:11],
                plotSubtitle="(<span style='color:orangered3;'>up-regulated</span><span style='color:white;'> l</span><span style='color:orangered3;'>maltases</span>; down-regulated<span style='color:white;'> l</span>genes expressed in hemocytes;<br> <span style='color:dodgerblue2;'>gene <i>Galphaf</i></span>; <span style='color:green3;'>gene <i>tobi</i></span>)")
ggsave(file="Cluster9Genes.pdf", p, width=6, height=4.5, units="in")

