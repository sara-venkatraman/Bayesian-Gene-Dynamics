Dynamic models of co-regulated gene expression
================

The code in this repository implements the Bayesian lead-lag
*R*<sup>2</sup> (LLR2) algorithm introduced in Venkatraman et al. (2021)
for biologically-informed clustering and network analysis of time-course
gene expression data. For details of how it works, see Section 3 of our
[paper](https://www.biorxiv.org/content/10.1101/2021.07.08.451684v1).
Below is a short demo.

This repository consists of the following R files:

-   `1_LLR2.R`: contains the function `LLR2`, which computes the
    Bayesian LLR2 similarity matrix for a given gene expression dataset
-   `2_PlottingFunctions.R`: contains functions for plotting temporal
    gene expression trajectories
-   `3_Results.R`: produces all figures and tables in our paper

### Getting started

We begin by loading a *N* × *n* time-course gene expression dataset,
consisting of *N* genes whose expressions are recorded at *n* time
points, and a *N* × *N* prior adjacency matrix. The latter contains
entries `0`, `1`, or `NA` to indicate that a biological relationship
between the corresponding two genes is unlikely, likely, or unknown
according to external sources.

``` r
# Read gene expression data and prior adjacency matrix
geneData <- read.csv("Data/GeneData.csv", row.names=1)
priorMatrix <- read.csv("Data/PriorMatrix.csv", row.names=1)

# Define hours corresponding to each time point
hours <- c(0, 1, 2, 4, 5, 6, 8, 10, 12, 14, 16, 20, 24, 30, 36, 42, 48)
```

This data comes from Schlamp et al. (2021), whose GitHub repository can
be found
[here](https://github.com/florschlamp/Drosophila_Immunity_TimeSeries),
and the external information comes from the
[STRING](https://string-db.org/) database (Szklarczyk et al. 2019).

Next we run the Bayesian LLR2 algorithm. For this demo, we’ll run it on
a random subset of 150 genes.

``` r
# Select random subset of genes
set.seed(12657)
indices <- sample(1:nrow(geneData), size=150)
geneData <- geneData[indices,]
priorMatrix <- priorMatrix[indices, indices]

# Run Bayesian LLR2 similarity matrix computations
source("1_LLR2.R")
bayesLLR2 <- LLR2(geneData, hours, bayes=TRUE, priorMatrix, writeToCSV=FALSE)
```

### Clustering with the Bayesian LLR2 metric

We now perform hierarchical clustering on the Bayesian LLR2 similarity
matrix.

``` r
# Turn Bayesian LLR2 similarity matrix into a distance matrix
distanceMatrix <- as.dist(1 - bayesLLR2$LLR2Mat)
hierClust <- hclust(distanceMatrix, method="ward.D")

# Divide into 12 clusters
geneClusters <- cutree(hierClust, k=12)
```

Next we plot the genes in each cluster:

``` r
# Load plotting functions (we use Plot.Genes() below)
source("2_PlottingFunctions.R")

# Initialize list of plots and choose plot colors
plotList <- list()
plotColors <- c(brewer.pal(8, "Dark2"), brewer.pal(4, "Set1"))

# Plot expression trajectories of genes in each cluster
for(i in 1:12) {
  numGenes <- table(geneClusters)[i]
  plotList[[i]] <- Plot.Genes(geneData[geneClusters == i,], hours, points=FALSE, plotColors=rep(plotColors[i], numGenes), plotTitle=paste("Cluster ", i, " (", numGenes, " genes)", sep=""), plotLegend=FALSE, lineOpacity=0.6)
}
grid.arrange(grobs=plotList, ncol=4)
```

![](Demo_files/figure-gfm/unnamed-chunk-4-1.png)<!-- -->

### References

<div id="refs" class="references csl-bib-body hanging-indent">

<div id="ref-schlamp2021dense" class="csl-entry">

Schlamp, Florencia, Sofie YN Delbare, Angela M Early, Martin T Wells,
Sumanta Basu, and Andrew G Clark. 2021. “Dense Time-Course Gene
Expression Profiling of the Drosophila Melanogaster Innate Immune
Response.” *BMC Genomics* 22 (1): 1–22.

</div>

<div id="ref-szklarczyk2019string" class="csl-entry">

Szklarczyk, Damian, Annika L Gable, David Lyon, Alexander Junge, Stefan
Wyder, Jaime Huerta-Cepas, Milan Simonovic, et al. 2019. “STRING V11:
Protein–Protein Association Networks with Increased Coverage, Supporting
Functional Discovery in Genome-Wide Experimental Datasets.” *Nucleic
Acids Research* 47 (D1): D607–13.

</div>

<div id="ref-venkatraman2021empirical" class="csl-entry">

Venkatraman, Sara, Sumanta Basu, Andrew G Clark, Sofie Delbare, Myung
Hee Lee, and Martin T Wells. 2021. “An Empirical Bayes Approach to
Estimating Dynamic Models of Co-Regulated Gene Expression.” *bioRxiv*.

</div>

</div>
