# --- Script 3: Plotting time-course gene expression data ---

library(ggplot2)          # General plotting
library(plotly)           # Interactive plotting
library(igraph)           # Network plotting
library(RColorBrewer)     # Color palettes
library(scales)           # Color opacity
library(ggtext)           # Markdown parser
library(latex2exp)        # LaTeX parser
library(directlabels)     # Labels on line graphs
library(reshape2)         # Dataset reshaping

# Input: geneName, the name of a gene as a string
# Output: an natural cubic spline interpolant of the gene's expression data;
# the interpolant is a function that can be evaluated at any real number.
Expression.Profile.Interpolant <- function(geneName) {
  profile <- geneData[geneName,]
  numPoints <- length(profile)
  splinefun(x=hours, y=profile, method="natural")
}  

# Plots the time profiles of specified genes as smooth lines.
# Input:
# - genesToPlot: vector of gene names as strings
# - plotColors: optional vector of colors for each line
# - plotTitle: string that will be printed at the top of the plot. Default none.
# - plotSubtitle: string that will be printed underneath the title. Default none.
# - points: logical indicating whether the observed data points should be drawn as points. Default TRUE.
# - plotLegend: logical indicating whether a legend should be displayed. Default TRUE.
# - plotGrid: logical indicating whether grid lines should be drawn. Default TRUE.
# - lineLabels: logical indicating whether gene names should be printed at the end of the corresponding lines. Default FALSE.
# - pointSize: numeric size of observed data points, if points is TRUE. Default 1.5.
# - lineOpacity: opacity of plotted lines. Default 1 (opaque); decrease for transparency.
# - legendSize: font size of legend text. Default 10.
# - legendPos: position of legend ("top", "bottom", "left", "right"). Default "bottom". 
# - titleSize: font size of title text. Default 12.
# - plotTimes: length of time over which gene expression should be plotted. Default hours (full time period).
# Output: ggplot object.
Plot.Gene.Group <- function(genesToPlot, plotColors, plotTitle="", plotSubtitle="", points=TRUE, plotLegend=TRUE,
                            plotGrid=TRUE, lineLabels=FALSE, pointSize=1.5, lineOpacity=1, legendSize=10, legendPos="bottom",
                            titleSize=12, plotTimes=hours) {

  # Set the colors for each line if none are provided
  if(missing(plotColors)) {
    plotColors <- c(brewer.pal(n=8, name="Dark2"), brewer.pal(n=8, name="Set2"), brewer.pal(n=12, name="Paired"))
    if(length(genesToPlot) <= length(plotColors))
      plotColors <- plotColors[1:length(genesToPlot)]
    else
      plotColors <- rep(plotColors, 20)[1:length(genesToPlot)]
  }
  
  # Reduce opacity of each line if needed
  if(lineOpacity != 1)
    plotColors <- alpha(plotColors, lineOpacity)
  
 # Create a dataframe of the expression data for each gene
  exprData <- data.frame(plotTimes, t(geneData[genesToPlot, 1:length(plotTimes)]))
  colnames(exprData) <- c("time", genesToPlot)
  
  # Get the timepoints at each gene's spline interpolant will be plotted
  interpTimes <- seq(from=plotTimes[1], to=plotTimes[length(plotTimes)], length.out=500)
  
  # Evaluate each gene's spline interpolant at each time point (this is the
  # data that will be plotted)
  interpData <- data.frame(interpTimes, matrix(0, nrow=length(interpTimes), ncol=length(genesToPlot)))  
  for(i in 1:length(genesToPlot))
    interpData[,i+1] <- Expression.Profile.Interpolant(genesToPlot[i])(interpTimes)
  colnames(interpData) <- c("time", genesToPlot)
  
  # Initial plot: smooth lines for each gene
  p <- ggplot(melt(interpData, id.var="time"), aes(x=time, y=value, col=variable)) +
    geom_line() + scale_color_manual(values=plotColors) + theme_bw() + 
    labs(x="Time (hours)", y="Expression") 
  
  # If desired, draw points of specified size at each observed timepoint
  if(points == TRUE)
    p <- p + geom_point(data=melt(exprData, id.var="time"), mapping=aes(x=time, y=value, col=variable), size=pointSize) 
  
  # If plot title is supplied, center it and parse as markdown. If plot title is 
  # desired but not supplied, set it to be a comma-separated list of the gene names.
  if(plotTitle != "") {
    if(plotTitle == TRUE)
      plotTitle <- paste(genesToPlot, collapse=", ")
    p <- p + labs(title=plotTitle) + theme(plot.title=element_text(size=titleSize, hjust=0.5)) + theme(plot.title=element_markdown(lineheight=1.1))
  }
  
  # If a plot subtitle is supplied, center it and parse as markdown
  if(plotSubtitle != "")
    p <- p + labs(subtitle=plotSubtitle) + theme(plot.subtitle=element_text(hjust=0.5)) + theme(plot.subtitle=element_markdown(lineheight=1.1))
  
  # If not desired, remove grid lines
  if(plotGrid == FALSE)
    p <- p + theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank())
  
  # If desired, add gene names as labels next to each line
  if(lineLabels == TRUE)
    p <- p + geom_dl(aes(label=variable), method=list(dl.trans(x=x+0.2), "last.points", cex=0.6)) + expand_limits(x=tail(plotTimes,1)+5)
  
  # If desired, add a legend of the specified size at the specified position
  if(plotLegend == TRUE)
    p <- p + theme(legend.title=element_blank(), legend.position=legendPos, legend.text=element_text(size=legendSize), legend.background=element_rect(size=0.1, linetype="solid", color="black"))
  else
    p <- p + theme(legend.position="none")
  
  p
}

# Turns a similarity matrix into pairwise similarity list. This is a 
# helper function for Draw.R2.Scatterplot.
# Input: a similarity matrix with labeled rows and columns
# Output: a dataframe with one row per pair of items. Each row is of the
# form [item1Name   item2Name   similarity].
Vectorize.Similarity.Matrix <- function(similarityMatrix) {
  namePairs <- t(combn(colnames(similarityMatrix), 2))
  data.frame(namePairs, value=similarityMatrix[namePairs])
}

# Draws a scatterplot of lead-lag R^2 values for all gene pairs within a 
# specified set of genes (shown in Figure 9 in our paper). X axis displays
# the LLR2_other value and Y axis displays LLR2 - LLR2_own; see section 2.2.
# Input:
# - LLR2, LLR2other, and LLR2own: similarity matrices produced by the
#   script "2_LLR2Calculatioons.R".
# - priorMatrix: matrix of prior information, loaded by "1_DatasetLoader.R"
# - geneSubset: vector of gene names as strings
# - interactive: logical indicating whether or not the scatterplot should be
#   interactive; if TRUE (default), hovering over the points shows the gene names.
# - plotTitle: optional string that will be printed at the top of the plot.
# Output: a ggplot object if interactive is FALSE, or a plotly object if TRUE.
Draw.R2.Scatterplot <- function(LLR2, LLR2other, LLR2own, priorMatrix, geneSubset, interactive=TRUE, plotTitle="") {
  # Turn R^2 similarity matrices into pairwise distance lists for both axes
  xAxis <- Vectorize.Similarity.Matrix(LLR2other[geneSubset, geneSubset])
  yAxis <- Vectorize.Similarity.Matrix(LLR2[geneSubset, geneSubset] - LLR2own[geneSubset, geneSubset])
  
  # Turn priorMatrix into pairwise prior list, treated as factor
  priorVec <- Vectorize.Similarity.Matrix(priorMatrix[geneSubset, geneSubset])
  priorVec$value[is.na(priorVec$value)] <- "NA"
  priorVec$value <- factor(priorVec$value, levels=c("0", "NA", "1"))
  
  # Arrange scatterplot data into a dataframe
  plotData <- data.frame(xAxis=round(xAxis$value, 3), yAxis=round(yAxis$value, 3), prior=priorVec$value)
  
  # Re-order points so that points with prior=0 are drawn first, and points
  # with prior=1 are drawn last.
  plotData <- plotData[order(plotData$prior),]
  
  # Create gene name labels for interactive scatterplot
  pointLabels <- paste(xAxis$X1, xAxis$X2, sep=", ")
  
  # Initial scatterplot
  p <- ggplot(plotData, aes(x=xAxis, y=yAxis, color=prior, text=pointLabels)) + 
    geom_point(size=0.9) + theme_light() + 
    scale_color_manual(name="Prior", labels=c("0", "NA", "1"), values=c(alpha("gray67", 0.9), alpha("navy", 0.65), alpha("orangered3",0.8))) + 
    theme(plot.title=element_text(hjust=0.5))

  # Further settings for interactive scatterplot
  if(interactive == TRUE) {
    p <- p + labs(title=ifelse(plotTitle == "", "Bayesian lead-lag R-squared values", plotTitle), x="LLR2_other", y="LLR2 - LLR2_own")
    ggplotly(p) %>% layout(legend=list(title=list(text="Prior"), orientation="h", xanchor="center", x=0.5, y=-0.2))
  }
  # Further settings for non-interactive scatterplot
  else {
    p + labs(title=ifelse(class(plotTitle) != "expression" && plotTitle == "", latex2exp::TeX("Bayesian lead-lag $R^2$ values"), plotTitle), x=latex2exp::TeX("LLR^2_{other}"), y=latex2exp::TeX("LLR^2 - LLR^2_{own}")) + 
      theme(legend.position="bottom", legend.background=element_rect(size=0.1, linetype="solid", color="black")) + guides(color=guide_legend(override.aes=list(size=3)))
  }
}


