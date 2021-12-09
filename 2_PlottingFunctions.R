# ======= Script 2: Plotting time-course gene expression data =======

library(ggplot2)          # General plotting
library(plotly)           # Interactive plotting
library(igraph)           # Network plotting
library(RColorBrewer)     # Color palettes
library(scales)           # Color opacity
library(ggtext)           # Markdown parser
library(directlabels)     # Labels on line graphs
library(reshape2)         # Dataset reshaping
library(gridExtra)        # Arrange plots on grids

 # ========================================================================

#' @description 
#' Plot the temporal expression trajectories of one or more genes as 
#' smooth lines. 
#'
#' @param geneData N-by-n numeric dataframe containing the N temporal 
#' expression trajectories of length-n to be plotted. Row names should be
#' the names of the genes.
#' @param timePoints Length-n numeric vector of time points corresponding
#' to the times at which the columns of geneData were recorded.
#' @param plotColors Optional vector of colors to use for each line (gene),
#' of the same number of rows as \code{geneData}.
#' @param plotTitle Optional title to place at the top of the plot, 
#' centered, or TRUE if the desired title is a list of the gene names. 
#' Can include Markdown/HTML tags.
#' @param plotSubtitle An optional title to place under the plot title,
#' centered. Can include Markdown/HTML tags.
#' @param titleSize Numeric font size of the title text. Default 12.
#' @param points Logical indicating whether the observed data points should
#' be plotted atop their corresponding spline interpolants. Default TRUE.
#' @param pointSize Numeric size of observed data points, if \code{points}
#' is TRUE. Default 1.5.
#' @param plotLegend Logical indicating whether a legend should be displayed.
#' Default TRUE.
#' @param legendSize Numeric font size of the legend text. Default 10.
#' @param legendPos Position of legend. Possible values are "top", "bottom",
#' "left", "right". Default "bottom".
#' @param plotGrid Logical indicating whether background grid lines 
#' should be drawn. Default TRUE.
#' @param lineLabels Logical indicating whether or not gene names should
#' be printed at the end (right-hand side) of their plotted trajectories. 
#' Default FALSE.
#' @param axisLabels Length-2 list with items "x" and "y" specifying x- and
#' y-axis labels. Can include Markdown/HTML tags. Default 
#' \code{list(x="Time (hours)", y="Expression")}.
#' @param lineOpacity Numeric opacity of plotted lines. Default 1 (opaque);
#' decrease for transparency.
#'
#' @return A ggplot2 object.
#' @export
Plot.Genes <- function(geneData, timePoints, plotColors, plotTitle="",
                       plotSubtitle="", titleSize=12, points=TRUE, pointSize=1.5,
                       plotLegend=TRUE, legendSize=10, legendPos="bottom",
                       plotGrid=TRUE, lineLabels=FALSE, lineOpacity=1,
                       axisLabels=list(x="Time", y="Expression")) {
  
  # Set the colors for each line if none are provided
  if(missing(plotColors)) {
    plotColors <- c(brewer.pal(n=8, name="Dark2"), brewer.pal(n=8, name="Set2"), brewer.pal(n=12, name="Paired"))
    plotColors <- rep(plotColors, length.out=nrow(geneData))
  }
  
  # Adjust opacity of colors
  plotColors <- alpha(plotColors, lineOpacity)
  
  # Define many time points at which to evaluate each gene's spline interpolant
  interpTimes <- seq(from=timePoints[1], to=tail(timePoints, 1), length.out=500)
  
  # Evaluate each gene's spline interpolant at each time point. This is the
  # data that will be plotted.
  interpData <- data.frame(interpTimes, matrix(0, nrow=500, ncol=nrow(geneData)))
  interpData[,-1] <- apply(geneData, 1, function(y) splinefun(x=timePoints, y=y, method="natural")(interpTimes))
  colnames(interpData) <- c("time", rownames(geneData))
  
  # Reshape interpolation data into long format for use with ggplot
  meltInterpData <- melt(interpData, id.var="time")
  
  # Construct initial plot: smooth lines for each gene
  p <- ggplot(meltInterpData, aes(x=time, y=value, col=variable)) +
    geom_line() + scale_color_manual(values=plotColors) + theme_bw() + 
    labs(x=axisLabels$x, y=axisLabels$y) + # "Expression (log<sub>2</sub>-fold change)", "")) +
    theme(axis.title.x=element_markdown(), axis.title.y=element_markdown())
  
  # If desired, draw points of specified size at each observed time
  if(points == TRUE) {
    pointData <- data.frame(timePoints, t(geneData))
    colnames(pointData) <- c("time", rownames(geneData))
    meltPointData <- melt(pointData, id.var="time")
    p <- p + geom_point(data=meltPointData, mapping=aes(x=time, y=value, col=variable), size=pointSize) 
  }
  
  # If plot title is desired but not supplied, set it to be a comma-separated
  # list of the gene names
  if(plotTitle == TRUE)
    plotTitle <- paste(rownames(geneData), collapse=", ")
  
  # Add the plot title, centered and with line height 1.1, with parsing
  p <- p + labs(title=plotTitle) + 
    theme(plot.title=element_text(size=titleSize, hjust=0.5)) + 
    theme(plot.title=element_markdown(lineheight=1.1))
  
  # Add the plot subtitle, centered and with line height 1.1, with parsing
  p <- p + labs(subtitle=plotSubtitle) +
    theme(plot.subtitle=element_text(hjust=0.5)) + 
    theme(plot.subtitle=element_markdown(lineheight=1.1))
  
  # Remove grid lines if not desired
  if(plotGrid == FALSE)
    p <- p + theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank())
  
  # If desired, add gene names as labels next to each line
  if(lineLabels == TRUE)
    p <- p + geom_dl(aes(label=variable), method=list(dl.trans(x=x+0.2), "last.points", cex=0.6)) +
    expand_limits(x=tail(timePoints, 1) + 5)
  
  # If desired, add a legend with specified font size and position
  if(plotLegend == TRUE)
    p <- p + theme(legend.title=element_blank(), legend.position=legendPos, 
                   legend.text=element_text(size=legendSize), 
                   legend.background=element_rect(size=0.1, linetype="solid", color="black"))
  else
    p <- p + theme(legend.position="none")
  
  # Return the ggplot object
  p
}

# ========================================================================

#' @description
#' Draw a scatterplot of lead-lag R^2 (LLR2) values, with LLR2_other on the 
#' horizontal axis and LLR2 - LLR2_own on the vertical axis. Examples can
#' be found in Figure 8 of Venkatraman et al. 2021.
#'
#' @param LLR2Mat Length-3 list of N-by-N matrices named \code{LLR2Mat},
#' \code{LLR2Mat.own}, and \code{LLR2Mat.other}. This list can be obtained
#' by running the function \code{LLR2}.
#' @param priorMatrix N-by-N prior adjacency matrix.
#' @param geneSubset Length-k list of gene names, where k \eqn{leq} N, 
#' whose pairwise LLR2 values are to be included in the scatterplot.
#' @param interactive Logical indicating whether the plot should be made
#' interactive, i.e. with labels displayed upon hovering.
#' @param plotTitle Optional title to place at the top of the plot, 
#' centered. Can include Markdown/HTML tags.
#'
#' @return A ggplot2 object.
#' @export
Plot.LLR2.Scatterplot <- function(LLR2Matrices, priorMatrix, geneSubset,
                                  interactive=TRUE,
                                  plotTitle="Lead-lag R<sup>2</sup> values") {
  # Define a function for vectorizing a symmetric matrix
  Vectorize.Symm.Mat <- function(symmMat) {
    namePairs <- t(combn(colnames(symmMat), 2))
    data.frame(namePairs, value=symmMat[namePairs])
  }
  
  # Turn LLR2 matrices into pairwise distance lists for both axes
  xAxis <- Vectorize.Symm.Mat(LLR2Matrices$LLR2Mat.other[geneSubset, geneSubset])
  yAxis <- Vectorize.Symm.Mat(LLR2Matrices$LLR2Mat[geneSubset, geneSubset] - 
                                         LLR2Matrices$LLR2Mat.own[geneSubset, geneSubset])
  
  # Turn priorMatrix into pairwise prior list, treated as factor
  priorVec <- Vectorize.Symm.Mat(priorMatrix[geneSubset, geneSubset])
  priorVec$value[is.na(priorVec$value)] <- "NA"
  priorVec$value <- factor(priorVec$value, levels=c("0", "NA", "1"))
  
  # Arrange scatterplot data into a dataframe
  plotData <- data.frame(xAxis=round(xAxis$value, 3), 
                         yAxis=round(yAxis$value, 3), 
                         prior=priorVec$value)
  
  # Re-order data so that points with prior=0 are drawn first, and points
  # with prior=1 are drawn last.
  plotData <- plotData[order(plotData$prior),]
  
  # Create gene name labels for interactive scatterplot
  pointLabels <- paste(xAxis$X1, xAxis$X2, sep=", ")
  
  # Define gray, blue, and red colors for the scatterplot
  g <- alpha("gray67", 0.9)
  b <- alpha("navy", 0.65)
  r <- alpha("orangered3",0.8)
  
  # Construct initial plot
  p <- ggplot(plotData, aes(x=xAxis, y=yAxis, color=prior, text=pointLabels)) + 
    geom_point(size=0.9) + theme_light() + 
    scale_color_manual(name="Prior", labels=c("0", "NA", "1"), values=c(g, b, r)) + 
    theme(plot.title=element_text(hjust=0.5)) + 
    labs(title=plotTitle, x="LLR<sup>2</sup><sub>other</sub>",
         y="LLR<sup>2</sup> - LLR<sup>2</sup><sub>own</sub>")
  
  # Further settings for interactive scatterplot
  if(interactive == TRUE) {
    ggplotly(p) %>% layout(legend=list(title=list(text="Prior"), orientation="h",
                                       xanchor="center", x=0.5, y=-0.2))
  }
  # Further settings for non-interactive scatterplot
  else {
    p + theme(plot.title=element_markdown(), axis.title.x = element_markdown(), axis.title.y = element_markdown()) +
      theme(legend.position="bottom", legend.background=element_rect(size=0.1, linetype="solid", color="black")) +
      guides(color=guide_legend(override.aes=list(size=3)))
  }
}