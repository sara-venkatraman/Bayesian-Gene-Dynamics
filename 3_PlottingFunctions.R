library(igraph)
library(RColorBrewer)
library(ggplot2)
library(plotly)

Expression.Profile.Interpolant <- function(geneName) {
  profile <- geneData[geneName,]
  numPoints <- length(profile)
  interpolant <- splinefun(x=hours, y=profile, method="natural")
  return(interpolant)
}

Time.Profile.Extrema <- function(genesToPlot) {
  minima <- c();  maxima <- c()
  for(i in 1:length(genesToPlot)) {
    profile <- geneData[genesToPlot[i],]
    minima[i] <- min(profile);  maxima[i] <- max(profile)
  }
  return(list(min=min(minima), max=max(maxima)))
}

Plot.Gene.Group <- function(genesToPlot) {
  plotColors <- c(brewer.pal(n=8, name="Dark2"), brewer.pal(n=8, name="Set2"), brewer.pal(n=12, name="Paired"))[1:length(genesToPlot)]
  plotTitle <- strwrap(paste(genesToPlot, collapse=", "))
  plotExtrema <- Time.Profile.Extrema(genesToPlot)
  interp <- Expression.Profile.Interpolant(genesToPlot[1])
  profile <- geneData[genesToPlot[1],]
  curve(interp, from=0, to=hours[length(hours)], col=plotColors[1], xlab="Time", ylab="Expression (log-fold)", ylim=c(plotExtrema$min, plotExtrema$max), lwd=1.5, main=plotTitle, cex.main=.5)
  points(hours, profile, pch=20, col=plotColors[1])
  for(i in 2:length(genesToPlot)) {
    interp <- Expression.Profile.Interpolant(genesToPlot[i])
    profile <- geneData[genesToPlot[i],]
    curve(interp, from=0, to=hours[length(hours)], col=plotColors[i], add=T, xlab="Time", ylab="Expression (log-fold)", lwd=1.5)
    points(hours, profile, pch=20, col=plotColors[i])
  }
  # legend("bottomright", legend=genesToPlot, col=plotColors, fill=plotColors, cex=.6) 
}

Vectorize.Labeled.Square.Matrix <- function(labeledMatrix) {
  matLabels <- rownames(labeledMatrix)
  size <- nrow(labeledMatrix)
  matVector <- matrix(0, nrow=sum(1:(size-1)))
  pairLabels <- c();  k <- 1
  for(i in 1:(size-1)) {
    for(j in (i+1):size) {
      matVector[k] <- labeledMatrix[i,j]
      pairLabels[k] <-  paste(matLabels[i], ", ", matLabels[j], sep="")
      k <- k+1 } }
  rownames(matVector) <- pairLabels
  matVector
}

Draw.R2.Scatterplot <- function(matrix1, matrix2, priorMatrix, geneSubset) {
  xAxis <- Vectorize.Labeled.Square.Matrix(matrix1[geneSubset, geneSubset])
  yAxis <- Vectorize.Labeled.Square.Matrix(matrix2[geneSubset, geneSubset])
  if(missing(priorMatrix))
    prior <- matrix(0L, nrow=nrow(xAxis), ncol=1)
  else {
    prior <- Vectorize.Labeled.Square.Matrix(priorMatrix[geneSubset, geneSubset])
    prior[is.na(prior)] <- 0
  }
  plotData <- as.data.frame(cbind(round(xAxis, 3), round(yAxis, 3), prior))
  colnames(plotData) <- c("x.axis", "y.axis", "prior")
  plotData$prior <- as.factor(plotData$prior)
  hoverText <- row.names(plotData)
  p <- ggplot(plotData, aes(x=x.axis, y=y.axis, color=prior, text=hoverText)) + 
    geom_point(size=0.8) + xlab('LLR2_other') + ylab('LLR2 - LLR2_own') + 
    ggtitle("Lead-lag R^2 Values") + theme_light() + theme(legend.position="none") +
    scale_color_manual(values=c("navy", "orangered3"))
  ggplotly(p)
}
