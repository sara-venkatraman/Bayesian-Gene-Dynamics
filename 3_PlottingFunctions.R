library(igraph)
library(RColorBrewer)
library(ggplot2)
library(plotly)
library(latex2exp)
library(scales)
library(gplots)
library(reshape2)
library(gridExtra)
library(directlabels)

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

Plot.Gene.Group <- function(genesToPlot, monochrome=F, points=T, 
                            add=F, plotLegend=F, plotGrid=F, gg=F, 
                            lineLabels=F, plotColors, 
                            legendSize, legendPos, plotTitle,
                            titleSize, genesForExtrema) {
  if(monochrome == TRUE) {
    if(missing(plotColors))
      plotColors <- rep(alpha("blue", 0.18), length(genesToPlot))
    else
      plotColors <- rep(alpha(plotColors, 0.2), length(genesToPlot))
  } else {
    if(missing(plotColors)) {
      plotColors <- c(brewer.pal(n=8, name="Dark2"), 
                      brewer.pal(n=8, name="Set2"), 
                      brewer.pal(n=12, name="Paired"),
                      brewer.pal(n=9, name="Pastel1"))[1:length(genesToPlot)] 
    }
  }
  if(missing(plotTitle))
    plotTitle = ""
  else if(plotTitle == TRUE)
    plotTitle <- strwrap(paste(genesToPlot, collapse=", "))
  
  if(missing(genesForExtrema)) 
    plotExtrema <- Time.Profile.Extrema(genesToPlot)
  else
    plotExtrema <- Time.Profile.Extrema(genesForExtrema)

  if(gg == FALSE) {
    if(missing(titleSize))
      titleSize <- 1
    if(missing(legendSize))
      legendSize <- 1
    if(missing(legendPos))
      legendPos <- "bottomright"
    interp <- Expression.Profile.Interpolant(genesToPlot[1])
    profile <- geneData[genesToPlot[1],]
    curve(interp, from=0, to=hours[length(hours)], col=plotColors[1], xlab="Time", 
            ylab="Expression (log-fold)", ylim=c(plotExtrema$min, plotExtrema$max), 
            lwd=1.5, main=plotTitle, cex.main=titleSize, add=add)
    if(plotGrid == T) {
      grid(col="gray88", lwd=1.2, lty=1)
      curve(interp, from=0, to=hours[length(hours)], col=plotColors[1], xlab="Time", 
            ylab="Expression (log-fold)", ylim=c(plotExtrema$min, plotExtrema$max), 
            lwd=1.5, main=plotTitle, cex.main=titleSize, add=T)
    }
    if(points == T) 
      points(hours, profile, pch=19, col=plotColors[1])
    for(i in 2:length(genesToPlot)) {
      interp <- Expression.Profile.Interpolant(genesToPlot[i])
      profile <- geneData[genesToPlot[i],]
      curve(interp, from=0, to=hours[length(hours)], col=plotColors[i], add=T, xlab="Time", ylab="Expression (log-fold)", lwd=1.5)
      if(points == T)
        points(hours, profile, pch=19, col=plotColors[i])
    }
    if(plotLegend == T)
      legend(legendPos, legend=genesToPlot, col=plotColors, fill=plotColors, cex=legendSize)
  } else {
    if(missing(titleSize))
      titleSize <- 15
    if(missing(legendSize))
      legendSize <- 10
    if(missing(legendPos)) {
      legendPos <- "bottom"
    }
    interpTimes <- seq(from=hours[1], to=hours[length(hours)], length.out=500)
    interpData <- data.frame(interpTimes, matrix(0, nrow=length(interpTimes), ncol=length(genesToPlot)))
    exprData <- data.frame(hours, t(geneData[genesToPlot,]))
    for(i in 1:length(genesToPlot))
      interpData[,i+1] <- Expression.Profile.Interpolant(genesToPlot[i])(interpTimes)
    colnames(interpData) <- c("time", genesToPlot)
    colnames(exprData) <- c("time", genesToPlot)
    p <- ggplot(melt(interpData, id.var="time"), aes(x=time, y=value, col=variable)) + geom_line() + scale_color_manual(values=plotColors)
    if(points == TRUE)
      p <- p + geom_point(data=melt(exprData, id.var="time"), mapping=aes(x=time, y=value, col=variable))
    p <- p + ggtitle(plotTitle) + theme_bw() + xlab("Hours after infection") + ylab(TeX("Expression ($\\log_2$-fold change)")) +
      theme(plot.title=element_text(size=titleSize)) + theme(plot.title = element_text(hjust=0.5))
    if(plotGrid == FALSE)
      p <- p + theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank())
    if(lineLabels == TRUE) {
      p <- p + geom_dl(aes(label=variable), method=list(dl.trans(x=x+0.2), "last.points", cex=0.6)) +
        expand_limits(x=53)
      }
    if(plotLegend == TRUE) {
      p <- p + theme(legend.title=element_blank()) + theme(legend.position=legendPos) +
        theme(legend.text=element_text(size=legendSize)) + 
        theme(legend.background=element_rect(size=0.1, linetype="solid", color="black"))
    } else
      p <- p + theme(legend.position = "none")
    return(p)
  }
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

Draw.R2.Scatterplot <- function(matrix1, matrix2, priorMatrix, geneSubset, interactive=T) {
  xAxis <- Vectorize.Labeled.Square.Matrix(matrix1[geneSubset, geneSubset])
  yAxis <- Vectorize.Labeled.Square.Matrix(matrix2[geneSubset, geneSubset])
  if(missing(priorMatrix))
    prior <- matrix(0L, nrow=nrow(xAxis), ncol=1)
  else {
    prior <- Vectorize.Labeled.Square.Matrix(priorMatrix[geneSubset, geneSubset])
    prior[is.na(prior)] <- 0
  }
  if(interactive) {
    plotData <- as.data.frame(cbind(round(xAxis, 3), round(yAxis, 3), prior))
    colnames(plotData) <- c("x.axis", "y.axis", "prior")
    plotData$prior <- as.factor(plotData$prior)
    hoverText <- row.names(plotData)
    p <- ggplot(plotData, aes(x=x.axis, y=y.axis, color=prior, text=hoverText)) + 
      geom_point(size=0.8) + xlab('LLR2_other') + ylab('LLR2 - LLR2_own') + 
      ggtitle("Lead-lag R^2 Values") + theme_light() + theme(legend.position="none") +
      scale_color_manual(values=c("navy", "orangered3"))
    ggplotly(p)
  } else {
    plot(xAxis[prior == 0], yAxis[prior == 0], pch=20, col="navy", xlab=TeX("$LLR2_{other}$"), ylab=TeX("LLR2 - LLR2_{own}$"))
    points(xAxis[prior == 1], yAxis[prior == 1], pch=20, col="orangered3", xlab=TeX("$LLR2_{other}$"), ylab=TeX("LLR2 - LLR2_{own}$"))
    # plot(xAxis, yAxis, pch=20, col=ifelse(prior == 1, "orangered3", "navy"),
    #      xlab=TeX("$LLR2_{other}$"), ylab=TeX("LLR2 - LLR2_{own}$"))
  }
}

