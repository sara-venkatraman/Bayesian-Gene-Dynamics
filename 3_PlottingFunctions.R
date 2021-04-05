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
library(ggtext)

Expression.Profile.Interpolant <- function(geneName, method="spline") {
  profile <- geneData[geneName,]
  numPoints <- length(profile)
  if(method == "spline")
    interpolant <- splinefun(x=hours, y=profile, method="natural")
  else if(method == "linear")
    interpolant <- approxfun(x=hours, y=profile, method="linear")
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
                            add=F, plotLegend=F, plotGrid=F, 
                            lineLabels=F, pointSize=1.5, lineOpacity=1,
                            useSplines=T, fitSplines=T, plotColors, legendSize, legendPos, 
                            plotTitle, titleSize, genesForExtrema, plotTimes, subtitle) {
  if(monochrome == TRUE) {
    if(missing(plotColors))
      plotColors <- rep(alpha("blue", 0.18), length(genesToPlot))
    else
      plotColors <- rep(alpha(plotColors, 0.2), length(genesToPlot))
  } else {
    if(missing(plotColors)) {
      plotColors <- c(brewer.pal(n=8, name="Dark2"), 
                      brewer.pal(n=8, name="Set2"), 
                      brewer.pal(n=12, name="Paired"))
      if(length(genesToPlot) <= length(plotColors))
        plotColors <- plotColors[1:length(genesToPlot)]
      else
        plotColors <- rep(plotColors, 20)[1:length(genesToPlot)]
    }
  }
  if(lineOpacity != 1)
    plotColors <- alpha(plotColors, lineOpacity)
  if(missing(plotTitle))
    plotTitle = ""
  else if(class(plotTitle) != "expression" && plotTitle == TRUE)
    plotTitle <- strwrap(paste(genesToPlot, collapse=", "))
  
  if(missing(genesForExtrema)) 
    plotExtrema <- Time.Profile.Extrema(genesToPlot)
  else
    plotExtrema <- Time.Profile.Extrema(genesForExtrema)

  if(missing(titleSize))
    titleSize <- 15
  
  if(missing(legendSize))
    legendSize <- 10
  
  if(missing(legendPos))
    legendPos <- "bottom"
  
  if(missing(plotTimes))
    plotTimes <- hours
  
  exprData <- data.frame(plotTimes, t(geneData[genesToPlot, 1:length(plotTimes)]))
  colnames(exprData) <- c("time", genesToPlot)
  
  if(useSplines == T) {
    interpTimes <- seq(from=plotTimes[1], to=plotTimes[length(plotTimes)], length.out=500)
    interpData <- data.frame(interpTimes, matrix(0, nrow=length(interpTimes), ncol=length(genesToPlot)))  
    for(i in 1:length(genesToPlot))
      interpData[,i+1] <- Expression.Profile.Interpolant(genesToPlot[i], method=ifelse(fitSplines == T, "spline", "linear"))(interpTimes)
    colnames(interpData) <- c("time", genesToPlot)
    
    p <- ggplot(melt(interpData, id.var="time"), aes(x=time, y=value, col=variable)) + geom_line() + scale_color_manual(values=plotColors)
  } else {
    p <- ggplot(melt(exprData, id.var="time"), aes(x=time, y=value, col=variable)) + geom_line() + scale_color_manual(values=plotColors)
  }
  if(points == TRUE)
    p <- p + geom_point(data=melt(exprData, id.var="time"), mapping=aes(x=time, y=value, col=variable), size=pointSize) 
  
  p <- p + theme_bw() + labs(title=plotTitle, x="Hours after infection", y=TeX("Expression ($\\log_2$-fold change)")) +
    theme(plot.title=element_text(size=titleSize)) + theme(plot.title=element_text(hjust=0.5)) + theme(plot.title=element_markdown(lineheight=1.1))
  
  if(!missing(subtitle))
    p <- p + labs(subtitle=subtitle) + theme(plot.subtitle=element_text(hjust=0.5)) + theme(plot.subtitle=element_markdown(lineheight=1.1))
  
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

Draw.R2.Scatterplot <- function(matrix1, matrix2, priorMatrix, geneSubset, bayes, interactive=T, plotLegend=T, legendPos="bottom") {
  xAxis <- Vectorize.Labeled.Square.Matrix(matrix1[geneSubset, geneSubset])
  yAxis <- Vectorize.Labeled.Square.Matrix(matrix2[geneSubset, geneSubset])
  if(missing(priorMatrix))
    prior <- matrix(0L, nrow=nrow(xAxis), ncol=1)
  else {
    prior <- Vectorize.Labeled.Square.Matrix(2*priorMatrix[geneSubset, geneSubset])
    prior[is.na(prior)] <- 1
  }
  if(bayes == TRUE)
    plotTitle <- ifelse(interactive == T, "Bayesian lead-lag R^2 values", TeX("(b) Bayesian lead-lag $R^2$ values"))
  else
    plotTitle <- ifelse(interactive == T, "Non-Bayesian lead-lag R^2 values", TeX("(a) Non-Bayesian lead-lag $R^2$ values"))
  plotData <- as.data.frame(cbind(round(xAxis, 3), round(yAxis, 3), prior))
  colnames(plotData) <- c("x.axis", "y.axis", "prior")
  plotData$prior <- as.factor(plotData$prior)
  plotData <- plotData[order(plotData$prior),]
  hoverText <- row.names(plotData)
  
  p <- ggplot(plotData, aes(x=x.axis, y=y.axis, color=prior, text=hoverText)) + 
    geom_point(size=0.9) + theme_light() + ggtitle(plotTitle) + 
    scale_color_manual(name="Prior", labels=c("0", "NA", "1"), values=c(alpha("gray67", 0.9), alpha("navy", 0.65), alpha("orangered3",0.8))) + 
    theme(plot.title = element_text(hjust=0.5))
  
  if(plotLegend == TRUE) {
    p <- p + theme(legend.position=legendPos) + guides(color=guide_legend(override.aes=list(size=3))) +
      theme(legend.background=element_rect(size=0.1, linetype="solid", color="black"))
  } else
    p <- p + theme(legend.position = "none")
  
  if(interactive) {
    p <- p + xlab("LLR2_other") + ylab("LLR2 - LLR2_own")
    return(ggplotly(p))
  } else {
    p <- p + xlab(TeX("LL$R^2_{other}$")) + ylab(TeX("LL$R^2$ - LL$R^2_{own}$"))
    return(p)
  }
}


