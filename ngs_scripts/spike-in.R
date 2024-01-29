# import your avgprof.RData first from ngs.plot.r output
# import function
library(caTools)
progpath <- 'ngs_scripts/'
source(file.path(progpath, 'parse.args.r'))
source(file.path(progpath, 'genedb.r'))
source(file.path(progpath, 'plotlib.r'))
source(file.path(progpath, 'coverage.r'))
plotmat <- function(regcovMat, title2plot, plot.colors, bam.pair, xticks, 
                    pts, m.pts, f.pts, pint, shade.alp=0, confiMat=NULL, mw=1, 
                    misc.options=list(yscale='auto', legend=T, box=T, vline=T, 
                                      xylab=T, line.wd=3)) {
  # Plot avg. profiles and standard errors around them.
  # Args:
  #   regcovMat: matrix for avg. profiles.
  #   title2plot: profile names, will be shown in figure legend.
  #   plot.colors: vector of color specifications for all curves.
  #   bam.pair: boolean for bam-pair data.
  #   xticks: X-axis ticks.
  #   pts: data points
  #   m.pts: middle part data points
  #   f.pts: flanking part data points
  #   pint: tag for point interval
  #   shade.alp: shading area alpha
  #   confiMat: matrix for standard errors.
  #   mw: moving window size for smoothing function.
  #   misc.options: list of misc. options - y-axis scale, legend, box around plot, 
  #       verticle lines, X- and Y-axis labels, line width.
  
  # Smooth avg. profiles if specified.
  if(mw > 1){
    regcovMat <- as.matrix(runmean(regcovMat, k=mw, alg='C', 
                                   endrule='mean'))
  }
  
  # Choose colors.
  if(any(is.na(plot.colors))) {
    ncurve <- ncol(regcovMat)
    if(ncurve <= 8) {
      suppressMessages(require(RColorBrewer, warn.conflicts=F))
      col2use <- brewer.pal(ifelse(ncurve >= 3, ncurve, 3), 'Dark2')
      col2use <- col2use[1:ncurve]
    } else {
      col2use <- rainbow(ncurve)
    }
  } else {
    col2use <- plot.colors
  }
  col2use <- col2alpha(col2use, 0.8)
  
  # Plot profiles.
  ytext <- ifelse(bam.pair, "log2(Fold change vs. control)", 
                  "Read count Per Million mapped reads") # Read count Per Million mapped reads.   Sense/Anti-sense ratio
  xrange <- 0:pts
  y.lim <- NULL
  if(length(misc.options$yscale) == 2) {
    y.lim <- misc.options$yscale
  }
  matplot(xrange, regcovMat, 
          xaxt='n', type="l", col=col2use, ylim=y.lim,
          lty="solid", lwd=misc.options$line.wd, frame.plot=F, ann=F, cex = 2)
  if(misc.options$xylab) {
    title(xlab="Genomic Region (5' -> 3')", ylab=ytext)
  }
  axis(1, at=xticks$pos, labels=xticks$lab, lwd=3, lwd.ticks=3)
  if(misc.options$box) {
    # box around plot.
    box()
  }
  
  # Add shade area.
  if(shade.alp > 0){
    for(i in 1:ncol(regcovMat)){
      v.x <- c(xrange[1], xrange, xrange[length(xrange)])
      v.y <- regcovMat[, i]
      v.y <- c(0, v.y, 0)
      col.rgb <- col2rgb(col2use[i])
      p.col <- rgb(col.rgb[1, 1], col.rgb[2, 1], col.rgb[3, 1], 
                   alpha=shade.alp * 255, maxColorValue=255)
      polygon(v.x, v.y, density=-1, border=NA, col=p.col)
    }
  }
  
  # Add standard errors.
  if(!is.null(confiMat)){
    v.x <- c(xrange, rev(xrange))
    for(i in 1:ncol(confiMat)){
      v.y <- c(regcovMat[, i] + confiMat[, i], 
               rev(regcovMat[, i] - confiMat[, i]))
      col.rgb <- col2rgb(col2use[i])
      p.col <- rgb(col.rgb[1, 1], col.rgb[2, 1], col.rgb[3, 1], 
                   alpha=0.2 * 255, maxColorValue=255)
      polygon(v.x, v.y, density=-1, border=NA, col=p.col)
    }
  }
  
  if(misc.options$vline) {
    # Add gray lines indicating feature boundaries.
    if(pint) {
      abline(v=f.pts, col="gray", lwd=2)
    } else {
      abline(v=f.pts - 1, col="gray", lwd=2)
      abline(v=f.pts + m.pts - 2, col="gray", lwd=2)
    }
  }
  
  if(misc.options$legend) {
    # Legend.
    legend("topright", title2plot, text.col=col2use)
  }
}

#CPM return to count 
#regcovMat2<- t(t(regcovMat)*v.lib.size)

# your scale factor 
#scale.factor <- c(1.000,1.077,0.892,0.814,0.834 ,0.787 ,0.828 ,0.797 ,0.769 ,
#                  0.751 ,0.869 ,0.867 ,0.859 ,0.909 ,1.010 ,0.816 ,0.961,
#                  1.195 ,0.895 ,1.129 ,0.738 ,1.321 ,1.188 ,1.251 ,2.092,0.880 ,1.269 ,0.908)

v.lib.size <- c(0.293501438, 0.701632865,0.651298554,1.269100272,0.888765669,1,0.495115869,
                0.553770819, 1.70207875,1.490058412,0.996417176,0.626633071,2.127043951,0.791760634,1.513450775,
                0.93275325, 1.234578182)

#v.lib.size <- c(1,1.70207875,0.93275325,1.490058412,0.495115869,1.234578182,0.293501438,0.996417176)


scale.factor <- c(0.747037946, 0.524566014, 0.528231178, 1.075928653, 1.390177927,1,1.132923739,
                  0.917350636,0.796787753,0.617696787,0.719544608,0.545587757,0.885339643,0.758129339,0.778126961,
                  1.215715968, 0.792122022)

# multiply by scale factor
regcovMat<- t(t(regcovMat)*v.lib.size)

regcovMat2<- t(t(regcovMat)*scale.factor)
confiMat2 <- confiMat

# import the antisense
regcovMat<- -t(t(regcovMat)*v.lib.size)
regcovMat<- t(t(regcovMat)*scale.factor)

regcovComb <- cbind(regcovMat2, regcovMat) 
confiMatComb <- cbind(confiMat2, confiMat)

#ratio <- regcovMat2/-regcovMat


#regcovMat2<- regcovMat2/1e7
out.plot <- 'figures/sc_sense_anti.pdf'
pdf(out.plot, width=plot.width, height=plot.height, pointsize=font.size)

c("#D1E5F0","#92C5DE", "#B2182B","#D6604D",
  "#F4A582","#FDDBC7" ,"#4393C3","#2166AC")

t <- c('#999999','#E69F00','#56B4E9', 
'#800080','#ADD8E6','#FFC0CB', "#009999")

prof.misc$legend <- FALSE
plotmat(regcovMat2[,c(1:7,16)], ctg.tbl$title[c(1:7,16)], plot.colors = c('#999999',"#FFABC7", "#FFABC7","#D6604D","#D6604D",
                                                                                   "#B2182B","#B2182B","#92C5DE"), bam.pair, xticks, pts, 
        m.pts, f.pts, pint, shade.alp, confiMat = confiMat2[,c(1:7,16)], mw = 3, prof.misc)



plotmat(regcovComb[,c(1:7,16,18:24,33)], ctg.tbl$title[c(1:7,16)], 
        plot.colors = c('#999999',"#FFABC7", "#FFABC7","#D6604D","#D6604D",
                                                   "#B2182B","#B2182B","#92C5DE"), bam.pair, xticks, pts, 
        m.pts, f.pts, pint, shade.alp, confiMat = confiMatComb[,c(1:7,16,18:24,33)], mw = 3, prof.misc)

out.dev <- dev.off()


plotmat(regcovMat2[,c(8:15,17)], ctg.tbl$title[c(8:15,17)], plot.colors = c("#009999","#009999","#C1E999","#C1E999",'#E99F00','#E99F00',
                                                                                     '#800080','#800080','#B29AE6'), bam.pair, xticks, pts, 
        m.pts, f.pts, pint, shade.alp, confiMat = confiMat2[,c(8:15,17)], mw = 3, prof.misc)

plotmat(log10(ratio_EU), ctg.tbl$title, 
        plot.colors = c("#009999","#C1E999","#B2182B",'#E99F00',
                        '#800080','#B29AE6'), bam.pair, xticks, pts, 
        m.pts, f.pts, pint, shade.alp, confiMat = NULL, mw = 3, prof.misc)

plotmat(log10(ratio[,c(1:7,16)]), ctg.tbl$title[c(1:7,16)], 
        plot.colors = c('#999999',"#FFABC7", "#FFABC7","#D6604D","#D6604D",
                        "#B2182B","#B2182B","#92C5DE"), bam.pair, xticks, pts, 
        m.pts, f.pts, pint, shade.alp, confiMat = NULL, mw = 3, prof.misc)

plotmat(log10(ratio[,c(8:15,17)]), ctg.tbl$title[c(8:15,17)], plot.colors = c("#009999","#009999","#C1E999","#C1E999",'#E99F00','#E99F00',
                                                                            '#800080','#800080','#B29AE6'), bam.pair, xticks, pts, 
        m.pts, f.pts, pint, shade.alp, confiMat = NULL, mw = 3, prof.misc)

plotmat(regcovComb[,c(8:15,17,25:32,34)], ctg.tbl$title[c(8:15,17)], plot.colors = c("#009999","#009999","#C1E999","#C1E999",'#E99F00','#E99F00',
                                                   '#800080','#800080','#B29AE6'), bam.pair, xticks, pts, 
        m.pts, f.pts, pint, shade.alp, confiMat = confiMatComb[,c(8:15,17,25:32,34)], mw = 3, prof.misc)

plotmat(regcovComb[,c(8:15,17,25:32,34)], ctg.tbl$title[c(8:15,17)], plot.colors = c("#009999","#FFABC7","#C1E999","#D6604D",'#92C5DE','#E99F00',
                                                                                     '#800080','#B2182B','#B29AE6'), bam.pair, xticks, pts, 
        m.pts, f.pts, pint, shade.alp, confiMat = confiMatComb[,c(8:15,17,25:32,34)], mw = 3, prof.misc)


plotmat(recovRatio, ctg.tbl$title, plot.colors = c('#999999','#999999', "#D6604D","#D6604D",
                                                  "#009999","#009999","#B2182B","#B2182B"), bam.pair, xticks, pts, 
        m.pts, f.pts, pint, shade.alp, confiMat = NULL, mw = 3, prof.misc)

# plot sense antisense simutaneously:
out.plot <- 'figures/bivalent_sense_anti.pdf'
pdf(out.plot, width=plot.width, height=plot.height, pointsize=font.size)
regcovMat2<- t(t(regcovMat)*scale.factor)
confiMat2 <- confiMat
regcovMat<- t(t(regcovMat)*scale.factor)
regcovMat <- -regcovMat
regcovComb <- cbind(regcovMat2, regcovMat) 
confiMatComb <- cbind(confiMat2, confiMat)
plotmat(regcovComb, ctg.tbl$title, plot.colors = c('#999999','#999999', "#D6604D","#D6604D",
                                                   "#009999","#009999","#B2182B","#B2182B"), bam.pair, xticks, pts, 
        m.pts, f.pts, pint, shade.alp, confiMat = confiMatComb, mw = 3, prof.misc)
abline(h = 0, col="gray", lwd=2)

#no legend
plotmat(regcovComb, title2plot = NULL, plot.colors = c('#999999','#999999', "#D6604D","#D6604D",
                                                       "#009999","#009999","#B2182B","#B2182B"), bam.pair, xticks, pts, 
        m.pts, f.pts, pint, shade.alp, confiMat = confiMatComb, mw = 3, prof.misc)


# plotmat(regcovMat2, ctg.tbl$title, ctg.tbl$color, bam.pair, xticks, pts, 
#         m.pts, f.pts, pint, shade.alp, confiMat, mw, prof.misc)
abline(v = 40, col="red", lwd=2)
abline(h = 0.075, col="red", lwd=2)

out.dev <- dev.off()


 plotmat <- function(regcovMat, title2plot, plot.colors, bam.pair, xticks, 
                    pts, m.pts, f.pts, pint, shade.alp=0, confiMat=NULL, mw=1, 
                    misc.options=list(yscale='auto', legend=T, box=T, vline=T, 
                                      xylab=T, line.wd=3)) {
  # Plot avg. profiles and standard errors around them.
  # Args:
  #   regcovMat: matrix for avg. profiles.
  #   title2plot: profile names, will be shown in figure legend.
  #   plot.colors: vector of color specifications for all curves.
  #   bam.pair: boolean for bam-pair data.
  #   xticks: X-axis ticks.
  #   pts: data points
  #   m.pts: middle part data points
  #   f.pts: flanking part data points
  #   pint: tag for point interval
  #   shade.alp: shading area alpha
  #   confiMat: matrix for standard errors.
  #   mw: moving window size for smoothing function.
  #   misc.options: list of misc. options - y-axis scale, legend, box around plot, 
  #       verticle lines, X- and Y-axis labels, line width.
  
  # Smooth avg. profiles if specified.
  if(mw > 1){
    regcovMat <- as.matrix(runmean(regcovMat, k=mw, alg='C', 
                                   endrule='mean'))
  }
  
  # Choose colors.
  if(any(is.na(plot.colors))) {
    ncurve <- ncol(regcovMat)
    if(ncurve <= 8) {
      suppressMessages(require(RColorBrewer, warn.conflicts=F))
      col2use <- brewer.pal(ifelse(ncurve >= 3, ncurve, 3), 'Dark2')
      col2use <- col2use[1:ncurve]
    } else {
      col2use <- rainbow(ncurve)
    }
  } else {
    col2use <- plot.colors
  }
  col2use <- col2alpha(col2use, 0.8)
  
  # Plot profiles.
  ytext <- ifelse(bam.pair, "log2(Fold change vs. control)", 
                  "Read count Per Million mapped reads")
  xrange <- 0:pts
  y.lim <- NULL
  if(length(misc.options$yscale) == 2) {
    y.lim <- misc.options$yscale
  }
  matplot(xrange, regcovMat, 
          xaxt='n', type="l", col=col2use, ylim=y.lim,
          lty="solid", lwd=misc.options$line.wd, frame.plot=F, ann=F)
  if(misc.options$xylab) {
    title(xlab="Genomic Region (5' -> 3')", ylab=ytext)
  }
  axis(1, at=xticks$pos, labels=xticks$lab, lwd=3, lwd.ticks=3)
  if(misc.options$box) {
    # box around plot.
    box()
  }
  
  # Add shade area.
  if(shade.alp > 0){
    for(i in 1:ncol(regcovMat)){
      v.x <- c(xrange[1], xrange, xrange[length(xrange)])
      v.y <- regcovMat[, i]
      v.y <- c(0, v.y, 0)
      col.rgb <- col2rgb(col2use[i])
      p.col <- rgb(col.rgb[1, 1], col.rgb[2, 1], col.rgb[3, 1], 
                   alpha=shade.alp * 255, maxColorValue=255)
      polygon(v.x, v.y, density=-1, border=NA, col=p.col)
    }
  }
  
  # Add standard errors.
  if(!is.null(confiMat)){
    v.x <- c(xrange, rev(xrange))
    for(i in 1:ncol(confiMat)){
      v.y <- c(regcovMat[, i] + confiMat[, i], 
               rev(regcovMat[, i] - confiMat[, i]))
      col.rgb <- col2rgb(col2use[i])
      p.col <- rgb(col.rgb[1, 1], col.rgb[2, 1], col.rgb[3, 1], 
                   alpha=0.2 * 255, maxColorValue=255)
      polygon(v.x, v.y, density=-1, border=NA, col=p.col)
    }
  }
  
  if(misc.options$vline) {
    # Add gray lines indicating feature boundaries.
    if(pint) {
      abline(v=f.pts, col="gray", lwd=2)
    } else {
      abline(v=f.pts - 1, col="gray", lwd=2)
      abline(v=f.pts + m.pts - 2, col="gray", lwd=2)
    }
  }
  
  if(misc.options$legend) {
    # Legend.
    legend("topright", title2plot, text.col=col2use)
  }
}



