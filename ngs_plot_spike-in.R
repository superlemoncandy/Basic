# Modified from ngs.plot.r, spike-in.R(MSQ)
# Purpose: Plot ngs.plot in a prettier way and normalize to spike-in

# Jan 2024.

# Run ngs.plot.r in the shell (for all samples) first and import avgprof.RData
# What is the avgprof.RData normalized to? (RPKM)

# content of Avg. profile R data:
#   plot.width, plot.height
#   regcovMat
#   ctg.tbl
#   bam.pair: boolen for bam-pair
#   xticks <- genXticks(reg2plot, pint, lgint, pts, flanksize, flankfactor, Labs)
#   pts: data points for avg. profile and standard errors
#   m.pts, f.pts
#   pint: point interval
#   shade.alp
#   confiMat
#   mw
#   prof.misc
#   se
#   v.lib.size
#   font.size,file=prof.dat

# color
library(ggsci)
npg <- pal_npg("nrc",alpha=0.7)(10)
npg <- c("#E64B35B2", "#4DBBD5B2", "#00A087B2", "#3C5488B2", "#F39B7FB2", "#8491B4B2", "#91D1C2B2",  
            "#DC0000B2", "#7E6148B2", "#B09C85B2")  



progpath <- "ngs_scripts"
source(file.path(progpath, 'parse.args.r'))
source(file.path(progpath, 'genedb.r'))
source(file.path(progpath, 'plotlib.r'))
source(file.path(progpath, 'coverage.r'))

# Matrix for avg. profiles: regcovMat.
# Matrix for standard errors: confiMat.
# Library size for normalization, calculate in QC_summary file: v.lib.size
factors <- read.csv("facto.csv")
v.lib.size <- factors$F1.factor
spike_in.factor <- factors$bw_scaleFactor

# import avg. Rdata
# Same
regcovMat1 <- regcovMat
confiMat1 <- confiMat

# Opp
regcovMat2 <- regcovMat
confiMat2 <- confiMat

# check if the two Matrics are different
identical(regcovMat1, regcovMat2)
identical(confiMat2, confiMat2) # True?

# normalize to spike-in
regcovMat3 <- t(t(regcovMat1)*v.lib.size)
regcovMat3 <- t(t(regcovMat3)*spike_in.factor)
# 注意此处并非矩阵乘法，而是元素级乘法

regcovMat4 <- t(t(regcovMat2)*v.lib.size)
regcovMat4 <- t(t(regcovMat4)*v.lib.size)

# Combine
regcovMat_Comb = cbind(regcovMat3, regcovMat4)
confiMat_Comb = cbind(confiMat1, confiMat2)

# Create image file and plot data into it
out.plot <- 'test.pdf'
pdf(out.plot, width=plot.width, height=plot.height, pointsize=font.size)
plotmat(regcovMat=regcovMat_Comb[,c(1:6)], ctg.tbl$title[c(1:6)], plot.colors = npg[1:6], bam.pair, xticks, pts, m.pts, f.pts, pint, shade.alp, confiMat=confiMat_Comb[,c(1:6)], mw, prof.misc)
out.dev <- dev.off()
