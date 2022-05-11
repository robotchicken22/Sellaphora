###############################################################################
# R script accompanying 
#
# M. Kloster, G. Kauer, B. Beszteri: SHERPA: An image 
# segmentation and outline feature extraction tool for 
# diatoms and other objects
#
# submitted to BMC Bioinformatics
#
# B. Beszteri 2013
#
# this example analysis used images accompanying
# Mann DG, McDonald SM, Bayer MM, Droop SJM, Chepurnov VA, Loke RE, 
# Ciobanu A & DuBuf JMH (2004). The Sellaphora pupula complex (Bacillariophyceae): 
# morphometric analysis, ultrastructure and mating data provide evidence for five 
# new species. Phycologia 43: 459-482.
###############################################################################

require(RColorBrewer)
require(rgl)

###############################################################################
# plots from Mann et al. 2004:
# (Fig. 14)

# load csv file exported by SHERPA:
sherpa.tab <- read.csv("EFA/Sellaphora with RATS with contour optimization.new rectangularity.csv", sep=";", dec=",")

sherpa.tab <- sherpa.tab[-344,] # this is a small centric valve, not a Sellaphora

im.map <- sherpa.tab[,c("Source.Image","Contour.Image")]
src.im <- regmatches(im.map[[1]], regexpr("\\d+ca", im.map[[1]], perl=T))
src.im <- sub("0+", "", src.im)
src.im <- sub("ca$", "", src.im)
co.im <- regmatches(im.map[[2]], regexpr("optimization\\.\\d+", im.map[[2]], perl=T))
co.im <- sub("optimization", "sellaphora", co.im)
names(src.im) <- co.im

# read demes table:
demes.tab <- read.table("sellaphora-demes.txt", sep="\t", header=T, row.names=1)
sherpa.tab$Deme <- factor(demes.tab[src.im,1])
sherpa.tab$DemeNr <- factor(demes.tab[src.im,2])
sherpa.tab <- sherpa.tab[sherpa.tab$Ranking.Index < 3,]

demes <- c("rectangular", "capitate", "small", "lanceolate", "obese", "neat")

# setup palette for plotting:
pal <- brewer.pal(6, "Dark2")


# re-scale units to µms:
# make sure you only exeute the following three lines
# once!
scalefact <- 0.55
sherpa.tab$Width <- sherpa.tab$Width * scalefact
sherpa.tab$Height <- sherpa.tab$Height * scalefact

# plot:
par(mfrow=c(2,2))
# Fig 5-6:
boxplot(Width ~ Deme, data=sherpa.tab, col=pal[1:6], main="Fig. 5 (length by deme)")
boxplot(Height ~ Deme, data=sherpa.tab, col=pal[1:6], main="Fig. 6 (width by deme)")
# Fig. 10:
plot(sherpa.tab$Width, sherpa.tab$Height, col=pal[sherpa.tab$Deme], pch=16, xlab="length (µm)", ylab="width (µm)", main="Fig. 10")
legend("bottomright", pch=16, col=pal[1:6], legend=levels(sherpa.tab$Deme))
# Fig. 14:
plot(sherpa.tab$Rectangularity, sherpa.tab$Height, col=pal[sherpa.tab$Deme], pch=16, ylab="width (µm)", xlab="rectangularity", main="Fig. 14")
legend("topleft", pch=16, col=pal[1:6], legend=levels(sherpa.tab$Deme))


###############################################################################
# PCA of EFDs (Fig. 15):
source("R/SHERPA-2-R.R")
indir <- "EFA"
efds <- load.efds(indir, max.rank=2)
efds <- efds[-344,]  # remove centric valve
efd.pc <- prcomp(efds, center=T, scale=F)
windows()
plot(efd.pc$x[, c(1,2)], pch=18, col=pal[sherpa.tab$DemeNr])
legend("topleft", legend=demes, pt.bg=pal, cex=1, pch=rep(22,6), bty="n")

# PCs 1-3:
# rectangular and capitate are quite well separated on PC3
plot3d(efd.pc$x[,1:3], col=pal[sherpa.tab$DemeNr], size=10)

