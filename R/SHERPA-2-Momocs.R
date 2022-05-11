# R script accompanying 
#
# M. Kloster, G. Kauer, B. Beszteri: SHERPA: An image 
# segmentation and outline feature extraction tool for 
# diatoms and other objects
#
# submitted to BMC Bioinformatics
#
# B. Beszteri 2013

##########################################################################
# for using the R package Momocs for further outline analyses:
# load outline points from a set of SHERPA output files
# return as a Coo object

require(Momocs)

load.outlines <- function(indir) {
   outline.files <- dir(path=indir, pattern="*.XY_EFA.csv", full.names=T)
   ou.list <- list()
   for(f in outline.files) {
      outline <- read.csv(f, sep=";", dec=",", na.strings="n. def.")
      ou.mat <- as.matrix(outline)
      f1 <- sub(".*/", "", f)
      f1 <- sub("\\.XY_EFA.csv", "", f1)
      ou.list[[f1]] <- ou.mat
   }
   ou.coo <- Coo(ou.list)
   return(ou.coo)
}

