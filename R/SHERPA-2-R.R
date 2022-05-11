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
# template function to load EFD coefficients from SHERPA output files:
# returns a matrix with objects in rows, EFD coefficients in columns

load.efds <- function(indir, max.rank=2) {
   # files w/ outline coordinates used for EFA:
   efou.files <- dir(path=indir, pattern="*.EFDs.csv", full.names=T)
   # read settings for SHERPA run:
   settings.file <- dir(path=indir, pattern="*.settings.csv", full.names=T)
   settings <- read.csv(settings.file, sep=";", dec=",", quote="\"")
   n.coef <- settings$NumberOfFrequenciesForEFD
   # read results table for rank filtering:
   results.file <- sub("\\.settings", "", settings.file)
   results.tab <- read.csv(results.file, sep=";", dec=",", quote="\"")
   included <- which(results.tab$Ranking.Index <= max.rank)
   # use only objects passing rank filter:
   efou.files <- efou.files[included]
   # set up EFD matrix:
   n.obj <- length(included)
   efd.mat <- matrix(NA, nrow=n.obj, ncol=4*n.coef)
   fnames <- character(n.obj)
   cc <- 1
   for(f in efou.files) {
	  efd.tab <- read.csv(f, sep=";", dec=",", na.strings="n. def.", row.names=1)
      efd.mat[cc,] <- as.vector(as.matrix(efd.tab[2:(n.coef+1),5:8]))
	  fn <- strsplit(f, "/", fixed=T)
	  fnames[cc] <- sub(".EFDs.csv$", "", fn[[1]][length(fn[[1]])])
      cc <- cc + 1
   }
   # convert EFD matrix to data frame:
   efds <- data.frame(efd.mat)
   row.names(efds) <- fnames
   names(efds) <- c(paste("A", 1:n.coef, sep=""),
                       paste("B", 1:n.coef, sep=""),
                       paste("C", 1:n.coef, sep=""),
                       paste("D", 1:n.coef, sep=""))
   return(efds)
}

