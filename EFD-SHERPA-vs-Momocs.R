source("R/SHERPA-2-Momocs.R")
source("R/SHERPA-2-R.R")

indir <- "EFA"
EFDs <- load.efds(indir)

outlines <- load.outlines(indir)

# note that the below needs a fix of a Momocs function
# due to the fact that the efourier function resamples the
# outline points even if this is unnecessary, leading to loss of information
# this should lead to an error message:

efds <- eFourier(outlines, norm=T, start=F, nb.h=30, smooth.it=0)

# to fix this, load a modified version of the Momocs function coo.sample, do:
source("R/eFourier-fix.R")

efds <- eFourier(outlines, norm=T, start=F, nb.h=30, smooth.it=0)

# check that EFDs from SHERPA and from Momocs are virtually identical:
range(efds@coe - EFDs)
