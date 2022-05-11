# this file contains code from the Momcs R package by Bonhomme et al.
# see the description file of the package below
# the function coo.sample has been modified to correctly
# handle the case when a opened set of points (i.e., start
# and end point are not identical)needs to be resampled 
# to the identical number of points
# the other two functions / setMethod calls are identically 
# copied from Momocs and are only included to make sure that
# the modified version of coo.sample is called when using
# efourier or eFourier
#
# B. Beszteri 
# Alfred Wegener Institute Helmholtz 
# Centre for Polar and Marine Research
# 2013

#Package: Momocs
#Type: Package
#Title: Shape Analysis of Outlines
#Version: 0.2-04
#Date: 2013-05-06
#Author: Vincent Bonhomme, Sandrine Picq, Julien Claude with
#        contributions from Dan Dkin, Cedric Gaucherel, Ricardo Kriebel,
#        Neus Martinez, Marcelo Reginato, Norbert Telmon, Asher
#        Wishkerman.
#Maintainer: Vincent Bonhomme <bonhomme.vincent@gmail.com>
#Description: Momocs is intended to ease and popularize shape analysis
#        of outlines (especially using elliptical Fourier analysis). It
#        mostly hinges on the functions developed in Morphometrics with
#        R (Claude, 2008). From outline extraction of images and
#        elliptical Fourier calculation to multivariate analysis and the
#        visualization of transformations within the morphological
#        space, Momocs provides a complete and convenient toolkit to
#        specialists within every field that are, or may be, interested
#        in morphological comparisons of outlines.
#License: GPL (>= 2)
#URL: http://www.vincentbonhomme.fr/Momocs
#Depends: methods, sp, jpeg, ade4, spdep, shapes, ape, R(>= 2.15.0)
#Collate: global.R Coe.R Coo.R
#Packaged: 2013-05-08 11:57:47 UTC; vincent
#NeedsCompilation: no
#Repository: CRAN
#Date/Publication: 2013-05-08 15:47:08


coo.sample     <- function (coo, n) {
  if (is.matrix(coo)){
    if (n == nrow(coo) && any (coo[1,] != coo[nrow(coo),])) return(coo);
	sampled <- round(seq(1, nrow(coo), len = n + 1)[-(n + 1)])
    return(coo[sampled, ])
  } else if (is.list(coo)) {
    if (n == length(coo$x) && (coo$x[1] != coo$x[length(coo$x)] || coo$y[1] != coo$y[length(coo$y)])) return(coo);
    sampled <- round(seq(1, length(coo$x), len = n + 1)[-(n + 1)])
    return(list(x=coo$x[sampled], y=coo$y[sampled]))  
  } else stop("A list of a coordinate matrix must be provided to coo.sample")}


efourier  <- function (coo, nb.h = 32, smooth.it = 0, silent = FALSE) {
  if (is.matrix(coo)) coo <- m2l(coo)
  if (is.closed(coo)) coo <- coo.unclose(coo)
  if (missing(nb.h))  {
    nb.h <- length(coo$x)/2 - 1 # should not be 1
    warning(paste("nb.h not provided and set to", nb.h))}
  if(nb.h * 2 > length(coo$x)) {
    nb.h = floor(length(coo$x)/2)-1 # should not be -1
    if (!silent){
    warning("The number of harmonics to calculate should be lower than half the number of points. 
    The number of harmonics used has been set to: ", nb.h)}}
  if (nb.h == -1) {
    nb.h = floor(length(coo$x)/2)-1 # should not be -1
    if (!silent){
    cat("The number of harmonics used has been set to: ", nb.h)}}
  if (smooth.it!=0) { coo <- coo.smooth(coo, smooth.it)}
  coo <- coo.sample(coo, nb.h * 2)
  p <- length(coo$x)
  Dx <- coo$x - coo$x[c(p, (1:p - 1))]
  Dy <- coo$y - coo$y[c(p, (1:p - 1))]
  Dt <- sqrt(Dx^2 + Dy^2)
  t1 <- cumsum(Dt)
  t1m1 <- c(0, t1[-p])
  T <- sum(Dt)
  an <- bn <- cn <- dn <- numeric(nb.h)
  for (i in 1:nb.h) {
    Ti <- (T/(2 * pi^2 * i^2))
    r <- 2 * i * pi
    an[i] <- Ti * sum((Dx/Dt) * (cos(r * t1/T) - cos(r * t1m1/T)))
    bn[i] <- Ti * sum((Dx/Dt) * (sin(r * t1/T) - sin(r * t1m1/T)))
    cn[i] <- Ti * sum((Dy/Dt) * (cos(r * t1/T) - cos(r * t1m1/T)))
    dn[i] <- Ti * sum((Dy/Dt) * (sin(r * t1/T) - sin(r * t1m1/T)))
  }
  ao <- 2 * sum(coo$x * Dt/T)
  co <- 2 * sum(coo$y * Dt/T)
  return(list(an = an, bn = bn, cn = cn, dn = dn, ao = ao, co = co))}
 
 
setMethod(f="eFourier", signature="Coo", definition=
  function(Coo,
           nb.h      = 32,
           smooth.it = 0,
           norm       = TRUE,
           start     = FALSE) {
    q <- floor(min(Coo@coo.len)/2) - 1
    if (missing(nb.h))  {
      nb.h <- if (q >= 32) { 32 } else { q }
      cat(paste("  * nb.h not provided and set to", nb.h))}
    if(nb.h  > (q+1)*2) {
      nb.h <- q # should not be 1
      warning("At least one outline has ", (q+1)*2, " coordinates. The number of harmonics has been set to: ", q)}
    coo <- Coo@coo
    col.n <- paste0(rep(LETTERS[1:4], each = nb.h), rep(1:nb.h, times = 4))
    coe <- matrix(ncol = 4 * nb.h, nrow = length(coo), dimnames = list(names(coo), col.n))
    for (i in seq(along = coo)) {
      ef <- efourier(coo[[i]], nb.h = nb.h, smooth.it = smooth.it, silent = FALSE)
      if (norm) {
        ef <- efourier.norm(ef, start=start)
        if (ef$A[1] < 0) {
          ef$A <- (-ef$A)
          ef$B <- (-ef$B)
          ef$C <- (-ef$C)
          ef$D <- (-ef$D)
          ef$lnef <- (-ef$lnef)}
        coe[i, ] <- c(ef$A, ef$B, ef$C, ef$D)
      } else {
        coe[i, ] <- c(ef$an, ef$bn, ef$cn, ef$dn)}}
    return(Coe(coe, fac=Coo@fac, method="eFourier"))})
