This folder contains R code and data accompanying
M. Kloster, G. Kauer, B. Beszteri: SHERPA: An image segmentation and 
  outline feature extraction tool for diatoms and other objects

submitted to BMC Bioinformatics and illustrate the Sellaphora analysis from the paper
which were using images obtained by D.G. Mann et al. for their 2004 paper


Mann DG, McDonald SM, Bayer MM, Droop SJM, Chepurnov VA, Loke RE, 
  Ciobanu A & DuBuf JMH (2004). The Sellaphora pupula complex (Bacillariophyceae): 
  morphometric analysis, ultrastructure and mating data provide evidence for five 
  new species. Phycologia 43: 459-482.

For reproduing the analyses presented in the paper, start R, and change working directory
into this folder. Then open the script file Sellaphora.R and execute it line by line using Ctrl+R.
The required packages can be installed by issuing 

> install.packages(c("rgl", "RColorBrewer"))

on the R prompt.

Besides these packages, the script uses a template function showing how elliptic Fourier 
descriptor invariants written by SHERPA can be read into R for further analyses. This function 
is stored in the file SHERPA-2-R.R in the R subfolder and can be executed similarly to the above
(after changing working directory into this folder).


Another script (EFD-SHERPA-vs-Momocs.R) illustrates how outline points exported from SHERPA
can be read into a Coo data structure of the Momocs R package by V. Bonhomme & al. The core 
data import function can be found in SHERPA-2-Momocs.R. For this to work, you need to have 
Momocs installed; if this is not the case, do

> install.packages("Momocs")

Please note that the script needs a fix of a Momocs function due to the fact that the efourier
function resamples the outline points even if this is unnecessary, leading to loss of 
information, which causes an error message. 
To fix this, please follow the instructions in the script file.

