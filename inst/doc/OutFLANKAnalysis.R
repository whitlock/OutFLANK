## ------------------------------------------------------------------------
if (!("devtools" %in% installed.packages())){install.packages(devtools)}
library(devtools)

if (!("qvalue" %in% installed.packages())){TODO}
if (!("vcfR" %in% installed.packages())){install.packages("vcfR")} 

#if (!("OutFLANK" %in% installed.packages())){devtools::install_github("whitlock/OutFLANK")} ## will need to UNCOMMENT this link for final .Rmd!

devtools::install_github("whitlock/OutFLANK", ref="development", force=TRUE) # will need to DELETE this link for final .Rmd!

library(OutFLANK)  # outflank package
library(vcfR)

## ------------------------------------------------------------------------
data("sim1a")
str(sim1a)

# Sample sizes of individuals within populations
table(sim1a$pop)

## ------------------------------------------------------------------------
my_fst <- MakeDiploidFSTMat(t(sim1a$G), locusNames = sim1a$position, popNames = sim1a$pop)
head(my_fst)

## ---- fig.width=6--------------------------------------------------------
plot(my_fst$He, my_fst$FST)

## ---- fig.width=6--------------------------------------------------------
plot(my_fst$FST, my_fst$FSTNoCorr)
abline(0,1)

## ------------------------------------------------------------------------
data("which_pruned")
head(which_pruned)

## ------------------------------------------------------------------------
#### Evaluating OutFLANK with trimmed SNPs ####
out_trim <- OutFLANK(my_fst[which_pruned,], NumberOfSamples=39, qthreshold = 0.05, Hmin = 0.1)
str(out_trim)
head(out_trim$results)

## ---- fig.width=6--------------------------------------------------------
OutFLANKResultsPlotter(out_trim, withOutliers = TRUE,
                       NoCorr = TRUE, Hmin = 0.1, binwidth = 0.001, Zoom =
                         FALSE, RightZoomFraction = 0.05, titletext = NULL)

## Zoom in on right tail
OutFLANKResultsPlotter(out_trim , withOutliers = TRUE,
                       NoCorr = TRUE, Hmin = 0.1, binwidth = 0.001, Zoom =
                         TRUE, RightZoomFraction = 0.15, titletext = NULL)

## ---- fig.width=6--------------------------------------------------------
hist(out_trim$results$pvaluesRightTail)

## ---- fig.width=6--------------------------------------------------------
P1 <- pOutlierFinderChiSqNoCorr(my_fst, Fstbar = out_trim$FSTNoCorrbar, 
                                   dfInferred = out_trim$dfInferred, qthreshold = 0.05, Hmin=0.1)
head(P1)
tail(P1)
# notice how the output is ordered differently

my_out <- P1$OutlierFlag==TRUE
plot(P1$He, P1$FST, pch=19, col=rgb(0,0,0,0.1))
points(P1$He[my_out], P1$FST[my_out], col="blue")

hist(P1$pvaluesRightTail)
# check the P-value histogram for the full set of data
# if there are outliers, it should look flat with an inflation near 0

## ---- fig.width=7--------------------------------------------------------
plot(P1$LocusName[P1$He>0.1], P1$FST[P1$He>0.1],
     xlab="Position", ylab="FST", col=rgb(0,0,0,0.2))
  points(P1$LocusName[my_out], P1$FST[my_out], col="magenta", pch=20)  

## ------------------------------------------------------------------------
data("muts")
muts[muts$prop>0.1,]

