## ------------------------------------------------------------------------
if (!("devtools" %in% installed.packages())){install.packages(devtools)}
library(devtools)

if (!("qvalue" %in% installed.packages())){TODO}
if (!("vcfR" %in% installed.packages())){install.packages("vcfR")} 

devtools::install_github("whitlock/OutFLANK")

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

