## ------------------------------------------------------------------------
if (!("devtools" %in% installed.packages())){install.packages(devtools)}
library(devtools)
if (!("bigsnpr" %in% installed.packages())){devtools::install_github("privefl/bigsnpr")}
if (!("bigstatsr" %in% installed.packages())){devtools::install_github("privefl/bigsnpr")}
if (!("qvalue" %in% installed.packages())){TODO}
if (!("vcfR" %in% installed.packages())){install.packages("vcfR")} 
#if (!("OutFLANK" %in% installed.packages())){devtools::install_github("whitlock/OutFLANK")} 

#devtools::install_github("whitlock/OutFLANK", ref="development", force=TRUE) # will need to DELETE this link for final .Rmd!
library(OutFLANK)  # outflank package
library(vcfR)
library(bigsnpr)   # package for LD pruning
library(bigstatsr) # package for LD pruning 

## ------------------------------------------------------------------------
data("sim1a")
str(sim1a)
data(muts)

# Sample sizes of individuals within populations
table(sim1a$pop)

## ------------------------------------------------------------------------
my_fst <- MakeDiploidFSTMat(t(sim1a$G), locusNames = sim1a$position, popNames = sim1a$pop)

## ---- fig.width=6--------------------------------------------------------
plot(my_fst$He, my_fst$FST)

## ---- fig.width=6--------------------------------------------------------
plot(my_fst$FST, my_fst$FSTNoCorr)

## ----LD pruning----------------------------------------------------------
#### LD Pruning ####
G<-add_code256(big_copy(t(sim1a$G),type="raw"),code=bigsnpr:::CODE_012)
newpc<-snp_autoSVD(G=G,infos.chr =sim1a$chromosome,infos.pos = sim1a$position)
which_pruned <- attr(newpc, which="subset") # Indexes of remaining SNPS after pruning
length(which_pruned)

## ------------------------------------------------------------------------
#### Evaluating OutFLANK with pruned data ####
out_trim <- OutFLANK(my_fst[which_pruned,], NumberOfSamples=39, qthreshold = 0.05)
str(out_trim)
#head(out_trim$results)

## ---- fig.width=6--------------------------------------------------------
OutFLANKResultsPlotter(out_trim, withOutliers = TRUE,
                       NoCorr = TRUE, Hmin = 0.1, binwidth = 0.001, Zoom =
                         FALSE, RightZoomFraction = 0.05, titletext = NULL)

## Zoom in on right tail
OutFLANKResultsPlotter(out_trim , withOutliers = TRUE,
                       NoCorr = TRUE, Hmin = 0.1, binwidth = 0.001, Zoom =
                         TRUE, RightZoomFraction = 0.15, titletext = NULL)

## ------------------------------------------------------------------------
hist(out_trim$results$pvaluesRightTail)

## ---- fig.width=6--------------------------------------------------------
P_all <- pChiSqNoCorr(my_fst, Fstbar = out_trim$FSTNoCorrbar, 
                                   dfInferred = out_trim$dfInferred)

head(P_all)
sum(P_all$FSTNoCorr>0.05)
plot(P_all$FSTNoCorr, P_all$Pval, xlim=c(0,0.2))
hist(P_all$FSTNoCorr[P_all$He>0.1], xlim=c(0,0.2), breaks=100)

## Check the P-value histogram
  hist(P_all$Pval, breaks=50)

## Control for false discovery rate
  q <- qvalue(P_all$Pval, fdr.level = 0.05)
  plot(P_all$FST, q$qvalues)
  P_all$q <-  q$qvalues
  
## My outliers with a q-value of less than 0.01
  my_out <- which(P_all$q < 0.01)

## ---- fig.width=7--------------------------------------------------------
plot(P_all$LocusName[P_all$He>0.1], P_all$FST[P_all$He>0.1],
     xlab="Position", ylab="FST", col=rgb(0,0,0,0.2))
  points(P_all$LocusName[my_out], P_all$FST[my_out], col="magenta", pch=20)  

## ------------------------------------------------------------------------
data(muts)
muts[muts$prop>0.1,]

## ------------------------------------------------------------------------
obj.vcfR <- read.vcfR("../data/sim1a.vcf.gz")

geno <- extract.gt(obj.vcfR) # Character matrix containing the genotypes
position <- getPOS(obj.vcfR) # Positions in bp
chromosome <- getCHROM(obj.vcfR) # Chromosome information

G <- matrix(NA, nrow = nrow(geno), ncol = ncol(geno))

G[geno %in% c("0/0", "0|0")] <- 0
G[geno  %in% c("0/1", "1/0", "1|0", "0|1")] <- 1
G[geno %in% c("1/1", "1|1")] <- 2

table(as.vector(G))

