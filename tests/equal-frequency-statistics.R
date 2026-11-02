library(OutFLANK)

# Each population has two individuals of each genotype (0, 1 and 2).
# Equal population frequencies must not become eight zero statistics.
counts <- matrix(2, nrow = 2, ncol = 3)
result <- WC_FST_Diploids_2Alleles(counts)
expected <- list(
  He = 1/2, FST = -2/15, T1 = -1/30, T2 = 1/4,
  FSTNoCorr = 0, T1NoCorr = 0, T2NoCorr = 1/4,
  meanAlleleFreq = 1/2
)
stopifnot(isTRUE(all.equal(result, expected, tolerance = 1e-12)))

# Check propagation through the public multi-locus interface.
genotypes <- matrix(rep(c(0, 1, 2), 4), ncol = 1)
populations <- rep(c("A", "B"), each = 6)
observed <- MakeDiploidFSTMat(genotypes, "equal_frequencies", populations)
stopifnot(nrow(observed) == 1L,
          identical(as.character(observed$LocusName), "equal_frequencies"),
          isTRUE(all.equal(as.numeric(observed[1, -1]),
                           unname(unlist(expected)), tolerance = 1e-12)))

# Unequal frequencies: these expected values are unchanged by the fix.
# Population counts are (3,2,1) and (1,2,3), with six individuals each.
ordinary <- WC_FST_Diploids_2Alleles(rbind(c(3, 2, 1), c(1, 2, 3)))
expected_ordinary <- list(
  He = 1/2, FST = 1/10, T1 = 1/36, T2 = 5/18,
  FSTNoCorr = 1/5, T1NoCorr = 1/18, T2NoCorr = 5/18,
  meanAlleleFreq = 1/2
)
stopifnot(isTRUE(all.equal(ordinary, expected_ordinary, tolerance = 1e-12)))
