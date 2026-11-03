# Regression test for issue #29 ----------------------------------------------
# Untested loci (negative FST or low heterozygosity) carry NA flags by design.
# They must not cause `if (any(...))` or `if (all(...))` to fail.

if (requireNamespace("qvalue", quietly = TRUE)) {
  source("R/FST functions.R")
  source("R/Likelihood functions for OutFLANK.R")
  source("R/OutFLANK.R")

  # Keep this regression test focused on NA-flag handling, not qvalue fitting.
  qvalue <- function(p, fdr.level = 0.05, pi0.method = "bootstrap") {
    list(qvalues = rep(1, length(p)), significant = rep(FALSE, length(p)))
  }

  set.seed(29)
  n <- 100L
  fst <- c(-0.02, -0.01, runif(n - 2L, 0.05, 0.25))
  x <- data.frame(
    LocusName = paste0("locus-", seq_len(n)),
    He = runif(n, 0.15, 0.45),
    FST = fst,
    T1 = pmax(fst, 0) * 0.25,
    T2 = runif(n, 0.1, 0.3),
    FSTNoCorr = fst,
    T1NoCorr = pmax(fst, 0) * 0.25,
    T2NoCorr = runif(n, 0.1, 0.3),
    meanAlleleFreq = runif(n, 0.1, 0.9)
  )

  result <- OutFLANK(x, NumberOfSamples = 2L)
  stopifnot(is.list(result), nrow(result$results) == n)
}
