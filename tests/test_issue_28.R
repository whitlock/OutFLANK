# Regression test for issue #28 ----------------------------------------------
# `meanAlleleFreq` should use the documented focal (genotype-2) allele.

source("R/Fst Diploids.R")

counts <- matrix(c(701, 1595, 2704,
                   4524, 180, 296),
                 nrow = 2L, byrow = TRUE,
                 dimnames = list(NULL, c("0", "1", "2")))

result <- WC_FST_Diploids_2Alleles(counts)
expected <- ((2704 + 1595 / 2) + (296 + 180 / 2)) / 10000
stopifnot(abs(result$meanAlleleFreq - expected) < 1e-12)
