# Regression tests for issue #20 ---------------------------------------------
# Loci observed in fewer than two populations cannot have an FST estimate.

source("R/Fst Diploids.R")

populations <- rep(c("A", "B", "C"), each = 5L)
one_population <- c(rep(0, 5L), rep(9, 10L))
all_missing <- rep(9, 15L)

one_result <- getFSTs_diploids(populations, one_population)
all_result <- getFSTs_diploids(populations, all_missing)

stopifnot(is.na(one_result$FST), is.na(one_result$FSTNoCorr))
stopifnot(is.na(all_result$FST), is.na(all_result$FSTNoCorr))
