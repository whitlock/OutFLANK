library(OutFLANK)

# Individuals are rows; loci are columns. Neither population has heterozygotes.
x <- matrix(c(0, 0, 2, 0, 2, 2), ncol = 1)
pops <- rep(c("A", "B"), each = 3)
result <- MakeDiploidFSTMat(x, "locus1", pops)
stopifnot(nrow(result) == 1L, is.finite(result$FSTNoCorr))
stopifnot(identical(result,
  MakeDiploidFSTMat(as.data.frame(x), "locus1", pops)))

# Other valid subsets, including the documented missing-data code.
for (values in list(c(0, 0, 1, 0, 1, 1), c(1, 1, 2, 1, 2, 2),
                    c(0, 0, 2, 0, 2, 9), c(0, 1, 2, 1, 2, 2))) {
  ans <- MakeDiploidFSTMat(matrix(values, ncol = 1), "locus1", pops)
  stopifnot(nrow(ans) == 1L, is.finite(ans$FSTNoCorr))
}

expect_error_text <- function(expr, text) {
  err <- tryCatch({ force(expr); NULL }, error = identity)
  stopifnot(inherits(err, "error"),
            grepl(text, conditionMessage(err), fixed = TRUE))
}

for (bad in c(NA_real_, NaN, Inf, -1, 3, 0.5)) {
  invalid <- x
  invalid[1, 1] <- bad
  expect_error_text(MakeDiploidFSTMat(invalid, "locus1", pops),
                    "0, 1, 2 or 9 (missing)")
}
expect_error_text(MakeDiploidFSTMat(matrix("0", 6, 1), "locus1", pops),
                  "numeric matrix")
expect_error_text(MakeDiploidFSTMat(as.vector(x), "locus1", pops),
                  "numeric matrix")
expect_error_text(MakeDiploidFSTMat(matrix(numeric(), 0, 0), NULL, NULL),
                  "nonempty")
expect_error_text(MakeDiploidFSTMat(x, "locus1", pops[-1]),
                  "one entry per row (individual)")
expect_error_text(MakeDiploidFSTMat(x, c("locus1", "locus2"), pops),
                  "one entry per column (locus)")
