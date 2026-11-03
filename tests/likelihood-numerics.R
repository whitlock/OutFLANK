library(OutFLANK)

# Internal helper: the public API is unchanged by this proposal.
negll <- getFromNamespace("negLLdfFstTrim", "OutFLANK")

# Original expression, retained here only as a regression reference.
original_negll <- function(Fst, df, Fstbar, low, high) {
  upper_gamma <- function(a, z) pgamma(z, a, lower.tail = FALSE) * gamma(a)
  (df * Fst + df * Fstbar * log(2) - df * Fstbar * log(df) -
     (df - 2) * Fstbar * log(Fst) + df * Fstbar * log(Fstbar) +
     2 * Fstbar * log(upper_gamma(df/2, df * low/(2 * Fstbar)) -
                       upper_gamma(df/2, df * high/(2 * Fstbar)))) /
    (2 * Fstbar)
}

# Agreement with the original expression in its stable numerical range.
for (df in c(2, 3, 10, 50, 100)) {
  observed <- negll(c(.06, .1, .18), df, .1, .05, .2)
  expected <- original_negll(c(.06, .1, .18), df, .1, .05, .2)
  stopifnot(isTRUE(all.equal(observed, expected, tolerance = 1e-10)))
}

# Reproduce gamma overflow, then compare the fix with an independent formula.
stopifnot(!is.finite(suppressWarnings(
  original_negll(.1, 400, .1, .05, .2))))
expected <- -dchisq(400, 400, log = TRUE) - log(400/.1) +
  log(pchisq(800, 400) - pchisq(200, 400))
observed <- negll(.1, 400, .1, .05, .2)
stopifnot(is.finite(observed),
          isTRUE(all.equal(observed, expected, tolerance = 1e-10)))

# A correctly normalized truncated density integrates to one.
# Exercise lower-tail, central and upper-tail intervals.
for (bounds in list(c(.001, .002), c(.05, .2), c(.4, .5))) {
  mass <- integrate(function(x) exp(-negll(x, 10, .1,
                                          bounds[1], bounds[2])),
                    lower = bounds[1], upper = bounds[2])$value
  stopifnot(abs(mass - 1) < 1e-7)
}
