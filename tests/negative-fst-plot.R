library(OutFLANK)

# Capture plotting inputs without requiring a graphics device.
# The calculations and public wrapper remain the package's own functions.
plotter <- getFromNamespace("FstDistPlotter", "OutFLANK")
capture <- new.env(parent = environment(plotter))
capture$hist <- function(x, breaks, ...) {
  capture$histogram <- graphics::hist(x, breaks = breaks, plot = FALSE)
  invisible(capture$histogram)
}
capture$lines <- function(x, y, ...) {
  capture$curve_x <- x
  capture$curve_y <- y
  invisible(NULL)
}
environment(plotter) <- capture
capture$FstDistPlotter <- plotter
wrapper <- OutFLANKResultsPlotter
environment(wrapper) <- capture

# Public NoCorr=FALSE path: all observations, including negative estimates,
# must be counted. The old zero-based breaks produce a histogram error.
values <- c(-.02, .01, .04)
output <- list(dfInferred = 2, FSTbar = .01,
  results = data.frame(He = rep(.5, 3), FST = values,
                       OutlierFlag = rep(FALSE, 3)))
wrapper(output, NoCorr = FALSE)
stopifnot(sum(capture$histogram$counts) == length(values),
          min(capture$histogram$breaks) <= min(values),
          max(capture$histogram$breaks) >= max(values))

# The overlaid curve must use the adjusted x coordinates.
x <- capture$curve_x
expected <- length(values) *
  (pchisq((x + .0025)/.01 * 2, df = 2) -
     pchisq((x - .0025)/.01 * 2, df = 2))
stopifnot(isTRUE(all.equal(capture$curve_y, expected, tolerance = 1e-12)))

# Nonnegative inputs retain their original breaks and curve values.
values <- c(0, .01, .04)
plotter(2, values, .01)
breaks <- seq(0, ceiling(max(values) * 100)/100 + .005, by = .005)
i <- seq_along(breaks)
original_curve <- length(values) *
  (pchisq(((i - .5) * .005)/.01 * 2, df = 2) -
     pchisq(((i - 1.5) * .005)/.01 * 2, df = 2))
stopifnot(sum(capture$histogram$counts) == length(values),
          isTRUE(all.equal(capture$histogram$breaks, breaks)),
          isTRUE(all.equal(capture$curve_y, original_curve,
                           tolerance = 1e-12)))
