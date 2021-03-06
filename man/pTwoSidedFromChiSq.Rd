% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/OutFLANK.R
\name{pTwoSidedFromChiSq}
\alias{pTwoSidedFromChiSq}
\title{P-values for test of neutrality}
\usage{
pTwoSidedFromChiSq(x, df)
}
\arguments{
\item{DataList}{A data frame with a row for each locus, that includes at least a column for $FSTNoCorr and $He.}

\item{Fstbar}{Mean Fst (without sample size correction) as inferred from neutral loci or OutFLank}

\item{dfInferred}{The inferred degrees of freedom of the chi-square distribution describing neutral Fst values.}

\item{Hmin}{Minimum heterozygosity (H) to exclude low H alleles}
}
\value{
Returns a data frame with the original data, and two new columns appended:
\itemize{
\item $pvalues the p-value for a locus, with extremely large values of FST near 0
\item $pvaluesRightTail the one-sided (right tail) p-value for a locus
}

#@export
}
\description{
Calculates P-values for test of neutrality for a list of loci, using input of an inferred degrees of freedom for the chi-square and mean Neutral FST
}
