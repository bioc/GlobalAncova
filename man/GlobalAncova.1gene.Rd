\name{GlobalAncova.1gene}
\alias{GlobalAncova.1gene}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Test for differential gene expression for just one gene}
\description{
 A squared t-statistic which is equivalent to a F-statistic is computed.
 Corresponding p-values and permutation p-values are provided.
}
\usage{
GlobalAncova.1gene(xx, group, perm = 10000)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{xx}{A vector of gene expression data. The data should be properly normalized beforehand
             (and log- or otherwise transformed). Missing values are not allowed.}
  \item{group}{A vector with the group membership information. In the given version group must be coded as \code{0-1}.}
  \item{perm}{The number of permutations to be used. The default is 10,000.}
}
\value{
  A vector giving the F value, permutation p--value and theoretical p--value.
}
\author{
  Reinhard Meister \email{meister@tfh-berlin.de}\cr
  Ulrich Mansmann \email{mansmann@ibe.med.uni-muenchen.de}\cr
  Manuela Hummel \email{hummel@ibe.med.uni-muenchen.de}
  }
\note{For internal use.}
\seealso{\code{\link{GlobalAncova}}, \code{\link{GlobalAncova.closed}}, \code{\link{Plot.genes}}, \code{\link{Plot.subjects}}}
\keyword{internal}% at least one, from doc/KEYWORDS

