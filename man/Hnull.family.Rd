\name{Hnull.family}
\alias{Hnull.family}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Family of null hypotheses closed under intersection}
\description{
  Creates new groups of genes by building unions of the given pathways.
  Each group of genes corresponds to the null hypothesis that those genes are
  not differentially expressed. All groups together build a family
  of null hypotheses that is closed under intersection.
}
\usage{
Hnull.family(test.genes)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{test.genes}{A list of pathways that shall be tested, each containing gene names.}
}
\value{A list of groups of genes each corresponding to one hypothesis.}
\references{
  Marcus, R., Peritz, E. and Gabriel, K. R., 1976, On closed testing procedures
  with special reference to ordered analysis of variance, \emph{Biometrika} 63 (3): 655--660.
  }
\author{
  Reinhard Meister \email{meister@tfh-berlin.de}\cr
  Ulrich Mansmann \email{mansmann@ibe.med.uni-muenchen.de}\cr
  Manuela Hummel \email{hummel@ibe.med.uni-muenchen.de}
  }
\note{For internal use.}
\seealso{\code{\link{GlobalAncova}}, \code{\link{GlobalAncova.closed}}, \code{\link{Plot.genes}}, \code{\link{Plot.subjects}}}
\keyword{internal}% at least one, from doc/KEYWORDS

