\name{Plot.genes}
\alias{Plot.genes}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Genes Plot for Global Ancova}
\description{
  Produces a plot to show the influence of individual genes on the test result produced by 'GlobalAncova'.
}
\usage{
Plot.genes(xx, group, covars = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{xx}{A matrix of gene expression data, where columns correspond to samples
              and rows to genes. The data should be properly normalized beforehand
             (and log- or otherwise transformed). Missing values are not allowed. Gene
             and sample names can be included as the row and column names of xx.}
  \item{group}{A vector with the group membership information. In the given version group must be coded as \code{0-1}.}
  \item{covars}{A vector or matrix which contains for each sample the covariate information.}
}
\author{
  Reinhard Meister \email{meister@tfh-berlin.de}\cr
  Ulrich Mansmann \email{mansmann@ibe.med.uni-muenchen.de}
  }
\note{This work was supported by the NGFN project 01 GR 0459, BMBF, Germany.}
\seealso{\code{\link{GlobalAncova}}, \code{\link{Plot.subjects}}, \code{\link{GlobalAncova.closed}}}
\examples{
data(p53.signalling)
data(group.info)

Plot.genes(xx=p53.signalling, group=group.info, covars=NULL)
}
\keyword{ hplot }
