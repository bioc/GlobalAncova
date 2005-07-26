\name{Plot.subjects}
\alias{Plot.subjects}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Subjects Plot for Global Ancova  }
\description{
  Produces a plot to show the influence of the samples on the test result produced by 'GlobalAncova'.
  }
\usage{
Plot.subjects(xx, group, covars = NULL, sort = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{xx}{A matrix of gene expression data, where columns correspond to samples
              and rows to genes. The data should be properly normalized beforehand
                  (and log- or otherwise transformed). Missing values are not allowed.
                  Gene and sample names can be included as the row and column names of xx.}
  \item{group}{A vector with the group membership information. In the given version group must be coded as \code{0-1}.}
  \item{covars}{A vector or matrix which contains for each sample the covariate information.}
  \item{sort}{Should the samples be ordered by group?}
  }
\author{
  Reinhard Meister \email{meister@tfh-berlin.de}\cr
  Ulrich Mansmann \email{mansmann@ibe.med.uni-muenchen.de}
  }
\note{This work was supported by the NGFN project 01 GR 0459, BMBF, Germany.}
\seealso{\code{\link{GlobalAncova}}, \code{\link{Plot.genes}}, \code{\link{GlobalAncova.closed}}}
\examples{
data(p53.signalling)
data(group.info)

Plot.subjects(xx=p53.signalling, group=group.info, covars=NULL, sort=FALSE)
}
\keyword{ hplot }
