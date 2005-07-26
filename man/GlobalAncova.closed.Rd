\name{GlobalAncova.closed}
\alias{GlobalAncova.closed}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Closed testing procedure for testing several pathways using GlobalAncova  }
\description{
  Computation of a closed testing procedure for several pathways as an alternative of
  correcting for multiple testing. Starting from the pathways of interest a family of
  null hypotheses is created that is closed under intersection. Each null hypothesis
  can be rejected at a given level if it is rejected along with all hypotheses included in it.
  For details see Marcus, R., Peritz, E. and Gabriel, K. R., 1976, On closed testing procedures
  with special reference to ordered analysis of variance, \emph{Biometrika} 63 (3): 655--660.
  }
\usage{
GlobalAncova.closed(xx, group, covars = NULL, test.genes, previous.test = NULL, level = 0.05, perm = 10000)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{xx}{A matrix of gene expression data, where columns correspond to samples
              and rows to genes. The data should be properly normalized beforehand (and log- or
              otherwise transformed). Missing values are not allowed. Gene and sample names can
              be included as the row and column names of xx.}
  \item{group}{A vector with the group membership information. In the given version group must be coded as \code{0-1}.}
  \item{covars}{A vector or matrix which contains for each sample the covariate information.}
  \item{test.genes}{A list of pathways that shall be tested, each containing gene names.}
  \item{previous.test}{The output of a call to 'GlobalAncova' with specified option
                           'test.genes' according to the pathways of interest.}
  \item{level}{The level of significance of the testing procedure.}
  \item{perm}{The number of permutations to be used. The default is 10,000.}
  }
\value{
  A list with components
  \item{new.data}{Family of null hypotheses (vectors of genes to be tested simultaneously with 'GlobalAncova').}
  \item{test.results}{Test results for each pathway of interest and all hypotheses included in it.}
  \item{significant}{Names of the significant pathways.}
  \item{not.significant}{Names of the not significant pathways.}
  }
\author{
  Reinhard Meister \email{meister@tfh-berlin.de}\cr
  Ulrich Mansmann \email{mansmann@ibe.med.uni-muenchen.de}\cr
  Manuela Hummel \email{hummel@ibe.med.uni-muenchen.de}
  }
\note{This work was supported by the NGFN project 01 GR 0459, BMBF, Germany.}
\seealso{\code{\link{GlobalAncova}}, \code{\link{Plot.genes}}, \code{\link{Plot.subjects}}}
\keyword{ models }% at least one, from doc/KEYWORDS

