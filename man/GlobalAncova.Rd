\name{GlobalAncova}
\alias{GlobalAncova}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Global test for differential gene expression}
\description{
 Computation of the sum of squares decomposition of gene-expression values according to a two way layout
 with gene and group as factors. Gene specific covariate effects are allowed. A MC approximation to a permutation-test
 for group main effect and gene: group interaction is provided.
 If just one gene should be tested a squared t-statistic which is equivalent to a F-statistic is computed.
 Corresponding p-values and permutation p-values are provided.
 }
\usage{
GlobalAncova(xx, group, covars = NULL, perm = 10000, test.genes = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{xx}{A matrix of gene expression data, where columns correspond to samples
             and rows to genes. The data should be properly normalized beforehand
             (and log- or otherwise transformed). Missing values are not allowed.
             Gene and sample names can be included as the row and column
             names of \code{xx}.}
  \item{group}{A vector with the group membership information. In the given version group must be coded as \code{0-1}.}
  \item{covars}{A vector or matrix which contains for each sample the covariate information.}
  \item{perm}{The number of permutations to be used. The default is 10,000.}
  \item{test.genes}{Vector of genes that shall be tested or list of pathways, each containing gene names.}
  }
\value{
  An ANOVA table
}
\references{
Mansmann, U. and Meister, R., 2005, Testing differential gene expression
in functional groups, \emph{Methods Inf Med} 44 (3).
}
\author{
  Reinhard Meister \email{meister@tfh-berlin.de}\cr
  Ulrich Mansmann \email{mansmann@ibe.med.uni-muenchen.de}
  }
\note{This work was supported by the NGFN project 01 GR 0459, BMBF, Germany.}
\seealso{\code{\link{Plot.genes}}, \code{\link{Plot.subjects}}, \code{\link{GlobalAncova.closed}}}
\examples{
set.seed(123)
data(p53.signalling)
data(cov.info)
data(group.info)
table.1   <- GlobalAncova(xx=p53.signalling, group=group.info, covars=NULL, perm=10000, test.genes = NULL)
table.2   <- GlobalAncova(xx=p53.signalling, group=group.info, covars=cov.info, perm=10000, test.genes = NULL)
table.sex <- GlobalAncova(xx=p53.signalling, group=cov.info[,1], covars=NULL, perm=10000, test.genes = NULL)
table.loc <- GlobalAncova(xx=p53.signalling, group=cov.info[,2], covars=NULL, perm=10000, test.genes = NULL)
}
\keyword{ models }
