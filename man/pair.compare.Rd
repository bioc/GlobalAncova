\name{pair.compare}
\alias{pair.compare}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Pairwise comparisons of factor levels within GlobalAncova}
\description{Pairwise comparisons of gene expression in different levels of a factor by GlobalAncova tests.
  The method uses the reduction in residual sum of squares obtained when two respective factor levels are set to the same level.
  Holm-adjusted permutation-based p-values are given.
}


\usage{
pair.compare(xx, formula, group, model.dat = NULL, test.genes = NULL, perm = 10000)}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{xx}{Matrix of gene expression data, where columns correspond to samples
             and rows to genes. The data should be properly normalized beforehand
             (and log- or otherwise transformed). Missing values are not allowed.
             Gene and sample names can be included as the row and column
             names of \code{xx}.}
  \item{formula}{Model formula for the linear model.}
  \item{group}{Factor for which pairwise comparisons shall be calculated.}
  \item{model.dat}{Data frame that contains all the variable information for each sample.}
  \item{test.genes}{Vector of gene names or a list where each element is a vector of gene names.}
  \item{perm}{Number of permutations to be used for the permutation approach. The default is 10,000.}
}
\value{
  An ANOVA table, or list of ANOVA tables for each gene set, for the pairwise comparisons.
}
%\references{!!!}
\author{Ramona Scheufele \email{ramona.scheufele@charite.de} \cr
  Reinhard Meister \email{meister@tfh-berlin.de}\cr
  Manuela Hummel \email{hummel@ibe.med.uni-muenchen.de} \cr
  Urlich Mansmann \email{mansmann@ibe.med.uni-muenchen.de}}
\note{This work was supported by the NGFN project 01 GR 0459, BMBF, Germany.}
\seealso{\code{\link{GlobalAncova}}, \code{\link{GlobalAncova.decomp}}}
\examples{
data(vantVeer)
data(phenodata)
data(pathways)

pair.compare(xx = vantVeer, formula = ~ grade, group = "grade", model.dat = phenodata, test.genes = pathways[1:3], perm = 100)
}
\keyword{ models }% at least one, from doc/KEYWORDS

