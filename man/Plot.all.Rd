\name{Plot.all}
\alias{Plot.all}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Combined visualization of sequential decomposition and influence of single genes on the GlobalAncova statistic}
\description{Plot that combines \code{\link{Plot.genes}} and \code{\link{Plot.sequential}} into one graphic.
}

\usage{
Plot.all(xx, formula, model.dat = NULL, test.genes = NULL, name.geneset = "")}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{xx}{Matrix of gene expression data, where columns correspond to samples
             and rows to genes. The data should be properly normalized beforehand
             (and log- or otherwise transformed). Missing values are not allowed.
             Gene and sample names can be included as the row and column
             names of \code{xx}.}
  \item{formula}{Model formula for the linear model.}
  \item{model.dat}{Data frame that contains all the variable information for each sample.}
  \item{test.genes}{Vector of gene names or gene indices specifying a gene set.}
  \item{name.geneset}{Name of the plotted geneset.}
}
\author{Ramona Scheufele \email{ramona.scheufele@charite.de} \cr
  Reinhard Meister \email{meister@tfh-berlin.de}\cr
  Manuela Hummel \email{hummel@ibe.med.uni-muenchen.de} \cr
  Urlich Mansmann \email{mansmann@ibe.med.uni-muenchen.de}}
\note{This work was supported by the NGFN project 01 GR 0459, BMBF, Germany.}
\seealso{\code{\link{Plot.genes}}, \code{\link{Plot.sequential}}, \code{\link{GlobalAncova.decomp}}, \code{\link{GlobalAncova}}}
\examples{
data(vantVeer)
data(phenodata)
data(pathways)

Plot.all(vantVeer, formula = ~ ERstatus + metastases + grade, model.dat = phenodata, test.genes = pathways[[3]], name.geneset = "cell cycle pathway")
}
\keyword{ hplot }% at least one, from doc/KEYWORDS

