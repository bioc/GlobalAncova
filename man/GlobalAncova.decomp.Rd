\name{GlobalAncova.decomp}
\alias{GlobalAncova.decomp}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{GlobalAncova with sequential and type III sum of squares decomposition and adjustment for global covariates}
\description{Computation of a F-test for the association between expression values
  and clinical entities. 
  The test is carried out by comparison of corresponding linear models
  via the extra sum of squares principle.
  In models with various influencing factors extra sums of squares can be treated with sequential and type III
  decomposition. Adjustment for global covariates, e.g. gene expression values in normal tissue as compared
  to tumour tissue, can be applied.  
  Given theoretical p-values may not be appropriate due to correlations and non-normality. The functions
  are hence seen more as a descriptive tool.
}


\usage{
GlobalAncova.decomp(xx, formula, model.dat = NULL, method = c("sequential", "type3", "all"),  test.genes = NULL, genewise = FALSE, zz = NULL, zz.per.gene = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{xx}{Matrix of gene expression data, where columns correspond to samples
             and rows to genes. The data should be properly normalized beforehand
             (and log- or otherwise transformed). Missing values are not allowed.
             Gene and sample names can be included as the row and column
             names of \code{xx}.}
  \item{formula}{Model formula for the linear model.}
  \item{model.dat}{Data frame that contains all the variable information for each sample.}
  \item{method}{Whether sequential or type III decomposition or both should be calculated.}
  \item{test.genes}{Vector of gene names or a list where each element is a vector of gene names.}
  \item{genewise}{Shall the sequential decomposition be displayed for each single gene in a (small) gene set?}
  \item{zz}{Global covariate, i.e. matrix of same dimensions as \code{xx}.}
  \item{zz.per.gene}{If set to \code{TRUE} the adjustment for the global covariate is applied on a gene-wise basis.}
}
\value{
%  If \code{test.genes = NULL} a list with components
%  \item{sequential}{ANOVA table for the sequential decomposition.}
%  \item{typeIII}{ANOVA table for the type III decomposition.}
%  \item{adjustment}{Sum of squares and degrees of freedom of the adjustment for a global covariate.}
  
%  If a collection of gene sets is provided in \code{test.genes} 
  
%  If parameter \code{genewise} is set to \code{TRUE} a list with the names of tested terms, sums of squares, degrees of freedom,
%  mean squares, F-statistics and theoretical p-values for each gene in the gene set is returned.
  Depending on parameters \code{test.genes}, \code{method} and \code{genewise} ANOVA tables, or lists of ANOVA tables for each 
  decomposition and/or gene set, or lists with components of ANOVA tables for each gene are returned.
}
%\references{!!!}
\author{Ramona Scheufele \email{ramona.scheufele@charite.de} \cr
  Reinhard Meister \email{meister@tfh-berlin.de}\cr
  Manuela Hummel \email{hummel@ibe.med.uni-muenchen.de} \cr
  Urlich Mansmann \email{mansmann@ibe.med.uni-muenchen.de}}
\note{This work was supported by the NGFN project 01 GR 0459, BMBF, Germany.}
\seealso{\code{\link{Plot.sequential}}, \code{\link{pair.compare}}, \code{\link{GlobalAncova}}}
\examples{
data(vantVeer)
data(phenodata)
data(pathways)

# sequential or type III decomposition
GlobalAncova.decomp(xx = vantVeer, formula = ~ grade + metastases + ERstatus, model.dat = phenodata, method = "sequential", test.genes = pathways[1:3])
GlobalAncova.decomp(xx = vantVeer, formula = ~ grade + metastases + ERstatus, model.dat = phenodata, method = "type3", test.genes = pathways[1:3]) 

# adjustment for global covariate
data(colon.tumour)
data(colon.normal)
data(colon.pheno)
GlobalAncova.decomp(xx = colon.tumour, formula = ~ UICC.stage + sex + location, model.dat = colon.pheno, method = "all", zz = colon.normal)
}
\keyword{ models }% at least one, from doc/KEYWORDS

