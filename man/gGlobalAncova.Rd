\name{gGlobalAncova}
\alias{gGlobalAncova}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Generalized GlobalAncova group test}
\description{Computation of a permutation test for the association between sets of variables (e.g. genes, SNPs, ...) and clinical entities. The variables can be continuous, binary, categorical, ordinal, or of mixed types.
  The test is carried out by comparing the deviances of the full generalized linear model and the reduced model lacking the design parameters of interest. The variable-wise models are summarized to a global test statistic for the complete set.
}


\usage{
gGlobalAncova(data, formula.full, formula.red=~1, model.dat, Sets, sumstat=sum, perm=10000)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{\code{data.frame} of variables to be tested in sets (columns=variables); (multi-) categorical variables should be factors, ordinal variables should be ordered factors}
  \item{formula.full}{model formula for the full model}
  \item{formula.red}{model formula for the reduced model (that does not contain the terms of interest)}
  \item{model.dat}{\code{data.frame} of regressors, containing variables specified in \code{formula.full} and \code{formula.red}}
  \item{Sets}{vector of names or indices of variables or list of those, defining sets of variables}
  \item{sumstat}{function for summarizing univariate test statistics; default is \code{sum}}
  \item{perm}{number of permutations}
}
\value{A data.frame with test statistic and p-value for each tested set.}
\author{Reinhard Meister \email{meister@beuth-hochschule.de}\cr
  Manuela Hummel \email{m.hummel@dkfz.de}}
  
\note{The test is fast for categorical data and categorical design variable. For other types of variables and more complex designs it is rather slow.

This work was supported by BMBF grant 01ZX1309B, Germany.}
\seealso{\code{\link{GlobalAncova}}}
\examples{
data(bindata)
gGlobalAncova(bindata[,-1], formula.full = ~group, model.dat = bindata, perm = 1000)
}

\keyword{ models }% at least one, from doc/KEYWORDS

