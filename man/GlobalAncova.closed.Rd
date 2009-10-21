\name{GlobalAncova.closed}
\alias{GlobalAncova.closed}
\alias{GlobalAncova.closed,matrix,list,formula,formula,ANY,missing,missing,missing-method}
\alias{GlobalAncova.closed,matrix,list,formula,missing,ANY,missing,missing,character-method}
\alias{GlobalAncova.closed,matrix,list,missing,missing,missing,ANY,ANY,missing-method}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Closed testing procedure for testing several groups of genes using GlobalAncova }
\description{
  Computation of a closed testing procedure for several groups of genes, e.g. pathways, as an alternative of
  correcting for multiple testing. Starting from the pathways of interest a family of
  null hypotheses is created that is closed under intersection. Each null hypothesis
  can be rejected at a given level if it is rejected along with all hypotheses included in it.

  There are three possible ways of using \code{GlobalAncova}.
  Also \code{GlobalAncova.closed} can be invoked with these three alternatives.
}

\section{Methods}{
\describe{
\item{xx = "matrix", test.genes="list", formula.full = "formula", formula.red = "formula", model.dat = "ANY", group = "missing", covars = "missing", test.terms = "missing"}{In
this method, besides the expression matrix \code{xx} and the list of gene groups
\code{test.genes}, model formulas for the full
and reduced model and a data frame \code{model.dat} specifying corresponding model
terms have to be given. Terms that are included in the full but not in the reduced
model are those whose association with differential expression will be tested.
The arguments \code{group}, \code{covars} and \code{test.terms} are '"missing"'
since they are not needed for this method.}

\item{xx = "matrix", test.genes="list", formula.full = "formula", formula.red = "missing", model.dat = "ANY", group = "missing", covars = "missing", test.terms = "character"}{In
 this method, besides the expression matrix \code{xx} and the list of gene groups
\code{test.genes}, a model formula for the full
model and a data frame \code{model.dat} specifying corresponding model
terms are required. The character argument \code{test.terms} names the terms of interest
whose association with differential expression will be tested.
The arguments \code{formula.red}, \code{group} and \code{covars} are '"missing"'
since they are not needed for this method.}

\item{xx = "matrix", test.genes="list", formula.full = "missing", formula.red = "missing", model.dat = "missing", group = "ANY", covars = "ANY", test.terms = "missing"}{Besides
 the expression matrix \code{xx} and the list of gene groups
\code{test.genes} a clinical variable \code{group} is
required. Covariate adjustment is possible via the argument \code{covars} but
more complex models have to be specified with the methods described above.
This method emulates the function call in the first version of the package.
The arguments \code{formula.full}, \code{formula.red}, \code{model.dat} and
\code{test.terms} are '"missing"' since they are not needed for this method.}
}}


\usage{
\S4method{GlobalAncova.closed}{matrix,list,formula,formula,ANY,missing,missing,missing}(xx, test.genes, 
          formula.full, formula.red, model.dat, previous.test, level, method = c("permutation","approx"), perm = 10000, 
          max.group.size = 2500, eps = 1e-16, acc = 50)

\S4method{GlobalAncova.closed}{matrix,list,formula,missing,ANY,missing,missing,character}(xx, test.genes, 
          formula.full, test.terms, model.dat, previous.test, level, method = c("permutation","approx"), perm = 10000, 
          max.group.size = 2500, eps = 1e-16, acc = 50)

\S4method{GlobalAncova.closed}{matrix,list,missing,missing,missing,ANY,ANY,missing}(xx, test.genes, 
          group, covars = NULL, previous.test, level, method = c("permutation","approx"), perm = 10000, 
          max.group.size = 2500, eps = 1e-16, acc = 50)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{xx}{Matrix of gene expression data, where columns correspond to samples
             and rows to genes. The data should be properly normalized beforehand
             (and log- or otherwise transformed). Missing values are not allowed.
             Gene and sample names can be included as the row and column
             names of \code{xx}.}
  \item{test.genes}{A list of named pathways that shall be tested, each containing vectors of gene names.}
  \item{previous.test}{The output of a call to \code{GlobalAncova} with specified option
                           \code{test.genes} according to the pathways of interest (optional).}
  \item{level}{The global level of significance of the testing procedure.}
  \item{formula.full}{Model formula for the full model.}
  \item{formula.red}{Model formula for the reduced model (that does not contain the terms of interest).}
  \item{model.dat}{Data frame that contains all the variable information for each sample.}
  \item{group}{Vector with the group membership information.}
  \item{covars}{Vector or matrix which contains the covariate information for each sample.}
  \item{test.terms}{Character vector that contains names of the terms of interest.}
  \item{method}{Raw p-values can be calculated permutation-based (\code{"permutation"}) or by means of an approximation (\code{"approx"}).}
  \item{perm}{Number of permutations to be used for the permutation approach. The default is 10,000.}
  \item{max.group.size}{Maximum size of a gene set for which the asymptotic p-value is calculated. 
                        For bigger gene sets the permutation approach is used.}
  \item{eps}{Resolution of the asymptotic p-value.}
  \item{acc}{Accuracy parameter needed for the approximation. Higher values indicate higher accuracy.}
}
\value{
  A list with components
  \item{new.data}{Family of null hypotheses (vectors of genes to be tested simultaneously with \code{GlobalAncova}).}
  \item{test.results}{Test results for each pathway of interest and all hypotheses included in it.}
  \item{significant}{Names of the significant pathways.}
  \item{not.significant}{Names of the non significant pathways.}
}
\references{Marcus, R., Peritz, E. and Gabriel, K.R., 1976, On closed testing procedures
  with special reference to ordered analysis of variance, \emph{Biometrika} 63 (3): 655--660.}
\author{  Reinhard Meister \email{meister@tfh-berlin.de}\cr
  Ulrich Mansmann \email{mansmann@ibe.med.uni-muenchen.de}\cr
  Manuela Hummel \email{hummel@ibe.med.uni-muenchen.de}
}
\note{This work was supported by the NGFN project 01 GR 0459, BMBF, Germany.}
\seealso{\code{\link{GlobalAncova}}, \code{\link{Plot.genes}}, \code{\link{Plot.subjects}}}
\keyword{ models }% at least one, from doc/KEYWORDS

