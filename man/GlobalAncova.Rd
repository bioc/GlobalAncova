\name{GlobalAncova}
\alias{GlobalAncova}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Global test for differential gene expression}
\description{Computation of a F-test for the association between expression values
  and clinical entities. In many cases a two way layout
  with gene and a dichotomous group as factors will be considered.
  However, adjustment for other covariates and the analysis of arbitrary clinical variables,
  interactions, gene co-expression, time series data and so on is also possible.
  The test is carried out by comparison of corresponding linear models
  via the extra sum of squares principle.
  Corresponding p-values, permutation p-values and/or asymptotic p-values are given.

  There are three possible ways of using \code{GlobalAncova}. The general way is to define
formulas for the full and reduced model, respectively, where the formula terms correspond
to variables in \code{model.dat}.
An alternative is to specify the full model and the name of the model terms
that shall be tested regarding differential expression.
In order to make this layout compatible with the
function call in the first version of the package there is also a method where simply
a group variable (and possibly covariate information) has to be given. This is maybe
the easiest usage in cases where no 'special' effects like e.g. interactions are
of interest.
}

\section{Methods}{
\describe{
\item{xx = "matrix", formula.full = "formula", formula.red = "formula", model.dat = "ANY", group = "missing", 
      covars = "missing", test.terms = "missing"}{In this method, besides the expression matrix \code{xx}, model formulas for the full
and reduced model and a data frame \code{model.dat} specifying corresponding model
terms have to be given. Terms that are included in the full but not in the reduced
model are those whose association with differential expression will be tested.
The arguments \code{group}, \code{covars} and \code{test.terms} are '"missing"'
since they are not needed for this method.}

\item{xx = "matrix", formula.full = "formula", formula.red = "missing", model.dat = "ANY", group = "missing", 
      covars = "missing", test.terms = "character"}{In this method, besides the expression matrix \code{xx}, a model formula for the full
model and a data frame \code{model.dat} specifying corresponding model
terms are required. The character argument \code{test.terms} names the terms of interest
whose association with differential expression will be tested. The basic idea
behind this method is that one can select single terms, possibly from the list of
terms provided by previous \code{GlobalAncova} output, and test them without having
to specify each time a model formula for the reduced model.
The arguments \code{formula.red}, \code{group} and \code{covars} are '"missing"'
since they are not needed for this method.}


\item{xx = "matrix", formula.full = "missing", formula.red = "missing", model.dat = "missing", 
      group = "ANY", covars = "ANY", test.terms = "missing"}{Besides the expression matrix \code{xx} a clinical variable \code{group} is
required. Covariate adjustment is possible via the argument \code{covars} but
more complex models have to be specified with the methods described above.
This method emulates the function call in the first version of the package.
The arguments \code{formula.full}, \code{formula.red}, \code{model.dat} and
\code{test.terms} are '"missing"' since they are not needed for this method.}
}}


\usage{
\S4method{GlobalAncova}{matrix,formula,formula,ANY,missing,missing,missing}(xx, formula.full, formula.red, model.dat, 
          test.genes = NULL, method = c("permutation","approx","both","Fstat"), perm = 10000, max.group.size = 2500, eps = 1e-16, acc = 50)

\S4method{GlobalAncova}{matrix,formula,missing,ANY,missing,missing,character}(xx, formula.full, test.terms, model.dat, 
          test.genes = NULL, method = c("permutation","approx","both","Fstat"), perm = 10000, max.group.size = 2500, eps = 1e-16, acc = 50)

\S4method{GlobalAncova}{matrix,missing,missing,missing,ANY,ANY,missing}(xx, group, covars = NULL,   
          test.genes = NULL, method = c("permutation","approx","both","Fstat"), perm = 10000, max.group.size = 2500, eps = 1e-16, acc = 50)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{xx}{Matrix of gene expression data, where columns correspond to samples
             and rows to genes. The data should be properly normalized beforehand
             (and log- or otherwise transformed). Missing values are not allowed.
             Gene and sample names can be included as the row and column
             names of \code{xx}.}
  \item{formula.full}{Model formula for the full model.}
  \item{formula.red}{Model formula for the reduced model (that does not contain the terms of interest.)}
  \item{model.dat}{Data frame that contains all the variable information for each sample.}
  \item{group}{Vector with the group membership information.}
  \item{covars}{Vector or matrix which contains the covariate information for each sample.}
  \item{test.terms}{Character vector that contains names of the terms of interest.}
  \item{test.genes}{Vector of gene names or a list where each element is a vector of gene names.}
  \item{method}{p-values can be calculated permutation-based (\code{"permutation"}) or by means of an approximation for a mixture of chi-square 
                distributions (\code{"approx"}).
                Both p-values are provided when specifying \code{method = "both"}. With option \code{"Fstat"} only the global F-statistics are
                returned without p-values or further information.}
  \item{perm}{Number of permutations to be used for the permutation approach. The default is 10,000.}
  \item{max.group.size}{Maximum size of a gene set for which the asymptotic p-value is calculated. 
                        For bigger gene sets the permutation approach is used.}
  \item{eps}{Resolution of the asymptotic p-value.}
  \item{acc}{Accuracy parameter needed for the approximation. Higher values indicate higher accuracy.}
}
\value{
  If \code{test.genes = NULL} a list with components
  \item{effect}{Name(s) of the tested effect(s)}
  \item{ANOVA}{ANOVA table}
  \item{test.result}{F-value, theoretical p-value, permutation-based and/or asymptotic p-value}
  \item{terms}{Names of all model terms}
  
  If a collection of gene sets is provided in \code{test.genes} a matrix is returned whose columns show the number of genes, value of the 
  F-statistic, theoretical p-value, permutation-based and/or asymptotic p-value for each of the gene sets.
}
\references{Mansmann, U. and Meister, R., 2005, Testing differential gene expression
in functional groups, \emph{Methods Inf Med} 44 (3).}
\author{Reinhard Meister \email{meister@tfh-berlin.de}\cr
  Ulrich Mansmann \email{mansmann@ibe.med.uni-muenchen.de}\cr
  Manuela Hummel \email{hummel@ibe.med.uni-muenchen.de}}
\note{This work was supported by the NGFN project 01 GR 0459, BMBF, Germany.}
\seealso{\code{\link{Plot.genes}}, \code{\link{Plot.subjects}}, \code{\link{GlobalAncova.closed}}}
\examples{
data(vantVeer)
data(phenodata)
data(pathways)

set.seed(123)
GlobalAncova(xx = vantVeer, formula.full = ~metastases + ERstatus, formula.red = ~ERstatus, model.dat = phenodata, test.genes=pathways[1], method="both", perm = 100)
set.seed(123)
GlobalAncova(xx = vantVeer, formula.full = ~metastases + ERstatus, test.terms = "metastases", model.dat = phenodata, test.genes=pathways[1], method="both", perm = 100)
set.seed(123)
GlobalAncova(xx = vantVeer, group = phenodata$metastases, covars = phenodata$ERstatus, test.genes=pathways[1], method="both", perm = 100)
}
\keyword{ models }% at least one, from doc/KEYWORDS

