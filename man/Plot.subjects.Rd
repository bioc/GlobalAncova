\name{Plot.subjects}
\alias{Plot.subjects}
\alias{Plot.subjects,matrix,formula,formula,ANY,missing,missing,missing-method}
\alias{Plot.subjects,matrix,formula,missing,ANY,missing,missing,character-method}
\alias{Plot.subjects,matrix,missing,missing,missing,ANY,ANY,missing-method}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Subjects Plot for GlobalAncova}
\description{
  Produces a plot to show the influence of the samples on the test result produced by \code{\link{GlobalAncova}}.

  There are three possible ways of using \code{GlobalAncova}.
  Also \code{Plot.subjects} can be invoked with these three alternatives.
  }

\section{Methods}{
\describe{
\item{xx = "matrix", formula.full = "formula", formula.red = "formula", model.dat = "ANY", group = "missing", covars = "missing", test.terms = "missing"}{In this method, besides the expression matrix \code{xx}, model formulas for the full
and reduced model and a data frame \code{model.dat} specifying corresponding model
terms have to be given. Terms that are included in the full but not in the reduced
model are those whose association with differential expression will be tested.
The arguments \code{group}, \code{covars} and \code{test.terms} are '"missing"'
since they are not needed for this method.}

\item{xx = "matrix", formula.full = "formula", formula.red = "missing", model.dat = "ANY", group = "missing", covars = "missing", test.terms = "character"}{In this method, besides the expression matrix \code{xx}, a model formula for the full
model and a data frame \code{model.dat} specifying corresponding model
terms are required. The character argument \code{test.terms} names the terms of interest
whose association with differential expression will be tested.
The arguments \code{formula.red}, \code{group} and \code{covars} are '"missing"'
since they are not needed for this method.}

\item{xx = "matrix", formula.full = "missing", formula.red = "missing", model.dat = "missing", group = "ANY", covars = "ANY", test.terms = "missing"}{Besides the expression matrix \code{xx} a clinical variable \code{group} is
required. Covariate adjustment is possible via the argument \code{covars} but
more complex models have to be specified with the methods described above.
This method emulates the function call in the first version of the package.
The arguments \code{formula.full}, \code{formula.red}, \code{model.dat} and
\code{test.terms} are '"missing"' since they are not needed for this method.}
}}

\usage{
\S4method{Plot.subjects}{matrix,formula,formula,ANY,missing,missing,missing}(xx, formula.full, formula.red, model.dat, group,covars = NULL, test.terms,test.genes, Colorgroup = NULL, sort = FALSE, legendpos = "topright", returnValues = FALSE, bar.names, ...)

\S4method{Plot.subjects}{matrix,formula,missing,ANY,missing,missing,character}(xx, formula.full,formula.red, model.dat, group,covars = NULL, test.terms,test.genes, Colorgroup = NULL, sort = FALSE, legendpos = "topright", returnValues = FALSE, bar.names, ...)

\S4method{Plot.subjects}{matrix,missing,missing,missing,ANY,ANY,missing}(xx, formula.full, formula.red, model.dat, group,covars = NULL, test.terms,test.genes, Colorgroup = NULL, sort = FALSE, legendpos = "topright", returnValues = FALSE, bar.names, ...)
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
  \item{test.genes}{Vector of gene names or gene indices specifying the gene set. If missing, the plot refers to all genes in \code{xx}.}
  \item{Colorgroup}{Character variable giving the group that specifies coloring.
        If the function is called using the argument \code{group} then this variable is assumed to
        be relevant for coloring.}
  \item{sort}{Should the samples be ordered by \code{colorgroup}?}
  \item{legendpos}{Position of the legend (a single
     keyword from the list '"bottomright"', '"bottom"', '"bottomleft"',
     '"left"', '"topleft"', '"top"', '"topright"', '"right"' and
     '"center"').}
  \item{returnValues}{Shall bar heights (subject-wise reduction in sum of squares) be returned?}
  \item{bar.names}{Vector of bar labels. If missing, column names of \code{xx} are taken.}
  \item{\dots}{Graphical parameters for specifying colors, titles etc.}
  }
\author{
  Reinhard Meister \email{meister@tfh-berlin.de}\cr
  Ulrich Mansmann \email{mansmann@ibe.med.uni-muenchen.de}\cr
  Manuela Hummel \email{m.hummel@dkfz.de}
  }
\note{This work was supported by the NGFN project 01 GR 0459, BMBF, Germany.}
\seealso{\code{\link{GlobalAncova}}, \code{\link{Plot.genes}}, \code{\link{Plot.sequential}}}
\examples{
data(vantVeer)
data(phenodata)
data(pathways)

Plot.subjects(xx = vantVeer, formula.full = ~metastases + ERstatus, formula.red = ~ERstatus, model.dat = phenodata, test.genes = pathways[[3]], colorgroup = "metastases")
Plot.subjects(xx = vantVeer, formula.full = ~metastases + ERstatus, test.terms = "metastases", model.dat = phenodata, test.genes = pathways[[3]], colorgroup = "metastases")
Plot.subjects(xx = vantVeer, group = phenodata$metastases, covars = phenodata$ERstatus, test.genes = pathways[[3]])
}
\keyword{ hplot }
