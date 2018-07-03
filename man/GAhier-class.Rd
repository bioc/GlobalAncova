\name{GAhier class}
\Rdversion{1.1}
\docType{class}
\alias{GAhier}
\alias{GAhier-class}
\alias{Plot.hierarchy}
\alias{Plot.hierarchy,GAhier-method}
\alias{results}
\alias{results,GAhier-method}
\alias{show,GAhier-method}
\alias{sigEndnodes}
\alias{sigEndnodes,GAhier-method}

\title{Class \code{"GAhier"}}
\description{Class for storing results of hierarchical testing procedure performed by \code{\link{gGlobalAncova.hierarchical}}}
%\section{Objects from the Class}{
%Objects can be created by calls of the form \code{new("GAhier", ...)}.
%%  ~~ describe objects here ~~ 
%}

\section{Slots}{
  \describe{
    \item{\code{clustervariables}:}{Object of class \code{"list"} containing names of variables in each tested cluster}
    \item{\code{p.values}:}{Object of class \code{"list"} containing p-values for each tested cluster}
    \item{\code{alpha}:}{Object of class \code{"numeric"}; chosen global significance level; if \code{K} had been specified, this additionally contains the adjusted significance levels for the \code{K} sub-hierarchies}
    \item{\code{n.variables}:}{Object of class \code{"numeric"}; number of variables in total}
    \item{\code{permstats}:}{Object of class \code{"matrixOrNULL"}; if \code{returnPermstats} had been set to \code{TRUE}, this is a matrix containing individual statistics for all variables for all permutations, otherwise \code{NULL}}
  }
}

\section{Methods}{
  \describe{
    \item{show}{\code{signature(object = "GAhier")}: Show general information and significant end nodes}
    \item{results}{\code{signature(object = "GAhier")}: Get a \code{data.frame} with significant end nodes, number and names of variables included in each node and corresponding p-value}
    \item{sigEndnodes}{\code{signature(object = "GAhier")}: Get names of signficant end nodes}
    \item{Plot.hierarchy}{\code{signature(object = "GAhier")}: Plot hierarchy dendrogram, where significant nodes (and branches to those nodes) are highlighted}
	 }
}

\usage{
\S4method{show}{GAhier}(object)
\S4method{results}{GAhier}(object)
\S4method{sigEndnodes}{GAhier}(object, onlySingleton=FALSE)
\S4method{Plot.hierarchy}{GAhier}(object, dend, col=1:2, lwd=1:2, collab, returndend=FALSE, cex.labels=1.5, ...)
}

\arguments{
  \item{object}{object of class \code{GAhier}}
  \item{onlySingleton}{if \code{TRUE}, only names of singleton variables within the significant nodes are returned}
  \item{dend}{\code{\link{dendrogram}} object specifying the hierarchy of the variables}
  \item{col}{colors for significant and non-significant nodes and branches, respectively}
  \item{lwd}{line width for branches to non-significant and significant nodes, respectively}
  \item{collab}{vector of colors for coloring dendrogram leave labels (can be independent of significant/non-significant nodes); has to be named according to variable names; if missing, significant and non-significant variables are colored using colors defined in \code{col}}
  \item{returndend}{if \code{TRUE}, updated \code{\link{dendrogram}} object is returned (e.g. for use in further plots)}
  \item{cex.labels}{size of leave labels}
  \item{\dots}{further graphical parameters, passed to \code{\link{plot.dendrogram}}}
}

%\references{%%  ~~put references to the literature/web site here~~}
\author{Manuela Hummel \email{m.hummel@dkfz.de}}
\note{Coloring the dendrogram in \code{Plot.hierarchy} is based on functionality from the \pkg{globaltest} package}

\seealso{\code{\link{gGlobalAncova.hierarchical}}}

\examples{
showClass("GAhier")

# see examples in documentation of gGlobalAncova.hierarchical
}
\keyword{classes}
