\name{Plot.subjects}
\alias{Plot.subjects-methods}
\alias{Plot.subjects,matrix,formula,formula,ANY,missing,missing,missing-method}
\alias{Plot.subjects,matrix,formula,missing,ANY,missing,missing,character-method}
\alias{Plot.subjects,matrix,missing,missing,missing,numeric,ANY,missing-method}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Methods for Function Plot.subjects}
\description{
  There are three possible ways of using \code{GlobalAncova}, use \code{methods ? GlobalAncova}
  for getting more information.
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

\item{xx = "matrix", formula.full = "missing", formula.red = "missing", model.dat = "missing", group = "numeric", covars = "ANY", test.terms = "missing"}{Besides the expression matrix \code{xx} a clinical variable \code{group} is
required. Covariate adjustment is possible via the argument \code{covars} but
more complex models have to be specified with the methods described above.
This method emulates the function call in the first version of the package.
The arguments \code{formula.full}, \code{formula.red}, \code{model.dat} and
\code{test.terms} are '"missing"' since they are not needed for this method.}
}}

\keyword{ hplot }
