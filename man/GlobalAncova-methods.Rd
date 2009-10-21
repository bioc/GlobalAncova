\name{GlobalAncova-methods}
\alias{GlobalAncova-methods}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Methods for Function GlobalAncova}
\description{There are three possible ways of using \code{GlobalAncova}. The general way is to define
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
\item{xx = "matrix", formula.full = "formula", formula.red = "formula", model.dat = "ANY", group = "missing", covars = "missing", test.terms = "missing"}{In this method, besides the expression matrix \code{xx}, model formulas for the full
and reduced model and a data frame \code{model.dat} specifying corresponding model
terms have to be given. Terms that are included in the full but not in the reduced
model are those whose association with differential expression will be tested.
The arguments \code{group}, \code{covars} and \code{test.terms} are '"missing"'
since they are not needed for this method.}

\item{xx = "matrix", formula.full = "formula", formula.red = "missing", model.dat = "ANY", group = "missing", covars = "missing", test.terms = "character"}{In this method, besides the expression matrix \code{xx}, a model formula for the full
model and a data frame \code{model.dat} specifying corresponding model
terms are required. The character argument \code{test.terms} names the terms of interest
whose association with differential expression will be tested. The basic idea
behind this method is that one can select single terms, possibly from the list of
terms provided by previous \code{GlobalAncova} output, and test them without having
to specify each time a model formula for the reduced model.
The arguments \code{formula.red}, \code{group} and \code{covars} are '"missing"'
since they are not needed for this method.}


\item{xx = "matrix", formula.full = "missing", formula.red = "missing", model.dat = "missing", group = "ANY", covars = "ANY", test.terms = "missing"}{Besides the expression matrix \code{xx} a clinical variable \code{group} is
required. Covariate adjustment is possible via the argument \code{covars} but
more complex models have to be specified with the methods described above.
This method emulates the function call in the first version of the package.
The arguments \code{formula.full}, \code{formula.red}, \code{model.dat} and
\code{test.terms} are '"missing"' since they are not needed for this method.}
}}

\keyword{ methods }% at least one, from doc/KEYWORDS

