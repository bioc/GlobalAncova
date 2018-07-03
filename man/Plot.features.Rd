\name{Plot.features}
\alias{Plot.features}
\title{Features Plot for generalized Global Ancova}
\description{Produces a plot to show the influence of individual variables on the test result produced by \code{\link{gGlobalAncova}}. The variables can be continuous, binary, categorical, ordinal, or of mixed types.}
\usage{
Plot.features(data, formula.full, formula.red = ~1, model.dat, Set, returnValues = FALSE, ...)
}

\arguments{
  \item{data}{\code{data.frame} of variables to be tested in sets (columns=variables); (multi-) categorical variables should be factors, ordinal variables should be ordered factors}
  \item{formula.full}{model formula for the full model}
  \item{formula.red}{model formula for the reduced model (that does not contain the terms of interest)}
  \item{model.dat}{\code{data.frame} of regressors, containing variables specified in \code{formula.full} and \code{formula.red}}
  \item{Set}{optional vector of names or indices of variables, defining the set of variables to plot; if missing, all variables in \code{data} are shown}
  \item{returnValues}{shall variable-wise statistics = bar heights be returned?}
  \item{\dots}{graphical parameters passed to \code{barplot}}
}

%\details{}

\value{
If \code{returnValues = TRUE}, a vector with the bar heights is returned.
}
%\references{}
\author{Manuela Hummel \email{m.hummel@dkfz.de}}
%\note{}


\seealso{\code{\link{gGlobalAncova}}, \code{\link{Plot.genes}}}

\examples{
data(bindata)
Plot.features(bindata[,-1], formula.full = ~group, model.dat = bindata)
}
\keyword{ hplot }% use one of  RShowDoc("KEYWORDS")
