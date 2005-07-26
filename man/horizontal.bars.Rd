\name{horizontal.bars}
\alias{horizontal.bars}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Plotting horizontal bars}
\description{
  Function for plotting horizontal bars which is needed in 'Plot.genes' and 'Plot.subjects'.
}
\usage{
horizontal.bars(x, xlabel = "", ylabel = "", color = NULL, labelsize = 0.75, bar.names = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{A vector of values to be plot.}
  \item{xlabel}{A title for the x axis.}
  \item{ylabel}{A title for the y axis.}
  \item{color}{Colours of the bars.}
  \item{labelsize}{The magnification to be used for axis annotation
          relative to the current.}
  \item{bar.names}{Labels of the bars.}
}
\author{Reinhard Meister \email{meister@tfh-berlin.de}\cr
  Ulrich Mansmann \email{mansmann@ibe.med.uni-muenchen.de}}
\note{For internal use.}
\seealso{\code{\link{Plot.genes}}, \code{\link{Plot.subjects}}, \code{\link{GlobalAncova}}, \code{\link{GlobalAncova.closed}}}
\keyword{internal}% __ONLY ONE__ keyword per line
