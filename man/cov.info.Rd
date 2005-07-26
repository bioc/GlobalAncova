\name{cov.info}
\alias{cov.info}
\docType{data}
\title{Data of the covariates}
\description{
  Covariate data for the p53 example:
  Date of covariates sex (\code{0}-male/\code{1}-female) and tumour location (\code{1}-distal/\code{0}-proximal)  }
\usage{data(cov.info)}
\format{
  The format is:
 num [1:36, 1:2] 1 1 1 1 0 0 1 1 0 1 ...
 - attr(*, "dimnames")=List of 2
  ..$ : NULL
  ..$ : chr [1:2] "sex" "loc"
  }
\examples{
data(cov.info)
#str(cov.info)
#table(cov.info[,"sex"])
}
\keyword{datasets}
