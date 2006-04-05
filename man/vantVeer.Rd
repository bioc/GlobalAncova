\name{vantVeer}
\alias{vantVeer}
\docType{data}
\title{ Gene expression data }
\description{
Normalized gene expression data for the van t'Veer example:
A subset of 96 samples without BRCA1 or BRCA2 mutations and 1113 genes associated with nine
cancer related pathways (see also \code{?pathways}) was chosen.
}
\usage{data(vantVeer)}
\format{
 The format is: \cr
\code{num [1:1113, 1:96] 0.13 0.936 -0.087 0.118 0.168 -0.081 0.023 -0.086 -0.154 0.025 ...} \cr
 \code{- attr(*, "dimnames")=List of 2}                                      \cr
  \code{..$ : chr [1:1113] "AW025529" "NM_001648" "NM_001753" "NM_003298" ...}  \cr
  \code{..$ : chr [1:96] "1" "2" "3" "4" ...}
  }
\examples{
data(vantVeer)
#str(vantVeer)
}
\keyword{datasets}
