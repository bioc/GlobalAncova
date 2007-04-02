\name{colon.normal}
\alias{colon.normal}
\docType{data}
\title{ Gene expression data }
\description{
Normalized gene expression data of 12 patients with colorectal cancer.
Samples are taken from inside the tumours. Additionally, from same patients samples are
taken from normal tissue, see \code{\link{colon.normal}}. The expression matrix is only an
exemplary subset of 1747 probe sets associated with cell proliferation.
}
\usage{data(colon.normal)}
\format{
 The format is: \cr
\code{num [1:1747, 1:12]  8.74 10.53  8.48 12.69  8.55 ...} \cr
 \code{- attr(*, "dimnames")=List of 2}                                      \cr
  \code{..$ : chr [1:1747] "200808_s_at" "215706_x_at" "217185_s_at" "202136_at" ...}  \cr
  \code{..$ : chr [1:12] "Co10.N.E.84.F.CEL" "Co14.N.E.89.F.CEL" "Co17.N.E.1037.F.CEL" "Co1.N.E.31.F.CEL" ...}
  }
\references{Groene, J. et al., 2006, Transcriptional census of 36 microdissected colorectal cancers yields a gene signature 
to distinguish UICC II and III, \emph{Int J Cancer} 119(8):1829--36.}
\examples{
data(colon.normal)
#str(colon.normal)
}
\keyword{datasets}
