\name{colon.tumour}
\alias{colon.tumour}
\docType{data}
\title{ Gene expression data }
\description{
Normalized gene expression data of 12 patients with colorectal cancer.
Samples are taken from inside the tumours. Additionally, from same patients samples are
taken from normal tissue, see \code{\link{colon.normal}}. The expression matrix is only an
exemplary subset of 1747 probe sets associated with cell proliferation.
}
\usage{data(colon.tumour)}
\format{
 The format is: \cr
\code{num [1:1747, 1:12]  8.77 10.40  8.52 12.86  8.28 ...} \cr
 \code{- attr(*, "dimnames")=List of 2}                                      \cr
  \code{..$ : chr [1:1747] "200808_s_at" "215706_x_at" "217185_s_at" "202136_at" ...}  \cr
  \code{..$ : chr [1:12] "Co10.T.IT.83.F.CEL" "Co14.T.IT.88.F.CEL" "Co17.T.IT.563.F.CEL" "Co1.T.IT.30.F.CEL" ...}
  }
\references{Groene, J. et al., 2006, Transcriptional census of 36 microdissected colorectal cancers yields a gene signature 
to distinguish UICC II and III, \emph{Int J Cancer} 119(8):1829--36.}
\examples{
data(colon.tumour)
#str(colon.tumour)
}
\keyword{datasets}
