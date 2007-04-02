\name{colon.pheno}
\alias{colon.pheno}
\docType{data}
\title{Covariate information for the colon data}
\description{
Covariate data for the colon data example:
\item{sex}{Sex of the patient.}
\item{age}{Age of the patient.}
\item{location}{Location of the tumour.}
\item{grade}{Histologic tumour grade.}
\item{UICC.stage}{UICC stage of colorectal carcinoma.}
}
\usage{data(colon.pheno)}
\format{
 The format is:    \cr
 \code{`data.frame':   12 obs. of  5 variables:}  \cr
 \code{$ sex       : Factor w/ 2 levels "0","1": 2 2 1 2 2 1 2 1 2 1 ...}  \cr
 \code{$ age       : int  71 76 63 73 58 66 60 66 86 76 ...}   \cr
 \code{$ location  : Factor w/ 2 levels "distal","proximal": 1 1 1 1 1 1 1 1 2 1 ...} \cr
 \code{$ grade     : Factor w/ 2 levels "2","3": 1 1 2 2 1 2 1 2 2 2 ...}   \cr
 \code{$ UICC.stage: Factor w/ 2 levels "2","3": 2 1 2 1 2 1 1 1 2 1 ...}
  }
\references{Groene, J. et al., 2006, Transcriptional census of 36 microdissected colorectal cancers yields a gene signature 
to distinguish UICC II and III, \emph{Int J Cancer} 119(8):1829--36.}
\examples{
data(colon.pheno)
#str(colon.pheno)
}
\keyword{datasets}
