\name{phenodata}
\alias{phenodata}
\docType{data}
\title{Covariate information for the van t'Veer data}
\description{
Covariate data for the van t'Veer example:
\describe{
  \item{Sample}{Sample number.}
  \item{metastases}{Development of distant metastases within five years (\code{0}-no/\code{1}-yes).}
  \item{grade}{Tumor grade (three ordere levels).}
  \item{ERstatus}{Estrogen receptor status (\code{pos}-positive/\code{neg}-negative).}
}
}
\usage{data(phenodata)}
\format{
 The format is:    
 \describe{
   \item{\code{'data.frame'}:}{96 obs. of  4 variables:}  
   \item{\code{$Sample}:}{int  1 2 3 4 5 6 7 8 9 10 ...}  
   \item{\code{$metastases}:}{int  0 0 0 0 0 0 0 0 0 0 ...}   
   \item{\code{$grade}:}{int  2 1 3 3 3 2 1 3 3 2 ...}   
   \item{\code{$ERstatus}:}{Factor w/ 2 levels "neg","pos": 2 2 1 2 2 2 2 1 2
                           2 ...} 
  }
}
\examples{
data(phenodata)
#str(phenodata)
}
\keyword{datasets}
