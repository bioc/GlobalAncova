\name{p53.signalling}
\alias{p53.signalling}
\docType{data}
\title{ Gene expression Data }
\description{
  RMA normalized gene expression data of 45 probesets on the HGU133A GeneChip which belong to the p53 signalling pathway.
  For details see Mansmann, U. and Meister, R., 2005, Testing differential gene expression
  in functional groups, \emph{Methods Inf Med} 44 (3).
  }
\usage{data(p53.signalling)}
\format{
  The format is:
 num [1:45, 1:36] 8.14 7.38 9.46 8.18 7.15 ...
 - attr(*, "dimnames")=List of 2
  ..$ : chr [1:45] "204859\\_s\\_at" "211553\\_x\\_at" "211554\\_s\\_at" "212672\\_at" ...
  ..$ : chr [1:36] "UICC.2" "UICC.3" "UICC.3" "UICC.2" ...
  }
\examples{
data(p53.signalling)
#boxplot(p53.signalling[group.info==0,1],p53.signalling[group.info==1,1])
}
\keyword{datasets}
