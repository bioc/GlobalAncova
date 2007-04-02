\name{Multiple testing on the GO graph}

  \alias{makeGOstructure}
  \alias{getFocus}
  \alias{GAGO}


\title{Multiple testing on the GO graph}

\description{Three functions adapted from package \pkg{globaltest} to test (part of) the GO graph for
association of the gene expression profile of GO terms with a
response variable. Used together, these functions return
multiplicity-adjusted p-values calculated using the Focus Level
procedure that preserves the structure of the GO graph.}


  \usage{
  GAGO(..., GO, focus, maxalpha = 0.05, stopafter = 100, verbose = FALSE)
  makeGOstructure(data, annotation, top, only.ids, ontology = c("BP", "CC", "MF"), unreliable)
  getFocus(GOstructure, maxatoms = 10)}

\arguments{
    \item{...}{Arguments describing the tests to be performed are passed on to \code{\link{GlobalAncova}}. 
                Note that only the approximative version of \code{\link{GlobalAncova}} is used here and hence the parameter 
                \code{method} is not available. Even though the number of permutations (\code{perm}) may be specified since
                very large gene sets (with more genes than \code{max.group.size}) are treated with the permutation test.}
    \item{GO}{An object of class \code{\link[globaltest]{GOstructure}} describing the structure of the GO
    graph. This object should be created using \code{makeGOstructure}.}
    \item{focus}{A vector of GO ids to describe te focus level. Typically made using \code{getFocus}.}
    \item{maxalpha}{The maximum multiplicity-adjusted p-value. The algorithm will stop when this value is exceeded.}
    \item{stopafter}{The maximum number of significant GO terms to be found. The algorithm will stop when this value is exceeded.}
    \item{verbose}{If set to \code{TRUE}, prints much more extensive progress information.}
    \item{data}{The data set for which the \code{\link[globaltest]{GOstructure}} object is to be made. Can be either an
    \code{\link[Biobase:exprSet-class]{exprSet}} or a vector of probe indicators.}
    \item{annotation}{The name of the metaData annotation package to be used, or a vector of EntrezGene ids with the same
    length as the number of genes in \code{data}.}
    \item{top}{The node to be used as top node of the \code{\link[globaltest]{GOstructure}}. Defaults
    to one of "GO:0008150", "GO:0005575", "GO:0003674", depending on he ontology chosen. The
    \code{\link[globaltest]{GOstructure}} object will only contain offspring of the chosen top node.}
    \item{only.ids}{A vector of GO ids. If this is supplied, GO ids not appearing in this list
    will not appear in the \code{\link[globaltest]{GOstructure}} object.}
    \item{ontology}{The ontology to be used. One of Biological Process (BP), Cellular
    Component (CC) or Molecular Function (MF).}
    \item{unreliable}{Can be used to designate one or more of the GO evidence codes as unreliable.
    Value must be a vector containing one or more from
    "IC", "IDA", "IEA", "IEP", "IGI", "IMP", "IPI", "ISS", "NAS", "ND", "RCA", "TAS", "NR".}
    \item{GOstructure}{An object of class \code{\link[globaltest]{GOstructure}}.}
    \item{maxatoms}{A tuning parameter that governs the choice of the focus level. \code{getFocus} will
    choose those GO terms in the focus level in such a way that all offspring gene sets of each focus level
    term can be constructed as unions of no more than \code{maxatoms} atom gene sets. Default value of
    \code{maxatoms} is 10. Higher values quickly lead to slower computation. Lower values typically lead to reduced power.}
}

\details{These functions should be used in the following order.
First use \code{makeGOstructure} to make a GO graph tailored to a
specific data set. Then \code{getFocus} can be used to choose a
focus level. Finally \code{GAGO} performs the focus level procedure.}

\value{The function returns a named vector of multiplicity-adjusted
p-values. Adjusted p-values of GO terms not appearing in this vector
are larger than the chosen value of \code{maxalpha}.}

\note{Function \code{GAGO} corresponds to function \code{gtGO} in package \pkg{globaltest}. The
    difference is in the use of the \code{GlobalAncova} test instead of \code{globaltest} within the focus level procedure.}

\references{Goeman, J.J. and Mansmann, U., Family-wise error rate on the directed acyclic 
    graph of Gene Ontology, submitted.} 

\author{Jelle Goeman: \email{j.j.goeman@lumc.nl}; Jan Oosting; Manuela Hummel}

\seealso{\code{\link[globaltest]{gtGO}}, \code{\link{GlobalAncova}}, \code{\link[globaltest]{globaltest}}, \code{\link[globaltest]{GOstructure}}}

\examples{
  # see vignettes of packages GlobalAncova and globaltest and help of gtGO
}


\keyword{htest}
