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
  }

\arguments{
    \item{...}{Arguments describing the tests to be performed are passed on to \code{\link{GlobalAncova}}. 
                Note that only the approximative version of \code{\link{GlobalAncova}} is used here and hence the parameter 
                \code{method} is not available. Even though the number of permutations (\code{perm}) may be specified since
                very large gene sets (with more genes than \code{max.group.size}) are treated with the permutation test.}
    \item{GO}{An object of class
      \code{\link[globaltest:v4_GOstructure]{GOstructure}} describing
      the structure of the GO 
    graph. This object should be created using \code{makeGOstructure}.}
    \item{focus}{A vector of GO ids to describe te focus level. Typically made using \code{getFocus}.}
    \item{maxalpha}{The maximum multiplicity-adjusted p-value. The algorithm will stop when this value is exceeded.}
    \item{stopafter}{The maximum number of significant GO terms to be found. The algorithm will stop when this value is exceeded.}
    \item{verbose}{If set to \code{TRUE}, prints much more extensive progress information.}
}

\details{Previous to a call to \code{GAGO}, first use
  \code{\link[globaltest:v4_getFocus]{makeGOstructure}} to make a GO
  graph tailored to a 
specific data set. Then \code{\link[globaltest:v4_getFocus]{getFocus}}
can be used to choose a 
focus level. Finally \code{GAGO} performs the focus level procedure.}

\value{The function returns a named vector of multiplicity-adjusted
p-values. Adjusted p-values of GO terms not appearing in this vector
are larger than the chosen value of \code{maxalpha}.}

\note{Function \code{GAGO} corresponds to function \code{gtGO} in package \pkg{globaltest}. The
    difference is in the use of the \code{GlobalAncova} test instead of \code{globaltest} within the focus level procedure.}

\references{Goeman, J.J. and Mansmann, U., Family-wise error rate on the directed acyclic 
    graph of Gene Ontology, submitted.} 

\author{Jelle Goeman: \email{j.j.goeman@lumc.nl}; Jan Oosting; Manuela Hummel}

\seealso{\code{\link[globaltest:genesettesting]{gtGO}},
  \code{\link{GlobalAncova}}, 
  \code{\link[globaltest:v4_globaltest]{globaltest}},  
        \code{\link[globaltest:v4_GOstructure]{GOstructure}},
	\code{\link[globaltest:v4_getFocus]{makeGOstructure}},
	\code{\link[globaltest:v4_getFouse]{getFocus}}} 

\examples{
  # see vignettes of packages GlobalAncova and globaltest and help of gtGO
}


\keyword{htest}
