\name{GlobalAncova gene set testing methods}

\alias{GAGO}
\alias{GABroad}


\title{Gene set testing of gene set databases using GlobalAncova}

\description{Three functions adapted from package \pkg{globaltest} to test gene sets from databases for
association of the gene expression profile with a
response variable. Three function are provided for Gene Ontology and for the Broad Institute's gene sets.}


  \usage{
GAGO (xx, ..., id, annotation, probe2entrez,
           ontology = c("BP", "CC", "MF"), minsize=1, maxsize=Inf,
           multtest = c("holm", "focuslevel", "BH", "BY"),
           focuslevel = 10, sort = TRUE)

GABroad (xx, ..., id, annotation, probe2entrez, collection,
           category = c("c1", "c2", "c3", "c4", "c5"),
           multtest = c("holm", "BH", "BY"), sort = TRUE)
  }

\arguments{
    \item{xx}{Matrix of gene expression data, where columns correspond to samples
             and rows to genes. Gene names have to be included as the row names of \code{xx}}
    \item{...}{Arguments describing the tests to be performed are passed on to \code{\link{GlobalAncova}}. 
                Note that only the approximative version of \code{\link{GlobalAncova}} is used here and
                hence the parameter \code{method} is not available. Even though the number of permutations
                (\code{perm}) may be specified since very large gene sets (with more genes than
                \code{max.group.size}) are treated with the permutation test.}
    \item{id}{The identifier(s) of gene sets to be tested (character vector). If omitted, tests all
              gene sets in the database.}
    \item{annotation}{The name of the probe annotation package for the microarray that was used, or
                      the name of the genome wide annotation package for the species
                      (e.g. org.Hs.eg.db for human). If an organism package is given, the argument
                      \code{probe2entrez} must be supplied.}
    \item{probe2entrez}{Use only if no probe annotation package is available. A mapping from probe
                      identifiers to entrez gene ids. May be an environment, named list or named vector.}
    \item{multtest}{The method of multiple testing correction. Choose from: Benjamini and Hochberg FDR
                    control (BH); Benjamini and Yekutieli FDR control (BY) or Holm familywise error
                    control (holm). For \code{GAGO} also the focus level method is available.
                    See \code{\link[globaltest:multtest]{focusLevel}}.}
    \item{sort}{If \code{TRUE}, sorts the results to increasing p-values.}
    \item{ontology}{The ontology or ontologies to be used. Default is to use all three ontologies.}
    \item{minsize}{The minimum number of probes that may be annotated to a gene set. Gene sets with
                  fewer annotated probes are discarded.}
    \item{maxsize}{The maximum number of probes that may be annotated to a gene set. Gene sets with
                  more annotated probes are discarded.}
    \item{focuslevel}{The focus level to be used for the focus level method. Either a vector of gene
                  set ids, or a numerical level. In the latter case, \code{\link[globaltest:multtest]{findFocus}}
                  is called with \code{maxsize} at the specified level to find a focus level.}
    \item{collection}{The Broad gene set collection, created by a call to
                    \code{\link[GSEABase:getObjects]{getBroadSets}}.}
    \item{category}{The subcategory of the Broad collection to be tested. The default is to test all sets.}
}

\details{These are utility functions to make it easier to do gene set testing of gene sets available
        in gene set databases. The functions automatically retrieve the gene sets, preprocess and
        select them, perform global test, do multiple testing correction, and sort the results on
        the basis of their p-values.
        All functions require that \code{annotate} and the appropriate annotation packages are installed.
        \code{GAGO} requires the
        \code{GO.db} package; \code{GABroad} requires the user to download the XML file "msigdb_v2.5.xml"
        from \\ \code{http://www.broad.mit.edu/gsea/downloads.jsp}, and to preprocess that file using
        the \code{\link[GSEABase:getObjects]{getBroadSets}} function.
}


\value{The function returns a data frame with raw and multiplicity-adjusted
p-values for each gene set.}

\note{Functions \code{GAGO}and \code{GABroad} correspond to functions \code{gtGO},
      and \code{gtBroad} in package \pkg{globaltest}. The
    difference is in the use of the \code{GlobalAncova} test instead of \code{gt} within the procedures.}

\references{Goeman, J.J. and Mansmann, U., Multiple testing on the directed acyclic graph of Gene Ontology.
          Bioinformatics 2008; 24(4): 537-44.}

\author{Jelle Goeman: \email{j.j.goeman@lumc.nl}; Jan Oosting; Manuela Hummel}

\seealso{\code{\link[globaltest:genesettesting]{gtGO}},
  \code{\link[globaltest:genesettesting]{gtKEGG}},
  \code{\link[globaltest:genesettesting]{gtBroad}},
  \code{\link{GlobalAncova}}, 
  \code{\link[globaltest:gt]{gt}},
 }
 

\examples{
  # see vignettes of packages GlobalAncova and globaltest and help of gtGO
}


\keyword{htest}
