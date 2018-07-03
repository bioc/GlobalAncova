\name{gGlobalAncova.hierarchical}
\alias{gGlobalAncova.hierarchical}
\title{Hierarchical testing procedure using generalized GlobalAncova}

\description{Hierarchical testing procedure according to Meinshausen (2008) screening for groups of related variables within a hierarchy instead of screening individual variables independently. Groups are tested by the generalized GlobalAncova approach. The family-wise error rate is simultaneously controlled over all levels of the hierarchy. In order to reduce computational complexity for large hierarchies, a "short cut" is implemented, where the testing procedure is applied separately to K sub-hierarchies. The p-values are adjusted such that they are identical to the ones obtained when testing the complete hierarchy.}

\usage{
gGlobalAncova.hierarchical(data, H, formula.full, formula.red=~1, model.dat, sumstat=sum, 
                                       alpha=0.05, K, perm=10000, returnPermstats=FALSE, permstats)
}
\arguments{
  \item{data}{\code{data.frame} of variables (columns=variables) to be tested hierarchically; (multi-) categorical variables should be factors, ordinal variables should be ordered factors}
  \item{H}{dendrogram object specifying the hierarchy of the variables; \code{labels(H)} has to coincide with \code{colnames(data)}}
  \item{formula.full}{model formula for the full model}
  \item{formula.red}{model formula for the reduced model (that does not contain the terms of interest)}
  \item{model.dat}{\code{data.frame} of regressors, containing variables specified in \code{formula.full} and \code{formula.red}}
  \item{sumstat}{function for summarizing univariate test statistics; default is \code{sum}}
  \item{alpha}{global significance level}
  \item{K}{optional integer; if this is specified, "short cut" on hierarchical testing will be applied separately to \code{K} sub-hierarchies}
  \item{perm}{number of permutations}
  \item{returnPermstats}{if \code{TRUE}, the variable-wise statistics for all permutations are returned}
  \item{permstats}{if variable-wise permutation statistics were calculated previously, they can be provided in order not to repeat permutation testing (but only the hierarchical prodcedure); useful e.g. if procedure is run again with different \code{alpha} and/or hierarchy \code{H}; NOTE: \code{data}, \code{formula.full} and \code{formula.red} must be identical to the previous call}
}

\details{The hierarchical procedure starts with testing the global null hypothesis that all variables are not associated with the design of interest, and then moves down the given hierarchy testing subclusters of variables. A subcluster is only tested if the null hypothesis corresponding to its ancestor cluster could be rejected. The p-values are adjusted for multiple testing according to cluster size \eqn{p_{C,adj} = p_C m/|C|}, where \eqn{m} is the total number of variables and \eqn{|C|} is the number of variables in cluster \eqn{C}.

If \code{K} is specified and the procedure is split to \code{K} sub-hierarchies containing \eqn{m_1, \ldots, m_K} variables, p-values are additionally adjusted by \eqn{\tau = m / m_k, k=1, \ldots, K}, such that resulting p-values are identical to the ones obtained when testing the complete hierarchy
\deqn{p_{C,adj,k} \cdot \tau = p_C m_k/|C| \cdot m/m_k = p_{C,adj}}
}

\value{an object of class \code{\link{GAhier}}}

\references{Meinshausen N, 2008. Hierarchical testing of variable importance. Biometrika, 95(2):265}
\author{Manuela Hummel \email{m.hummel@dkfz.de}}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{\code{\link{gGlobalAncova}}, \code{\link{GAhier}}, \code{\link{Plot.hierarchy}}}

\examples{
data(bindata)
X <- as.matrix(bindata[,-1])

# get a hierarchy for variables
dend <- as.dendrogram(hclust(dist(t(X))))

# hierarchical test
set.seed(555)
res <- gGlobalAncova.hierarchical(X, H = dend, formula.full = ~group, model.dat = bindata, alpha = 0.05, perm = 1000)
res
results(res)

# get names of significant clusters
sigEndnodes(res)

# visualize results
Plot.hierarchy(res, dend)

# starting with 3 sub-hierarchies
set.seed(555)
res2 <- gGlobalAncova.hierarchical(X, H = dend, K = 3, formula.full = ~group, model.dat = bindata, alpha = 0.05, perm = 1000)

results(res2)
}
\keyword{ models }
