# focus-level procedure from Jelle Goeman (combination of closed testing and Holm correction)
#   for multiple testing correction of global tests on the GO,
#   modified for the use of 'GlobalAncova' instead of 'globaltest'

GAGO <- function(..., GO, focus, maxalpha = 0.05, stopafter = 100, verbose = FALSE){
# ...: parameters for function 'GlobalAncova' (parameter 'method' is not available since only approximative test is performed) 
# GO: object of class 'GOstructure'
# focus: vector of GO ids to describe te focus level
# maxalpha: maximum multiplicity-adjusted p-value
# stopafter: maximum number of significant GO terms to be found
# verbose: should extensive progress information be printed

  # check input
  if (missing(GO) || !(is(GO, "GOstructure"))) stop("Provide GO as a GOstructure object")
  if (!all(focus %in% GO@ids)) stop("Not all focus terms found in the GOstructure object")

  # Get the raw p-values
  pFocus <- pGAapprox(..., test.genes=genesets(GO)[focus])    
  
  # Prepare the objects that store the test results
  sigFocus <- logical(length(focus))                # Keeps track of focus nodes which are declared significant
  names(sigFocus) <- focus
  emptyFocus <- logical(length(focus))              # Keeps track of subgraphs which are completely significant
  names(emptyFocus) <- focus
  atoms <- vector("list", length(focus))            # Stores the atoms of offspring sets of each focus node used to make unions
  names(atoms) <- focus
  unions <- vector("list", length(focus))           # Keeps track of all unions that may be tested
  names(unions) <- focus
  pUnions <- vector("list", length(focus))          # Stores the p-values of unions that may be tested
  names(pUnions) <- focus

  sigUnions <- vector("list", length(focus))        # Keeps track of unions that have so far been called significant
  names(sigUnions) <- focus
  offspring <- lapply(as.list(focus), function(term) { # Stores the offspring nodes as unions of atoms
    offnames <- unlist(GO@offspring[term])
    off <- matrix(,length(offnames),0)
    rownames(off) <- offnames
    off
  })
  names(offspring) <- focus
  sigOffspring <- vector("list", length(focus))     # Stores which offspring has already been declared significant
  names(sigOffspring) <- focus
  adjustedP <- rep(NA, length(GO))             # Prepare the vector of corrected p-values
  names(adjustedP) <- GO@ids

  # Find all GO terms above the focus level
  forefathers <- unique(unlist(GO@ancestors[focus]))
  forefathers <- setdiff(forefathers, focus)
  forefathers <- setdiff(forefathers, unique(unlist(GO@offspring[focus])))

  # Initialize
  alpha <- 0
  holm <- length(focus)
  indirectAffected <- list()
  ready <- FALSE
  change <- TRUE

  while (! ready) {

    # Find the focus terms where the action is
    affected1 <- focus[(!sigFocus) & (pFocus * holm <= alpha)]    # Find newly significant focus level terms
    opensubgraphs <- focus[sigFocus & (!emptyFocus)]              # Find the subgraphs in which something may happen
    minP <- sapply(opensubgraphs, function(ff) min(pUnions[[ff]][!sigUnions[[ff]]]))  # Find the smallest p-value not yet declared significant
    if (length(minP) > 0) {
      names(minP) <- opensubgraphs
      affected2 <- names(minP * holm <= alpha)                      # Find those subgraphs for which that minimum is significant
    } else {
      affected2 <- character(0)
    }
    affected <- unique(c(affected1, affected2, names(indirectAffected)))  # Include those significant through upward propagation

    newSigGO <- character(0)
    change <- FALSE

    for (term in affected) {

      # If a focus term was declared significant. Start a new subgraph
      if (!sigFocus[term]) {
        sigFocus[term] <- TRUE
        newSigGO <- c(newSigGO, term)                       # A GO term has been declared significant
        if (verbose) cat("Significant:", term, "\n")
        alloffspring <- GO@offspring[[term]]
        if (length(alloffspring) > 0) {
          offspringSets <- GO@genesets[alloffspring]
          TermAtoms <- globaltest:::getAtoms(offspringSets)                                  # Get the set of atoms
          atoms[[term]] <- TermAtoms
          unions[[term]] <- matrix(rep(TRUE, length(TermAtoms)), 1, length(TermAtoms))  # Prepare the matrix of unions of atoms
          sigUnions[[term]] <- FALSE
          pUnions[[term]] <- pGAapprox(..., test.genes=list(unlist(TermAtoms)))
          offspring[[term]] <- matrix(sapply(TermAtoms, function(y)            # Write all offspring as unions of atoms
            sapply(offspringSets, function(x) all(y %in% x))), ncol = length(TermAtoms))
          rownames(offspring[[term]]) <- names(offspringSets)
          sigOffspring[[term]] <- logical(length(offspringSets))
          names(sigOffspring[[term]]) <- names(offspringSets)
        } else {                                                                # An empty version for when the focus term is an end node
          sigUnions[[term]] <- logical(0)
          atoms[[term]] <- list()
          unions[[term]] <- matrix(,0,0)
          offspring[[term]] <- matrix(,0,0)
          sigOffspring[[term]] <- logical(0)
        }
        change <- TRUE
      }
      # Propagate significance from offspring that was declared significant in another subgraph
      propagate <- indirectAffected[[term]]
      propagate <- propagate[!(propagate %in% term)]
      propagate <- propagate[!sigOffspring[[term]][propagate]]
      if (length(propagate) > 0) {
        pPatterns <- matrix(,0,length(atoms[[term]]))   # Get all superset patterns of a propagated GO term
        for (i in 1:length(propagate)) {
          basePattern <- offspring[[term]][propagate[[i]],] # The propagated GO term as unions of atoms
          ancestorPatterns <- fillPattern(basePattern)      # All superset patterns as unions of atoms
          ancestorPatterns <- ancestorPatterns[!globaltest:::intersectPatterns(ancestorPatterns, pPatterns),,drop=FALSE]    # Remove duplicates
          pPatterns <- rbind(pPatterns, ancestorPatterns)
        }
        matchedPatterns1 <- globaltest:::intersectPatterns(pPatterns, unions[[term]])  # These new patterns were already testable
        matchedPatterns2 <- globaltest:::intersectPatterns(unions[[term]], pPatterns)  # These already testable patterns are now called significant
        newPatterns <- pPatterns[!matchedPatterns1,, drop=FALSE]          # These are new patterns
        if (verbose) cat("\tSignificance of", nrow(newPatterns), "genesets from", propagate, " propagated to", term, "\n")
        unions[[term]] <- rbind(unions[[term]], newPatterns)              # Make all supersets testable and give them raw p 0
        pUnions[[term]][which(matchedPatterns2)] <- 0
        pUnions[[term]] <- c(pUnions[[term]], rep(0, nrow(newPatterns)))
        sigUnions[[term]] <- c(sigUnions[[term]], rep(FALSE, nrow(newPatterns)))    # These patterns will be called significant later
        change <- TRUE
      }

      # Find the testable non-significant unions with small p-values
      newsigs <- (pUnions[[term]] * holm <= alpha) & !sigUnions[[term]]
      newsigs <- which(newsigs)

      # Find whether there are GO terms among the newly significant unions
      sigPatterns <- unions[[term]][newsigs,, drop=FALSE]
      offspringPatterns <- offspring[[term]][!sigOffspring[[term]],,drop=FALSE]     # Patterns of not already significant offspring
      matched <- globaltest:::intersectPatterns(offspringPatterns, sigPatterns)                  # Any new offspring terms among the significant patterns?
      newSigOffspring <- rownames(offspringPatterns)[matched]
      if (verbose && (length(setdiff(newSigOffspring, propagate))>0)) {
        cat("Significant:", setdiff(newSigOffspring, propagate), "\n")
      }
      newSigGO <- c(newSigGO, newSigOffspring)
      sigOffspring[[term]][newSigOffspring] <- TRUE

      # Expand the current subgraph by finding new unions that may now be tested, and find the p-values
      newpatterns <- matrix(, 0, length(atoms[[term]]))
      for (ix in newsigs) {
        sigUnionsMatrix <- unions[[term]][sigUnions[[term]], ,drop = FALSE]         # The unions so far significant for this term
        pattern <- unions[[term]][ix,]      # the significant pattern
        if (sum(pattern) > 1) {
          for (iy in (1:length(pattern))[pattern]) { # For-loop over the TRUEs
            newpattern <- pattern
            newpattern[iy] <- FALSE # All direct subsets have one extra FALSE
            reallyNew <- !any(globaltest:::intersectPatterns(matrix(newpattern,nrow=1), unions[[term]])) # Is the pattern really new?
            if (reallyNew) {      # Are the other parents of this pattern also present?
              parentspresent <- all(sapply((1:length(newpattern))[!pattern], function(iz) { # Loop over the FALSEs of pattern
                newpatternparent <- newpattern
                newpatternparent[iz] <- TRUE                        # The parents of newpattern have one extra TRUE
                any(apply(sigUnionsMatrix, 1, function(pttn) all(pttn == newpatternparent))) # Is this superset present among the significant tests?
              }))
              if (parentspresent) {
                newpatterns <- rbind(newpatterns, newpattern)
                rownames(newpatterns) <- NULL
              }
            }
          }
        }
        sigUnions[[term]][ix] <- TRUE # Only now call the term itself significant (prevents duplicate patterns)
      }
      unions[[term]] <- rbind(unions[[term]], newpatterns)
      sigUnions[[term]] <- c(sigUnions[[term]], rep(FALSE, nrow(newpatterns)))
      if (length(newpatterns) > 0) {                        # Calculate and store their p-values
        if (verbose) cat("\ttests:", nrow(newpatterns), "in", term, "\n")
        newsets <- lapply(as.list(1:nrow(newpatterns)), function(i) { # Assemble gene sets from the atoms
          unlist(atoms[[term]][newpatterns[i,]])
        })
        newpvalues <- pGAapprox(..., test.genes=newsets)        
        pUnions[[term]] <- c(pUnions[[term]], newpvalues)
        change <- TRUE
      }

      # Is the subgraph emptied?
      empty <- (length(atoms[[term]]) == 0) || (all(sigOffspring[[term]]))
      if (empty) {
        emptyFocus[term] <- TRUE
      }
    }

    if (length(newSigGO)>0) {
      change <- TRUE                # Don't change alpha yet; there may be propagation to be done
      adjustedP[newSigGO] <- alpha  # Adjusted P for newly affected GO terms
    }

    # Find all forefather terms that are now significant
    affectedGOforefathers <- unique(unlist(GO@ancestors[newSigGO])) # Find GO terms above the focus level significant through propagation
    affectedGOforefathers <- intersect(affectedGOforefathers, forefathers)
    newlyAffectedGOforefathers <- intersect(affectedGOforefathers, names(adjustedP[is.na(adjustedP)]))
    adjustedP[newlyAffectedGOforefathers] <- alpha
    if (verbose && (length(newlyAffectedGOforefathers) > 0))
      cat("Significant through upward propagation:", newlyAffectedGOforefathers, "\n")

    # For all significant GO terms, see if they are also present in other subtrees
    affectedGOwithParents <- lapply(as.list(newSigGO), function(nsg) {
      present1 <- sapply(offspring, function(off) nsg %in% rownames(off))
      present2 <- focus %in% nsg
      focus[present1 | present2]
    })      # Creates a list of new significant GO terms, listing the subtrees they appear in
    names(affectedGOwithParents) <- newSigGO
    indirectAffected <- globaltest:::turnListAround(affectedGOwithParents)  # Reverses the list to a list of subtrees, listing the significant GO offspring terms

    # Output progress information
    if (!verbose) {
      cat(paste(rep("\b", 48), collapse=""))
      cat("Alpha = ", format(alpha, digits=3, scientific=T, width= 10), ". Significant GO terms: ", format(sum(!is.na(adjustedP)), width=5), ".", sep="")
      flush.console()
    }

    # If nothing happened, increase alpha
    if (!change) {
      allPs <- sort(c(pFocus, unlist(pUnions)))
      allPs <- allPs[allPs * holm > alpha]
      alpha <- allPs[1] * holm
      names(alpha) <- NULL
      if (verbose) cat("Alpha =", format(alpha, digits=3, scientific=T, width= 10), "\n")
    } else {
      if (verbose && (holm > sum(!emptyFocus))) cat("Holm's factor:", sum(!emptyFocus), "\n")
      holm <- sum(!emptyFocus)
    }
    ready <- (holm == 0) || (alpha > maxalpha) || (sum(!is.na(adjustedP)) >= stopafter)
  }
  if (!verbose) cat("\n")
  adjustedP[!is.na(adjustedP)]
}


################################################################################
################################################################################

# simpler GlobalAncova-function that only returns p-values
setGeneric("pGAapprox", function(xx,formula.full,formula.red,model.dat,group,covars=NULL,
                            test.terms,test.genes=NULL,max.group.size=2500,perm=10000,eps=1e-16,acc=50)
           standardGeneric("pGAapprox"))
# xx: expression matrix (rows=genes, columns=subjects)
# formula.full: model formula for the full model
# formula.red: model formula for the reduced model
# model.dat: data frame that contains the group and covariable information
# test.genes: vector of probeset names or a list where each element is a vector of probeset names 
# max.group.size: if a gene set is larger than max.group.size the permutation test is applied even if method="approx"
# perm: number of permutations 
# eps: resuolution of approximation
# acc: additional accuracy parameter for approximation (the higher the value the higher the accuracy) 


################################# general function #############################

setMethod("pGAapprox", signature(xx="matrix",formula.full="formula",formula.red="formula",
                           model.dat="ANY",group="missing",covars="missing",test.terms="missing"),
          definition = function(xx,formula.full,formula.red,model.dat,test.genes=NULL,
                           max.group.size=2500,perm=10000,eps=1e-16,acc=50){
  # test for model.dat
  if(!is.data.frame(model.dat))
    stop("'model.dat' has to be a data frame")

  .pGAapprox(xx=xx,formula.full=formula.full,formula.red=formula.red,model.dat=model.dat,test.genes=test.genes,
            max.group.size=max.group.size,perm=perm,eps=eps,acc=acc)
})



########################## function for 2 groups ###############################

setMethod("pGAapprox", signature(xx="matrix",formula.full="missing",formula.red="missing",
                           model.dat="missing",group="ANY",covars="ANY",test.terms="missing"),
          definition = function(xx,group,covars=NULL,test.genes=NULL,
                           max.group.size=2500,perm=10000,eps=1e-16,acc=50){

  # parameter names
  group.name   <- deparse(substitute(group))

  if(is.null(dim(covars)))
    covar.names <- deparse(substitute(covars))
  else
    covar.names <- colnames(covars)

  # get formulas and 'model.dat' out of 'group' and 'covars'
  res          <- group2formula(group=group, group.name=group.name, covars=covars, covar.names)
  formula.full <- res$formula.full
  formula.red  <- res$formula.red
  model.dat    <- res$model.dat

  # then apply the usual function
  .pGAapprox(xx=xx,formula.full=formula.full,formula.red=formula.red,
                model.dat=model.dat,test.genes=test.genes,max.group.size=max.group.size,perm=perm,eps=eps,acc=acc)
})


############################# with 'test.terms' ################################

setMethod("pGAapprox", signature(xx="matrix",formula.full="formula",formula.red="missing",
                           model.dat="ANY",group="missing",covars="missing",test.terms="character"),
          definition = function(xx,formula.full,test.terms,model.dat,test.genes=NULL,
                           max.group.size=2500,perm=10000,eps=1e-16,acc=50){

  # test for model.dat
  if(!is.data.frame(model.dat))
    stop("'model.dat' has to be a data frame")

  # test for 'test.terms'
  terms.all <- test.terms
  D.full    <- model.matrix(formula.full, model.dat)
  terms.all <- colnames(D.full)

  # are all terms variables compatible with 'model.dat'?
  if(!all(test.terms %in% terms.all))
    stop("'test.terms' are not compatible with the specified models")

  D.full <- model.matrix(formula.full, data=model.dat)
  D.red  <- D.full[,!(colnames(D.full) %in% test.terms), drop=FALSE]

  # then apply the usual function
  .pGAapprox(xx=xx,formula.full=formula.full,D.red=D.red,
                model.dat=model.dat,test.genes=test.genes,max.group.size=max.group.size,perm=perm,eps=eps,acc=acc)
})


################################################################################

# main function 
.pGAapprox <- function(xx,formula.full,formula.red=NULL,D.red=NULL,model.dat,test.genes=NULL,
                     max.group.size=2500,perm=10000,eps=1e-16,acc=50){

    if(!is.list(test.genes))
        test.genes <- list(test.genes)

    xx2        <- xx[unique(unlist(test.genes)),,drop=F]
    N.Genes    <- sapply(test.genes, length)
    N.Subjects <- ncol(xx2)
    N.tests    <- length(test.genes)

    # design matrices
    D.full     <- model.matrix(formula.full, data=model.dat)
    if(is.null(D.red))
        D.red  <- model.matrix(formula.red,  data=model.dat)
    N.par.full <- ncol(D.full)
    N.par.red  <- ncol(D.red)

    # checking for NA's
    if(nrow(D.full) < N.Subjects)
        stop("Missing values in the model variables")

     # extra sum of squares
     genewiseSS <- genewiseGA(xx2, D.full, D.red=D.red)
     SS.extra   <- sapply(test.genes, function(x) sum(genewiseSS[x,"nominator"]))

     # get eigen values
     ew.H.nom   <- eigen(hat.matrix(D.full) - hat.matrix(D.red))$values

     # remove gene sets that are bigger than 'max.group.size' -> treat those with permutation test
     toobig <- sum(N.Genes > max.group.size) > 0

     w <- which(N.Genes <= max.group.size)
     test.genes.red <- test.genes[w]

     ew.cov     <- sapply(test.genes.red, function(y) eigen(cov.shrink(t(xx2[y,,drop=FALSE]), verbose=FALSE), only.values = TRUE)$values)
     if(!is.list(ew.cov))    # if only one gene group is tested
         ew.cov <- list(ew.cov)
     # all pairwise products of eigen values
     ew.nom   <- sapply(ew.cov, function(x) as.vector(outer(x, ew.H.nom, "*")))
     if(!is.list(ew.nom))    # if only one gene group is tested
         ew.nom <- list(ew.nom)

     # compute approximate p-values
     p.approx <- rep(NA, N.tests)
     if(length(w) > 0)
       for(i in 1:length(w))
           p.approx[w[i]] <- .pAsymptotic(x = SS.extra[w[i]], lams = ew.nom[[i]], eps = eps, acc = acc)

     # for large gene sets use permutation test
     if(toobig) {
        w <- which(N.Genes > max.group.size)
        test.genes.red <- test.genes[w]
        xx3 <- xx2[unique(unlist(test.genes.red)),,drop=F]
        DF.full  <- N.Genes[w] * (N.Subjects - N.par.full)
        DF.extra <- N.Genes[w] * (N.par.full - N.par.red)
        SS.full    <- sapply(test.genes.red, function(x) sum(genewiseSS[x,"denominator"]))
        SS.extra   <- SS.extra[w]
        MS.full    <- SS.full / DF.full
        MS.extra   <- SS.extra / DF.extra
        F.value <- MS.extra / MS.full

        p.approx[w] <- resampleGA(xx3, D.full, D.red, perm, test.genes.red, F.value, DF.full, DF.extra)
     }

    return(p.approx)
}





