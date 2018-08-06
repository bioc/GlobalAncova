# geneset testing function from Jelle Goeman
#   for multiple testing correction of global tests on the GO (choose between focus level, BH, BY and Holm method),
#   modified for the use of 'GlobalAncova' instead of 'globaltest'

GAGO <- function(xx, ..., id, annotation, probe2entrez, ontology = c("BP", "CC", "MF"),
                  minsize=1, maxsize=Inf,
                  multtest = c("holm", "focuslevel", "BH", "BY"), focuslevel = 10,
                  sort = TRUE) {
# xx: expression matrix (rows=genes, columns=subjects); rownames have to be gene identifiers
# ...: parameters for function 'GlobalAncova' (parameter 'method' is not available since only approximative test is performed)
# id: gene set identifiers. If omitted, tests all gene sets in the database
# annotation: name of the probe annotation package or the name of the genome wide annotation package for the species
# probe2entrez: Use only if no probe annotation package is available. A mapping from probe identifiers to entrez gene ids. May be an environment, named list or named vector.
# ontology: choose one or more ontologies. Default is to use all three ontologies
# minsize, maxsize: size restriction for GO terms
# multtest: method for multiple testing correction
# focuslevel: The focus level to be used for the focus level method. Either a vector of gene set ids, or a numerical level. In the latter case, 'findFocus' is called with 'maxsize' at the specified level to find a focus level
# sort: if TRUE, sorts the results to increasing p-values

  # get the right annotation package
  if (substr(annotation, nchar(annotation)-2, nchar(annotation)) == ".db")
    annotation <- substr(annotation, 1, nchar(annotation)-3)
  package <- paste(annotation, ".db", sep="")
  #require(package, character.only=TRUE) || stop("package ", package, " is not available")
  requireNamespace(package) || stop("package ", package, " is not available")
  requireNamespace("GO.db") || stop("package GO.db is not available")
  requireNamespace("annotate") || stop("package annotate is not available")

  # check whether "org" package is given
  if (substr(annotation,1,4) == "org.")
    extension <- "GO2ALLEGS"
  else
    extension <- "GO2ALLPROBES"
  GOOBJECT <- eval(as.name(paste(annotation, extension, sep="")))

  # reduce the terms
  if (missing(id))
    id <- AnnotationDbi::mappedkeys(GOOBJECT)
  myGOTERM <- annotate::lookUp(id, "GO", "TERM")

  # get the right ontology/ies
  if (!missing(ontology)) {
    choose <- sapply(myGOTERM, Ontology) %in% ontology
    id <- id[choose]
    myGOTERM <- myGOTERM[choose]
  }

  # retrieve sets
  sets <- annotate::lookUp(id, annotation, extension)
  sets <- lapply(sets, function(st) if (all(is.na(st))) character(0) else st)

  # map back
  if (!missing(probe2entrez)) {
    if (is.environment(probe2entrez) || is(probe2entrez, "AnnDbBimap")) probe2entrez <- as.list(probe2entrez)
    if (is.list(probe2entrez)) probe2entrez <- unlist(probe2entrez)
    sets <- lapply(sets, function(set) {
      names(probe2entrez)[probe2entrez %in% set]
    })
  }

  # reduce to genes available in xx
  probes <- rownames(xx)
  sets <- lapply(sets, function(set) intersect(set, probes))

  # size restrictions
  size <- sapply(sets, length)
  choose <- size >= minsize & size <= maxsize
  sets <- sets[choose]
  myGOTERM <- myGOTERM[choose]
  if (length(sets) == 0) stop("No GO terms with size between \"minsize\" and \"maxsize\"")

  # perform tests and do multiple testing
  if (length(sets) > 1) {
    multtest <- match.arg(multtest)
    if (multtest == "focuslevel") {
      ancestors <- lapply(as.list(ontology), function(ont) {
        ext <- paste(ont, "ANCESTOR", sep="")
        GOOBJ <- eval(as.name(paste("GO", ont, "ANCESTOR", sep="")))
        ontid <- intersect(AnnotationDbi::keys(GOOBJ), id)
        if (length(ontid) > 0) annotate::lookUp(ontid, "GO", ext) else list()
      })
      ancestors <- do.call(c, ancestors)
      offspring <- lapply(as.list(ontology), function(ont) {
        ext <- paste(ont, "OFFSPRING", sep="")
        GOOBJ <- eval(as.name(paste("GO", ont, "OFFSPRING", sep="")))
        ontid <- intersect(AnnotationDbi::keys(GOOBJ), id)
        if (length(ontid) > 0) annotate::lookUp(ontid, "GO", ext) else list()
      })
      offspring <- do.call(c, offspring)
      offspring <- sapply(offspring, function(os) if (all(is.na(os))) character(0) else os)
      if (is.numeric(focuslevel))
        focuslevel <- globaltest::findFocus(sets, ancestors = ancestors, offspring = offspring, maxsize = focuslevel)

      # function that only takes gene names of a set and only returns p-value
      helpGA <- function(set){
        pGAapprox(xx=xx, ..., test.genes=set)
      }

      res <- globaltest::focusLevel(helpGA, sets=sets, focus=focuslevel, ancestors = ancestors, offspring = offspring)
    } else {
      raw.p <- pGAapprox(xx=xx, ..., test.genes=sets)
      adj.p <- p.adjust(raw.p, method=multtest)
      res <- data.frame(raw.p, adj.p)
      dimnames(res) <- list(names(sets), c("raw.p", multtest))
    }
    # add term names
    res$Term <- sapply(myGOTERM, function(mgt) if (is(mgt, "GOTerms")) Term(mgt) else "")
    res$Term[is.na(res$Term)] <- ""

  } else {
    res <- pGAapprox(xx=xx, ..., test.genes=sets)
    names(res) <- names(sets)
  }

  if (sort & !is.null(dim(res)))
    res <- res[order(res[,2]),]

  res
}


# !! TODO:
# functions to get the (significant) leaf nodes and to visualize the significant subgraph



################################################################################

# function to test KEGG pathways

GAKEGG <- function(xx, ..., id, annotation, probe2entrez,
                  multtest = c("holm", "BH", "BY"), sort = TRUE) {
# xx: expression matrix (rows=genes, columns=subjects); rownames have to be gene identifiers
# ...: parameters for function 'GlobalAncova' (parameter 'method' is not available since only approximative test is performed)
# id: gene set identifiers. If omitted, tests all gene sets in the database
# annotation: name of the probe annotation package or the name of the genome wide annotation package for the species
# probe2entrez: Use only if no probe annotation package is available. A mapping from probe identifiers to entrez gene ids. May be an environment, named list or named vector.
# multtest: method for multiple testing correction
# sort: if TRUE, sorts the results to increasing p-values

  # get the right annotation package
  if (substr(annotation, nchar(annotation)-2, nchar(annotation)) == ".db")
    annotation <- substr(annotation, 1, nchar(annotation)-3)
  package <- paste(annotation, ".db", sep="")
  #require(package, character.only=TRUE) || stop("package ", package, " is not available")
  requireNamespace(package) || stop("package ", package, " is not available")
  requireNamespace("KEGG.db") || stop("package KEGG.db is not available")
  requireNamespace("annotate") || stop("package annotate is not available")

  # check whether "org" package is given
  if (substr(annotation,1,4) == "org.")
    extension <- "PATH2EG"
  else
    extension <- "PATH2PROBE"
  KEGGOBJECT <- eval(as.name(paste(annotation, extension, sep="")))

  # default terms
  if (missing(id))
    id <- AnnotationDbi::mappedkeys(KEGGOBJECT)

  # retrieve sets
  sets <- annotate::lookUp(id, annotation, extension)

  # map back
  if (!missing(probe2entrez)) {
    if (is.environment(probe2entrez) || is(probe2entrez, "AnnDbBimap")) probe2entrez <- as.list(probe2entrez)
    if (is.list(probe2entrez)) probe2entrez <- unlist(probe2entrez)
    sets <- lapply(sets, function(set) {
      names(probe2entrez)[probe2entrez %in% set]
    })
  }

  # reduce to genes available in xx
  probes <- rownames(xx)
  sets <- lapply(sets, function(set) intersect(set, probes))

  # remove sets without any gene in xx
  size <- sapply(sets, length)
  sets <- sets[size > 0]

  # perform tests and do multiple testing
  if (length(sets) > 1) {
    multtest <- match.arg(multtest)
    raw.p <- pGAapprox(xx=xx, ..., test.genes=sets)
    adj.p <- p.adjust(raw.p, method=multtest)
    res <- data.frame(raw.p, adj.p)
    dimnames(res) <- list(names(sets), c("raw.p", multtest))

    # add pathway names
    res$pathway <- unlist(annotate::lookUp(rownames(res), "KEGG", "PATHID2NAME"))
    res$pathway[is.na(res$pathway)] <- ""

  } else {
    res <- pGAapprox(xx, ..., test.genes=sets)
    names(res) <- names(sets)
  }

  if (sort & !is.null(dim(res)))
    res <- res[order(res[,2]),]

  res
}



################################################################################

# function to test Broad gene sets

GABroad <- function(xx, ..., id, annotation, probe2entrez, collection,
                  category = c("c1", "c2", "c3", "c4", "c5"),
                  multtest = c("holm", "BH", "BY"), sort = TRUE) {
# xx: expression matrix (rows=genes, columns=subjects); rownames have to be gene identifiers
# ...: parameters for function 'GlobalAncova' (parameter 'method' is not available since only approximative test is performed)
# id: gene set identifiers. If omitted, tests all gene sets in the database
# annotation: name of the probe annotation package or the name of the genome wide annotation package for the species
# probe2entrez: Use only if no probe annotation package is available. A mapping from probe identifiers to entrez gene ids. May be an environment, named list or named vector.
# collection: collection of Broad gene sets (provided by 'getBroadSets()' from 'GSEABase' package)
# category: which category/ies of Broad gene sets. By default all categories are tested
# multtest: method for multiple testing correction
# sort: if TRUE, sorts the results to increasing p-values

  # get the right annotation package
  if (substr(annotation, nchar(annotation)-2, nchar(annotation)) == ".db")
    annotation <- substr(annotation, 1, nchar(annotation)-3)
  package <- paste(annotation, ".db", sep="")
  #require(package, character.only=TRUE) || stop("package ", package, " is not available")
  requireNamespace(package) || stop("package ", package, " is not available")
  requireNamespace("GSEABase") || stop("package GSEABase is not available")

  # read the file
  if (missing(collection) || ! is(collection, "GeneSetCollection"))
    stop("Please specify a \"collection\", created with the getBroadSets function")

  # Get the right categories
  if (!missing(category)) {
    pw.cat <- sapply(sapply(collection, collectionType), GSEABase::bcCategory)
    collection <- collection[pw.cat %in% category]
  }

  # Get the right sets
  if (!missing(id)) {
    collection <- collection[id]
  }

  # Map symbol identifiers to anotation-specific identifiers
  collection <- GSEABase::mapIdentifiers(collection, GSEABase::AnnotationIdentifier(annotation))
  sets <- lapply(collection, geneIds)
  names(sets) <- names(collection)

  # map to probe identifiers
  if (!missing(probe2entrez)) {
    if (is.environment(probe2entrez) || is(probe2entrez, "AnnDbBimap")) probe2entrez <- as.list(probe2entrez)
    if (is.list(probe2entrez)) probe2entrez <- unlist(probe2entrez)
    sets <- lapply(sets, function(st) {
      names(probe2entrez)[probe2entrez %in% st]
    })
  }

  # reduce to genes available in xx
  probes <- rownames(xx)
  sets <- lapply(sets, function(set) intersect(set, probes))

  # remove sets without any gene in xx
  size <- sapply(sets, length)
  sets <- sets[size > 0]

  # perform tests and do multiple testing
  if (length(sets) > 1) {
    multtest <- match.arg(multtest)
    raw.p <- pGAapprox(xx=xx, ..., test.genes=sets)
    adj.p <- p.adjust(raw.p, method=multtest)
    res <- data.frame(raw.p, adj.p)
    dimnames(res) <- list(names(sets), c("raw.p", multtest))
  } else {
    res <- pGAapprox(xx, ..., test.genes=sets)
    names(res) <- names(sets)
  }

  if (sort & !is.null(dim(res)))
    res <- res[order(res[,2]),]

  res
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
          definition = function(xx,formula.full,model.dat,test.terms,test.genes=NULL,
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

     w <- which(N.Genes <= max.group.size)
     test.genes.red <- test.genes[w]

     ew.cov <- lapply(test.genes.red, function(y) eigen(cov.shrink(t(xx2[y,,drop=FALSE]), verbose=FALSE), only.values = TRUE)$values)
     # all pairwise products of eigen values
     ew.nom   <- lapply(ew.cov, function(x) as.vector(outer(x, ew.H.nom, "*")))

     # compute approximate p-values
     p.approx <- rep(NA, N.tests)
     if(length(w) > 0)
       for(i in 1:length(w))
           p.approx[w[i]] <- .pAsymptotic(x = SS.extra[w[i]], lams = ew.nom[[i]], eps = eps, acc = acc)

     # for large gene sets use permutation test
      if(any(N.Genes > max.group.size)) {
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

        p.approx[w] <- resampleGA(xx=xx3, formula.full=formula.full, D.full=D.full, D.red=D.red,
                                  model.dat=model.dat, perm=perm, test.genes=test.genes.red,
                                  F.value=F.value, DF.full=DF.full, DF.extra=DF.extra)
     }

    return(p.approx)
}






