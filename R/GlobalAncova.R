# !! test.genes nicht auf NULL
setGeneric("GlobalAncova", function(xx,formula.full,formula.red,model.dat,group,covars=NULL,
                            test.terms,test.genes,method=c("permutation","approx","both","Fstat"),perm=10000,max.group.size=2500,eps=1e-16,acc=50)
           standardGeneric("GlobalAncova"))
# !!
# this function computes a Global Ancova which tests for differential gene expression
# xx: expression matrix (rows=genes, columns=subjects)
# formula.full: model formula for the full model
# formula.red: model formula for the reduced model
# model.dat: data frame that contains the group and covariable information
# test.genes: vector of probeset names or a list where each element is a vector of
#             probeset names (e.g. for testing many pathways simultaneously)
# method: "permutation"/"approx"/"both"/"Fstat" - permutation test/approximative test/both/only calculation of F-statistics
# perm: number of permutations
# max.group.size: if a gene set is larger than max.group.size the permutation test is applied even if method="approx"
# eps: resuolution of approximation
# acc: additional accuracy parameter for approximation (the higher the value the higher the accuracy) 


################################# general function #############################

setMethod("GlobalAncova", signature(xx="matrix",formula.full="formula",formula.red="formula",
                           model.dat="ANY",group="missing",covars="missing",test.terms="missing"),
# !!
          definition = function(xx,formula.full,formula.red,model.dat,test.genes,
                           method=c("permutation","approx","both","Fstat"),perm=10000,max.group.size=2500,eps=1e-16,acc=50){
# !!

  # test for model.dat
  if(!is.data.frame(model.dat))
    stop("'model.dat' has to be a data frame")

  expr.test(xx=xx,formula.full=formula.full,formula.red=formula.red,model.dat=model.dat,
            test.genes=test.genes,method=method,perm=perm,max.group.size=max.group.size,eps=eps,acc=acc)
}
)


########################## function for 2 groups ###############################

setMethod("GlobalAncova", signature(xx="matrix",formula.full="missing",formula.red="missing",
                           model.dat="missing",group="ANY",covars="ANY",test.terms="missing"),
# !!
          definition = function(xx,group,covars=NULL,test.genes,
                           method=c("permutation","approx","both","Fstat"),perm=10000,max.group.size=2500,eps=1e-16,acc=50){                           
# !!

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
  expr.test(xx=xx,formula.full=formula.full,formula.red=formula.red,
                model.dat=model.dat,test.genes=test.genes,method=method,perm=perm,max.group.size=max.group.size,eps=eps,acc=acc)
}
)


############################# with 'test.terms' ################################

setMethod("GlobalAncova", signature(xx="matrix",formula.full="formula",formula.red="missing",
                           model.dat="ANY",group="missing",covars="missing",test.terms="character"),
# !!  Achtung: test.terms muß NACH model.dat kommen, so wie in d. Signatur !!
          definition = function(xx,formula.full,model.dat,test.terms,test.genes,
                           method=c("permutation","approx","both","Fstat"),perm=10000,max.group.size=2500,eps=1e-16,acc=50){
# !!

  # test for model.dat
  if(!is.data.frame(model.dat))
    stop("'model.dat' has to be a data frame")

  # test for 'test.terms'
  D.full    <- model.matrix(formula.full, data=model.dat)
  terms.all <- colnames(D.full)

  # are all terms variables compatible with 'model.dat'?
  if(!all(test.terms %in% terms.all))
    stop("'test.terms' is not compatible with the specified models")

  D.red  <- D.full[,!(colnames(D.full) %in% test.terms), drop=FALSE]

  # then apply the usual function
  expr.test(xx=xx,formula.full=formula.full,D.red=D.red,
                model.dat=model.dat,test.genes=test.genes,method=method,perm=perm,max.group.size=max.group.size,eps=eps,acc=acc)
}
)


################################################################################
################################################################################

# main function of GlobalAncova
# !! test.genes nicht mehr auf NULL setzen zwecks Fehlerabfrage
expr.test <- function(xx,formula.full,formula.red=NULL,D.red=NULL,model.dat,test.genes,
                     method=c("permutation","approx","both","Fstat"),perm=10000,max.group.size=2500,eps=1e-16,acc=50){ 
                           
    # check if formula variables correspond to columns of 'model.dat'
    design.terms <- as.character(terms(formula.full)@variables)[-1]
    if(!all(design.terms %in% names(model.dat)))
      stop("terms in 'formula.full' do not match variables in 'model.dat'")  

    # check if dimensions of expression and phenotype data match
    if(ncol(xx) != nrow(model.dat))
      stop("number of samples in expression matrix 'xx' differs from number of samples in 'model.dat'")

    # check 'test.genes'
    if(!(missing(test.genes) || is.vector(test.genes)))
      #stop("'test.genes' does not define valid gene sets")
      stop("'test.genes' should be a vector or list")
# !!                  
                              
    # if just one gene should be tested
    if(is.vector(xx))
        xx <- t(as.matrix(xx))

    if(is.null(rownames(xx)))
        rownames(xx) <- 1:nrow(xx)
# !!
    if(missing(test.genes))
        test.genes <- list(rownames(xx))
# !!
    if(!is.list(test.genes))
        test.genes <- list(test.genes)
    if(is.numeric(unlist(test.genes)))
        test.genes <- lapply(test.genes, as.character)
# !!
    if(!all(unlist(test.genes) %in% rownames(xx)))
      stop("gene names in 'test.genes' do not correspond to gene names in 'xx'")
# !!

    xx2        <- xx[unique(unlist(test.genes)),,drop=F]
    N.Genes    <- sapply(test.genes, length)
    N.Subjects <- ncol(xx2)
    N.tests    <- length(test.genes)

    # design matrices
    D.full     <- model.matrix(formula.full, data=model.dat)
    if(is.null(D.red))
        D.red  <- model.matrix(formula.red,  data=model.dat)
        
# !! 
    # check if reduced model is included in full model
    terms.full <- colnames(D.full)
    terms.red <- colnames(D.red)
    if(!all(terms.red %in% terms.full))
      #stop("the reduced model is not part of the full model")
      stop("full model and reduced model are not nested")
# !!

    N.par.full <- ncol(D.full)
    N.par.red  <- ncol(D.red)

    # checking for NA's
    if(nrow(D.full) < N.Subjects)
        stop("Missing values in the model variables")

    # degrees of freedom
    DF.full  <- N.Genes * (N.Subjects - N.par.full)
    DF.extra <- N.Genes * (N.par.full - N.par.red)

    # sums of squares
    genewiseSS <- genewiseGA(xx2, D.full, D.red=D.red)
    SS.full    <- sapply(test.genes, function(x) sum(genewiseSS[x,"denominator"]))
    SS.extra   <- sapply(test.genes, function(x) sum(genewiseSS[x,"nominator"]))
    MS.full    <- SS.full / DF.full
    MS.extra   <- SS.extra / DF.extra

    # F statistic and theoretical p-value for each gene group
    method  <- match.arg(method)
    F.value <- MS.extra / MS.full
    if(method == "Fstat") {
      test.result <- cbind(genes=N.Genes, F.value=F.value)
      return(test.result)
    }

    # permutation p-values
    p.value <- p.perm <- NULL
    # !!
    if(method == "permutation" || method == "both") {     # einfaches | könnte zu logischem Vektor führen -> Warnung
    # !!
# !!
        p.perm <- resampleGA(xx2, formula.full, D.full, D.red, model.dat, perm, test.genes, F.value, DF.full, DF.extra)
# !!
        p.value <- cbind(p.value, p.perm)
    }

    # asymptotic p-values
    # !!
    if(method == "approx" || method == "both"){
    # !!
        # compute eigen values of (H.full-H.red) and XX'
        require(corpcor)
        ew.H.nom <- eigen(hat.matrix(D.full) - hat.matrix(D.red))$values

        # remove gene sets that are bigger than 'max.group.size' -> treat those with permutation test
        toobig <- sum(N.Genes > max.group.size) > 0
        if(toobig)
          warning("p-values of gene sets bigger than 'max.group.size' are calculated permutation-based")
          
        w <- which(N.Genes <= max.group.size)
        test.genes.red <- test.genes[w]    
        
# !!
        ew.cov <- lapply(test.genes.red, function(y) eigen(cov.shrink(t(xx2[y,,drop=FALSE]), verbose=FALSE), only.values = TRUE)$values)
        # all pairwise products of eigen values
        ew.nom   <- lapply(ew.cov, function(x) as.vector(outer(x, ew.H.nom, "*")))
# !!

        # compute approximate p-values
        p.approx <- rep(NA, N.tests)
        if(length(w) > 0)
          for(i in 1:length(w))
              p.approx[w[i]] <- .pAsymptotic(x = SS.extra[w[i]], lams = ew.nom[[i]], eps = eps, acc = acc)
        p.value <- cbind(p.value, p.approx)
        
        # make permutation test for large gene sets
        # !!
        if(is.null(p.perm) && toobig) {        # einfaches & könnte zu logischem Vektor führen -> Warnung
        # !!
           w <- which(N.Genes > max.group.size)
           test.genes.red <- test.genes[w]
           p.perm <- rep(NA, N.tests)
# !!
           p.perm[w] <- resampleGA(xx2, formula.full, D.full, D.red, model.dat, perm, test.genes.red, F.value[w], DF.full[w], DF.extra[w])
# !!
           p.value <- cbind(p.value, p.perm)
        }
    }
    # results
    test.result <- cbind(F.value, p.value)

    if(N.tests == 1){
        term.names          <- effectnames(D.full, D.red)
        ANOVA.tab           <- cbind(c(SS.extra,SS.full), c(DF.extra,DF.full), c(MS.extra,MS.full))
        dimnames(ANOVA.tab) <- list(c("Effect", "Error"), c("SSQ", "DF", "MS"))
        result <- list("effect"=term.names,"ANOVA"=ANOVA.tab,"test.result"=t(test.result),"terms"=colnames(D.full))
    }

    else
        result <- cbind("genes"=N.Genes, test.result)

    return(result)
}


################################################################################

# builds the hat matrix
hat.matrix <- function(x)
  x %*% solve(t(x) %*% x) %*% t(x)

# computes residuals for given response 'xx' and design matrix D
row.orth2d <- function(xx, D)
  xx %*% (diag(dim(D)[1]) - hat.matrix(D))


################################################################################

# computes nominator and denominator of GlobalAncova statistic for each single gene
genewiseGA <- function(xx, D.full, D.red=NULL, SS.red.i=NULL){  
    R.full <- row.orth2d(xx, D.full)
    if(is.null(SS.red.i)){
      R.red  <- row.orth2d(xx, D.red) 
      SS.red.i <- rowSums(R.red*R.red)
    }

    # denominator: residual sum of squares in the full model
    SS.full.i <- rowSums(R.full*R.full)

    # nominator: extra residual sum of squares
    SS.extra.i <- SS.red.i - SS.full.i
    
    return(cbind(nominator=SS.extra.i, denominator=SS.full.i))
}


################################################################################

# !!
# Funktion f. d. Permutationstest jetzt in permutation.R
# !!

################################################################################

# !!
# Funktionen f. d. asymptotischen p-Werte jetzt in 'approximation.R'
# !!

################################################################################

# extracts the name of the tested effect
effectnames <- function(D.full, D.red){
# (because there can be problems with interactions:
#  e.g '~group*covar' yields 'group,covar,group:covar' but ~covar*group yields '..,covar:group')
 names.all  <- union(colnames(D.full), colnames(D.red))
 no.effect  <- intersect(colnames(D.full), colnames(D.red))

 # are there interaction effects?
 interact   <- grep(":", names.all)
 if(length(interact) > 0)
 {
   # split the interaction terms and look if components are the same
   split.names <- strsplit(names.all[interact], ":")
   id          <- sapply(split.names, sort)
   for(i in 1:length(split.names))
   {
     for(j in 1:length(split.names))
     {
       # if there are identical columns, these represent the same factor and hence
       #  are added to the 'non-effect-terms'
       # !!
       if(identical(id[,i],id[,j]) && i!=j)
       # !!
         no.effect <- c(no.effect, names.all[interact[c(i,j)]])
     }
   }
 }
 # the remaining terms are the tested effects
 effect.names <- names.all[!(names.all %in% no.effect)]
 return(effect.names)
}


################################################################################

# derives 'formula.full', 'formula.red' and 'model.dat' out of 'group' and 'covars'
group2formula <- function(group, group.name, covars, covar.names){
# group: group variable
# group.name: character: name of group variable
# covars: covariate information
# covar.names: character: names of covraiates

  # model matrix
  if(is.null(covars)){
    model.dat <- data.frame(group)
    names(model.dat)[1] <- group.name
  }
  else{
    model.dat <- data.frame(group, covars)
    names(model.dat) <- c(group.name, covar.names)
  }

  # model formulas
  formula.full <- paste("~", group.name)

  # if there are no covariates
  if(is.null(covars)){
    formula.full <- as.formula(formula.full)
    formula.red  <- ~ 1
  }
  else{
    formula.full <- as.formula(paste(formula.full, "+", paste(covar.names, collapse="+")))
    formula.red  <- as.formula(paste("~", paste(covar.names, collapse="+")))
  }
  
  return(list(formula.full=formula.full, formula.red=formula.red, model.dat=model.dat))
}



