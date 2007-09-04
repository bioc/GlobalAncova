setGeneric("GlobalAncova", function(xx,formula.full,formula.red,model.dat,group,covars=NULL,
                            test.terms,test.genes=NULL,method=c("permutation","approx","both","Fstat"),perm=10000,max.group.size=2500,eps=1e-16,acc=50)
           standardGeneric("GlobalAncova"))
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
          definition = function(xx,formula.full,formula.red,model.dat,test.genes=NULL,
                           method=c("permutation","approx","both","Fstat"),perm=10000,max.group.size=2500,eps=1e-16,acc=50){
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
          definition = function(xx,group,covars=NULL,test.genes=NULL,
                           method=c("permutation","approx","both","Fstat"),perm=10000,max.group.size=2500,eps=1e-16,acc=50){                           
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
          definition = function(xx,formula.full,test.terms,model.dat,test.genes=NULL,
                           method=c("permutation","approx","both","Fstat"),perm=10000,max.group.size=2500,eps=1e-16,acc=50)
{
  # test for model.dat
  if(!is.data.frame(model.dat))
    stop("'model.dat' has to be a data frame")

  # test for 'test.terms'
  terms.all <- test.terms
  D.full    <- model.matrix(formula.full, model.dat)
  terms.all <- colnames(D.full)

  # are all terms variables compatible with 'model.dat'?
  if(!all(test.terms %in% terms.all))
    stop("'test.terms' is not compatible with the specified models")

  D.full <- model.matrix(formula.full, data=model.dat)
  D.red  <- D.full[,!(colnames(D.full) %in% test.terms), drop=FALSE]

  # then apply the usual function
  expr.test(xx=xx,formula.full=formula.full,D.red=D.red,
                model.dat=model.dat,test.genes=test.genes,method=method,perm=perm,max.group.size=max.group.size,eps=eps,acc=acc)
}
)


################################################################################
################################################################################

# main function of GlobalAncova
expr.test <- function(xx,formula.full,formula.red=NULL,D.red=NULL,model.dat,test.genes=NULL,
                     method=c("permutation","approx","both","Fstat"),perm=10000,max.group.size=2500,eps=1e-16,acc=50){             
    # if just one gene should be tested
    if(is.vector(xx))
        xx <- t(as.matrix(xx))

    if(is.null(rownames(xx)))
        rownames(xx) <- 1:nrow(xx)
    if(is.null(test.genes))
        test.genes <- list(rownames(xx))
    if(!is.list(test.genes))
        test.genes <- list(test.genes)
    if(is.numeric(unlist(test.genes)))
        test.genes <- lapply(test.genes, as.character)

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
    #p.value <- pf(F.value, DF.extra, DF.full, lower.tail=F)

    # permutation p-values
    p.value <- p.perm <- NULL
    if(method == "permutation" | method == "both") {
        p.perm <- resampleGA(xx2, D.full, D.red, perm, test.genes, F.value, DF.full, DF.extra)
        p.value <- cbind(p.value, p.perm)
    }

    # asymptotic p-values
    if(method == "approx" | method == "both"){
        # compute eigen values of (H.full-H.red) and XX'
        require(corpcor)
        ew.H.nom <- eigen(hat.matrix(D.full) - hat.matrix(D.red))$values

        # remove gene sets that are bigger than 'max.group.size' -> treat those with permutation test
        toobig <- sum(N.Genes > max.group.size) > 0
        if(toobig)
          warning("p-values of gene sets bigger than 'max.group.size' are calculated permutation-based")
          
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
        p.value <- cbind(p.value, p.approx)
        
        # make permutation test for large gene sets
        if(is.null(p.perm) & toobig) {
           w <- which(N.Genes > max.group.size)
           test.genes.red <- test.genes[w]
           p.perm <- rep(NA, N.tests)
           p.perm[w] <- resampleGA(xx2, D.full, D.red, perm, test.genes.red, F.value[w], DF.full[w], DF.extra[w])
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

# conducts the permutation test
resampleGA <- function(xx, D.full, D.red, perm, test.genes, F.value, DF.full, DF.extra){
    N.Subjects  <- ncol(xx)
    rr     <- row.orth2d(xx, D.red)

    # sum of squares in reduced model do not have to be re-calculated in each permutation
    genewSS  <- genewiseGA(xx, D.full, D.red=D.red)
    SS.red.i <- genewSS[,1] + genewSS[,2]   # SS.red = SS.extra + SS.full

    D.full.perm <- D.full
    test.col <- !colnames(D.full) %in% colnames(D.red)
    count <- numeric(length(test.genes))
    for(i in 1:perm) {
      # permute only values of interesting variables
      ord <- sample(N.Subjects)
  		D.full.perm[,test.col] <- D.full[ord, test.col]	 		 
      genewSS.perm <- genewiseGA(rr, D.full.perm, SS.red.i=SS.red.i) 
      F.perm       <- sapply(test.genes, function(x) sum(genewSS.perm[x,1]) / sum(genewSS.perm[x,2])) / (DF.extra / DF.full)
      count        <- count + (F.perm > F.value)
    }
    
    return(count / perm)
}


################################################################################

# calculates the asymptotic p-value using the method of Robbins and Pitman (1949)
.pAsymptotic <- function(x, lams, eps, acc) {
# x: quantile
# lams: vector of eigenvalues
# eps: accuracy
# acc: accuracy for removing small eigenvalues

  lams <- .weed(lams = lams, accuracy = acc)
  lams <- sort(lams, decreasing=TRUE)
  m <- length(lams)
  if (m == 0) 
    p <- NA
  else {
    bet <- min(lams)

    # get an upper bound to the number of iterations needed
    Q2 <- qnorm(eps)^2
    maxiter <- trunc(0.5 * (x/bet + Q2 + sqrt(2*Q2*x/bet + Q2*Q2) - m))

    # starting values
    d <- numeric(maxiter)
    c <- numeric(maxiter + 1)
    c[1] <- prod(sqrt(bet / lams))
    restc <- 1 - c[1]
    chi <- pchisq(x / bet, df = m, lower.tail = FALSE)
    partialsum <- c[1] * chi
    dbase <- (1 - bet / lams)
    ready <- FALSE
    ix <- 1

    # iterate!
    while (!ready) {
      d[ix] <- 0.5 * sum(dbase^ix)
      c[ix+1] <- mean(c[1:ix] * d[ix:1])
      if (restc > 100 * .Machine$double.neg.eps) {
        restc <- restc - c[ix+1]
      } else {
        restc <- c[ix+1] * lams[1] / lams[m]
      }
      chi <- pchisq(x / bet, df = m + 2 * ix, lower.tail = FALSE)
      partialsum <- partialsum + c[ix+1] * chi
      error <- restc * pchisq(x / bet, df = m + 2 * ix + 2, lower.tail = TRUE)
      ready <- (error < eps) || (ix == maxiter) || (error / partialsum < 10^-4)
      ix <- ix + 1
    }
    p <- partialsum + restc
    if (p < eps) { p <- 0 }
  }
  p
}

################################################################################

# Removes extremely small eigenvalues
.weed <- function(lams, accuracy) {
# lams: vector of eigenvalues
  if (missing(accuracy)) {
    thresh <- 1/50
  } else {
    thresh <- 1/accuracy
  }
  lams <- -sort(-lams)
  m <- length(lams)
  if (lams[1] == 0) {
    lams <- numeric(0)
    m <- 0
  } else {
    while (lams[m] / lams[1] < thresh) {
      q <- m-1
      r <- m-2
      lams[q] <- lams[q] + lams[m]
      while ((r > 0) && (lams[r] < lams[q])) {
        lams[r:q] <- mean(lams[r:q])
        r <- r - 1
      }
      m <- q
    }
    lams <- lams[1:m]
  }
  lams
}


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
       if(identical(id[,i],id[,j]) & i!=j)
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




