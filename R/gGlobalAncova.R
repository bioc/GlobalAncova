
######################################
### univariate test statistics
######################################

### categorical response 

## calculate deviances; taken from glm.fit

getdev <- function(x, y, weights=rep(1, nobs)){
  
  nobs <- NROW(y)
  #  if(missing(weights))
  #    weights <- rep(1, nobs)
  
  offset = rep(0, nobs)
  
  # settings for binomial family
  family <- binomial()
  dev.resids <- family$dev.resids
  variance <- family$variance
  linkinv <- family$linkinv
  #mustart
  eval(family$initialize)
  eta <- family$linkfun(mustart)
  mu.eta <- family$mu.eta
  
  # further settings
  start <- NULL
  coefold <- start
  
  # deviance
  mu <- linkinv(eta)
  devold <- sum(dev.resids(y, mu, weights))
  
  unless.null <- function(x, if.null) if (is.null(x)) 
    if.null else x
  valideta <- unless.null(family$valideta, function(eta) TRUE)
  validmu <- unless.null(family$validmu, function(mu) TRUE)
  
  # iterations
  control <- do.call("glm.control", list())
  for (iter in 1L:control$maxit) {
    good <- weights > 0
    varmu <- variance(mu)[good]
    if (anyNA(varmu)) 
      stop("NAs in V(mu)")
    if (any(varmu == 0)) 
      stop("0s in V(mu)")
    mu.eta.val <- mu.eta(eta)
    if (any(is.na(mu.eta.val[good]))) 
      stop("NAs in d(mu)/d(eta)")
    good <- (weights > 0) & (mu.eta.val != 0)
    if (all(!good)) {
      conv <- FALSE
      warning(gettextf("no observations informative at iteration %d", iter), domain = NA)
      break
    }
    z <- (eta - offset)[good] + (y - mu)[good]/mu.eta.val[good]
    w <- sqrt((weights[good] * mu.eta.val[good]^2)/variance(mu)[good])
    fit <- .Call(stats:::C_Cdqrls, x[good, , drop = FALSE] * w, z * w, min(1e-07, control$epsilon/1000), check = FALSE)
    if (any(!is.finite(fit$coefficients))) {
      conv <- FALSE
      warning(gettextf("non-finite coefficients at iteration %d", iter), domain = NA)
      break
    }
    if (nobs < fit$rank) 
      stop(sprintf(ngettext(nobs, "X matrix has rank %d, but only %d observation", 
                            "X matrix has rank %d, but only %d observations"), fit$rank, nobs), domain = NA)
    start[fit$pivot] <- fit$coefficients
    eta <- drop(x %*% start)
    mu <- linkinv(eta <- eta + offset)
    dev <- sum(dev.resids(y, mu, weights))
    if (control$trace) 
      cat("Deviance = ", dev, " Iterations - ", iter, "\n", sep = "")
    boundary <- FALSE
    if (!is.finite(dev)) {
      if (is.null(coefold)) 
        stop("no valid set of coefficients has been found: please supply starting values", call. = FALSE)
      warning("step size truncated due to divergence", call. = FALSE)
      ii <- 1
      while (!is.finite(dev)) {
        if (ii > control$maxit) 
          stop("inner loop 1; cannot correct step size", call. = FALSE)
        ii <- ii + 1
        start <- (start + coefold)/2
        eta <- drop(x %*% start)
        mu <- linkinv(eta <- eta + offset)
        dev <- sum(dev.resids(y, mu, weights))
      }
      boundary <- TRUE
      if (control$trace) 
        cat("Step halved: new deviance = ", dev, "\n", sep = "")
    }
    if (!(valideta(eta) && validmu(mu))) {
      if (is.null(coefold)) 
        stop("no valid set of coefficients has been found: please supply starting values", call. = FALSE)
      warning("step size truncated: out of bounds", call. = FALSE)
      ii <- 1
      while (!(valideta(eta) && validmu(mu))) {
        if (ii > control$maxit) 
          stop("inner loop 2; cannot correct step size", call. = FALSE)
        ii <- ii + 1
        start <- (start + coefold)/2
        eta <- drop(x %*% start)
        mu <- linkinv(eta <- eta + offset)
      }
      boundary <- TRUE
      dev <- sum(dev.resids(y, mu, weights))
      if (control$trace) 
        cat("Step halved: new deviance = ", dev, "\n", sep = "")
    }
    if (abs(dev - devold)/(0.1 + abs(dev)) < control$epsilon) {
      conv <- TRUE
      coef <- start
      break
    }
    else {
      devold <- dev
      coef <- coefold <- start
    }
  }
  
  # null deviance
  #wtdmu <- sum(weights * y) / sum(weights)
  #nulldev <- sum(dev.resids(y, wtdmu, weights))
  
  #return(list(dev=dev, nulldev=nulldev))
  return(dev)
}



### categorical/ordinal response - general deviance test

devtest <- function(y, D.full, D.red, family=c("binomial", "multinomial", "propodds")){
  require(VGAM)
  OK <- !is.na(y)
  y <- y[OK]
  D.full <- D.full[OK,,drop=FALSE]
  D.red <- D.red[OK,,drop=FALSE]
  
  # binomial
  family <- match.arg(family)
  if(family == "binomial"){
    dev.full <- getdev(x=D.full, y=y)
    dev.red  <- getdev(x=D.red, y=y)
    statistic <- dev.red - dev.full
  }
  
  # multinomial or ordinal
  else{
    # get model formulas
    int <- "(Intercept)" %in% colnames(D.full)
    # with intercept
    if(int){
      covars.full <- paste(colnames(D.full)[-1], collapse="+")
      covars.red  <- paste(c("1", colnames(D.red)[-1]), collapse="+")
    }
    # without
    else{
      covars.full <- paste(c("-1", colnames(D.full)), collapse="+")
      covars.red  <- paste(c("-1", colnames(D.red)), collapse="+")
    }
    
    if("y" %in% covars.full)
      covars.full[covars.full == "y"] <- "y.1"
    
    data <- data.frame(y=y, D.full)
    form.full  <- formula(paste("y ~", covars.full))
    form.red   <- formula(paste("y ~", covars.red))
    model.full <- vglm(form.full, family=family, data=data)
    model.red  <- vglm(form.red, family=family, data=data)
    statistic  <- deviance(model.red) - deviance(model.full)
    
    # multinomial: get corresponding statistic for chi^2 distribution w. G-1 df (with G = number of categories of tested variable)
    # (relevant for globaltest w. mixed data, but also in case not all y's have observations in all categories...)
    if(family == "multinomial"){
      K <- length(na.omit(unique(y)))
      df <- model.red@df.residual - model.full@df.residual
      p <- pchisq(statistic, df=df)
      statistic <- qchisq(p, df=df / (K-1))
    }
  }
  
  return(statistic)
}  



## equivalent alternative when both regressor and dependent variables are categorical: G^2 statistic
# + transformation to chi^2 distribution w. G-1 df (with G = number of categories of tested variable)

Gsquared <- function(y, D.full, D.red){
  OK <- !is.na(y)
  y <- factor(y[OK])
  
  # group variable
  testcol <- setdiff(colnames(D.full), colnames(D.red))
  # ..binary
  if(length(testcol) == 1)
    x <- D.full[OK, testcol]
  # ..multicategorical
  else
    x <- colSums(t(D.full[OK, testcol]) * 1:length(testcol))
  
  # contingency table  
  tab <- table(y, x)
  K <- nrow(tab)
  L <- ncol(tab)
  
  # if there are 0s in the table, use glm (otherwise G2 will be NaN)
  if(any(tab == 0)){
    if(K == 2){
      test <- glm(y ~ x, family="binomial")
      G2 <- test$null.deviance - test$deviance
    }
    
    # in case x is multinomial
    else if(K > 2){
      require(VGAM)
      model.full <- vglm(y ~ x, family=multinomial)
      model.red <- vglm(y ~ 1, family=multinomial)
      G2 <- deviance(model.red) - deviance(model.full)
    }
  }
  
  # else, calculate it by standard formula
  else{
    n <- sum(tab)
    sr <- rowSums(tab)
    sc <- colSums(tab)
    E <- outer(sr, sc, "*")/n
    G2 <- 2 * sum(tab * log(tab / E))
  }
  
  # transform to chi^2 distribution w. G-1 df (with G = number of categories of tested variable)
  # (relevant for globaltest w. mixed data, but also in case not all y's have observations in all categories...)
  if(K > 2){
    df <- (K-1) * (L-1)
    p <- pchisq(G2, df=df)
    G2 <- qchisq(p, df=df / (K-1))
  }
  
  return(G2)
}



### continuous response - reduction in sum of squares (GlobalAncova nominator)

hat.matrix <- function(x)
  x %*% solve(t(x) %*% x) %*% t(x)

row.orth2d <- function(xx, D)
  xx %*% (diag(dim(D)[1]) - hat.matrix(D))


SSred <- function(y, D.full, D.red){
  OK <- !is.na(y)
  y <- y[OK]
  D.full <- D.full[OK,,drop=FALSE]
  D.red <- D.red[OK,,drop=FALSE]
  
  R.full <- row.orth2d(y, D.full)
  R.red  <- row.orth2d(y, D.red) 
  SS.red <- sum(R.red * R.red)
  SS.full <- sum(R.full * R.full)
  
  statistic <- SS.red - SS.full
  return(statistic)
}



#########################################
### global permutation test
#########################################

## all permutation statistics for all variables - when regressor and outcome are both categorical (for general kxl tables)
# + transformation to chi^2 distribution w. G-1 df (with G = number of categories of tested variable)
Tperm.catcat <- function(x, Perms){
  # x: data.frame (or matrix) of categorical variables
  # Perms: matrix of permutations of group variable
  
  # code x and Perms numeric from 0 to k-1/l-1
  x     <- t(apply(x, 2, function(b) as.numeric(factor(b))))
  Perms <- apply(Perms, 2, function(b) as.numeric(factor(b)))
  
  # total counts
  n <- ncol(x)
  
  # numbers of categories -1
  K <- apply(x, 1, max)  # outcome variables could have different numbers of categories
  k <- max(K)     
  l <- max(Perms) # just one group variable
  
  # marginal counts for group variable
  n.cols <- table(Perms[,1])
  
  G2 <- matrix(0, nrow=nrow(x), ncol=ncol(Perms))
  for(i in 1:k){
    xi <- x == i
    
    # marginal count for outcome variable, category i
    n.rows <- rowSums(xi)
    
    # remove variables that do not have category i
    zero <- n.rows == 0
    xi <- xi[!zero,]
    
    for(j in 1:l){
      Permsj <- Perms == j
      
      # counts for crosstable entry i.j for all variables, for all permutations
      n.ij <- xi %*% Permsj
      
      # expected counts
      e.ij <- n.cols[j] * n.rows[!zero] / n
      
      # G^2 summand
      G2[!zero,] <- G2[!zero,] + n.ij * log(n.ij / e.ij) 
    }
  }
  
  G2 <- 2 * G2
  
  # set permutations corresponding to crosstables with 0 entries (for any variable) to NA
  na <- apply(G2, 2, function(x) any(is.na(x)))
  #G2 <- G2[, !na, drop=FALSE]
  G2[, na] <- NA
  
  # transform to chi^2 distribution w. G-1 df (with G = number of categories of tested variable)
  # (relevant for globaltest w. mixed data, but also in case not all x's have observations in all categories...)
  if(k > 2){
    for(i in 1:nrow(x)){
      if(K[i] > 2){
        df <- (K[i]-1) * (l-1)
        p <- pchisq(G2[i,], df=df)
        G2[i,] <- qchisq(p, df=df / (K[i]-1))
      }
    }
  }
  
  rownames(G2) <- rownames(x)
  return(G2)
}



## univariate observed and permutation test statistics, also for mixed data

gGAteststats <- function(data, formula.full, formula.red=~1, model.dat, perm=10000){
  # data: data.frame (columns=variables)
  # formula.full, formula.red: models to be compared
  # model.dat: data.frame of regressors, containing variables specified in formula.full and formula.red
  # perm: number of permutations; if set to 0, only observed statistics are computed
  
  if(is.matrix(data)){
    data0 <- data
    if(is.null(colnames(data0)))
      colnames(data0) <- 1:ncol(data0)
    data <- data.frame(data)
    names(data) <- colnames(data0)  # prevent adding 'X' at variable names in case input data had no colnames
  }
  
  # remove missing values in regressors
  design.terms <- as.character(attr(terms(formula.full), "variables"))[-1]
  test.dat <- model.dat[,design.terms]
  OK <- complete.cases(test.dat)
  data <- data[OK,,drop=FALSE]
  model.dat <- model.dat[OK,,drop=FALSE]
  test.dat  <- model.dat[,design.terms]
  n <- nrow(model.dat)
  
  # remove outcome variables with only one observed level
  nlevels <- apply(data, 2, function(x) length(na.omit(unique(x))))
  if(any(nlevels < 2)){
    w <- which(nlevels < 2)
    data <- data[,-w]
    warning(paste(length(w), "variables removed with only one level"))
  }
  
  # model matrices
  D.full <- model.matrix(formula.full, data=model.dat)
  D.red  <- model.matrix(formula.red, data=model.dat)
  
  # model columns to be tested
  testcol <- which(!(colnames(D.full) %in% colnames(D.red)))
  
  # permutations
  if(perm > 0){
    nperm <- .nPerms(D.full, model.dat, formula.full)
    faccounts <- nperm$counts
    nperm <- nperm$nPerms
  
    # all permutations if possible
    if(nperm < perm){
      print(paste("enumerating all", nperm, "permutations"))
      if(is.null(faccounts))  # design with continuous covariates
        Perms <- .allperms(1:n)
      else
        Perms <- .allpermsG(faccounts, faccounts) # factorial design
      Perms <- apply(Perms, 2, order)  # get possible orderings
    }
  
    # random permutations otherwise
    else
      Perms <- replicate(perm, sample(n))
  }
  
  # group variables according to data type
  dc <- sapply(data, data.class)
  
  comments <- length(unique(dc)) > 1 & perm > 0
  
  # get binary variables
  bin <- sapply(data, function(x) length(na.omit(unique(x))) == 2)
  dc[bin] <- "binary"
  
  Tobs.bin <- Tobs.cat  <- Tobs.ord  <- Tobs.quant  <- NULL
  if(perm > 0)
    Tperm.bin <- Tperm.cat <- Tperm.ord <- Tperm.quant <- NULL
  
  ## continuous variables
  if(any(dc == "numeric")){
    if(comments)
      print("testing quantitative variables")
    xx <- data[,dc == "numeric", drop=FALSE]
    
    # observed univariate statistics
    Tobs.quant <- sapply(xx, SSred, D.full, D.red)
    
    # univariate permutation statistics
    if(perm > 0)
      Tperm.quant <- apply(Perms, 2, function(x) 
        sapply(xx, SSred, D.full=cbind(D.full[,-testcol, drop=FALSE], D.full[x, testcol, drop=FALSE]), D.red=D.red))
  }
  
  ## ordinal variables
  if(any(dc == "ordered")){
    if(comments)
      print("testing ordinal variables")
    xx <- data[,dc == "ordered", drop=FALSE]
    
    # observed univariate statistics
    Tobs.ord <- sapply(xx, devtest, D.full, D.red, family="propodds")
    
    # univariate permutation statistics
    if(perm > 0)
      Tperm.ord <- apply(Perms, 2, function(x) 
        sapply(xx, devtest, D.full=cbind(D.full[,-testcol, drop=FALSE], D.full[x, testcol, drop=FALSE]), D.red=D.red, family="propodds"))
  }
  
  # check if tested term is single categorical variable
  cat.group <- length(design.terms) == 1 & (data.class(test.dat) != "numeric" | length(unique(test.dat)) == 2)
  
  ## binary variables
  if(any(dc == "binary")){
    if(comments)
      print("testing binary variables")
    xx <- data[,dc == "binary", drop=FALSE]
    
    # if tested term is single categorical variable, use G^2 statistic
    if(cat.group){
      # observed univariate statistics
      Tobs.bin <- sapply(xx, Gsquared, D.full, D.red)
      
      # univariate permutation statistics
      if(length(testcol) == 1)
        group <- D.full[,testcol]
      else
        group <- colSums(t(D.full[,testcol]) * 1:length(testcol))
      
      if(perm > 0){
        Perms2 <- apply(Perms, 2, function(x) group[x])
        Tperm.bin <- Tperm.catcat(xx, Perms2)
      }
    }
    
    # in case of covariables etc., use deviance test
    else{
      Tobs.bin <- sapply(xx, devtest, D.full, D.red, family="binomial")
      if(perm > 0)
        Tperm.bin <- apply(Perms, 2, function(x) 
          sapply(xx, devtest, D.full=cbind(D.full[,-testcol, drop=FALSE], D.full[x, testcol, drop=FALSE]), D.red=D.red, family="binomial"))
    }
  }
  
  ## multi-categorical variables
  if(any(dc == "factor")){
    if(comments)
      print("testing multi-categorical variables")
    xx <- data[,dc == "factor", drop=FALSE]
    
    # if tested term is single categorical variable, use G^2 statistic
    if(cat.group){
      # observed univariate statistics
      Tobs.cat <- sapply(xx, Gsquared, D.full, D.red)
      
      # univariate permutation statistics
      if(length(testcol) == 1)
        group <- D.full[,testcol]
      else
        group <- colSums(t(D.full[,testcol]) * 1:length(testcol))
      
      if(perm > 0){
        Perms2 <- apply(Perms, 2, function(x) group[x])
        Tperm.cat <- Tperm.catcat(xx, Perms2)
      }
    }
    
    # in case of covariables etc., use deviance test
    else{
      Tobs.bin <- sapply(xx, devtest, D.full, D.red, family="multinomial")
      if(perm > 0)
        Tperm.bin <- apply(Perms, 2, function(x) 
          sapply(xx, devtest, D.full=cbind(D.full[,-testcol, drop=FALSE], D.full[x, testcol, drop=FALSE]), D.red=D.red, family="multinomial"))
    }
  }
  
  # put together
  Tobs <- c(Tobs.quant, Tobs.ord, Tobs.bin, Tobs.cat)
  if(perm > 0){
    Tperm <- rbind(Tperm.quant, Tperm.ord, Tperm.bin, Tperm.cat)
    rownames(Tperm) <- sub("Tperm.", "", rownames(Tperm))  # in case there is only 1 variable for a data type, 'Tperm.' will be added to the variable name...
  }
  
  # bring in original order
  Tobs <- Tobs[names(data)]
  if(perm > 0)
    Tperm <- Tperm[names(data),]
  
  if(perm > 0)
    return(list(Tobs=Tobs, Tperm=Tperm))
  else
    return(Tobs)
}



## 'generalized GlobalAncova' group test, also for mixed data

gGlobalAncova <- function(data, formula.full, formula.red=~1, model.dat, Sets, sumstat=sum, perm=10000){
  # data: data.frame of variables to be tested in sets (columns=variables)
  # formula.full, formula.red: models to be compared
  # model.dat: data.frame of regressors, containing variables specified in formula.full and formula.red
  # Sets: vector of variable names or indices or list of those, defining sets of variables
  # sumstat: function for summarizing univariate test statistics
  # perm: number of permutations
  
  if(!missing(Sets)){
    vars <- unique(unlist(Sets))
    if(!is.character(vars)){
      vars <- colnames(data)[vars]
      Sets <- lapply(Sets, function(x) colnames(data)[x])
    }
    data <- data[,vars]
  }
  else
    Sets <- colnames(data)
  
  if(!is.list(Sets))
    Sets <- list(Sets)
  
  unistats <- gGAteststats(data=data, formula.full=formula.full, formula.red=formula.red, model.dat=model.dat, perm=perm)
  Tobs  <- unistats$Tobs
  Tperm <- unistats$Tperm
  
  # observed global statistics 
  Sobs <- sapply(Sets, function(x) sumstat(Tobs[x]))
  
  # global permutation statistics 
  Sperm <- apply(Tperm, 2, function(x) sapply(Sets, function(y) sumstat(x[y])))
  if(is.null(dim(Sperm)))
    Sperm <- matrix(Sperm, ncol=ncol(Tperm))
  
  # mid p-values (see thesis of Monika Jelizarow, p. 44)
  p <- rowMeans(Sperm > Sobs, na.rm=TRUE) + 0.5 * rowMeans(Sperm == Sobs, na.rm=TRUE)
  
  return(data.frame(p.value=p, Statistic=Sobs))
}
