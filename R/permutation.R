
# main function for permutation test
resampleGA <- function(xx, formula.full, D.full, D.red, model.dat, perm, test.genes, F.value,
                         DF.full, DF.extra){
  N.Subjects  <- ncol(xx)
  rr     <- row.orth2d(xx, D.red)

  # sum of squares in reduced model do not have to be re-calculated in each permutation
  genewSS  <- genewiseGA(xx, D.full, D.red=D.red)
  SS.red.i <- genewSS[,1] + genewSS[,2]   # SS.red = SS.extra + SS.full

  D.full.perm <- D.full
  test.col <- !colnames(D.full) %in% colnames(D.red)
  count <- numeric(length(test.genes))

  test.genes0 <- test.genes
  if (length(test.genes) == 1) {
      test.genes <- unlist(test.genes)
      test.genes <- which( rownames(xx) %in% test.genes ) -1
      num.test.genes <- length(test.genes)
  } else {
      test.genes <- sapply(test.genes, function(x) which(rownames(xx) %in% x))
      num.test.genes <- sapply(test.genes, length)
      test.genes <- unlist(test.genes)-1
  }

  nperm <- .nPerms(D.full, model.dat, formula.full)
  faccounts <- nperm$counts
  nperm <- nperm$nPerms

  use.permMat <- 0
  permMat     <- 1
  nperm.used  <- perm
  if(nperm < perm){
    use.permMat <- 1
    nperm.used <- nperm

    print(paste("enumerating all", nperm, "permutations"))
    if(is.null(faccounts))  # design with continuous covariates
      permMat <- .allperms(1:N.Subjects)
    else
      permMat <- .allpermsG(faccounts, faccounts) # factorial design
    permMat <- apply(permMat, 2, order)  # get possible orderings
  }

  test.col0 <- which(test.col) -1 
  count <- .C ("permut", PACKAGE="GlobalAncova", as.double(D.full),as.integer(nrow(D.full)),
                   as.integer(ncol(D.full)) , as.double(D.full.perm) ,
                   as.double (D.red) , as.integer(nrow(D.red)),as.integer(ncol(D.red)),
                   as.integer(N.Subjects) ,
                   as.double(rr),as.integer(nrow(rr)), as.integer(ncol(rr)) ,
                   as.double(SS.red.i) , perm=as.integer(nperm.used) ,
                   as.integer(test.col0) , as.integer (length(test.col0)) , as.double(F.value),
                   as.double(DF.full), as.double(DF.extra) , as.integer(permMat-1) ,
                   as.integer(test.genes) , as.integer(num.test.genes) ,
                   as.integer(length(num.test.genes)) , count=integer(length(num.test.genes)) ,
                   n.singular=as.integer(0) , as.integer(use.permMat)
               )

  p.value <- count$count / (nperm.used - count$n.singular)
  return(p.value)
}


################################################################################

# functions for getting all possible permutations
# (adapted from package 'globaltest')


#==========================================================
# Calculates the number of permutations for multiple groups
#==========================================================
.mchoose <- function(counts) {
  out <- choose(sum(counts), counts[1])
  if (length(counts) > 2)
    out <- out * .mchoose(counts[-1])
  out
}

# number of permutations can be adjusted since e.g. factor combination 112233 yields the same as 332211...
.nPermsG <- function(counts, grouping) {
  total <- .mchoose(counts)
  if (any(!is.na(grouping))) {
    correction <- prod(factorial(sapply(unique(grouping[!is.na(grouping)]), function(cc) sum(grouping == cc, na.rm=TRUE))))
  } else {
    correction <- 1
  }
  total / correction
}



#==========================================================
# Calculates the number of permutations
#==========================================================
.nPerms <- function(D.full, model.dat, formula.full) {
# D.full: full design matrix
# model.dat: data frame with covariate information
# formula.full: model formula

  n <- nrow(D.full)

  # get all variables from 'model.dat' involved in the design
  design.terms <- as.character(terms(formula.full)@variables)[-1]
  design.terms <- intersect(design.terms, names(model.dat))

  # is there any continuous variable in the design -> then all n! permutations may be different
  continuous <- sapply(model.dat[,design.terms, drop=F], is.numeric)  # are there continuous covariates
  varlength  <- sapply(model.dat[,design.terms, drop=F], function(x) length(unique(x)))  # two-group variables may be 'numeric' and not 'factor'
  continuous <- ifelse(continuous & (varlength > 2), TRUE, FALSE)
  if(any(continuous)){
    out <- ifelse(n <= 100, factorial(n), Inf)
    counts <- NULL
  }

  else{
    # get all *different* rows in the design matrix
    unique.rows <- unique(D.full)
    Y <- numeric(nrow(unique.rows))
    for(i in 1:nrow(unique.rows)){
      equal.rows <- t(D.full) == unique.rows[i,]
      equal.rows <- apply(equal.rows, 2, all)
      Y[equal.rows] <- i
    }
    counts <- sapply(unique(Y), function(x) sum(Y == x))
    out <- .nPermsG(counts, counts)
  }
  return(list(nPerms=out, counts=counts))
}



#==========================================================
# Lists all permutations for the multiple-group case
#==========================================================
.allpermsG <- function(counts, grouping) {
  n <- sum(counts)
  if (n == 1) {
    app <- which.max(counts)
  } else {
    total <- .nPermsG(counts, grouping)
    app <- matrix(,n,total)
    choosable <- (counts > 0) & (is.na(grouping) | (1:length(counts) %in% match(unique(grouping[!is.na(grouping)]), grouping)))
    choosable <- (1:length(counts))[choosable]
    ix <- 0
    for (iy in choosable) {
      countstemp <- counts
      countstemp[iy] <- counts[iy] - 1
      groupingtemp <- grouping
      groupingtemp[iy] <- NA
      size <- .nPermsG(countstemp, groupingtemp)
      app[1,(ix+1):(ix+size)] <- iy
      app[2:n, (ix+1):(ix+size)] <- .allpermsG(countstemp, groupingtemp)
      ix <- ix + size
    }
  }
  app
}


#==========================================================
# Lists all permutations for the continuos case
#==========================================================
.allperms <- function(nums) {
# nums: indices to be permuted

  n <- length(nums)
  if (n == 1) {
    app <- nums
  } else {
    app <- matrix(,n,factorial(n))
    for (ix in 1:length(nums)) {
      range <- 1:factorial(n-1) + (ix - 1) * factorial(n-1)
      app[1,range] <- nums[ix]
      app[2:n,range] <- .allperms(nums[-ix])
    }
  }
  app
}
