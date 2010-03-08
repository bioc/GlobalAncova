setGeneric("GlobalAncova.closed", function(xx,test.genes,formula.full,formula.red,model.dat,group,covars=NULL,
                                  test.terms,previous.test=NULL,level=0.05,method=c("permutation","approx"),perm=10000,max.group.size=2500,eps=1e-16,acc=50)
           standardGeneric("GlobalAncova.closed"))
# xx: expression matrix (rows=genes, columns=subjects)
# test.genes: list of pathways (each containing a vector of genes) that shall be tested
#             and adjusted using the closed testing procedure
# formula.full: model formula for the full model
# formula.red: model formula for the reduced model
# model.dat: data frame that contains the group and covariable information
# group: group variable
# covars: covariate information
# test.terms: character vector of terms of interest
# previous.test: result of a GlobalAncova with many pathways simultaneously
# level: alpha
# method: "permutation"/"approx" - permutation test / approximative test
# perm: number of permutations
# max.group.size: if a gene set is larger than max.group.size the permutation test is applied even if method="approx"
# eps: resuolution of approximation
# acc: additional accuracy parameter for approximation (the higher the value the higher the accuracy) 


############################## allgemeine Funktion #############################

setMethod("GlobalAncova.closed", signature(xx="matrix",test.genes="list",formula.full="formula",formula.red="formula",
                          model.dat="ANY",group="missing",covars="missing",test.terms="missing"),
          definition = function(xx,test.genes,formula.full,formula.red,model.dat,
                       previous.test=NULL,level=0.05,method=c("permutation","approx"),perm=10000,max.group.size=2500,eps=1e-16,acc=50){
  # test for model.dat
  if(!is.data.frame(model.dat))
    stop("'model.dat' has to be a data frame")

  GA.closed(xx=xx,test.genes=test.genes,formula.full=formula.full,formula.red=formula.red,model.dat=model.dat,
                       previous.test=previous.test,level=level,method=method,perm=perm,max.group.size=max.group.size,eps=eps,acc=acc)
}
)


########################## 'alte' Fkt. für 2 Gruppen ###########################


setMethod("GlobalAncova.closed", signature(xx="matrix",test.genes="list",formula.full="missing",formula.red="missing",
                          model.dat="missing",group="ANY",covars="ANY",test.terms="missing"),
          definition = function(xx,test.genes,group,covars=NULL,previous.test=NULL,level=0.05,
                          method=c("permutation","approx"),perm=10000,max.group.size=2500,eps=1e-16,acc=50){
  # parameter names
  group.name   <- deparse(substitute(group))

  if(is.null(dim(covars)))
    covar.names <- deparse(substitute(covars))
  else
    covar.names <- colnames(covars)

  # get formulas and 'model.dat' out of 'group' and 'covars'
  res          <- group2formula(group=group, group.name=group.name, covars=covars)
  formula.full <- res$formula.full
  formula.red  <- res$formula.red
  model.dat    <- res$model.dat

  GA.closed(xx=xx,test.genes=test.genes,formula.full=formula.full,formula.red=formula.red,model.dat=model.dat,
                       previous.test=previous.test,level=level,method=method,perm=perm,max.group.size=max.group.size,eps=eps,acc=acc)
}
)


################### allgemeine Funktion m. Angabe v. 'terms' ###################

setMethod("GlobalAncova.closed", signature(xx="matrix",test.genes="list",formula.full="formula",formula.red="missing",
                          model.dat="ANY",group="missing",covars="missing",test.terms="character"),
          definition = function(xx,test.genes,formula.full,model.dat,test.terms,previous.test=NULL,level=0.05,
                        method=c("permutation","approx"),perm=10000,max.group.size=2500,eps=1e-16,acc=50){
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

  GA.closed(xx=xx,test.genes=test.genes,formula.full=formula.full,D.red=D.red,model.dat=model.dat,
                       previous.test=previous.test,level=level,method=method,perm=perm,max.group.size=max.group.size,eps=eps,acc=acc)
}
)


################################################################################
################################################################################

# main function of GlobalAncova.closed
GA.closed <- function(xx,test.genes,formula.full,formula.red=NULL,D.red=NULL,model.dat,previous.test,level,
        method=c("permutation","approx"),perm=10000,max.group.size=2500,eps=1e-16,acc=50){

  new.data      <- Hnull.family(test.genes)
  result        <- list()
  endresult     <- list()
  just.tested   <- NULL
  sig           <- 1:length(test.genes)
  nsig          <- NULL

  for(i in 1:length(test.genes)){
    # hypotheses that must be tested "before" testing the end node
    intersection       <- lapply(new.data, function(x) match(test.genes[[i]], x))  # Durchschn. zw. Kn. i u. jeweils anderem
    intlength          <- lapply(intersection, function(x) sum(!is.na(x)))
    included           <- lapply(intlength, function(x) x==length(test.genes[[i]])) # Kn., in denen Kn. i enthalten ist
    related            <- which(unlist(included))

    result.i           <- matrix(NA, length(related), 3)
    method             <- match.arg(method)
    dimnames(result.i) <- switch(method, 
                          "permutation" = list(names(new.data[related]), c("genes","F.value","p.perm")),
                          "approx"      = list(names(new.data[related]), c("genes","F.value","p.approx")))

    for(j in 1:length(related)){
        name.j         <- names(new.data[related[j]])
        tested         <- name.j %in% just.tested
        prev.tested    <- name.j %in% rownames(previous.test)

        # if node has not been tested before
	      if(tested == FALSE && prev.tested==FALSE){      
          if(is.null(D.red))
            ga <- expr.test(xx=xx[new.data[[related[j]]], ],formula.full=formula.full,
                             formula.red=formula.red,model.dat=model.dat,method=method,perm=perm,max.group.size=max.group.size,eps=eps,acc=acc)
          else
            ga <- expr.test(xx=xx[new.data[[related[j]]], ],formula.full=formula.full,
                             D.red=D.red,model.dat=model.dat,method=method,perm=perm,max.group.size=max.group.size,eps=eps,acc=acc)
          result.i[j, 1]   <- length(new.data[[related[j]]])
          pvals            <- round(ga$test.result[1:2], 4)
          # if 'method=approx' is chosen but some gene groups are bigger than 'max.group.size' 
          if(is.na(pvals[2]))
            pvals[2] <- round(ga$test.result[3], 4)  # then take the permutation p-value
          result.i[j, 2:3] <- pvals
          just.tested      <- c(just.tested, name.j)
	     }

	     # nodes that have been tested before within this function
	     else if(tested == TRUE && prev.tested == FALSE){
          list.ind           <- lapply(res.nodes, function(x) name.j %in% x)
          # the correct result has to be taken (there might be also rows with only NA's)
          list.ind2          <- 1
          if(sum(unlist(list.ind)) > 1){
            res.lines        <- lapply(result[list.ind==TRUE], function(x) x[rownames(x) == name.j])
            noNAs            <- lapply(res.lines, function(x) !is.na(sum(x)))
            list.ind2        <- which(unlist(noNAs)==TRUE)
          }
          row.ind            <- which(rownames(result[list.ind==TRUE][list.ind2[1]][[1]]) %in% name.j)
          result.i[j, ]      <- result[list.ind==TRUE][list.ind2[1]][[1]][row.ind, ]
      }

       # nodes that have been tested before with GlobalAncova (with specified 'test.genes')
       else {
          # if 'method=approx' is chosen but some gene groups are bigger than 'max.group.size'
          ifelse(method == "permutation",
            result.i[j,] <- c(previous.test[name.j,1:2], p.perm=previous.test[name.j,"p.perm"]), 
            ifelse(is.na(previous.test[name.j,"p.approx"]), # if 'method=approx' is chosen but some gene groups are bigger than 'max.group.size' 
              result.i[j,] <- c(previous.test[name.j,1:2], p.approx=previous.test[name.j,"p.perm"]),  # -> take permutation p-value
              result.i[j,] <- c(previous.test[name.j,1:2], p.approx=previous.test[name.j,"p.approx"]))) 
       }

       # stop if one of the hypotheses can not be rejected
       if(result.i[j, 3] > level){
      	  sig          <- sig[sig != i]
      	  nsig         <- c(nsig, i)
      	  break
       }
    }
    result            <- c(result, list(result.i))
    res.nodes         <- lapply(result, function(x) rownames(x))
  }
  names(result)       <- names(new.data)[1:length(test.genes)]
  sig.nodes           <- names(new.data[sig])
  nsig.nodes          <- names(new.data[nsig])
  endresult           <- list(new.data, result, sig.nodes, nsig.nodes)
  names(endresult)    <- c("new.data","test.results","significant","not.significant")
  return(endresult)
}



################################################################################

# builds needed intersections of null hypotheses
Hnull.family <- function(test.genes){
# test.genes: a list of pathways

  new.nodes     <- list()
  name          <- list()
  new.data      <- test.genes

  for(i in 1:length(test.genes)){
    for(j in 1:length(new.data)){
        temp.nodes          <- unique(c(new.data[[i]], new.data[[j]]))
        if(!identical(temp.nodes, new.data[[i]])){
          # if node (=null hypothesis) does not exist yet
          samenode          <- lapply(new.nodes, function(x) identical(sort(temp.nodes), x))
          if(sum(as.logical(samenode)) == 0){
	           name      <- c(name, paste(names(new.data[i]), names(new.data[j]), sep="."))
	           new.nodes <- c(new.nodes, list(sort(temp.nodes)))
	           namekegg  <- c(names(new.data), paste(names(new.data[i]), names(new.data[j]), sep="."))
	           new.data  <- c(new.data, list(sort(temp.nodes)))
	           names(new.data) <- namekegg
	        }
	     }
    }
  }
  return(new.data)
}


