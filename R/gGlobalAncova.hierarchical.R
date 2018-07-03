
#############################
### GAhier class definition
#############################

setClassUnion("matrixOrNULL", c("matrix", "NULL"))

setClass("GAhier",
         slots = c(clustervariables = "list", 
                   p.values = "list", 
                   alpha = "numeric", 
                   n.variables = "numeric",
                   permstats = "matrixOrNULL")
)


#########################################
### hierarchical testing
#########################################

## workhorse
hiertest <- function(H, Tobs, Tperm, sumstat, alpha){
  
  # go through dendrogram
  p.values <- clustervariables <- out <- Levels <- NULL
  m <- length(Tobs)
  clusters.k <- rep(1, m)
  for(k in 2:m){
    # only move on when there are significant clusters left
    if(length(out) == m)
      next
    
    clusters <- dendextend::cutree(H, k=k)
    clusters <- clusters[names(Tobs)]  # cutting dendrogram in subhierarchies might have changed the order...
    
    # set variables from previously non-significant clusters to NA
    clusters[out] <- NA
    
    # get new clusters to be tested (at each level there are only 2 new clusters, the others shall not be processed again)
    tab <- table(clusters.k, clusters)
    tmp <- apply(tab, 1, function(x) sum(x > 0) > 1)
    
    # don't move on when the new partition concerns only previous non-significant clusters
    # = when there is no new partition among the variables still under consideration
    if(all(!tmp))
      next
    
    print(paste("  level", k))
    Levels <- c(Levels, k)
    
    tmp <- tab[tmp,] > 0
    newclusters <- colnames(tab)[tmp]
    clusters.k <- clusters
    clusters[!(clusters %in% newclusters)] <- NA
    
    # observed global statistic for each cluster
    Sobs.k <- tapply(Tobs, clusters, sumstat)
    
    # global permutation statistic for each cluster
    Sperm.k <- apply(Tperm, 2, function(x) tapply(x, clusters, sumstat))
    
    # mid p-values
    p.k <- rowMeans(apply(Sperm.k, 2, function(x) x > Sobs.k)) + 0.5 * rowMeans(apply(Sperm.k, 2, function(x) x == Sobs.k))
    
    
    # Shaffer improvement: clusters with leaf node as sibling have larger effective size -> recieve less penalty (Meinshausen, p.272) 
    n.vars <- table(clusters)
    shaffer <- rev(n.vars == 1)
    
    # Meinshausen adjustment
    p.k[names(n.vars)] <- p.k[names(n.vars)] * m / (n.vars + shaffer)
    p.k <- pmin(p.k, 1)
    p.values <- c(p.values, list(p.k))
    
    # for non-significant clusters, remove respective variables from further testing
    nonsig <- p.k > alpha
    if(any(nonsig))
      out <- c(out, na.omit(names(clusters)[clusters %in% names(nonsig)[nonsig]]))
    
    # variables per tested cluster
    clustervars.k <- tapply(names(clusters), clusters, list)
    clustervariables <- c(clustervariables, list(clustervars.k))
  }
  names(p.values) <- names(clustervariables) <- paste("level", Levels, sep="")
  
  return(list(p.values=p.values, clustervariables=clustervariables))
}



## main function
gGlobalAncova.hierarchical <- function(data, H, formula.full, formula.red=~1, model.dat, sumstat=sum, 
                                       alpha=0.05, K, perm=10000, returnPermstats=FALSE, permstats){
  # data: data.frame (columns=variables)
  # H: dendrogram object specifying the hierarchy of the variables 
  # formula.full, formula.red: models to be compared
  # model.dat: data.frame of regressors, containing variables specified in formula.full and formula.red
  # family: type of response; defines suitable univariate statistics; ignored if alternative 'unistat' is given or if 'data' contains variables of mixed types
  # unistat: optional function for calculating univariate test statistic - ignored if 'data' contains variables of mixed types
  # sumstat: function for summarizing univariate test statistics
  # alpha: global alpha level   
  # K: real number; if this is specified, "short cut" on hierarchical testing will be applied on K subhierarchies
  #     -> correction for multiplicity: replace alpha by alpha / tau with tau = m / m.k
  #         with m = total number of variables, m.k = number of variables in k'th subtree
  # returnPermstats: if TRUE, the variable-wise statistics for all permutations are returned
  # permstats: if variable-wise permutation statistics were calculated previously, they can be provided in order
  #         not to repeat permutation testing (but only the hierarchical prodcedure)
  #         -> useful if procedure is run again w. different alpha and/or hierarchy H
  #         NOTE: data, formula.full and formula.red have to be identical to the previous call!
  
  if(!identical(sort(labels(H)), sort(colnames(data))))
    stop("'data' and 'H' do not seem to contain the same variables ('colnames(data)' and 'labels(H)' must coincide)")
  
  if(!missing(permstats))
    perm <- 0

  n.variables <- ncol(data)
  
  # univariate observed/permutation statistics 
  cat("testing global hypothesis...\n")
  unistats <- gGAteststats(data=data, formula.full=formula.full, formula.red=formula.red, model.dat=model.dat, perm=perm)
  if(missing(permstats)){
    Tobs  <- unistats$Tobs
    Tperm <- unistats$Tperm
  }
  else{
    Tobs  <- unistats
    Tperm <- permstats
  }
  
  # global p-value
  Sobs <- sumstat(Tobs)
  Sperm <- apply(Tperm, 2, sumstat)
  p <- mean(Sperm > Sobs) + 0.5 * mean(Sperm == Sobs)
  print(paste("global p-value =", p))
  if(p > alpha){
    print(paste("p >", alpha, "- no hierarchical testing is done"))
    return(p)
  }
  
  cat("hierarchical testing...\n")
  # hierarchical testing on complete hierarchy
  if(missing(K))
    out <- hiertest(H=H, Tobs=Tobs, Tperm=Tperm, sumstat=sumstat, alpha=alpha)
  
  # short cut - hierarchical testing on K subhierarchies
  else{
    # subdendrograms
    d <- dendextend::get_nodes_attr(H, "height", include_leaves=FALSE, na.rm=TRUE)
    d <- rev(sort(d))
    Hsub <- cut(H, h=d[K])$lower
    subtreevars <- sapply(Hsub, labels, simplify=FALSE)
    n.vars <- sapply(subtreevars, length)
    
    # alpha levels for subtrees
    #subtrees <- dendextend:::cutree(H, k=K)
    tau <- n.variables / n.vars
    alphas <- alpha / tau
    names(alphas) <- paste("subtree", 1:K, sep="")
    
    out <- list(p.values=NULL, clustervariables=NULL)
    for(i in 1:K){
      print(paste("subtree", i))
      #ind <- which(subtrees == i)
      ind <- subtreevars[[i]]
      
      # global tests for subtree roots
      Sobs.i <- sumstat(Tobs[ind])
      Sperm.i <- apply(Tperm[ind, , drop=FALSE], 2, sumstat)
      p.i <- mean(Sperm.i > Sobs.i) + 0.5 * mean(Sperm.i == Sobs.i)
      
      # only move on with subtrees where root is 'significant' - corrected for starting with subtrees
      if(p.i < alphas[i]){
        out.k <- hiertest(H=Hsub[[i]], Tobs=Tobs[ind], Tperm=Tperm[ind, , drop=FALSE], sumstat=sumstat, alpha=alphas[i])
        
        # re-adjust p-values by p * tau (to correspond to the global alpha)
        out.k$p.values <- lapply(out.k$p.values, function(x) x * tau[i])  
        out.k$p.values <- lapply(out.k$p.values, function(x) pmin(x, 1))
        
        names(out.k$p.values) <- names(out.k$clustervariables) <- paste(paste("K", i, sep=""), names(out.k$p.values), sep="_")
        out$p.values <- c(out$p.values, out.k$p.values)
        out$clustervariables <- c(out$clustervariables, out.k$clustervariables)
      }
    }
  }
  
  if(missing(K))
    obj <- new("GAhier", clustervariables=out$clustervariables, p.values=out$p.values, alpha=alpha, n.variables=n.variables)
  else
    obj <- new("GAhier", clustervariables=out$clustervariables, p.values=out$p.values, alpha=c(global=alpha, alphas), n.variables=n.variables)
  
  if(returnPermstats)
    obj@permstats <- Tperm
  
  return(obj)
}



######################################
### get significant end/leave nodes 
######################################

setGeneric("sigEndnodes", function(object, ...) {
  standardGeneric("sigEndnodes")
})


setMethod("sigEndnodes", signature="GAhier", function(object, onlySingleton=FALSE) {
  # onlySingleton: shall only single variables (leave nodes) be returned? in this case, or if all 
  #         significant end nodes are singletons anyways, a vector of variable names is returned; otherwise a list
  
  clusters <- object@clustervariables
  p        <- object@p.values
  alpha    <- object@alpha[1]
  
  m <- length(p)
  signodes <- sigp <- candidates <- NULL
  for(i in m:1){
    # get respective alpha for subtrees in case that hierarchical procedure was done w. shortcut 
#    if(length(alpha) > 1){
#      isub <- as.numeric(substr(names(p)[i], 2, 2))
#      a <- alpha[isub + 1]
#    }
#    else
#      a <- alpha
    
#    sig <- p[[i]] < a
    sig <- p[[i]] < alpha
    if(any(sig)){
      sigclust <- clusters[[i]][sig]
      sigp.i <- p[[i]][sig]
      names(sigclust) <- names(sigp.i) <- paste(names(p)[i], names(sigclust), sep=".")
      
      # don't move on where there was already a significant subcluster
      sigoffspring <- sapply(sigclust, function(x) any(candidates %in% x))
      sigclust <- sigclust[!sigoffspring]
      sigp.i <- sigp.i[!sigoffspring]
      
      signodes <- c(signodes, sigclust)
      sigp <- c(sigp, sigp.i)
      candidates <- c(candidates, unlist(sigclust))
    }
  }
  
  # if there are only leave nodes, return a vector
  if(length(signodes) == length(unlist(signodes)))
    return(unlist(signodes))
  
  # if only singletons shall be returned
  else if(onlySingleton){
    singletons <- Biobase::listLen(signodes) == 1
    return(unlist(signodes[singletons]))
  }
  
  else
    return(signodes)
  #return(list(endnodes=signodes, p=sigp))
})



######################################
### show/get results in a table 
######################################

setGeneric("results", function(object, ...) {
  standardGeneric("results")
})


setMethod("results", signature="GAhier", function(object){
  
  # significant end nodes
  signodes <- sigEndnodes(object)
  
  # variables per node
  #nvars <- NULL
  if(is.list(signodes)){
    nvars <- Biobase::listLen(signodes)
    signodes <- sapply(signodes, function(x) ifelse(length(x) < 4, paste(x, collapse=";"), 
                                                    paste(paste(x[1:3], collapse=";"), "...", sep=";")))
  }
  else
    nvars <- rep(1, length(signodes))
  
  # p-values
  pvals <- unlist(object@p.values)[names(signodes)]
  
  data.frame(variables=signodes, n.variables=nvars, p.value=pvals)
})


setMethod("show", signature="GAhier", function(object){
  cat(paste("results of hierarchical testing procedure for", object@n.variables, "variables\n"))
  
  alpha <- object@alpha
  cat(paste("global alpha =", alpha[1], "\n\n"))
  K <- length(alpha) - 1
  if(K > 0){
    cat(paste("procedure was split up into", K, "sub-hierarchies at alpha levels:", paste(base::format.pval(alpha[-1], digits=2), collapse=", "), "\n\n"))
    cat("significant end nodes (p-values are re-adjusted to global alpha):\n")
  }
  else
    cat("significant end nodes:\n")
  print(results(object))
})
  


#############################
### visualization
#############################

## highlight significant end nodes

coldend <- function(tree, siglist, col=1:2, lwd=1:2, pch=20){
  # tree: dendrogram object
  # siglist: list where each element contains the names of members in selected nodes
  # col: colors of non-selected and selected nodes
  # lwd: line width of branches to non-selected and selected nodes
  
  # significant node?
  uit <- tree
  sig <- any(sapply(siglist, function(x) all(x %in% labels(uit))))
  
  # end node?
  end <- any(sapply(siglist, function(x) setequal(x, labels(uit))))
  
  # set colors
  attr(uit, "edgePar") <- list(col=ifelse(sig, col[2], col[1]), lwd= ifelse(sig, lwd[2], lwd[1]))
  if(sig && end)
    attr(uit, "nodePar") <- list(col=col[2], pch=pch)
  
  # continue with the child branches
  if (!is.leaf(tree)) {
    select.branch <- 1:length(tree)
    for (i in 1:length(select.branch)) 
      uit[[i]] <- Recall(tree[[select.branch[i]]], siglist, col, lwd, pch)
  }
  return(uit)
}


setGeneric("Plot.hierarchy", function(object, ...) {
  standardGeneric("Plot.hierarchy")
})


setMethod("Plot.hierarchy", signature="GAhier", function(object, dend, col=1:2, lwd=1:2, collab, returndend=FALSE, cex.labels=1.5, ...) {
  # dend: dendrogram object specifying the hierarchy of the variables 
  # col: colors for significant and non-significant nodes + branches
  # lwd: line width for branches to non-significant and significant nodes
  # collab: vector of colors for coloring labels (can be chosen differently from sigleaves...);
  #         has to be named according to variable names; if missing same coloring as 'col' is made
  # returndend: shall adjusted dendrogram be returned?
  # cex.labels: size of leave labels
  
  # significant end nodes
  siglist <- sigEndnodes(object)
  
  # color for nodes + branches
  dend <- coldend(dend, siglist, col, lwd)
  
  # colors for labels - if not specified, highlight significant leave nodes
  if(missing(collab))
    collab <- ifelse(labels(dend) %in% siglist, col[2], col[1])
  else   # make sure to have the same order as in dendrogram
    collab <- collab[labels(dend)]
  
  dend <- dend %>% set("labels_col", collab) %>% set("labels_cex", cex.labels) 
  plot(dend, ...)
  
  if(returndend)
    return(dend) 
})


