"GlobalAncova.closed" <-
function(xx,group,covars=NULL,test.genes,previous.test=NULL,level=0.05,perm=10000)
{
  # xx: gene expression matrix
  # group: vector containing the group information; must be coded as 0-1
  # covars: covariate information
  # test.genes: list of pathways (each containing a vector of genes)
  # previous.test: result of a GlobalAncova with many pathways simultaneously
  # level: alpha
  # perm: number of permutations
  
  new.data      <- Hnull.family(test.genes)
  result        <- list()
  endresult     <- list()
  just.tested   <- NULL
  sig           <- 1:length(test.genes)
  nsig          <- NULL

  for(i in 1:length(test.genes))
  {
    # hypotheses that must be tested "before" testing the end node
    related            <- grep(names(new.data)[i], names(new.data))
    len                <- lapply(names(new.data[related]), nchar)
    related            <- related[order(unlist(len))]
    
    result.i           <- matrix(NA, length(related), 4)
    dimnames(result.i) <- list(names(new.data[related]), c("genes","F.value","p.value.perm","p.value.theo"))

    for(j in 1:length(related))
      {
        name.j         <- names(new.data[related[j]])
        tested         <- name.j %in% just.tested
        prev.tested    <- name.j %in% rownames(previous.test)
        
        # if node has not been tested before
	if(tested == FALSE & prev.tested==FALSE)
	{
          # if there is only one gene in the pathway
	  if(length(new.data[[related[j]]]) == 1)
	  {
	    ga               <- GlobalAncova(xx[new.data[[related[j]]],],group=group,perm=perm)
      	    result.i[j, 1]   <- length(new.data[[related[j]]])
	    result.i[j, 2:4] <- round(ga, 4)
	    just.tested      <- c(just.tested, name.j)
	  }
	  
	  else
	  {
	    ga               <- GlobalAncova(xx[new.data[[related[j]]], ],group=group,covars=covars,perm=perm)
	    result.i[j, 1]   <- length(new.data[[related[j]]])
	    result.i[j, 2:4] <- round(ga[[2]][1:3], 4)
	    just.tested      <- c(just.tested, name.j)
	  }
	}
	
	# nodes that have been tested before within this function
	else if(tested == TRUE & prev.tested == FALSE)
	{
          list.ind           <- lapply(res.nodes, function(x) name.j %in% x)
          
          # the correct result has to be taken (there might be also rows with only NA's)
          list.ind2          <- 1
          if(sum(unlist(list.ind)) > 1)
          {
            res.lines        <- lapply(result[list.ind==T], function(x) x[rownames(x) == name.j])
            noNAs            <- lapply(res.lines, function(x) !is.na(sum(x)))
            list.ind2        <- which(unlist(noNAs)==TRUE)
          }
          row.ind            <- which(rownames(result[list.ind==T][list.ind2[1]][[1]]) %in% name.j)
          result.i[j, ]      <- result[list.ind==T][list.ind2[1]][[1]][row.ind, ]
        }
        
        # nodes that have been tested before with GlobalAncova (with specified 'test.genes')
        else
          result.i[j, ]      <- previous.test[name.j, ]

        # stop if one of the hypotheses can not be rejected
	if(result.i[j, 3] > level)
	{
	  sig          <- sig[sig != i]
	  nsig         <- c(nsig, i)
	  break
	}
      }
    result            <- c(result, list(result.i))
    res.nodes          <- lapply(result, function(x) rownames(x))
  }
  names(result)       <- names(new.data)[1:length(test.genes)]
  sig.nodes           <- names(new.data[sig])
  nsig.nodes          <- names(new.data[nsig])
  #hypotheses         <- names(new.data)
  #endresult          <- list(hypotheses, new.data, result, sig.nodes, nsig.nodes)
  endresult           <- list(new.data, result, sig.nodes, nsig.nodes)
  #names(endresult)   <- c("hypotheses","new.data","test.results","significant","not.significant")
  names(endresult)    <- c("new.data","test.results","significant","not.significant")
  return(endresult)
}

