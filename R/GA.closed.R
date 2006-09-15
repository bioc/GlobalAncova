"GA.closed" <-
function(xx,test.genes,formula.full,formula.red=NULL,D.red=NULL,model.dat,previous.test,level,perm)
{
# xx: gene expression matrix
# test.genes: list of pathways (each containing a vector of genes) that shall be tested
#             and adjusted using the closed testing procedure
# formula.full: model formula for the full model
# formula.red: model formula for the reduced model
# model.dat: data frame that contains the group and covariable information
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
    intersection       <- lapply(new.data, function(x) match(test.genes[[i]], x))  # Durchschn. zw. Kn. i u. jeweils anderem
    intlength          <- lapply(intersection, function(x) sum(!is.na(x)))
    included           <- lapply(intlength, function(x) x==length(test.genes[[i]])) # Kn., in denen Kn. i enthalten ist
    related            <- which(unlist(included))

    result.i           <- matrix(NA, length(related), 4)
    dimnames(result.i) <- list(names(new.data[related]), c("genes","F.value","p.value","p.perm"))

    for(j in 1:length(related))
      {
        name.j         <- names(new.data[related[j]])
        tested         <- name.j %in% just.tested
        prev.tested    <- name.j %in% rownames(previous.test)

        # if node has not been tested before
	if(tested == FALSE & prev.tested==FALSE)
	{
          if(is.null(D.red))
            ga <- expr.test(xx=xx[new.data[[related[j]]], ],formula.full=formula.full,
                             formula.red=formula.red,model.dat=model.dat,perm=perm)
          else
            ga <- expr.test(xx=xx[new.data[[related[j]]], ],formula.full=formula.full,
                             D.red=D.red,model.dat=model.dat,perm=perm)
          result.i[j, 1]   <- length(new.data[[related[j]]])
          result.i[j, 2:4] <- round(ga$test.result, 4)
          just.tested      <- c(just.tested, name.j)
	}

	# nodes that have been tested before within this function
	else if(tested == TRUE & prev.tested == FALSE)
	{
          list.ind           <- lapply(res.nodes, function(x) name.j %in% x)
          # the correct result has to be taken (there might be also rows with only NA's)
          list.ind2          <- 1
          if(sum(unlist(list.ind)) > 1)
          {
            res.lines        <- lapply(result[list.ind==TRUE], function(x) x[rownames(x) == name.j])
            noNAs            <- lapply(res.lines, function(x) !is.na(sum(x)))
            list.ind2        <- which(unlist(noNAs)==TRUE)
          }
          row.ind            <- which(rownames(result[list.ind==TRUE][list.ind2[1]][[1]]) %in% name.j)
          result.i[j, ]      <- result[list.ind==TRUE][list.ind2[1]][[1]][row.ind, ]
        }

        # nodes that have been tested before with GlobalAncova (with specified 'test.genes')
        else
          result.i[j, ]      <- previous.test[name.j, ]

        # stop if one of the hypotheses can not be rejected
	if(result.i[j, "p.perm"] > level)
	{
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

