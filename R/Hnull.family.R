"Hnull.family" <-
function(test.genes)
{
  # 'Hnull.family' builds needed intersections of null hypotheses
  # test.genes: a list of pathways

  new.nodes     <- list()
  name          <- list()
  new.data      <- test.genes

  for(i in 1:length(test.genes))
  {
    for(j in 1:length(new.data))
    {
        temp.nodes          <- unique(c(new.data[[i]], new.data[[j]]))
        if(!identical(temp.nodes, new.data[[i]]))
	{
          # if node (=null hypothesis) does not exist yet
          samenode          <- lapply(new.nodes, function(x) identical(sort(temp.nodes), x))
          if(sum(as.logical(samenode)) == 0)
	  {
	    name            <- c(name, paste(names(new.data[i]), names(new.data[j]), sep="."))
	    new.nodes       <- c(new.nodes, list(sort(temp.nodes)))
	    namekegg        <- c(names(new.data), paste(names(new.data[i]), names(new.data[j]), sep="."))
	    new.data        <- c(new.data, list(sort(temp.nodes)))
	    names(new.data) <- namekegg
	  }
	}
    }
  }
  return(new.data)
}

