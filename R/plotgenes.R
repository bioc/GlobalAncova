"plotgenes" <-
function(xx, model.dat, colorgroup, redu.SSQ.Genes, msE.genes, legendpos)
{
# colorgroup: character variable giving the group that specifies coloring

  if(!is.character(colorgroup) & !is.null(colorgroup))
    stop("'colorgroup' has to be a character")

  if(is.null(rownames(xx)))
    rownames(xx) <- 1:dim(xx)[1]
  N.Genes 	 <- dim(xx)[1]

  # color palette
 palette(c("#931638",rgb(1,.95,0.1),"lightblue","NavyBlue","#F7B50C","lightgreen","grey","mistyrose","#008751",rgb(1,.2,.2)))

  # if a colorgroup variable is given and if it is not continuous
  colorgroup.vector <- as.numeric(model.dat[,colorgroup])
  N.groups <- length(unique(colorgroup.vector))
  if(N.groups > 0 & N.groups <= 10)
  {
    # in which group has a gene the highest expression
    means <- NULL
    for(elt in unique(colorgroup.vector))
       means <- cbind(means, apply(xx, 1, function(x) mean(x[colorgroup.vector==elt])))
    up <- apply(means, 1, function(x) unique(colorgroup.vector)[which(x == max(x))])

    # colors and labels for the legend
    color        <- numeric(length(up))
    colind       <- numeric(0)
    label        <- numeric(0)
    for(i in 1:N.groups)
    {
      if(sort(unique(colorgroup.vector))[i] %in% up)
      {
        color[up == sort(unique(colorgroup.vector))[i]] <- i
        colind     <- c(colind, i)
        label <- c(label, paste("max. expression in",colorgroup,"=",sort(unique(model.dat[,colorgroup]))[i]))
      }
    }
  }

  # for a continuous group variable
  else
    color        <- 1

  # plotting results
  #bars
  horizontal.bars(
        x         = rev(redu.SSQ.Genes),
        xlabel    = "Reduction in Sum of Squares",
        ylabel    = "Genes",
        color     = rev(color),
        bar.names = rev(rownames(xx))
  )

  # MSE-line
  pp    <- sort(c(-.5+(1:N.Genes),.5+(1:N.Genes)))
  vv    <- rev(rep(msE.genes,rep(2,N.Genes)))
  lines(vv,pp,type="s",lwd=2)

  # legend
  if(N.groups > 0 & N.groups <= 10)
    legend(legendpos, label, col=colind, pch=15)

  palette("default")
}

