"plotgenes" <-
function(xx, model.dat, Colorgroup, redu.SSQ.Genes, msE.genes, legendpos, returnValues=FALSE, col, xlab, ylab, ...)
{
# Colorgroup: character variable giving the group that specifies coloring
# returnValues: shall gene-wise reductions in sum of squares = bar heights be returned?

  if(!is.character(Colorgroup) & !is.null(Colorgroup))
    stop("'Colorgroup' has to be a character")

  if(is.null(rownames(xx)))
    rownames(xx) <- 1:dim(xx)[1]
  N.Genes    <- dim(xx)[1]

  if(missing(col))
    # default color palette
    palette(c("#931638",rgb(1,.95,0.1),"lightblue","NavyBlue","#F7B50C","lightgreen","grey","mistyrose","#008751",rgb(1,.2,.2)))
  else if(is.numeric(col))
    palette(palette()[rep(col,2)])      # 'rep(col,2), da d. Palette nicht nur aus 1 Farbe bestehen kann (falls nur 1 angeg.)
  else  
    palette(rep(col,2))
 
  # if a Colorgroup variable is given and if it is not continuous
  colorgroup.vector <- as.numeric(model.dat[,Colorgroup])
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
        label <- c(label, paste("max. expression in",Colorgroup,"=",sort(unique(model.dat[,Colorgroup]))[i]))
      }
    }
  }

  # for a continuous group variable
  else
    color        <- 1

  # plotting results
  if(missing(xlab))
    xlab <- "Reduction in Sum of Squares"
  if(missing(ylab))
    ylab <- "Genes"

  #bars
  horizontal.bars(
        x         = rev(redu.SSQ.Genes),
        xlab      = xlab,
        ylab      = ylab,
        color     = rev(color),
        bar.names = rev(rownames(xx)), ...
  )

  # MSE-line
  pp    <- sort(c(-.5+(1:N.Genes),.5+(1:N.Genes)))
  vv    <- rev(rep(msE.genes,rep(2,N.Genes)))
  lines(vv,pp,type="s",lwd=2)

  # legend
  if(N.groups > 0 & N.groups <= 10)
    legend(legendpos, label, col=colind, pch=15)

  palette("default")
  
  # return bar heights
  if(returnValues){
    names(redu.SSQ.Genes) <- rownames(xx)
    return(redu.SSQ.Genes)
  }
}
