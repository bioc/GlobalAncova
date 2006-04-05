"plotsubjects" <-
function(xx, model.dat, colorgroup, redu.SSQ.Subjects, sort=FALSE, legendpos)
{
# colorgroup: character variable giving the group that specifies coloring
# sort: shall samples be sorted by 'colorgroup'?

  if(!is.character(colorgroup) & !is.null(colorgroup))
    stop("'colorgroup' has to be a character")

  if(is.null(colnames(xx)))
    colnames(xx) <- seq(1:dim(xx)[2])
  N.Genes <- dim(xx)[1]

  # color palette
  palette(c("#931638",rgb(1,.95,0.1),"lightblue","NavyBlue","#F7B50C","lightgreen","grey","mistyrose","#008751",rgb(1,.2,.2)))

  # if a group variable is given and if it is not continuous
  colorgroup.vector <- as.numeric(model.dat[,colorgroup])
  N.groups          <- length(unique(colorgroup.vector))
  if(N.groups > 0 & N.groups <= 10)
  {
    color        <- numeric(N.groups)
    label        <- numeric(0)

    if(sort == TRUE)
    {
      x          <- redu.SSQ.Subjects[order(colorgroup.vector)]
      barnames   <- colnames(xx)[order(colorgroup.vector)]
      gr.sort    <- sort(colorgroup.vector)
      for(i in 1:N.groups)
      {
        color[gr.sort == unique(gr.sort)[i]] <- i
        label    <- c(label, paste(colorgroup, "=", unique(gr.sort)[i]))
      }
    }

    else
    {
      for(i in 1:N.groups)
      {
        x          <- redu.SSQ.Subjects
        barnames   <- colnames(xx)
        color[colorgroup.vector == sort(unique(colorgroup.vector))[i]] <- i
        label      <- c(label, paste(colorgroup, "=", sort(unique(model.dat[,colorgroup]))[i]))
      }
    }
  }

  # if no group is given or for a continuous group variable
  else
  {
    x            <- redu.SSQ.Subjects
    barnames     <- colnames(xx)
    color        <- 1
  }

# plotting results
  #bars
  horizontal.bars(
        x         = rev(x),
        xlabel    = "Reduction in Sum of Squares",
        ylabel    = "Subjects",
        color     = rev(color),
        bar.names = rev(barnames)
  )

  # legend
  if(N.groups > 0 & N.groups <= 10)
    legend(legendpos, label, col=1:N.groups, pch=15)

  palette("default")
}

