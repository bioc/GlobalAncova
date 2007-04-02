setGeneric("Plot.genes", function(xx,formula.full,formula.red,model.dat,group,covars=NULL,
                         test.terms,test.genes=NULL,Colorgroup=NULL,legendpos="topright",returnValues=FALSE,...)
            standardGeneric("Plot.genes"))
# xx: expression matrix (rows=genes, columns=subjects)
# formula.full: model formula for the full model
# formula.red: model formula for the reduced model
# model.dat: data frame that contains the group and covariable information
# group: group variable
# covars: covariate information
# test.terms: character vector of terms of interest
# Colorgroup: character variable giving the group that specifies coloring
# legendpos: position of the legend


############################# general function #################################

setMethod("Plot.genes", signature(xx="matrix",formula.full="formula",formula.red="formula",
                          group="missing",covars="missing",test.terms="missing"),
          definition = function(xx,formula.full,formula.red,model.dat,test.genes=NULL,Colorgroup=NULL,legendpos="topright",returnValues=FALSE,...){
  # test for model.dat
  if(!is.data.frame(model.dat))
    stop("'model.dat' has to be a data frame")

  # test for test.genes (i.e. only one gene set can be given and not a list of gene sets)
  if(!is.null(test.genes) & !(data.class(test.genes) %in% c("numeric","character"))) 
    stop("'test.genes' has to be a vector of gene names or indices")
    
  # get gene set 
  if(!is.null(test.genes))
    xx <- xx[test.genes,] 

  # basic analysis
  res <- reduSQ(xx=xx,formula.full=formula.full,formula.red=formula.red,model.dat=model.dat)
  redu.SSQ.Genes    <- res$redu.genes
  msE.genes         <- res$mse

  # plot
  plotgenes(xx=xx,model.dat=model.dat,Colorgroup=Colorgroup,redu.SSQ.Genes=redu.SSQ.Genes,
            msE.genes=msE.genes,legendpos=legendpos,returnValues=returnValues,...)
}
)


########################## function for 2 groups ###############################

setMethod("Plot.genes", signature(xx="matrix",formula.full="missing",formula.red="missing",
                          model.dat="missing",group="ANY",test.terms="missing"),
          definition = function(xx,group,covars=NULL,test.genes=NULL,Colorgroup=NULL,legendpos="topright",returnValues=FALSE,...){
  # test for test.genes (i.e. only one gene set can be given and not a list of gene sets)
  if(!is.null(test.genes) & !(data.class(test.genes) %in% c("numeric","character"))) 
    stop("'test.genes' has to be a vector of gene names or indices")
    
  # get gene set 
  if(!is.null(test.genes))
    xx <- xx[test.genes,] 

  # 'group' is assumed to be the variable relevant for coloring
  if(is.null(Colorgroup))
    Colorgroup <- deparse(substitute(group))

  # group name
  group.name   <- deparse(substitute(group))

  if(is.null(dim(covars)))
    covar.names <- deparse(substitute(covars))
  else
    covar.names <- colnames(covars)

  # get formulas and 'model.dat' out of 'group' and 'covars'
  res          <- group2formula(group=group, group.name=group.name, covars=covars, covar.names=covar.names)
  formula.full <- res$formula.full
  formula.red  <- res$formula.red
  model.dat    <- res$model.dat

  # basic analysis
  res <- reduSQ(xx=xx,formula.full=formula.full,formula.red=formula.red,model.dat=model.dat)
  redu.SSQ.Genes <- res$redu.genes
  msE.genes      <- res$mse

  # plot
  plotgenes(xx=xx,model.dat=model.dat,Colorgroup=Colorgroup,redu.SSQ.Genes=redu.SSQ.Genes,
            msE.genes=msE.genes,legendpos=legendpos,returnValues=returnValues,...)
}
)


############################# with 'test.terms' ################################

setMethod("Plot.genes", signature(xx="matrix",formula.full="formula",formula.red="missing",
                          group="missing",covars="missing",test.terms="character"),
          definition = function(xx,formula.full,test.terms,model.dat,test.genes=NULL,Colorgroup=NULL,legendpos="topright",returnValues=FALSE,...){
  # test for model.dat
  if(!is.data.frame(model.dat))
    stop("'model.dat' has to be a data frame")

  # test for test.genes (i.e. only one gene set can be given and not a list of gene sets)
  if(!is.null(test.genes) & !(data.class(test.genes) %in% c("numeric","character"))) 
    stop("'test.genes' has to be a vector of gene names or indices")
    
  # get gene set 
  if(!is.null(test.genes))
    xx <- xx[test.genes,] 

  # test for 'test.terms'
  terms.all <- test.terms
  D.full    <- model.matrix(formula.full, model.dat)
  terms.all <- colnames(D.full)

  # are all terms variables compatible with 'model.dat'?
  if(!all(test.terms %in% terms.all))
    stop("'test.terms' are not compatible with the specified models")

  D.red  <- D.full[,!(colnames(D.full) %in% test.terms), drop=F]

  # basic analysis
  res <- reduSQ(xx=xx,formula.full=formula.full,D.red=D.red,model.dat=model.dat)
  redu.SSQ.Genes <- res$redu.genes
  msE.genes      <- res$mse

  # plot
  plotgenes(xx=xx,model.dat=model.dat,Colorgroup=Colorgroup,redu.SSQ.Genes=redu.SSQ.Genes,
            msE.genes=msE.genes,legendpos=legendpos,returnValues=returnValues,...)
}
)

################################################################################
################################################################################

# main function
plotgenes <- function(xx, model.dat, Colorgroup, redu.SSQ.Genes, msE.genes, legendpos, returnValues=FALSE, col, xlab, ylab, ...){
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
  if(N.groups > 0 & N.groups <= 10){
    # in which group has a gene the highest expression
    means <- NULL
    for(elt in unique(colorgroup.vector))
       means <- cbind(means, apply(xx, 1, function(x) mean(x[colorgroup.vector==elt])))
    up <- apply(means, 1, function(x) unique(colorgroup.vector)[which(x == max(x))])

    # colors and labels for the legend
    color        <- numeric(length(up))
    colind       <- numeric(0)
    label        <- numeric(0)
    for(i in 1:N.groups){
      if(sort(unique(colorgroup.vector))[i] %in% up){
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


################################################################################

# function for plotting horizontal bars with labels added at right margin
# bars are determined by value of x which is assumed to be a vector
# no formal check of variables performed
# setting the plot region
# (also used in Plot.subjects')
horizontal.bars <- function(x, labelsize=.75, bar.names=NULL, color, xlim,...){
        if(missing(xlim)){  
          xlim    <- 0.05*c(-1,1)*range(x)+c(min(x),max(x))
          xlim[1] <- min(0,xlim[1])
        }
        n       <- length(x)
        ylim    <- c(0,n+1)

        # enlarging right margin for bar.names

        if(!is.null(bar.names)&length(bar.names)==n){
        names   <- TRUE
                   plot.new()
        w       <- 1.5 * max(strwidth(bar.names, "inches", labelsize))
        oldmai  <- par("mai")
                   par(mai=c(oldmai[1:3],max(w,oldmai[4])), new=T)
        }

        # plotting bars with border=F nothing appears color is NULL
        plot(0,type="n", xlim=xlim, ylim=ylim, yaxt="n",...)
                rect(rep(0,n),(1:n)-.3,x,(1:n)+.3, col=color, border="white")
                box()

        # adding bar.names at right margin
        if(names){
                axis(4,at=1:n,bar.names,cex.axis=labelsize,las=2)
                par(mai=oldmai)
        }
}

################################################################################

# computes the reduction in sum of squares for genes and subjects and the MSE
# (also used in Plot.subjects' and 'Plot.sequential')
reduSQ <- function(xx, formula.full, formula.red=NULL, D.red=NULL, model.dat){
  N.Subjects     <- ncol(xx)

  # design matrices
  D.full     <- model.matrix(formula.full, data=model.dat)
  if(is.null(D.red))
    D.red    <- model.matrix(formula.red,  data=model.dat)

  N.par.full <- ncol(D.full)
  N.par.red  <- ncol(D.red)

  # residuals
  R.full     <- row.orth2d(xx,D.full)
  R.red      <- row.orth2d(xx,D.red)

  # reduction sum of squares
  redu.sq           <- R.red^2 - R.full^2
  redu.SSQ.Genes    <- rowSums(redu.sq) / (N.par.full - N.par.red) 

  redu.SSQ.Subjects <- colSums(redu.sq)

  # mean square error
  msE.genes      <- rowSums(R.full^2) / (N.Subjects - N.par.full)

  return(list(redu.genes=redu.SSQ.Genes,redu.subjects=redu.SSQ.Subjects,mse=msE.genes))
}
