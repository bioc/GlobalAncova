setGeneric("Plot.subjects", function(xx,formula.full,formula.red,model.dat,group,covars=NULL,
                            test.terms,test.genes=NULL,Colorgroup=NULL,sort=FALSE,legendpos="topright",returnValues=FALSE,bar.names,...)
           standardGeneric("Plot.subjects"))
# xx: expression matrix (rows=genes, columns=subjects)
# formula.full: model formula for the full model
# formula.red: model formula for the reduced model
# model.dat: data frame that contains the group and covariable information
# group: group variable
# covars: covariate information
# test.terms: character vector of terms of interest
# test.genes: may define the relevant gene set 
# Colorgroup: character variable giving the group that specifies coloring
# sort: shall samples be sorted by 'Colorgroup'?
# legendpos: position of the legend
# returnValues: shall subject-wise reduction in sum of squares = bar heights be returned?
# bar.names: user specified bar names; if missing names of 'test.genes' or row names of 'xx' are taken
# ...: additional graphical parameters


################################# general function #############################

setMethod("Plot.subjects", signature(xx="matrix",formula.full="formula",formula.red="formula",
                          model.dat="ANY",group="missing",covars="missing",test.terms="missing"),
          definition = function(xx,formula.full,formula.red,model.dat,test.genes=NULL,Colorgroup=NULL,
                                sort=FALSE,legendpos="topright",returnValues=FALSE,bar.names,...){
  # test for model.dat
  if(!is.data.frame(model.dat))
    stop("'model.dat' has to be a data frame")

  # test for test.genes (i.e. only one gene set can be given and not a list of gene sets)
  if(!missing(test.genes)){
    if(!(data.class(test.genes) %in% c("numeric","character"))) 
      stop("'test.genes' has to be a vector of gene names or indices")
  }
    
  # get gene set 
  if(!missing(test.genes))
    xx <- xx[test.genes,,drop=FALSE] 

  # bar names
  if(is.null(colnames(xx)))
    colnames(xx) <- 1:ncol(xx)
  if(missing(bar.names))
    bar.names <- colnames(xx) 
  if(length(bar.names) != ncol(xx))
    stop("length of 'bar.names' not equal to sample size") 

  # basic analysis
  res <- reduSQ(xx=xx,formula.full=formula.full,formula.red=formula.red,model.dat=model.dat)
  redu.SSQ.Subjects    <- res$redu.subjects

  # plot
  plotsubjects(xx=xx,model.dat=model.dat,Colorgroup=Colorgroup,redu.SSQ.Subjects=redu.SSQ.Subjects,
                sort=sort,legendpos=legendpos,returnValues=returnValues,bar.names=bar.names,...)
}
)



########################## function for 2 groups ###############################

setMethod("Plot.subjects", signature(xx="matrix",formula.full="missing",formula.red="missing",
                          model.dat="missing",group="ANY",test.terms="missing"),
          definition = function(xx,group,covars=NULL,test.genes=NULL,Colorgroup=NULL,sort=FALSE,
                                legendpos="topright",returnValues=FALSE,bar.names,...){
  # test for test.genes (i.e. only one gene set can be given and not a list of gene sets)
  if(!missing(test.genes)){
    if(!(data.class(test.genes) %in% c("numeric","character"))) 
      stop("'test.genes' has to be a vector of gene names or indices")
  }
    
  # get gene set 
  if(!missing(test.genes))
    xx <- xx[test.genes,,drop=FALSE] 

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

  # bar names
  if(is.null(colnames(xx)))
    colnames(xx) <- 1:ncol(xx)
  if(missing(bar.names))
    bar.names <- colnames(xx) 
  if(length(bar.names) != ncol(xx))
    stop("length of 'bar.names' not equal to sample size") 

  # basic analysis
  res <- reduSQ(xx=xx,formula.full=formula.full,formula.red=formula.red,model.dat=model.dat)
  redu.SSQ.Subjects <- res$redu.subjects

  # plot
  plotsubjects(xx=xx,model.dat=model.dat,Colorgroup=Colorgroup,redu.SSQ.Subjects=redu.SSQ.Subjects,
              sort=sort,legendpos=legendpos,returnValues=returnValues,bar.names=bar.names,...)
}
)



############################# with 'test.terms' ################################

setMethod("Plot.subjects", signature(xx="matrix",formula.full="formula",formula.red="missing",
                          model.dat="ANY",group="missing",covars="missing",test.terms="character"),
          definition = function(xx,formula.full,model.dat,test.terms,test.genes=NULL,Colorgroup=NULL,
                                sort=FALSE,legendpos="topright",returnValues=FALSE,bar.names,...){
  # test for model.dat
  if(!is.data.frame(model.dat))
    stop("'model.dat' has to be a data frame")

  # test for test.genes (i.e. only one gene set can be given and not a list of gene sets)
  if(!missing(test.genes)){
    if(!(data.class(test.genes) %in% c("numeric","character"))) 
      stop("'test.genes' has to be a vector of gene names or indices")
  }
    
  # get gene set 
  if(!missing(test.genes))
    xx <- xx[test.genes,,drop=FALSE] 

  # test for 'test.terms'
  terms.all <- test.terms
  D.full    <- model.matrix(formula.full, model.dat)
  terms.all <- colnames(D.full)

  # are all terms variables compatible with 'model.dat'?
  if(!all(test.terms %in% terms.all))
    stop("'test.terms' are not compatible with the specified models")

  D.red  <- D.full[,!(colnames(D.full) %in% test.terms), drop=F]

  # bar names
  if(is.null(colnames(xx)))
    colnames(xx) <- 1:ncol(xx)
  if(missing(bar.names))
    bar.names <- colnames(xx) 
  if(length(bar.names) != ncol(xx))
    stop("length of 'bar.names' not equal to sample size") 

  # basic analysis
  res <- reduSQ(xx=xx,formula.full=formula.full,D.red=D.red,model.dat=model.dat)
  redu.SSQ.Subjects <- res$redu.subjects

  # plot
  plotsubjects(xx=xx,model.dat=model.dat,Colorgroup=Colorgroup,redu.SSQ.Subjects=redu.SSQ.Subjects,
              sort=sort,legendpos=legendpos,returnValues=returnValues,bar.names=bar.names,...)
}
)


################################################################################
################################################################################

# main function
plotsubjects <- function(xx, model.dat, Colorgroup, redu.SSQ.Subjects, sort=FALSE, legendpos, returnValues=FALSE, bar.names, col, xlab, ylab, ...){
  if(!is.character(Colorgroup) && !is.null(Colorgroup))  
    stop("'Colorgroup' has to be a character")

  N.Genes <- dim(xx)[1]

  if(missing(col))
    # color palette
    palette(c("#931638",rgb(1,.95,0.1),"lightblue","NavyBlue","#F7B50C","lightgreen","grey","mistyrose","#008751",rgb(1,.2,.2)))
  else if(is.numeric(col))
    palette(palette()[rep(col,2)])
  else
    palette(rep(col,2))

  # if a group variable is given and if it is not continuous
  colorgroup.vector <- as.numeric(model.dat[,Colorgroup])
  N.groups          <- length(unique(colorgroup.vector))
  if(N.groups > 0 && N.groups <= 10){
    color        <- numeric(N.groups)
    label        <- numeric(0)

    if(sort == TRUE){
      x          <- redu.SSQ.Subjects[order(colorgroup.vector)]
      bar.names   <- bar.names[order(colorgroup.vector)]
      gr.sort    <- sort(colorgroup.vector)
      for(i in 1:N.groups){
        color[gr.sort == unique(gr.sort)[i]] <- i
        label    <- c(label, paste(Colorgroup, "=", unique(gr.sort)[i]))
      }
    }

    else{
      for(i in 1:N.groups){
        x          <- redu.SSQ.Subjects
        color[colorgroup.vector == sort(unique(colorgroup.vector))[i]] <- i
        label      <- c(label, paste(Colorgroup, "=", sort(unique(model.dat[,Colorgroup]))[i]))
      }
    }
  }

  # if no group is given or for a continuous group variable
  else{
    x            <- redu.SSQ.Subjects
    color        <- 1
  }

# plotting results
  if(missing(xlab))
    xlab <- "Reduction in Sum of Squares"
  if(missing(ylab))
    ylab <- "Subjects"

  #bars
  horizontal.bars(
        x         = rev(x),
        xlab      = xlab,
        ylab      = ylab,
        color     = rev(color),
        bar.names = rev(bar.names), ...
  )

  # legend
  if(N.groups > 0 && N.groups <= 10)
    legend(legendpos, label, col=1:N.groups, pch=15)

  palette("default")

 # return bar heights
  if(returnValues){
    names(x) <- bar.names
    return(x)
  }
}
