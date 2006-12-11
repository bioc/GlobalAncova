setGeneric("Plot.subjects", function(xx,formula.full,formula.red,model.dat,group,covars=NULL,
                            test.terms,Colorgroup=NULL,sort=FALSE,legendpos="topright",returnValues=FALSE,...)
           standardGeneric("Plot.subjects"))
# xx: expression matrix (rows=genes, columns=subjects)
# formula.full: model formula for the full model
# formula.red: model formula for the reduced model
# model.dat: data frame that contains the group and covariable information
# group: group variable
# covars: covariate information
# test.terms: character vector of terms of interest
# Colorgroup: character variable giving the group that specifies coloring
# sort: shall samples be sorted by 'Colorgroup'?
# legendpos: position of the legend
# returnValues: shall subject-wise reduction in sum of squares = bar heights be returned?


############################# allgemeine Funktion ##############################

setMethod("Plot.subjects", signature(xx="matrix",formula.full="formula",formula.red="formula",
                          group="missing",covars="missing",test.terms="missing"),
          definition = function(xx,formula.full,formula.red,model.dat,Colorgroup=NULL,sort=FALSE,legendpos="topright",returnValues=FALSE,...)
{
  # test for model.dat
  if(!is.data.frame(model.dat))
    stop("'model.dat' has to be a data frame")

  # basic analysis
  res <- reduSQ(xx=xx,formula.full=formula.full,formula.red=formula.red,model.dat=model.dat)
  redu.SSQ.Subjects    <- res$redu.subjects

  # plot
  plotsubjects(xx=xx,model.dat=model.dat,Colorgroup=Colorgroup,redu.SSQ.Subjects=redu.SSQ.Subjects,sort=sort,legendpos=legendpos,returnValues=returnValues,...)
}
)



########################## 'alte' Fkt. für 2 Gruppen ###########################

setMethod("Plot.subjects", signature(xx="matrix",formula.full="missing",formula.red="missing",
                          model.dat="missing",group="ANY",test.terms="missing"),
          definition = function(xx,group,covars=NULL,Colorgroup=NULL,sort=FALSE,legendpos="topright",returnValues=FALSE,...)
{
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
  redu.SSQ.Subjects <- res$redu.subjects

  # plot
  plotsubjects(xx=xx,model.dat=model.dat,Colorgroup=Colorgroup,redu.SSQ.Subjects=redu.SSQ.Subjects,sort=sort,legendpos=legendpos,returnValues=returnValues,...)
}
)



################### allgemeine Funktion m. Angabe v. 'terms' ###################

setMethod("Plot.subjects", signature(xx="matrix",formula.full="formula",formula.red="missing",
                          group="missing",covars="missing",test.terms="character"),
          definition = function(xx,formula.full,test.terms,model.dat,Colorgroup=NULL,sort=FALSE,legendpos="topright",returnValues=FALSE,...)
{
  # test for model.dat
  if(!is.data.frame(model.dat))
    stop("'model.dat' has to be a data frame")

  # test for 'test.terms' and derive 'formula.red'
  #formula.red <- terms2formula(formula.full=formula.full, test.terms=test.terms, model.dat=model.dat)

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
  redu.SSQ.Subjects <- res$redu.subjects

  # plot
  plotsubjects(xx=xx,model.dat=model.dat,Colorgroup=Colorgroup,redu.SSQ.Subjects=redu.SSQ.Subjects,sort=sort,legendpos=legendpos,returnValues=returnValues,...)
}
)
