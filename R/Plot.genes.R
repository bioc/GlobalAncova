setGeneric("Plot.genes", function(xx,formula.full,formula.red,model.dat,group,covars=NULL,
                         test.terms,colorgroup=NULL,legendpos="topright")
            standardGeneric("Plot.genes"))
# xx: expression matrix (rows=genes, columns=subjects)
# formula.full: model formula for the full model
# formula.red: model formula for the reduced model
# model.dat: data frame that contains the group and covariable information
# group: group variable
# covars: covariate information
# test.terms: character vector of terms of interest
# colorgroup: character variable giving the group that specifies coloring
# legendpos: position of the legend


################################# general function #############################

setMethod("Plot.genes", signature(xx="matrix",formula.full="formula",formula.red="formula",
                          group="missing",covars="missing",test.terms="missing"),
          definition = function(xx,formula.full,formula.red,model.dat,colorgroup=NULL,legendpos="topright")
{
  # test for model.dat
  if(!is.data.frame(model.dat))
    stop("'model.dat' has to be a data frame")

  # basic analysis
  res <- reduSQ(xx=xx,formula.full=formula.full,formula.red=formula.red,model.dat=model.dat)
  redu.SSQ.Genes    <- res$redu.genes
  msE.genes         <- res$mse

  # plot
  plotgenes(xx=xx,model.dat=model.dat,colorgroup=colorgroup,redu.SSQ.Genes=redu.SSQ.Genes,msE.genes=msE.genes,legendpos=legendpos)
}
)


########################## function for 2 groups ################################

setMethod("Plot.genes", signature(xx="matrix",formula.full="missing",formula.red="missing",
                          model.dat="missing",group="numeric",test.terms="missing"),
          definition = function(xx,group,covars=NULL,colorgroup=NULL,legendpos="topright")
{
  # 'group' is assumed to be the variable relevant for coloring
  if(is.null(colorgroup))
    colorgroup <- deparse(substitute(group))

  # group name
  group.name   <- deparse(substitute(group))

  if(is.null(dim(covars)))
    covar.names <- deparse(substitute(covars))
  else
    covar.names <- colnames(covars)

  # get formulas and 'model.dat' out of 'group' and 'covars'
  res          <- group2formula(group=group, group.name=group.name, covars=covars, covar.names)
  formula.full <- res$formula.full
  formula.red  <- res$formula.red
  model.dat    <- res$model.dat

  # basic analysis
  res <- reduSQ(xx=xx,formula.full=formula.full,formula.red=formula.red,model.dat=model.dat)
  redu.SSQ.Genes <- res$redu.genes
  msE.genes      <- res$mse

  # plot
  plotgenes(xx=xx,model.dat=model.dat,colorgroup=colorgroup,redu.SSQ.Genes=redu.SSQ.Genes,msE.genes=msE.genes,legendpos=legendpos)
}
)


############################# with 'test.terms' ################################

setMethod("Plot.genes", signature(xx="matrix",formula.full="formula",formula.red="missing",
                          group="missing",covars="missing",test.terms="character"),
          definition = function(xx,formula.full,test.terms,model.dat,colorgroup=NULL,legendpos="topright")
{
  # test for model.dat
  if(!is.data.frame(model.dat))
    stop("'model.dat' has to be a data frame")

  # test for 'test.terms'
  terms.all <- test.terms
  D.full    <- model.matrix(formula.full, model.dat)
  terms.all <- colnames(D.full)

  # are all terms variables compatible with 'model.dat'?
  if(!all(test.terms %in% terms.all))
    stop("'test.terms' are not compatible with the specified models")

  D.red  <- D.full[,!(colnames(D.full) %in% test.terms), drop=FALSE]

  # basic analysis
  res <- reduSQ(xx=xx,formula.full=formula.full,D.red=D.red,model.dat=model.dat)
  redu.SSQ.Genes <- res$redu.genes
  msE.genes      <- res$mse

  # plot
  plotgenes(xx=xx,model.dat=model.dat,colorgroup=colorgroup,redu.SSQ.Genes=redu.SSQ.Genes,msE.genes=msE.genes,legendpos=legendpos)
}
)

