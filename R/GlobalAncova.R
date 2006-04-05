setGeneric("GlobalAncova", function(xx,formula.full,formula.red,model.dat,group,covars=NULL,test.terms,perm=10000,test.genes=NULL)
           standardGeneric("GlobalAncova"))


################################# general function #############################

setMethod("GlobalAncova", signature(xx="matrix",formula.full="formula",formula.red="formula",
                           model.dat="ANY",group="missing",covars="missing",test.terms="missing"),
          definition = function(xx,formula.full,formula.red,model.dat,perm=10000,test.genes=NULL)
{
# this function computes a Global Ancova which tests for differential gene expression
# xx: expression matrix (rows=genes, columns=subjects)
# formula.full: model formula for the full model
# formula.red: model formula for the reduced model
# model.dat: data frame that contains the group and covariable information
# perm: number of permutations
# test.genes: vector of probeset names or a list where each element is a vector of
#             probeset names (e.g. for testing many pathways simultaneously)

  # test for model.dat
  if(!is.data.frame(model.dat))
    stop("'model.dat' has to be a data frame")

  expr.test(xx=xx,formula.full=formula.full,formula.red=formula.red,model.dat=model.dat,perm=perm,test.genes=test.genes)
}
)


########################## function for 2 group ################################

setMethod("GlobalAncova", signature(xx="matrix",formula.full="missing",formula.red="missing",
                           model.dat="missing",group="numeric",covars="ANY",test.terms="missing"),
          definition = function(xx,group,covars=NULL,perm=10000,test.genes=NULL)
{
# xx: expression matrix (rows=genes, columns=subjects)
# group: group variable
# covars: covariate information
# perm: number of permutations
# test.genes: vector of probeset names or a list where each element is a vector of
#             probeset names (e.g. for testing many pathways simultaneously)

  # parameter names
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

  # then apply the usual function
  expr.test(xx=xx,formula.full=formula.full,formula.red=formula.red,
                model.dat=model.dat,perm=perm,test.genes=test.genes)
}
)


############################# with 'test.terms' ################################

setMethod("GlobalAncova", signature(xx="matrix",formula.full="formula",formula.red="missing",
                           model.dat="ANY",group="missing",covars="missing",test.terms="character"),
          definition = function(xx,formula.full,test.terms,model.dat,perm=10000,test.genes=NULL)
{
# this function computes a Global Ancova which tests for differential gene expression
# xx: expression matrix (rows=genes, columns=subjects)
# formula.full: model formula for the full model
# test.terms: character vector of terms of interest
# model.dat: data frame that contains the group and covariable information
# perm: number of permutations
# test.genes: vector of probeset names or a list where each element is a vector of
#             probeset names (e.g. for testing many pathways simultaneously)

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

  D.full <- model.matrix(formula.full, data=model.dat)
  D.red  <- D.full[,!(colnames(D.full) %in% test.terms), drop=F]

  # then apply the usual function
  expr.test(xx=xx,formula.full=formula.full,D.red=D.red,
                model.dat=model.dat,perm=perm,test.genes=test.genes)
}
)


