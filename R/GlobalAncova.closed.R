setGeneric("GlobalAncova.closed", function(xx,test.genes,formula.full,formula.red,model.dat,group,covars=NULL,
                                  test.terms,previous.test=NULL,level=0.05,perm=10000)
           standardGeneric("GlobalAncova.closed"))
# xx: expression matrix (rows=genes, columns=subjects)
# test.genes: list of pathways (each containing a vector of genes) that shall be tested
#             and adjusted using the closed testing procedure
# formula.full: model formula for the full model
# formula.red: model formula for the reduced model
# model.dat: data frame that contains the group and covariable information
# group: group variable
# covars: covariate information
# test.terms: character vector of terms of interest
# previous.test: result of a GlobalAncova with many pathways simultaneously
# level: alpha
# perm: number of permutations


############################## allgemeine Funktion #############################

setMethod("GlobalAncova.closed", signature(xx="matrix",test.genes="list",formula.full="formula",formula.red="formula",
                          group="missing",covars="missing",test.terms="missing"),
          definition = function(xx,test.genes,formula.full,formula.red,model.dat,
                       previous.test=NULL,level=0.05,perm=10000)
{
  # test for model.dat
  if(!is.data.frame(model.dat))
    stop("'model.dat' has to be a data frame")

  GA.closed(xx=xx,test.genes=test.genes,formula.full=formula.full,formula.red=formula.red,model.dat=model.dat,
                       previous.test=previous.test,level=level,perm=perm)
}
)


########################## 'alte' Fkt. für 2 Gruppen ###########################


setMethod("GlobalAncova.closed", signature(xx="matrix",test.genes="list",formula.full="missing",formula.red="missing",
                          model.dat="missing",group="ANY",test.terms="missing"),
          definition = function(xx,test.genes,group,covars=NULL,previous.test=NULL,level=0.05,perm=10000)
{
  # parameter names
  group.name   <- deparse(substitute(group))

  if(is.null(dim(covars)))
    covar.names <- deparse(substitute(covars))
  else
    covar.names <- colnames(covars)

  # get formulas and 'model.dat' out of 'group' and 'covars'
  res          <- group2formula(group=group, group.name=group.name, covars=covars)
  formula.full <- res$formula.full
  formula.red  <- res$formula.red
  model.dat    <- res$model.dat

  GA.closed(xx=xx,test.genes=test.genes,formula.full=formula.full,formula.red=formula.red,model.dat=model.dat,
                       previous.test=previous.test,level=level,perm=perm)
}
)


################### allgemeine Funktion m. Angabe v. 'terms' ###################

setMethod("GlobalAncova.closed", signature(xx="matrix",test.genes="list",formula.full="formula",formula.red="missing",
                          group="missing",covars="missing",test.terms="character"),
          definition = function(xx,test.genes,formula.full,test.terms,model.dat,previous.test=NULL,level=0.05,perm=10000)
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

  D.full <- model.matrix(formula.full, data=model.dat)
  D.red  <- D.full[,!(colnames(D.full) %in% test.terms), drop=FALSE]

  GA.closed(xx=xx,test.genes=test.genes,formula.full=formula.full,D.red=D.red,model.dat=model.dat,
                       previous.test=previous.test,level=level,perm=perm)
}
)



