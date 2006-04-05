"group2formula" <-
function(group, group.name, covars, covar.names)
{
# this function is used to derive 'formula.full', 'formula.red' and 'model.dat' out of 'group' and 'covars'
# group: group variable
# group.name: character: name of group variable
# covars: covariate information
# covar.names: character: names of covraiates

  # model matrix
  model.dat  <- data.frame(cbind(group, covars))
  if(is.null(covars))
    names(model.dat)[1] <- group.name
  else
    names(model.dat) <- c(group.name, covar.names)

  # model formulas
  formula.full <- paste("~", group.name)

  # if there are no covariates
  if(is.null(covars))
  {
    formula.full <- as.formula(formula.full)
    formula.red  <- ~ 1
  }
  else
  {
    formula.full <- as.formula(paste(formula.full, "+", paste(covar.names, collapse="+")))
    formula.red  <- as.formula(paste("~", paste(covar.names, collapse="+")))
  }
  return(list(formula.full=formula.full, formula.red=formula.red, model.dat=model.dat))
}

