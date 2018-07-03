
# plot contributions of individual variables to global test statistic, as Plot.genes but for gGlobalAncova

# - no coloring of bars (how could this be done in a senseful way...?)
# - no 'H0-line'

Plot.features <- function(data, formula.full, formula.red=~1, model.dat, Set, returnValues=FALSE, ...){
  # data: data.frame of variables to be tested in sets (columns=variables)
  # formula.full, formula.red: models to be compared
  # model.dat: data.frame of regressors, containing variables specified in formula.full and formula.red
  # Set: optional vector of variable names or indices, defining the set of variables to be plotted; if missing, all variables in 'data' are plotted
  # returnValues: shall variable-wise statistics = bar heights be returned?
  # ...: graphical parameters passed to 'barplot'
  
  if(missing(Set))
    Set <- 1:ncol(data)
  
  Tobs <- gGAteststats(data=data[, Set, drop=FALSE], formula.full=formula.full, formula.red=formula.red, model.dat=model.dat, perm=0)
  
  barplot(Tobs, horiz=TRUE, las=1, ...)
  
  if(returnValues)
    return(Tobs)
}
