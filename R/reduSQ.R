"reduSQ" <-
function(xx, formula.full, formula.red=NULL, D.red=NULL, model.dat)
{
# this function computes the reduction in sum of squares for genes and subjects and the MSE
# 'xx' is an expression matrix (rows=genes, columns=subjects)
# 'formula.full' is a model formula for the full model
# 'formula.red' is a model formula for the reduced model
# 'model.dat' is a data frame that contains the group and covariable information

  N.Subjects     <- ncol(xx)

  # Designmatrices
  D.full     <- model.matrix(formula.full, data=model.dat)
  if(is.null(D.red))
    D.red    <- model.matrix(formula.red,  data=model.dat)

  N.par.full <- ncol(D.full)
  N.par.red  <- ncol(D.red)

  # Residuals
  R.full     <- row.orth2d(xx,D.full)
  R.red      <- row.orth2d(xx,D.red)

  # reduction sum of squares
  redu.sq           <- R.red^2 - R.full^2
  redu.SSQ.Genes    <- rowSums(redu.sq) / (N.par.full - N.par.red) 

  redu.SSQ.Subjects <- colSums(redu.sq)

  # mean square error
  msE.genes      <- rowSums(R.full^2) / (N.Subjects - N.par.full)

  return(list(redu.genes=redu.SSQ.Genes,redu.subjects=redu.SSQ.Subjects,mse=msE.genes,redu.MS.genes=redu.SSQ.Genes))
}
