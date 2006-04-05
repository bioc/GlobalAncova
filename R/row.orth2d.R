"row.orth2d" <-
function(xx, D)
{
# this function computes residuals for given response 'xx' and design matrix D
  xx %*% (diag(dim(D)[1]) - D %*% solve(t(D) %*% D) %*% t(D))
}

