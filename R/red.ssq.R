"red.ssq" <-
function(xx, D.full, D.red)
{
# this function computes the reduction sum of squares between reduced and full model
# and the full sum of squares
  R.full   <- row.orth2d(xx,D.full)
  R.red    <- row.orth2d(xx,D.red)
  ssq.full <- sum(R.full^2)
  ssq.red  <- sum(R.red^2)
  red.ssq  <- ssq.red - ssq.full
  c(red.ssq, ssq.full)
}

