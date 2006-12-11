"genewiseGA" <- 
function(xx, D.full, D.red)
# this function computes nominator and denominator of GlobalAncova statistic for each single gene
{  
    R.full <- row.orth2d(xx, D.full)
    R.red  <- row.orth2d(xx, D.red)
 
    SS.red.i <- rowSums(R.red*R.red)

    # denominator: residual sum of squares in the full model
    SS.full.i <- rowSums(R.full*R.full)

    # nominator: extra residual sum of squares
    SS.extra.i <- SS.red.i - SS.full.i
    
    return(cbind(nominator=SS.extra.i, denominator=SS.full.i))
}
