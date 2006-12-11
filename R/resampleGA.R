"resampleGA" <-
function(xx, D.full, D.red, perm, test.genes, F.value, DF.full, DF.extra) 
# this function conducts the permutation test
{
    N.Subjects  <- ncol(xx)
    rr     <- row.orth2d(xx, D.red)
    rr.cov <- (diag(nrow(D.red)) - D.red %*% solve(t(D.red) %*% D.red) %*% t(D.red))
    rr     <- rr %*% diag(1 / sqrt(diag(rr.cov)))
    
    count  <- numeric(length(test.genes))
    for(i in 1:perm) {
      ord          <- sample(N.Subjects)
      genewSS.perm <- genewiseGA(rr[,ord,drop=F], D.full, D.red) 
      F.perm       <- sapply(test.genes, function(x) sum(genewSS.perm[x,1]) / sum(genewSS.perm[x,2])) / (DF.extra / DF.full)
      count        <- count + (F.perm > F.value)
    }
    
    return(count / perm)
}
