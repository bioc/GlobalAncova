"expr.test" <-
function(xx,formula.full,formula.red=NULL,D.red=NULL,model.dat,perm=10000,test.genes=NULL)
{
  # if just one gene should be tested
  if(is.vector(xx))
    xx <- t(as.matrix(xx))

  if(is.null(rownames(xx)))
    rownames(xx) <- 1:dim(xx)[1]
  if(is.null(test.genes))
    test.genes <- list(rownames(xx))
  if(!is.list(test.genes))
    test.genes <- list(test.genes)

  many.pathways           <- matrix(0, length(test.genes), 4)
  dimnames(many.pathways) <- list(names(test.genes), c("genes","F.value","p.value","p.perm"))

  for(j in 1:length(test.genes))
  {
    xx2 <- xx[test.genes[[j]], ]

    # if pathway contains just one gene
    if(is.vector(xx2))
      xx2 <- t(as.matrix(xx2))

    N.Genes 	<- dim(xx2)[1]
    N.Subjects 	<- dim(xx2)[2]

    # Design matrices
    D.full 	<- model.matrix(formula.full, data=model.dat)
    if(is.null(D.red))
      D.red 	<- model.matrix(formula.red,  data=model.dat)
    N.par.full	<- dim(D.full)[2]
    N.par.red	<- dim(D.red)[2]

    # extra sum of squares & full residual sum of squares
    extra.ssq	<- red.ssq(xx2, D.full, D.red)

    DF.ssq.res.full <- N.Genes * (N.Subjects - N.par.full)
    DF.ssq.res.red  <- N.Genes * (N.Subjects - N.par.red)
    DF.extra.ssq    <- DF.ssq.res.red - DF.ssq.res.full
    DF		    <- c(DF.extra.ssq, DF.ssq.res.full)

    # Test
    MS		<- extra.ssq /c(DF)
    F.value	<- MS[1] / MS[2]
    p.value	<- 1 - pf(F.value, DF[1], DF[2])

    # Resampling step
    rr	   <- row.orth2d(xx2, D.red)
    rr.cov <- (diag(dim(D.red)[1]) - D.red %*% solve(t(D.red) %*% D.red) %*% t(D.red))
    rr	   <- rr %*% diag(1 / sqrt(diag(rr.cov)))
    count  <- 0
    for(i in 1:perm)
    {
      ord         <- sample(N.Subjects)
      MS.resample <- red.ssq(rr[,ord], D.full, D.red) / DF
      count       <- count + ((MS.resample[1] / MS.resample[2]) > F.value)
    }
    p.perm <- count / perm

    # Results
    #effect.names          <- colnames(D.full)[!colnames(D.full) %in% colnames(D.red)]
    effect.names          <- effectnames(D.full, D.red)
    ANOVA.tab 	          <- cbind(extra.ssq, DF, MS)
    dimnames(ANOVA.tab)   <- list(c("Effect", "Error"), c("SSQ", "DF", "MS"))
    test.result   	  <- rbind(F.value, p.value, p.perm)
    dimnames(test.result) <- list(c("F.value", "p.value", "p.perm"), "")
    #ANOVA.list            <- list("effect"="","ANOVA"=ANOVA.tab,"test.result"=test.result,"terms"=colnames(D.full))
    ANOVA.list            <- list("effect"=effect.names,"ANOVA"=ANOVA.tab,"test.result"=test.result,"terms"=colnames(D.full))

    many.pathways[j, 1]   <- as.numeric(length(test.genes[[j]]))
    many.pathways[j, 2:4] <- round(test.result[1:3], 4)
  }

  if(length(test.genes) == 1)
    return(ANOVA.list)
  else
    return(many.pathways)
}

