"expr.test" <-
function(xx,formula.full,formula.red=NULL,D.red=NULL,model.dat,perm=10000,test.genes=NULL)
{
    # if just one gene should be tested
    if(is.vector(xx))
        xx <- t(as.matrix(xx))

    if(is.null(rownames(xx)))         
        rownames(xx) <- 1:nrow(xx)
    if(is.null(test.genes))
        test.genes <- list(rownames(xx))
    if(!is.list(test.genes))
        test.genes <- list(test.genes)
    if(is.numeric(unlist(test.genes)))
        test.genes <- lapply(test.genes, as.character)

    xx2        <- xx[unlist(unique(test.genes)),,drop=F]
    #N.Genes     <- nrow(xx2)
    N.Genes    <- sapply(test.genes, length)
    N.Subjects <- ncol(xx2)
  
    # design matrices
    D.full     <- model.matrix(formula.full, data=model.dat)
    if(is.null(D.red))
        D.red  <- model.matrix(formula.red,  data=model.dat)
    N.par.full <- ncol(D.full)
    N.par.red  <- ncol(D.red)

    # checking for NA's
    if(nrow(D.full) < N.Subjects)
        stop("Missing values in the model variables")
     
    # degrees of freedom
    DF.full  <- N.Genes * (N.Subjects - N.par.full)
    DF.extra <- N.Genes * (N.par.full - N.par.red)        
  
    # sums of squares
    genewiseSS <- genewiseGA(xx2, D.full, D.red) 
    SS.full    <- sapply(test.genes, function(x) sum(genewiseSS[x,"denominator"]))
    SS.extra   <- sapply(test.genes, function(x) sum(genewiseSS[x,"nominator"]))
    MS.full    <- SS.full / DF.full
    MS.extra   <- SS.extra / DF.extra
    
    # F statistic and theoretical p-value for each gene group
    F.value    <- MS.extra / MS.full
    p.value    <- pf(F.value, DF.extra, DF.full, lower.tail=F)

    # resampling step
    p.perm <- resampleGA(xx2, D.full, D.red, perm, test.genes, F.value, DF.full, DF.extra) 

    # results
    test.result           <- cbind(F.value, p.value, p.perm)
    colnames(test.result) <- c("F.value", "p.value", "p.perm")

    if(length(test.genes) == 1) {
        effect.names          <- effectnames(D.full, D.red)
        ANOVA.tab             <- cbind(c(SS.extra,SS.full), c(DF.extra,DF.full), c(MS.extra,MS.full))
        dimnames(ANOVA.tab)   <- list(c("Effect", "Error"), c("SSQ", "DF", "MS"))
        result <- list("effect"=effect.names,"ANOVA"=ANOVA.tab,"test.result"=t(test.result),"terms"=colnames(D.full))
    }

    else 
        result <- cbind("genes"=N.Genes, test.result)
    
    return(result)
}
