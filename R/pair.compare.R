
# function to compute pairwise comparisons for different factor levels
pair.compare<-function(xx,formula,group,model.dat=NULL,test.genes=NULL,perm=10000) # !! perm=10000 statt 1000
{
# !!
  # if just one gene should be tested
  if(is.vector(xx))
      xx <- t(as.matrix(xx))

  if(is.null(rownames(xx)))
    rownames(xx) <- 1:dim(xx)[1]
  if(is.null(test.genes))
    test.genes <- list(rownames(xx))
  if(!is.list(test.genes))
    test.genes <- list(test.genes)

  group.var<-factor(model.dat[,group])
  model.dat[,group]<-group.var
  model.work<-model.dat
  n.tests    <- length(test.genes)
# !!

  n.subjects <- dim(xx)[2]

  levels<-levels(group.var)
  nlevels<-length(levels)

  D.full<-model.matrix(formula,model.dat)

  # checking for NA's
  if(dim(D.full)[1] < n.subjects)
    stop("Missing values in the model variables")

# !!
  res <- list()
  for(j in 1:n.tests){
    # select test.genes
    xx2 <- xx[test.genes[[j]], ,drop=FALSE]
    n.genes <- nrow(xx2)
    df.full<-dim(D.full)[2]*n.genes
    ssq.full<-sum(row.orth2d(xx2,D.full)^2)
# !!

    red.ssq<-vector()
    df<-names<-vector()

    k<-1

    F.perm <- numeric()
    for (i in 1:(nlevels-1)){
     for(j in (i+1):nlevels){
      #set 2 levels to 1
      model.work[,group]<-set.pair(group.var,i,j)
      D.red<-model.matrix(formula,model.work)

      rr <- row.orth2d(xx2,D.red)                # !!  xx2
      red.ssq[k]<-sum(rr^2) - ssq.full
      names[k]<-paste(levels[i],":",levels[j])
      df[k]<-df.full-dim(D.red)[2]*n.genes

      F.perm.g <- numeric()
      for(g in 1:perm)
      {
        ord <- sample(n.subjects)
        ssq.full.perm <- sum(row.orth2d(rr[,ord], D.full)^2)
        red.ssq.perm <- sum(row.orth2d(rr[,ord], D.red)^2) - ssq.full.perm
        MS.full.perm <- ssq.full.perm / (n.subjects*n.genes-df.full)
        MS.red.perm <- red.ssq.perm / df[k]
        F.perm.g[g] <- MS.red.perm / MS.full.perm
      }
      F.perm <- rbind(F.perm, F.perm.g)

      model.work<-model.dat
      k<-k+1
    }
  }

  red.ssq<-c(red.ssq,ssq.full)
  names(red.ssq)<-c(names,"error")
  df<-c(df,n.subjects*n.genes-df.full)
  ms<-red.ssq/df
  f<-c(ms[1:(k-1)]/ms[k],NA)

  #f.count <- function(F.perm.g, f)
  #  sum(F.perm.g > f)

  p.perm.mat <- sweep(F.perm, 1, f[-k], ">")
  p.perm <- rowSums(p.perm.mat) / perm

# !!
  res <- c(res, list(cbind(SSQ=red.ssq,df=df,MS=ms,F=f,p.perm=c(p.perm,NA))))
  }
   if(n.tests == 1)
     res <- res[[1]]
   else
     names(res) <- names(test.genes)
  return (res)
# !!
}


set.pair<-function(var,i,j)
{
  levels<-levels(var)
  var[which(var==levels[j])]<-levels[i]
  levels.new<-levels[-j]
  return(factor(var,levels=levels.new))
}
