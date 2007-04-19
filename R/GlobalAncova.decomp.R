# this function computes the sequential and type III variance decomposition
# for a model specified in "formula". xx is the dataframe of the gene
# expression data and model.dat of the pheno data.
# It is possible to adjust the gene expressions by a global covariate
# specified in zz. This can be e.g. former gene expressions of the same
# subjects, or gene expressions from other parts of the body. By default
# the adjustment is performed with a common beta for all genes.
# Alternatively, an adjustment with a different beta for each gene can be
# used by setting zz.per.gene to TRUE.
# 'method': c("sequential","type3","all")

GlobalAncova.decomp <- function(xx, formula, model.dat=NULL, method=c("sequential","type3","all"),
              test.genes=NULL, genewise=FALSE, zz=NULL, zz.per.gene=FALSE){
  # if just one gene should be tested
  if(is.vector(xx))
      xx <- t(as.matrix(xx))

  if(is.null(rownames(xx)))
    rownames(xx) <- 1:dim(xx)[1]
  if(is.null(test.genes))
    test.genes <- list(rownames(xx))
  if(!is.list(test.genes))
    test.genes <- list(test.genes)

 # nicht-genweise
 if(!genewise){
   res <- list()
   for(j in 1:length(test.genes)){
      # select test.genes
      xx2 <- xx[test.genes[[j]], ,drop=FALSE]
      zz2 <- zz[test.genes[[j]], ,drop=FALSE]

      res <- c(res, list(decomp.ssq(xx=xx2, formula=formula, model.dat=model.dat,
                    zz=zz2, zz.per.gene=zz.per.gene, method=method)))
   }

   if(length(test.genes) == 1)
     res <- res[[1]]
   else
     names(res) <- names(test.genes)
 }

 # genweise
 else{
   if(length(test.genes) > 1)
     stop("genewise analysis only valid with one gene group")
   method <- match.arg(method)  
   if(method != "sequential")
     warning("genewise analysis yields only sequential decomposition")

   # select test.genes
   xx2 <- xx[test.genes[[1]], ,drop=FALSE]
   zz2 <- zz[test.genes[[1]], ,drop=FALSE]

   res <- decomp.ssq.genewise(xx=xx2, formula=formula, model.dat=model.dat) #,zz=zz2, zz.per.gene=zz.per.gene, method=method)
 }
 return(res)
}


################################################################################
################################################################################

# main function for the decomposition
decomp.ssq <- function(xx, formula, model.dat=NULL, method=c("sequential","type3","all"),  # !! 'method=...'
              zz=NULL, zz.per.gene=FALSE){

 #___1. adjustment for global covariates
 case <- 1
 ssq.total <- sum(xx^2)

 if (!is.null(zz))
  {
    if(!zz.per.gene){xx<-xx.adj(xx,zz);case<-2}
    else{xx<-xx.row.adj(xx,zz);case<-3}
  }

 # reduction ssq due to adjustment
 ssq.adj<-ssq.total-sum(xx^2)
 df.adj<- switch(case,0,1,dim(xx)[1])
 adjust<-matrix(c(ssq.adj,df.adj),1,2,dimnames=list("adjustment",c("ssq","df")))

 #___2. preparation
 # model matrix
  n.subjects<-dim(xx)[2]                     # number of observations
  p<-dim(xx)[1]                              # number of genes

  D<-model.matrix(formula,model.dat)

  # checking for NA's
  if(dim(D)[1] < n.subjects)
    stop("Missing values in the model variables")

  D.red<-matrix(1,n.subjects,1)             # design matrix for ~1

 #computing dfs using a "dummy anova"
  dummy.formula<-as.formula(paste("rep(1,n.subjects)~",as.character(formula[2])) )
  dummy.anova<-anova(lm(dummy.formula,model.dat))

  df<-dummy.anova$Df
  nt<- length(df)
  df<-df[-nt]
  terms<-rownames(dummy.anova)[-nt]

  if("(Intercept)" %in% colnames(D)){df<-c(1,df); terms<-c("Intercept",terms)}

  nt<-length(terms)                # number of model terms

  seq.ssq<-t3.ssq<-rep(0,nt+1)
  names(seq.ssq)<-names(t3.ssq)<-c(terms,"error")

 #SSmean
  seq.ssq[1]<-sum((xx%*%matrix(1/n.subjects,n.subjects,n.subjects))^2)

#___3. sequential decomposition

 #calculation of model sum of squares
  position<-cumsum(df)
  for (i in 2:nt)
  {
   D.full<-D[,1:position[i]]
   seq.ssq[i]<-red.ssq(xx,D.full,D.red)[1]
   D.red<-D.full
  }

 # SSresidual
  seq.ssq[nt+1]<-sum((row.orth2d(xx,D))^2)

 # df.residual
  df<-df*p
  df[nt+1]<-n.subjects*p-sum(df[1:nt])-adjust[2]

 # MS
  seq.ms<-seq.ssq/df

 # F-value
  seq.f<-c(seq.ms[1:nt]/seq.ms[nt+1],NA)

 # p-value
  seq.p <- pf(seq.f, df, df[nt+1], lower.tail=FALSE)

#___4. type III decomposition

  pos.start<-position
  pos.end <- c(1,(pos.start+1)[-nt])

  for (i in 1:nt)
  {
    test.cols<-pos.start[i]:pos.end[i]
    D.red<-D[,-test.cols,drop=F]
    t3.ssq[i]<-sum((row.orth2d(xx,D.red))^2)-seq.ssq[nt+1]
  }
 # SSresidual
    t3.ssq[nt+1]<-seq.ssq[nt+1]

 # MS
    t3.ms<-t3.ssq/df

 # F-value
    t3.f<-c(t3.ms[1:nt]/t3.ms[nt+1],NA)

 # p-value
    t3.p<-pf(t3.f,df,df[nt+1],lower.tail=F)

#___5. return ANOVA table

 sequential=cbind(seq.ssq,df,seq.ms,seq.f,seq.p)
 typeIII=cbind( t3.ssq,df, t3.ms, t3.f, t3.p)
 colnames(sequential)<-colnames(typeIII)<-c("SSQ","df","MS","F","p")

 res<-list (adjustment=adjust,
              sequential=sequential,
              typeIII=typeIII)

  method <- match.arg(method)        
  switch(method, sequential=res$sequential, type3=res$typeIII, all=res)
}



# Similarly to decomp.ssq, this function computes the decomposition of the
# given model, but on a genewise basis. The result is used by plot.ssq.genewise
decomp.ssq.genewise <- function(xx, formula, model.dat=NULL) {

  n.subjects<-dim(xx)[2]                     # number of observations
  p<-dim(xx)[1]                              # number of genes

  D<-model.matrix(formula,model.dat)

  # checking for NA's
  if(dim(D)[1] < n.subjects)
    stop("Missing values in the model variables")

  D.red<-matrix(1,n.subjects,1)             # design matrix for ~1

  dummy.formula<-as.formula(paste("rep(1,n.subjects)~",as.character(formula[2])) )
  dummy.anova<-anova(lm(dummy.formula,model.dat))

  df<-dummy.anova$Df
  nt<- length(df)
  df<-df[-nt]
  terms<-rownames(dummy.anova)[-nt]

  if("(Intercept)" %in% colnames(D)){df<-c(1,df); terms<-c("Intercept",terms)}

  nt<-length(terms)             # number of model terms

  ssq<-matrix(0,p,nt+1)
  rownames(ssq)<-rownames(xx)
  colnames(ssq)<-c(terms,"error")

  #SSmean
  ssq[,1]<-rowSums((xx%*%matrix(1/n.subjects,n.subjects,n.subjects))^2)

  #SSmodel
  position<-cumsum(df)
  for (i in 2:nt)
  {
    D.full<-D[,1:position[i]]
    ssq[,i] <- genewiseGA(xx, D.full, D.red=D.red)[,1]
    D.red<-D.full
  }

  #SSresidual
  ssq[,nt+1]<-rowSums((row.orth2d(xx,D))^2)

  all<-colSums(ssq)
  ssq<-rbind(ssq,all)

  #df.residual
  df[nt+1]<-n.subjects-sum(df[1:nt])
  df<-rbind(gene=df,all=df*p)

  colnames(df)<-colnames(ssq)

  #MS
  ms.gene<-t(t(ssq[1:p,,drop=F])/df["gene",])
  ms.all <-t(t(ssq["all",,drop=F])/df["all",]) #different dfs for row "all"
  ms<-rbind(ms.gene,ms.all)

  #F-values
  f<-ms[,1:nt]/ms[,nt+1]

  #p-values
  p.values<-t(pf(t(f),df[1:nt],df[nt+1],lower.tail=F))

  return(list(terms=colnames(ssq),SSQ=ssq,df=df,MS=ms,F=f,p=p.values))
}





# computes the I - Hat matrix for a given design matrix D
IminusHat <- function(D)
   diag(dim(D)[1]) - hat.matrix(D)


# computes the reduction sum of squares between reduced and full model and the full sum of squares
red.ssq <- function(xx, D.full, D.red){
  R.full   <- row.orth2d(xx,D.full)
  R.red    <- row.orth2d(xx,D.red)
  ssq.full <- sum(R.full^2)
  ssq.red  <- sum(R.red^2)
  red.ssq  <- ssq.red - ssq.full
  c(red.ssq, ssq.full)
}


# functions for the adjustment with a global covariate
xx.adj <- function(x,z){
  beta <- sum(x*z) / sum(z*z)
  x - beta * z
}

xx.row.adj <- function(xx, zz){
  beta <- rowSums(xx*zz) / rowSums(zz*zz)
  xx - sweep(zz, 1, beta, "*")
}







