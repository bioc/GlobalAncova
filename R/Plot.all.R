
# function that combines 'Plot.sequential' and 'Plot.genes'
Plot.all<-function(xx,formula,model.dat=NULL, test.genes=NULL, name.geneset="")  
{
 def.par <- par(no.readonly = TRUE) # save default, for resetting...

 if(!is.null(test.genes))
   xx <- xx[test.genes,]

  mat<-matrix(1:4,2,2)
  heights<-c(.9,.1)
  widths <-c(.5,.5)
  layout(mat,widths,heights)

#screen(1)
  plot.ssq.genewise2(xx,formula,model.dat,name.geneset)  

#screen(2)
  plot.ssq.all(xx,formula,model.dat)#,"(signature genes)")  

#screen(3)
  redu<-reduSQ(xx,formula,formula.red=~1,model.dat=model.dat)
  # center labels
  n.genes<-dim(xx)[1]
  all.name<-paste("all",n.genes,"genes")
  names.plus <- c(names(redu$redu.genes),all.name)
  names(redu$redu.genes) <- center.labels(names.plus)[1:n.genes]
  plotgenes2(redu$redu.genes, redu$mse)

#screen4
  red.genes <- sum(redu$redu.genes) / n.genes
  names(red.genes) <- all.name
  red.mse <- sum(redu$mse) / n.genes
  plotallgenes(red.genes, red.mse)
  
  par(def.par)#- reset to default
}


################################################################################
################################################################################


center.labels<-function(names)
{
 slength<-sapply(strsplit(names,NULL),length)
 maxlength<-max(slength)
 diff<-as.integer((maxlength-slength))#/2)
 new.names<-vector()
 for(i in 1:length(names))
 {
 spacec<-paste(rep(" ",diff[i]),collapse="")
 new.names[i]<-paste(spacec,names[i],spacec,sep="")
 }
 new.names
}


plot.ssq.genewise2<-function(xx,formula,model.dat=NULL, name.geneset="")  # !! 'name.geneset' statt 'name.pathw'
{
  # plots the genewise sequential sum of squares
  # adjusted by model SSQ (model - intercept =100%)

 ANOVA.g<-decomp.ssq.genewise(xx,formula,model.dat)
 ssq<-ANOVA.g$SSQ
 nt<-length(ANOVA.g$terms)-1   #number of terms
 p<-dim(ssq)[1]-1              #number of genes
 ssq<-ssq[p:1,]

 model.adj<-rowSums(ssq[,1:nt])-ssq[,"Intercept"]
 ssq<-ssq/model.adj
 rownames(ssq)<-NULL
 title<-paste("Sequential Sum of Squares -",name.geneset)
 if (nt<=5){color<-my.colors(5,reverse=T)[1:(nt-1)]}
 else      {color<-my.colors(nt-1,reverse=T)}

 op<-par(las=1,mar=c(3,0,4,0),cex.axis=.7)
 barplot(t(ssq[,2:nt]),xpd=F,width=0.715,space=.4,xlim=c(-.4,1),ylim=c(0,p) ,horiz=T,main=title,legend.text=F,axes=F,col=color)
 legend(x=c(-.38,-.02),y=c(p,p/3),ANOVA.g$terms[2:nt],fill=color)
 legend(x=c(-.38,-.02),y=c(p/3,0),c("MS model","MS error"),pch=c("-",NA),col=c("wheat",1),lty=c(NA,1),lwd=c(NA,2),pt.cex=c(6,1))
 axis(1,0:5/5)
 axis(4,1:p-.4,F)
 par(op)
}


plot.ssq.all<-function(xx,formula,model.dat=NULL)#,name.pathw="")   !!
{
  # plots the genewise sequential sum of squares
  # adjusted by model SSQ (model - intercept =100%)
  # only one bar for all genes

# !!
# ANOVA<-decomp.ssq(xx,,,formula,model.dat)
 ANOVA <- decomp.ssq(xx=xx, formula=formula, model.dat=model.dat)
# !!

 ssq<-ANOVA[,"SSQ"]
 nt<-length(rownames(ANOVA))-1   #number of terms
 p<-dim(xx)[1]                   #number of genes

 model.adj<-sum(ssq[1:nt])-ssq["Intercept"]
 ssq<-as.matrix(ssq[2:nt]/model.adj)
 colnames(ssq)<-NULL #paste("all",p,"genes")

 #title<-paste("Sequential Sum of Squares",name.pathw)  !!

 if (nt<=5){color<-my.colors(5,reverse=T)[1:(nt-1)]}
 else      {color<-my.colors(nt-1,reverse=T)}

 op<-par(las=1,mar=c(3,0,0,0),cex.axis=.7)
 barplot(ssq,xpd=F,width=0.715,space=.4,xlim=c(-.4,1) ,horiz=T,legend.text=F,axes=F,col=color)
 axis(1,0:5/5)
 par(op)
}


plotgenes2 <- function(redu.MS.Genes, msE.genes)
{
  N.Genes <- length(redu.MS.Genes)
  if(is.null(names(redu.MS.Genes)))
    names(redu.MS.Genes) <- 1:N.Genes

  #bars
  op<-par(las=1,mar=c(3,7,4,1),cex.axis=.7)
  barplot(rev(redu.MS.Genes),ylim=c(0,N.Genes),xlim=c(0,max(c(redu.MS.Genes,msE.genes))*1.05),
        xpd=F,width=0.715,space=.4,horiz=T,main="Mean Sum of Squares",col="wheat",border=NA)
  axis(2,1:N.Genes-.4,F)

  # MSE-line
  pp    <- sort(c(-.45+(1:N.Genes-.4),.55+(1:N.Genes-.4)))
  vv    <- rev(rep(msE.genes,rep(2,N.Genes)))
  lines(vv,pp,type="s",lwd=2)

  p<-par("usr")
  rect(0,min(pp),p[2],max(pp)+.1)

  par(op)
}


plotallgenes <- function(redu.MS.Genes, msE.genes)
{
  #bars
  op<-par(las=1,mar=c(3,7,0,1),cex.axis=.7)
  barplot(rev(redu.MS.Genes),xlim=c(0,max(c(redu.MS.Genes*1.05,msE.genes*1.05))),
        xpd=F,width=0.715,space=.4,horiz=T,col="wheat",border=NA)

  # MSE-line
  p<-par("usr")[3:4]
  r<-diff(p)*.01
  p<-p+c(r,-r)
  pp    <- sort(rep(p,2))
  vv    <- c(0,rep(msE.genes,2),0)
  lines(vv,pp,type="s",lwd=2)

  par(op)
}



