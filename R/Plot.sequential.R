# plots the genewise sequential sum of squares adjusted by model SSQ (model - intercept =100%)
Plot.sequential <- function(xx, formula, model.dat=NULL, test.genes=NULL, name.geneset="")  # !! 'test.genes', 'name.geneset' statt 'name.pathw'
{
# !!
 if(!is.null(test.genes))
   xx <- xx[test.genes,,drop=FALSE]
# !!
 ANOVA.g<-decomp.ssq.genewise(xx,formula,model.dat)
 ssq<-ANOVA.g$SSQ
 nt<-length(ANOVA.g$terms)-1   #number of terms
 p<-dim(ssq)[1]-1              #number of genes
 ssq<-ssq[p:1,]

 #adjustment
 model.adj<-rowSums(ssq[,1:nt])-ssq[,"Intercept"]
 ssq<-ssq/model.adj

 title<-paste("Sequential Sum of Squares",name.geneset)
 if (nt<=5){color<-my.colors(5,reverse=T)[1:(nt-1)]}
 else      {color<-my.colors(nt-1,reverse=T)}

 op<-par(las=1,mar=c(3,7,4,0),cex.axis=.7)
 barplot(t(ssq[,2:nt]),xpd=F,width=0.715,space=.4,xlim=c(0,1.4),ylim=c(0,p) ,horiz=T,main=title,legend.text=F,axes=F,col=color)
 legend(x=c(1.02,1.4),y=c(p,0),ANOVA.g$terms[2:nt],fill=color)
 axis(1,0:5/5)
 axis(2,1:p-.4,F)
 #op
}



my.colors<-function(nlevels=50,mat=NA, range=NA, reverse=F)
{
#produces a rainbow colourscale from darkblue to darkred
	ncolors<-nlevels
	if (is.numeric(mat))	{ncolors<-length(pretty(range(mat,na.rm=TRUE),nlevels))-1}
	if (is.numeric(range))	{ncolors<-length(pretty(range,nlevels))-1}

	nred<-as.integer(ncolors/10)
  if (nred==0) {col<-rainbow(ncolors-2*nred,start=0,end=4/6)}
	else
	{
		lcol<-length(col)
		middle<-ceiling(lcol/2)
		col<-rainbow(ncolors-2*nred+1,start=.95,end=4/6,gamma=.45)
		col<-col[-middle]
	}
	red<-rgb(seq(150,255,length=nred+1),rep(0,nred+1),rep(0,nred+1),maxColorValue=255)
	red<-red[-(nred+1)]
	blue<-rgb(rep(0,nred+1),rep(0,nred+1),seq(255,150,length=nred+1),maxColorValue=255)
	blue<-blue[-1]
	col<-c(red,col,blue)
	#col<-topo.colors(nlevels)
	if(!reverse)col<-col[(length(col)):1]
	col
}


