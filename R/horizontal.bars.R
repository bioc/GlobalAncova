"horizontal.bars" <-
function(x,xlabel="",ylabel="",color=NULL,labelsize=.75,bar.names=NULL)
{
        # function for plotting horizontal bars with labels added at right margin
        # bars are determined by value of x which is assumed to be a vector
        # no formal check of variables performed
        # setting the plot region

        xlim    <- 0.05*c(-1,1)*range(x)+c(min(x),max(x))
        xlim[1] <- min(0,xlim[1])
        n       <- length(x)
        ylim    <- c(0,n+1)

        # enlarging right margin for bar.names

        if(!is.null(bar.names)&length(bar.names)==n)
    {
        names   <- TRUE
                   plot.new()
        w       <- 1.5 * max(strwidth(bar.names, "inches", labelsize))
        oldmai  <- par("mai")
                   par(mai=c(oldmai[1:3],max(w,oldmai[4])), new=TRUE)
    }
        # plotting bars with border=FALSE nothing appears color is NULL


plot(0,type="n",xlim=xlim,ylim=ylim,yaxt="n",xlab=xlabel,ylab=ylabel)
                rect(rep(0,n),(1:n)-.3,x,(1:n)+.3,col=color,border="white")
                box()

        # adding bar.names at right margin

        if(names)
   {
                axis(4,at=1:n,bar.names,cex.axis=labelsize,las=2)
                par(mai=oldmai)
   }
}

