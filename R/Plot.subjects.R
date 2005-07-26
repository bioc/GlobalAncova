"Plot.subjects" <-
function (xx, group, covars = NULL, sort = FALSE)
{
# basic analysis

    N.subjects  <- dim(xx)[[2]]
    N.genes     <- dim(xx)[[1]]
    X.full      <- cbind(1, group, covars)
    X.redu      <- cbind(rep(1, N.subjects), covars)
    hat.matrix  <- function(x)
    {
        x %*% solve(t(x) %*% x) %*% t(x)
    }
    H.full      <- hat.matrix(X.full)
    H.redu      <- hat.matrix(X.redu)
    I           <- diag(N.subjects)
    X.addi      <- matrix(group %*% (I - H.redu), N.genes, N.subjects, byrow <- TRUE)
    project     <- function(v, w)
    {
        (sum(v * w)/sum(w * w)) * w
    }
    rr.full     <- xx %*% (I - H.full)
    rr.redu     <- xx %*% (I - H.redu)

#   modification if test for interaction only is intended
#   rr.addi      <- rr.redu - project(rr.redu, X.addi)

# reduction in sum of squares
    redu.sq           <- rr.redu^2-rr.full^2

    redu.SSQ.subjects <- apply(redu.sq,2,sum)
    redu.SSQ.genes    <- apply(redu.sq,1,sum)

# plotting results

    if(is.null(colnames(xx)))						
	colnames(xx) <- seq(1:dim(xx)[2])				

    if(sort == TRUE)							
    {									
	horizontal.bars(
		x         = rev(redu.SSQ.subjects[order(group)]),	
		xlabel    = "Reduction in Sum of Squares",		
    ylabel    = "Subjects",
		color     = 3-rev(sort(group)),				
    bar.names = rev(colnames(xx)[order(group)])
		)
    	legend("topright", c("group 0", "group 1"), col=c(3,2), pch=15)	
    }									
 
    else								
    horizontal.bars(
        x         = rev(redu.SSQ.subjects),
        xlabel    = "Reduction in Sum of Squares",
        ylabel    = "Subjects",
	      color     = 3-rev(group),
        bar.names = rev(colnames(xx))
    )
    legend("topright", c("group 0", "group 1"), col=c(3,2), pch=15) 	
}

