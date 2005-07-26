"Plot.genes" <-
function (xx, group, covars = NULL)
{
# basic analysis

   N.subjects   <- dim(xx)[[2]]
   N.genes      <- dim(xx)[[1]]
   X.full       <- cbind(1, group, covars)
   X.redu       <- cbind(rep(1, N.subjects), covars)
   hat.matrix   <- function(x)
   {
       x %*% solve(t(x) %*% x) %*% t(x)
   }
   H.full       <- hat.matrix(X.full)
   H.redu       <- hat.matrix(X.redu)
   I            <- diag(N.subjects)
   X.addi       <- matrix(group %*% (I - H.redu), N.genes, N.subjects, byrow <- TRUE)
   project      <- function(v, w)
   {
       (sum(v * w)/sum(w * w)) * w
   }
   rr.full      <- xx %*% (I - H.full)
   rr.redu      <- xx %*% (I - H.redu)

#   modification if test for interaction only is intended
#   rr.addi      <- rr.redu - project(rr.redu, X.addi)

# reduction sum of squares and mean square error

   redu.sq           <- rr.redu^2-rr.full^2

   redu.SSQ.subjects <- apply(redu.sq,2,sum)
   redu.SSQ.genes    <- apply(redu.sq,1,sum)
   msE.genes         <- rowSums(rr.full^2)/(N.subjects-dim(X.full)[2])

# determination of upregulation

   up         <- 0 < (apply(xx[,group==1],1,mean)-apply(xx[,group==0],1,mean))

# plotting results

   horizontal.bars(
    x         = rev(redu.SSQ.genes),
    xlabel    = "Reduction in Sum of Squares",
    ylabel    = "Genes",
    color     = 3-rev(up),
    bar.names = rev(rownames(xx))
   )

#pp       <- rev(sort(c(-.5+(1:N.genes),.5+(1:N.genes)))) <<<<DAS WAR LEIDER FALSCH !!!!!
pp        <- sort(c(-.5+(1:N.genes),.5+(1:N.genes)))
vv        <- rev(rep(msE.genes,rep(2,N.genes)))
lines(vv,pp,type="s",lwd=2)

legend("topright", c("higher expression in group 0", "higher expression in group 1"), col=c(3,2), pch=15) #!
}

