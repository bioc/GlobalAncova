"GlobalAncova" <-
function(xx,group,covars=NULL,perm=10000,test.genes=NULL)
{
if(is.null(test.genes))
  test.genes               <- list(rownames(xx))
if(!is.list(test.genes))
  test.genes               <- list(test.genes)

# if just one gene should be tested
if(is.vector(xx))
  ANOVA.list               <- GlobalAncova.1gene(xx, group, perm)

else
{
many.pathways              <- matrix(0, length(test.genes), 4)
dimnames(many.pathways)    <- list(names(test.genes), c("genes", "F.value", "p.value.perm", "p.value.theo"))

for(j in 1:length(test.genes))			
{						
xx2 <- xx[test.genes[[j]], ]

# dimensions
N.subjects                 <- dim(xx)[[2]]
N.genes                    <- dim(xx2)[[1]]

# if a pathway contains only one gene
if(is.vector(xx2))
{
        ANOVA.list            <- GlobalAncova.1gene(xx2, group, perm)
        many.pathways[j, 1]   <- 1
        many.pathways[j, 2:4] <- ANOVA.list
}

else
{

# design matrices
# Full Model
X.full      <- cbind(1,group,covars)

# Reduced Model gene and covars only
X.redu      <- cbind(rep(1,N.subjects),covars)

# function for hat matrix
        hat.matrix <- function(x){x%*%solve(t(x)%*%x)%*%t(x)}

# Hat matrices according to models
H.full      <- hat.matrix(X.full)
H.redu      <- hat.matrix(X.redu)
I           <- diag(N.subjects)

# Matrix with rows made up by complement of group orthogonal to H.redu
        #X.addi<-matrix(group%*%(I-H.redu),N.genes,N.subjects, byrow<-TRUE) 	

# Function for projection onto one variable
        #project<-function(v,w){(sum(v*w)/sum(w*w))*w}				

# Residuals for calculating sequential sum of squares
rr.full     <- xx2 %*% (I-H.full)
rr.redu     <- xx2 %*% (I-H.redu)
        #rr.addi<-rr.redu-project(rr.redu,X.addi)				

# Sum of squares
SS.total    <- sum((xx2-(sum(xx2)/(N.genes*N.subjects)))^2)
SS.redu     <- sum(rr.redu*rr.redu)
#SS.addi    <- sum(rr.addi*rr.addi)
SS.full     <- sum(rr.full*rr.full)

SS.adjst    <- SS.total-SS.redu
#SS.group   <- SS.redu -SS.addi
#SS.gr.ge   <- SS.addi -SS.full
SS.gr.x.ge  <- SS.redu-SS.full

# DF
        DF.adjst <- dim(X.redu)[[2]]*N.genes-1
        DF.group <- 1
        DF.gr.ge <- N.genes-1
        DF.resid <- N.subjects*N.genes-1-DF.adjst-DF.group-DF.gr.ge

# ANOVA
#SS        <- c(SS.total,SS.adjst,SS.group+SS.gr.ge,SS.full)
SS         <- c(SS.total,SS.adjst,SS.gr.x.ge,SS.full)
DF         <- c(N.genes*N.subjects-1,DF.adjst,DF.group+DF.gr.ge,DF.resid)
MS         <- SS/DF
#F.value   <- c(NA,NA,((SS.group+SS.gr.ge)/(DF.group+DF.gr.ge))/MS[4],NA)
F.value    <- c(NA, NA, MS[3]/MS[4], NA)
        
ANOVA.tab            <- cbind(SS, DF, MS)
dimnames(ANOVA.tab)  <- list(c("Total","Genes adjusted","GroupXGenes","Residual"), c("SS", "DF", "MS"))

# p-value by MC Approximation to Permutation Distribution
# if number of required permutations is small use enumeration 

# Check for enumeration or random selection of permutations
# for enumeration use Sandrin Dudoits mt.sample.label from multtest
perm.enum 	<- choose(N.subjects,sum(group==1))

if(perm < perm.enum)  
  method  	<- "mc"

else 
{ 
method	  	<- "exact"
require(multtest)
perm.mat  	<- mt.sample.label(group,B=perm.enum)
print(paste("enumerating all ",perm.enum," permutations"))
}

# Enable common code for both methods by func sample.use
sample.use 	<- function(i,method)
{
if(method=="mc") 
  sample(N.subjects)
else
  c((1:N.subjects)[perm.mat[i,]==0],(1:N.subjects)[perm.mat[i,]==1])
}

smaller.F 	<- 0

# Fix number of iterations
res.number      <- min(perm,perm.enum)
for(i in 1:res.number)
{

# use residuals under Ho: no group effects and interactions
order           <- sample.use(i,method)
xx.perm         <- rr.redu[,order]

# correct for permutation of covars if present
X.perm.full     <- X.full[order,]
X.perm.full[,2] <- group
H.perm.full     <- hat.matrix(X.perm.full)

rr.full.perm    <- xx.perm%*%(I-H.perm.full)
       
SS.resid        <- sum(rr.full.perm*rr.full.perm)

# use SS.redu from original fit
SS.gr.x.ge      <- SS.redu-SS.resid

smaller.F       <- smaller.F + ((SS.gr.x.ge/(DF.gr.ge+DF.group)/(SS.resid/DF.resid)) < F.value[3])
} 

p.value.perm    <- c(NA,NA,1-(smaller.F/res.number),NA)
p.value.theo    <- c(NA,NA,1-pf(F.value[3],(DF.group+DF.gr.ge),DF.resid),NA)

test.result           <- rbind(F.value[3], p.value.perm[3], p.value.theo[3])
dimnames(test.result) <- list(c("F.value", "p.value.perm", "p.value.theo"), "")
ANOVA.list            <- list("ANOVA.table"=ANOVA.tab,"test.result.GroupXGenes"=test.result)

many.pathways[j, 1]   <- as.numeric(length(test.genes[[j]]))
many.pathways[j, 2:4] <- round(test.result[1:3], 4)
}	
}
}

if(length(test.genes) == 1)
  return(ANOVA.list)	
else
  return(many.pathways)
}

