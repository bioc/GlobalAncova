"effectnames" <-
function(D.full, D.red)
{
# this function extracts the name of the tested effect
# (because there can be problems with interactions:
#  e.g '~group*covar' yields 'group,covar,group:covar' but ~covar*group yields '..,covar:group')
 names.all  <- union(colnames(D.full), colnames(D.red))
 no.effect  <- intersect(colnames(D.full), colnames(D.red))

 # are there interaction effects?
 interact   <- grep(":", names.all)
 if(length(interact) > 0)
 {
   # split the interaction terms and look if components are the same
   split.names <- strsplit(names.all[interact], ":")
   id          <- sapply(split.names, sort)
   for(i in 1:length(split.names))
   {
     for(j in 1:length(split.names))
     {
       # if there are identical columns, these represent the same factor and hence
       #  are added to the 'non-effect-terms'
       if(identical(id[,i],id[,j]) & i!=j)
         no.effect <- c(no.effect, names.all[interact[c(i,j)]])
     }
   }
 }
 # the remaining terms are the tested effects
 effect.names <- names.all[!(names.all %in% no.effect)]
 return(effect.names)
}




