
zeroAdjust3 <- function(x){
  N <- length(x)

 cual <- which(x==0) 
 for (jj in cual) {
   if(jj == N){
     x[jj] <- x[jj-1]
   } 
   if(jj == 1){
     x[jj] <- x[jj+1]
   } 
   if(jj != 1 & jj != N){
     x[jj] <- mean(c(x[jj-1], x[jj+1]))
   }
   
 }
 return(t(x))
}