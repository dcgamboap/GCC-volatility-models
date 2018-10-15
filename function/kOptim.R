
#----------------------------------------------------------------------------------
# K OPTIM
#----------------------------------------------------------------------------------
#
# name:         kOptim.R
# description:  search optim k to GGC measure using linear regressions  
# input: time series x and y   
#                
# output:	max value of k using BIC criteria, in the regressions: 
#                     x_0 ~ x1 +...+ xk + y0 +...+ yk
#                     y_0 ~ x0 +...+ xk + y1 +...+ yk

# 
# creadet by:   Carolina Gamboa
# v1: add a inf limit and calculate both regressions
#
#----------------------------------------------------------------------------------




kOptim <- function(x, y, k, kinf){
  
  
  N <- length(x)
  M_d <- matrix(nrow = N-k, ncol = 2*k+2)
  
  for(i in 1:(k+1)){
    M_d[, i]       <- x[i:(N-k+i-1)]
    M_d[, (i + k+ 1)] <- y[i:(N-k+i-1)]
  } 
  
  M_d <- data.frame( M_d)
  names(M_d) <- c(paste("x.lag", k:0, sep = ""), 
                  paste("y.lag", k:0, sep = ""))
  
  bic1 <- data.frame()

  for (jj in kinf:(k)){
    b1 <- BIC(model <- lm(x.lag0 ~ . , data = M_d[, c((k+1-jj):(k+1), 
                                                      (2*(k+1)-jj):(2*(k+1)))]))
    bic1 <- rbind(bic1, c(jj, b1))
  }
  
  names(bic1) <- c("k", "BIC")
  p1 <- which.min(bic1$BIC)
  
  if(p1 < k){
    bic2 <- data.frame()
    
    for (jj in kinf:(k)){
      b2 <- BIC(model <- lm(y.lag0 ~ . , data = M_d[, c((k+1-jj):(k+1), 
                                                        (2*(k+1)-jj):(2*(k+1)))]))
      bic2 <- rbind(bic2, c(jj, b2))
    }
    
    names(bic2) <- c("k", "BIC")
    p2 <- which.min(bic2$BIC)
    p  <- max(p1, p2)
  
  }
  
  if(!exists("p")) p <- p1
  
  return(p)
}

