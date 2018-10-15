

estimaP <- function(z, k){

  N <- length(z)
  M_d <- matrix(nrow = N-k, ncol = k+1)
  
  for(i in 1:(k+1)){
    M_d[, i]       <- z[i:(N-k+i -1)]
  }
  
  M_d <- data.frame( M_d)
  names(M_d) <- c(paste("lag", k:0, sep = ""))
  
  bic1 <- data.frame()
  for (jj in 1:(k)){
    
    #b <- BIC(model <- (lm(lag0 ~ . , data = M_d[, (k+1-jj):(k+1) ])))
    b <- BIC(arima(z, order = c(jj, 0, 0), method = "ML"))
    bic1 <- rbind(bic1, c(jj, b))
  }
  
  names(bic1) <- c("k", "BIC")
  p <- which.min(bic1$BIC)
  return(bic1)
}


