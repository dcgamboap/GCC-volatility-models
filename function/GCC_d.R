
GCC_d <- function(x, y, k){
  N <- length(x)
  
  if(k == 0){
    k <- floor(N/10)
  }
  # Falta entender cómo determinar el valor de k 
  
  M_xy <- matrix(nrow = N-k, ncol = 2*(k+1))
  
  for(i in 1:(k+1)){
    M_xy[, i]       <- x[i:(N-k+i-1)]
    M_xy[, i+(k+1)] <- y[i:(N-k+i-1)]
  }
  
  M_x <- M_xy[, 1:(k+1)]
  M_y <- M_xy[, (k+2):(2*(k+1))]
  
  # Sample corr
  
  R_xy <- cor(M_xy)
  R_x  <- cor(M_x)
  R_y  <- cor(M_y)
  
  GCC <- 1 - det(R_xy)^(1/ (1 * (k+1)))  / (det(R_x)^(1/ (1 * (k+1))) * det(R_y)^(1/ (2 * (k+1))))
  return(GCC)
  
}