
#----------------------------------------------------------------------------------
# ARCH(p) MODEL
#----------------------------------------------------------------------------------

# name:         SVMgen.R
# description:  generate a sample belongs ARCH model
# inputs:  
#               alpha <- SVM model parameters
#               error <- vector of noise in eq x = error*sigma
#               n     <- sample size, we use 100 iteration for chain burning
#                      this value must be grater than 100. Real size n - 100   
# output:       a sample of size n
#
# creadet by:   Carolina Gamboa


library(tseries)


genARCH <- function(alpha, error, n){
  
  p <- length(alpha)-1 # (p +1) is the number of parameters

  if (n <= 100){
    cat("error sample size!! \n")  # set 100 itera for burning the chain
  } 

  else {
    x     <- double(n)
    for(j in 1:p){
      x[j]  <- rnorm(1, sd = sqrt(alpha[1]/(1-sum(alpha[-1]))))
    }
    
    for(i in (p+1):n)  # Generate ARCH(1) process
    {
      sigma <- NULL
      for(j in 1:p){
        sigma <- sum(sigma, alpha[j+1]*x[i-j]^2)
      }
      sigma <- sqrt(alpha[1] + sigma)

      x[i] <- error[i]*sigma
    }
    
    x <- ts(x[101:n])
    return(x)
  }
  
}



# n = 460
# a2 <- c(0.0022, 0.322, 0.074, 0.093)
# a1 <- c(0.0126, 0.3526)

