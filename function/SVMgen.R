
#----------------------------------------------------------------------------------
# Stochastic model
#
# x_t = sigma_t * error
# log(sigma_t^2) = a[1] + a[2]*log(sigma_{t-1}^2) + .. + phi
#
#----------------------------------------------------------------------------------
#
# name:         SVMgen.R
# description:  generate a sample belongs stochastic volatility model
# inputs:  
#               a   <- SVM model parameters
#               phi <- vector of noise in eq for ln(sigma^2) from SVM
#               n   <- sample size, we use 100 iteration for chain burning
#                      this value must be grater than 100. Real size n - 100   
# output:       a sample of size n
#
# creadet by:   Carolina Gamboa
#
#----------------------------------------------------------------------------------

genSVM <- function(a, phi, n){
  # ARCH(2) coefficients
  if (n <= 100){
    cat("error sample size!! \n")
  } 
  else {
  error <- rnorm(n)
  sigma2_0 <- rlnorm(1, a[1]/(1-a[2]), sqrt(1/(1-a[2])))
  y <- NULL
  x <- NULL
  for(i in 1:n){
    y[i] <- exp(a[1] + a[2]*log(sigma2_0) + phi[i])
    sigma2_0 = y[i]
    x[i] <- error[i]*sqrt(y[i])
  }
  
  x <- ts(x[101:n])
  return(x)
  }
}

# n <- 90
# a <- c(0.831, 0.685)
# phi <- rnorm(n, 0, 0.265)
# x <- genSVM(a, phi, n)
# plot(x)

