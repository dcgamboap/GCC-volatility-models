
#----------------------------------------------------------------------------------
# simulations SV Model
#----------------------------------------------------------------------------------
#
# name:         simulationSVM.R
# description:  generate a three simple escenarios with a set of correlated 
#               time series
# functions:  
#                
# output:       
#
# creadet by:   Carolina Gamboa
#
#----------------------------------------------------------------------------------


#----------------------------------------------------------------------------------
# Paths
#----------------------------------------------------------------------------------

#funPath <- "C:\\Users\\dgamboa\\Documents\\programas\\function"
#funPath <-  "C:\\Users\\Carolina\\Documents\\2018\\TESIS\\programas\\function"

#----------------------------------------------------------------------------------
# Load functions
#----------------------------------------------------------------------------------
library(MASS)

setwd(funPath)
source("SVMgen.R")
source("tridiag.R")

#----------------------------------------------------------------------------------
# Main
#----------------------------------------------------------------------------------


simulationSVM <- function(n = n, nG1 = nG1, nG2 = nG2, type = type, a1, a2){
  
  nts <- nG1 + nG2 # number of time series 
  X1 <- NULL
  X2 <- NULL
  X3 <- NULL

  # first we generate a dependent phi noise
  Sigma <- matrix(1, nG1, nG1)
  for(j in 1:nG1){
    for(k in 1:nG1){
      if(j != k) Sigma[k, j] <- 0.9
    }
  }
  phi <- mvrnorm(n = n, rep(0, nG1), Sigma)
  
  
  for(i in 1:nts){
    #-----------------------------------------------------------------
    #first group: 10 full dependent time series with models 1 and 2
    # using phi noise.
    
    # model 1: 
    if(i %in% 1:5){ 
      x <- genSVM(a1, phi[,i], n)
      X1 <- cbind(X1, x)
    }  
    
    # model 2: 
      if(i %in% 6:nG1){ 
      x <- genSVM(a2, phi[,i], n)
      X2 <- cbind(X2, x)
    }
    
  }
    X12 <- cbind(X1, X2)

  #---------------------------------------------------------------------
  #first scenario: group 2 with 5  independent time series with model 1.
  #---------------------------------------------------------------------
  
  if (type == "scenario1"){
      
      # model 1: 
    for(i in 11:15){ 
      phi3 <- rnorm(n)
      x <- genSVM(a1, phi3, n) # a1 parameters model 1
      X3 <- cbind(X3, x)
    }
  }
      
  #--------------------------------------------------------------------------
  # second scenario: group 2 with 5 chain-dependent time series with model 1.
  #--------------------------------------------------------------------------
  
  if(type == "scenario2"){  
    X3 <- NULL
    
    # first we generate a chain dependent phi noise
    Sigma <- tridiag(rep(0.5, nG2-1), 
                     rep(0.5, nG2-1), 
                     rep(1, nG2)) 
    
    eigen(Sigma)
    phi <- mvrnorm(n = n, rep(0, nG2), Sigma)

    # generation of chain dependent time series using phi noise
    
        
    for(i in 1:nG2){
  
      # model 1: 
        x <- genSVM(a1, phi[, i], n) #a1 parameters model 1
        X3 <- cbind(X3, x)
    }
  }
  #--------------------------------------------------------------------------
  # third scenario: group 2 with 5 full-dependent time series with model 1.
  #--------------------------------------------------------------------------
  
  if(type == "scenario3"){    
    X3 <- NULL
    
    # first we generate a full dependent phi noise
    Sigma <- matrix(1, nG2, nG2)
    
    for(j in 1:nG2){
      for(k in 1:nG2){
        if(j != k) Sigma[k, j] <- 0.9
      }
    }
  
    phi <- mvrnorm(n = n, rep(0, nG2), Sigma)
    
    # full dependent series generation using phi
    
    for(i in 1:nG2){
      # model 1: 
      x <- genSVM(a1, phi[, i], n) # a1 parameters model 1
      X3 <- cbind(X3, x)
    }
  }  
  
  W <- cbind(X12, X3) # data scenario "type"
  colnames(W) <- paste(1:nts) # names columns

  return(W)
}




