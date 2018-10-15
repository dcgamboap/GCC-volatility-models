#----------------------------------------------------------------------------------
# clustering SV Model
#----------------------------------------------------------------------------------
#
# name:         GCCmatrix.R
# description:  built GCC similarity matrix based on dependences 
# functions:  
#                
# output:	- similarity matrix
        # - k
# creadet by:   Carolina Gamboa
#
#----------------------------------------------------------------------------------

#----------------------------------------------------------------------------------
# functions
#----------------------------------------------------------------------------------
# library

library('forecast')
library('tseries')

source("estimaARP.R")
source("kOptim.R")
source("GCC_d.R")


#............................................................................
# Main  
#............................................................................


GCCmatrix <- function(serie_0){
  # first we determine the upper limit for "k" using AR(P) model
  
  setAR <- function(x){
    aar <- ar(x, aic = TRUE)
    return(aar$order)
  }
  
  PP <- apply(serie_0, 2, setAR) # search max order in data set
  kSup <- max(PP) # upper limit for "k"
  
  # linear regression
  # we search optim k using  linear regressions between 0 and kSup
  kOp <- data.frame()
  kOp1 <- 1
  kinf <- 1
  for(jj in 1:(ncol(serie_0)-1)){
    for (ii in (jj+1):ncol(serie_0)){
      if(kOp1 < kSup){
        # kOp1 <- c(kOptim(serie_0[, jj], serie_0[, ii], kSup))
        kOp1 <- c(kOptim(serie_0[, jj], serie_0[, ii], kSup, kinf))
        
        kOp <- rbind(kOp, data.frame(i = jj, j = ii, kOp1))
        kinf <- kOp1
        # cat(jj, ii, ", ")
      }
    }
  }
  
  
  kDef <- max(kOp$kOp1)
  
  #------------------------------------------------------------------------------
  # similarity matrix
  
  nSerie <- ncol(serie_0)
  DM <- diag(x = 0, nrow = nSerie, ncol = nSerie)
  
  # construcción de la matriz de disimilaridad. 
  
  for(ii in 1:nSerie){
    for(jj in ii:nSerie){
      g <- GCC_d(serie_0[, ii], serie_0[, jj], kDef)
      DM[ii, jj] <- 1 - g
      DM[jj, ii] <- 1 - g
    }
    #cat(ii, ",", " ")
  }
  
  colnames(DM) <- colnames(serie_0)
  
  sale <- list(DM = DM,   k_GCC = kDef)
  return(sale)
}

