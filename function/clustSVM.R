#----------------------------------------------------------------------------------
# clustering SV Model
#----------------------------------------------------------------------------------
#
# name:         clustSVM.R
# description:  clustering time series based on dependences 
# functions:  
#                
# output:	- if if.plot is true a dendogram of mthod      
#		- RandAdjustIndex for the clustering and groups  	
# creadet by:   Carolina Gamboa
# v1: modify the funtion for select link method type. 
# v2: modify plot dendogram using jclust package
# v3: adding kinf for search the optim k to GGC measure
#
#----------------------------------------------------------------------------------

#----------------------------------------------------------------------------------
# functions
#----------------------------------------------------------------------------------
# library

library('ggplot2')
library('forecast')
library('tseries')
library(cluster)
library(ggplot2)
library(plyr)
library(zoo)
library(mclust)

source("estimaARP.R")
source("kOptim.R")
source("GCC_d.R")


#............................................................................
# Main  
#............................................................................


clustSVM <- function(serie_0, origClust1, if.plot, clus_method){
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
  
  # single linkage
  nClust <- length(unique(origClust1))
  Agnes.Euc.Single <- agnes(DM, diss = TRUE, method=clus_method)
  cl.Agnes.Euc.Single <- cutree(Agnes.Euc.Single, nClust)
  
  if(if.plot == TRUE){ 
    dagn  <- as.dendrogram(as.hclust(Agnes.Euc.Single))
#    op <- par(mar = par("mar") + c(0,0,0, 2)) # more space to the right
    plot(dagn, horiz = FALSE, center = TRUE,
         nodePar = list(lab.cex = 1, lab.col = "black", pch = NA))  
  }
  RAdjust <- adjustedRandIndex(origClust1, cl.Agnes.Euc.Single)
  sale <- list(groups = cl.Agnes.Euc.Single, RAdjust = RAdjust, k_GCC = kDef)
  return(sale)
}

