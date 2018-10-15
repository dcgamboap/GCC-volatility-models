#----------------------------------------------------------------------------------
# clustering SV Model
#----------------------------------------------------------------------------------
#
# name:         clust.R
# description:  clustering time series based on dependences 
# input: dissimilary matrix  
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


clust <- function(DM, origClust1, if.plot, clus_method){

  # built dendogram
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
  sale <- list(groups = cl.Agnes.Euc.Single, RAdjust = RAdjust)
  return(sale)
}

