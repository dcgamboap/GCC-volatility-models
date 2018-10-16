
#----------------------------------------------------------------------------------
# simulations ARCH Model
#----------------------------------------------------------------------------------
#
# name:         main.R
# description:  preliminar simulation study for ARCH 
# functions:  all methods use functions of function directory
#                
# output:    DATE_ARCH.RData   
#
# creadet by:   Carolina Gamboa
#
#----------------------------------------------------------------------------------
rm(list = ls())
#----------------------------------------------------------------------------------
# Paths
#----------------------------------------------------------------------------------

#funPath <- "C:\\Users\\dgamboa\\Documents\\programas\\function"
funPath <-  "C:\\Users\\Carolina\\Documents\\2018\\TESIS\\programas\\function"
datPath <-  "C:\\Users\\Carolina\\Documents\\2018\\TESIS\\data\\workingSpaceSave"
#----------------------------------------------------------------------------------
# function
#----------------------------------------------------------------------------------
setwd(funPath)
source("simulationARCH.R")
source("clust.R")
source("GCCmatrix.R")
library(fGarch)


#----------------------------------------------------------------------------------
# main
#----------------------------------------------------------------------------------

n <- 1100 # real size 900
nG1 <- 10
nG2 <- 5

# models parameters
a1 <- c(0.01, 0.3) # model 1
# a2 <- c(0.622, 0.0022, 0.074, 0.093) # model 2
a2 <- c(0.001, 0.6) # model 2

type1 <- "scenario1"
type2 <- "scenario2"
type3 <- "scenario3"

# original labels
origClust1 <- c(rep(1, 10), 2:6)
origClust2 <- c(rep(1, 10), rep(2, 5))
origClust3 <- origClust2


#-------------------------------------------------------------------------------
# simulation first scenario
#-------------------------------------------------------------------------------

ARI1 <- NULL
k1 <- NULL
ptm <- proc.time() 

for(i in 1:100){
  # i = 1
  X1 <- simulationARCH(n, nG1, nG2, type1, a1, a2)
  serie  <- X1
  serie2 <- X1^2 

  sigma_est <- apply(serie, 2, function(x) sqrt(garchFit(~garch(1, 0), data = x)@h.t) )
  cat(i, "\n")
  
  DM1  <- GCCmatrix(serie)
  DM2  <- GCCmatrix(serie2)
  DM3  <- GCCmatrix(sigma_est)
  

  r1  <- clust(DM1$DM, origClust1 = origClust1, F, "single")
  r2  <- clust(DM2$DM, origClust1 = origClust1, F,  "single")
  r3  <- clust(DM3$DM, origClust1 = origClust1, F,  "single")
  
  # r   <- c(r1$RAdjust, r2$RAdjust)
  # ARI1 <- rbind(ARI1, r)

  # # complete
  # 
  r1Com  <- clust(DM1$DM, origClust1 = origClust1, F, "complete")
  r2Com  <- clust(DM2$DM, origClust1 = origClust1, F, "complete")
  r3Com  <- clust(DM3$DM, origClust1 = origClust1, F, "complete")
  
  # # average
  # 
  r1Ave  <- clust(DM1$DM, origClust1 = origClust1, F, "average")
  r2Ave  <- clust(DM2$DM, origClust1 = origClust1, F, "average")
  r3Ave  <- clust(DM3$DM, origClust1 = origClust1, F, "average")
  
  #  r   <- c(r1Sin$RAdjust, r2Sin$RAdjust)
  r   <- c(r1$RAdjust, r2$RAdjust, r3$RAdjust, 
           r1Com$RAdjust, r2Com$RAdjust, r3Com$RAdjust, 
           r1Ave$RAdjust, r2Ave$RAdjust, r3Ave$RAdjust)
  
  k   <- c(DM1$k_GCC, DM2$k_GCC, DM3$k_GCC)
  
  k1 <- rbind(k1, k)
  ARI1 <- rbind(ARI1, r)
  }
ltm <- proc.time() - ptm

#-------------------------------------------------------------------------------
# simulation second scenario
#-------------------------------------------------------------------------------

ptm2 <- proc.time() 

ARI2 <- NULL
k2 <- NULL

for(i in 1:100){
  # i = 1
  X1 <- simulationARCH(n, nG1, nG2, type2, a1, a2) #simulacion de las series de tiempo
  
  serie  <- X1
  serie2 <- X1^2 
  # volatility estimation
  sigma_est <- apply(serie, 2, function(x) sqrt(garchFit(~garch(1, 0), data = x)@h.t) )
  
  cat(i, "\n")
  
  # Building simmilarity matrices using GCC distance
  DM1  <- GCCmatrix(serie) # original time series
  DM2  <- GCCmatrix(serie2) # squares time series
  DM3  <- GCCmatrix(sigma_est) # time series volatility
  
  
  # Clustering using single linkage
  r1  <- clust(DM1$DM, origClust1 = origClust2, F, "single")
  r2  <- clust(DM2$DM, origClust1 = origClust2, F, "single")
  r3  <- clust(DM3$DM, origClust1 = origClust2, F, "single")
  
  # # complete
  # 
  r1Com  <- clust(DM1$DM, origClust1 = origClust2, F, "complete")
  r2Com  <- clust(DM2$DM, origClust1 = origClust2, F, "complete")
  r3Com  <- clust(DM3$DM, origClust1 = origClust2, F, "complete")
  
  # # average
  # 
  r1Ave  <- clust(DM1$DM, origClust1 = origClust2, F, "average")
  r2Ave  <- clust(DM2$DM, origClust1 = origClust2, F, "average")
  r3Ave  <- clust(DM3$DM, origClust1 = origClust2, F, "average")
  

  r   <- c(r1$RAdjust, r2$RAdjust, r3$RAdjust, 
           r1Com$RAdjust, r2Com$RAdjust, r3Com$RAdjust, 
           r1Ave$RAdjust, r2Ave$RAdjust, r3Ave$RAdjust)
  
  k   <- c(DM1$k_GCC, DM2$k_GCC, DM3$k_GCC)

  ARI2 <- rbind(ARI2, r)
  k2 <- rbind(k2, k)
  
  
}

ltm2 <- proc.time() - ptm2

#-------------------------------------------------------------------------------
# simulation third scenario
#-------------------------------------------------------------------------------

ptm3 <- proc.time() 

ARI3 <- NULL
k3 <- NULL

for(i in 1:100){
  # i = 1
  X1 <- simulationARCH(n, nG1, nG2, type3, a1, a2)

  serie  <- X1
  serie2 <- X1^2 
  # volatility estimation
  sigma_est <- apply(serie, 2, function(x) sqrt(garchFit(~garch(1, 0), data = x)@h.t) )
  
  cat(i, "\n")
  
  # Building simmilarity matrices using GCC distance
  
  DM1  <- GCCmatrix(serie) #original time series
  DM2  <- GCCmatrix(serie2) # squares time series
  DM3  <- GCCmatrix(sigma_est) # volatility time series
  
  
  # # Clustering 
  # single
  r1  <- clust(DM1$DM, origClust1 = origClust3, F, "single")
  r2  <- clust(DM2$DM, origClust1 = origClust3, F, "single")
  r3  <- clust(DM3$DM, origClust1 = origClust3, F, "single")
  
  # # complete
  # 
  r1Com  <- clust(DM1$DM, origClust1 = origClust3, F, "complete")
  r2Com  <- clust(DM2$DM, origClust1 = origClust3, F, "complete")
  r3Com  <- clust(DM3$DM, origClust1 = origClust3, F, "complete")
  
  # # average
  # 
  r1Ave  <- clust(DM1$DM, origClust1 = origClust3, F, "average")
  r2Ave  <- clust(DM2$DM, origClust1 = origClust3, F, "average")
  r3Ave  <- clust(DM2$DM, origClust1 = origClust3, F, "average")
  
  r   <- c(r1$RAdjust, r2$RAdjust, r3$RAdjust, 
           r1Com$RAdjust, r2Com$RAdjust, r3Com$RAdjust, 
           r1Ave$RAdjust, r2Ave$RAdjust, r3Ave$RAdjust)
  
  k   <- c(DM1$k_GCC, DM2$k_GCC, DM3$k_GCC)
  
  k3 <- rbind(k3, k)
  ARI3 <- rbind(ARI3, r)
  

  }
ltm3 <- proc.time() - ptm3



#----------------------------------------------------------------------------------
# Consolidation 
#----------------------------------------------------------------------------------


resume <- cbind(colMeans(ARI1), colMeans(ARI2), colMeans(ARI3))
colnames(resume) <- c("Scenario 1", "Scenario 2", "Scenario 3")
rownames(resume) <- c("GGC single", "GGC2 single", "Volatility single",  
                      "GGC complete", "GGC2 complete", "Volatility complete", 
                      "GGC average", "GGC2 average", "Volatility average" )

library(xtable)
xtable(resume)
#-------------------------------------------------------------------------------
# outFiles
#-------------------------------------------------------------------------------


setwd(datPath)

# Saving R objects 
save(ARI1, ARI2, ARI3, k1, k2, k3,  file = "20181016_ARCH.RData")





