# rm(list = ls())   # Remove everything from the Global environment
# Ka=100
# method.in <- "unif"
# para1 <- 1
# para2 <- 20
# source("Simu1.R")

library(tidyverse)
library(GoFKernel) # For the inverse function of a CDF


# source("0327/EbayFunctions_Ver3.R")
source("Main.R")
source("distributions.R")
source("DataGeneration.R")
source("DisCompare.R")
source("SettingGeneration.R")
source("PT_MCMC.R")

#------------------------------------------------------------------------------------
## User Input: True lambda, Auction window, Number of independent auctions 'N.auction'
#------------------------------------------------------------------------------------

replication <- 100  # replication numbers of each such event i.e. selling an
# identical stuff on ebay with multiple independent auctions.

lambda0 <- 1
tau_a <- 100

#------------------------------------------------------------------------------------
## Choice of different true valuation distributions F, the corresponding choice of parameters,  
## and reserve prices are listed below:
#------------------------------------------------------------------------------------


distanceList <- list(
  KS.init = rep(0, replication),
  KS.mle = rep(0, replication),
  KS.pt = rep(0, replication),
  TV.init = rep(0, replication),
  TV.mle = rep(0, replication),
  TV.pt = rep(0, replication)
)


for(ii in 1:replication){
  
  set.seed(ii)
 
  SETTING_ALL=Setting_Gen_Tab1(Ka, method = method.in, para1, para2)
  x.med.vec <- SETTING_ALL$x.med.vec
  r_k <- SETTING_ALL$r_k
  reserve.price.Cutoff <- SETTING_ALL$reserve.price.Cutoff
  #------------------------------------------------------------------------------------
  ## Data generation, Initial Estimate, MLE, and HUlC
  #------------------------------------------------------------------------------------
  
  Raw_data_all_list <- Data_Gen_Raw(Ka, tau_a, lambda0,
                                    r_k, method = method.in, para1, para2)
  data <- Data_Gen_Observed(Raw_data_all_list)
  
  Init.Values <- Initialization_lambda_F(data, reserve.price.Cutoff)
  theta.init <- (1-Init.Values$F.y)/ c(1,(1-Init.Values$F.y)[-length(Init.Values$F.y)])
  MLE <- MLE_2ndprice(data, lambda.in = Init.Values$lambda, theta.in = theta.init)
  
  ## True F:
  trueF <- Setting_Gen_Tab1_F0(data, method = method.in, para1, para2)[[1]]
  
  
  ## Polya Tree method
  data_Polya <- Data_Gen_Unobserved_PolyaTree(Raw_data_all_list)
  Bayes_PT <- MCMC_PT(data_Polya, data$Z_iT_i[,1], x.max=max(x.med.vec))
  
  ## Kolmogorov-smirnov distance:
  distanceList$KS.init[ii] <- KS.d(Init.Values$F.y, trueF)
  distanceList$KS.mle[ii] <- KS.d(MLE$F.y, trueF)
  distanceList$KS.pt[ii] <- KS.d(Bayes_PT$F.y, trueF)
  
  
  e1.init <- AbscontDistribution(p = linear_interpolation_cdf(x = Init.Values$F.x, F.x = Init.Values$F.y) )
  e1.mle <- AbscontDistribution(p = linear_interpolation_cdf(x = MLE$F.x, F.x = MLE$F.y) )
  e1.pt <- AbscontDistribution(p = linear_interpolation_cdf(x = Bayes_PT$F.x, F.x = Bayes_PT$F.y) )
  
  distanceList$TV.init[ii] <- Setting_Gen_Tab1_TV(e1 = e1.init,data,trueF, method = method.in, para1, para2)[[1]]
  distanceList$TV.mle[ii] <- Setting_Gen_Tab1_TV(e1 = e1.mle,data,trueF, method = method.in, para1, para2)[[1]]
  distanceList$TV.pt[ii] <- Setting_Gen_Tab1_TV(e1 = e1.pt,data,trueF, method = method.in, para1, para2)[[1]]
  
  print(paste("ii =", ii, "init = ",distanceList$KS.init[ii],"mle = " , distanceList$KS.mle[ii], 
              "pt = ", distanceList$KS.pt[ii]))
}

KS.init.mean <- mean(distanceList$KS.init)
KS.mle.mean <- mean(distanceList$KS.mle)
KS.pt.mean <- mean(distanceList$KS.pt)
TV.init.mean <- mean(distanceList$TV.init)
TV.mle.mean <- mean(distanceList$TV.mle)
TV.pt.mean <- mean(distanceList$TV.pt)

print(paste("method.in =", method.in))
print(paste("K =",Ka))
print(paste("KS.mle.mean =", KS.mle.mean))
print(paste("KS.init.mean =", KS.init.mean))
print(paste("KS.pt.mean =", KS.pt.mean))
print(paste("TV.mle.mean =", TV.mle.mean))
print(paste("TV.init.mean =", TV.init.mean))
print(paste("TV.pt.mean =", TV.pt.mean))


save(list="distanceList", file = paste0("SIMU/Table1_", method.in, "_K=", Ka, ".RData"))
#######################################################################################
## Settings
#######################################################################################
## Choose the appropriate true F and comment others:

# method.in <- "unif"
# para1 <- 1
# para2 <- 20
# x.med.vec <- seq(0.1, 20, by= 0.995)
# r_k <- runif(Ka, 0.1, 3)
# reserve.price.Cutoff <- 1
# 

# method.in <- "piecewise_unif"
# para1 <- 2
# para2 <- 4
# x.med.vec <- seq(0, 4, by = 4/20)
# r_k <- runif(N.auction, 0.1, 1.5)
# reserve.price.Cutoff <- 1


# method.in <- "pareto"
# para1 <- 3
# para2 <- 100
# x.med.vec <- seq(0, 20, by= 1)
# r_k <- runif(N.auction, 0.001, 0.1)
# reserve.price.Cutoff <- 0.05


# method.in <- "gamma"
# para1 <- 10
# para2 <- 2
# x.med.vec <- seq(0, 12, by= 12/20)
# r_k <- runif(N.auction, 0.1, 3)
# reserve.price.Cutoff <- 1


# method.in <- "beta"
# para1 <- 2
# para2 <- 2
# x.med.vec <- seq(0, 1, by= 1/20)
# r_k <- runif(N.auction, 0.001, 0.2)
# reserve.price.Cutoff <- 0.05

