# rm(list = ls())   # Remove everything from the Global environment
# Ka=100
# method.in <- "unif"
# para1 <- 1
# para2 <- 20
# source("Simu1.R")

library(tidyverse)
library(GoFKernel) # For the inverse function of a CDF

library(msm) # for the truncated normal functions


# source("0327/EbayFunctions_Ver3.R")
source("Main.R")
source("distributions.R")
source("DataGeneration.R")
source("DisCompare.R")
source("SettingGeneration.R")
source("PT_MCMC.R")

source("Para_MCMC.R")
source("ProfitMetric.R")

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
load(file = paste0("SIMU/backup_S1/Table1_ProfitRes_", method.in, "_K=", Ka, ".RData"))
load(file = paste0("SIMU/backup_S1/Table1_", method.in, "_K=", Ka, ".RData"))


distanceList$KS.gamma = rep(NA, replication)
distanceList$KS.tNorm = rep(NA, replication)
distanceList$TV.gamma = rep(NA, replication)
distanceList$TV.tNorm = rep(NA, replication)


Profit_allList$max_price.gamma = rep(NA, replication)
Profit_allList$Profit.gamma = rep(NA, replication)
Profit_allList$est_Profit.gamma = rep(NA, replication)
Profit_allList$max_price.tNorm = rep(NA, replication)
Profit_allList$Profit.tNorm = rep(NA, replication)
Profit_allList$est_Profit.tNorm = rep(NA, replication)

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
  
  # Init.Values <- Initialization_lambda_F(data, reserve.price.Cutoff)
  # theta.init <- (1-Init.Values$F.y)/ c(1,(1-Init.Values$F.y)[-length(Init.Values$F.y)])
  # MLE <- MLE_2ndprice(data, lambda.in = Init.Values$lambda, theta.in = theta.init)
  
  ## True F:
  trueF <- Setting_Gen_Tab1_F0(data, method = method.in, para1, para2)[[1]]
  
  
  ## Polya Tree method
  data_Polya <- Data_Gen_Unobserved_PolyaTree(Raw_data_all_list)
  # Bayes_PT <- MCMC_PT(data_Polya, data$Z_iT_i[,1], x.max=max(x.med.vec))
  
  ## Gamma method
  Bayes_gamma <- MCMC_gamma(data_Polya, data$Z_iT_i[,1])
  ## Truncated Normal method
  ERROR1=try({Bayes_tNorm <- MCMC_tNorm(data_Polya, data$Z_iT_i[,1])}, silent=TRUE)
  if(class(ERROR1)=="try-error"){ ## strange error might occur in MCMC, rbinom
    print("Error in MCMC_tNorm")
    
    distanceList$Ks.gamma[ii] <- KS.d(Bayes_gamma$F.y, trueF)

    e1.gamma <- AbscontDistribution(p = linear_interpolation_cdf(x = Bayes_gamma$F.x, F.x = Bayes_gamma$F.y) )

    distanceList$TV.gamma[ii] <- Setting_Gen_Tab1_TV(e1 = e1.gamma,data,trueF, method = method.in, para1, para2)[[1]]

    Profit_gamma_list <- ProfitMetric(e1 = e1.gamma, method = method.in, para1, para2)

    Profit_allList$max_price.gamma[ii] <- Profit_gamma_list[[1]]

    Profit_allList$Profit.gamma[ii] <- Profit_gamma_list[[2]]

    Profit_allList$est_Profit.gamma[ii] <- Profit_gamma_list[[3]]

    
    print(paste("ii =", ii, "init = ",distanceList$KS.init[ii],"mle = " , distanceList$KS.mle[ii], 
                "pt = ", distanceList$KS.pt[ii], "gamma = ", distanceList$Ks.gamma[ii]))
    next
  }
  
  ## Kolmogorov-smirnov distance:
  # distanceList$KS.init[ii] <- KS.d(Init.Values$F.y, trueF)
  # distanceList$KS.mle[ii] <- KS.d(MLE$F.y, trueF)
  # distanceList$KS.pt[ii] <- KS.d(Bayes_PT$F.y, trueF)
  distanceList$Ks.gamma[ii] <- KS.d(Bayes_gamma$F.y, trueF)
  distanceList$KS.tNorm[ii] <- KS.d(Bayes_tNorm$F.y, trueF)
  
  
  # e1.init <- AbscontDistribution(p = linear_interpolation_cdf(x = Init.Values$F.x, F.x = Init.Values$F.y) )
  # e1.mle <- AbscontDistribution(p = linear_interpolation_cdf(x = MLE$F.x, F.x = MLE$F.y) )
  # e1.pt <- AbscontDistribution(p = linear_interpolation_cdf(x = Bayes_PT$F.x, F.x = Bayes_PT$F.y) )
  e1.gamma <- AbscontDistribution(p = linear_interpolation_cdf(x = Bayes_gamma$F.x, F.x = Bayes_gamma$F.y) )
  ERROR2 <- try({e1.tNorm <- AbscontDistribution(p = linear_interpolation_cdf(x = Bayes_tNorm$F.x, F.x = Bayes_tNorm$F.y) )
  }, silent=TRUE)
  if(class(ERROR2)=="try-error"){ ## strange error might occur in AbscontDistribution function due to randomness
    print("Error in e1.tNorm. Try again")
    e1.tNorm <- AbscontDistribution(p = linear_interpolation_cdf(x = Bayes_tNorm$F.x, F.x = Bayes_tNorm$F.y) )
  }
  
  
  # distanceList$TV.init[ii] <- Setting_Gen_Tab1_TV(e1 = e1.init,data,trueF, method = method.in, para1, para2)[[1]]
  # distanceList$TV.mle[ii] <- Setting_Gen_Tab1_TV(e1 = e1.mle,data,trueF, method = method.in, para1, para2)[[1]]
  # distanceList$TV.pt[ii] <- Setting_Gen_Tab1_TV(e1 = e1.pt,data,trueF, method = method.in, para1, para2)[[1]]
  distanceList$TV.gamma[ii] <- Setting_Gen_Tab1_TV(e1 = e1.gamma,data,trueF, method = method.in, para1, para2)[[1]]
  distanceList$TV.tNorm[ii] <- Setting_Gen_Tab1_TV(e1 = e1.tNorm,data,trueF, method = method.in, para1, para2)[[1]]
  
  ## Profit Metric:
  # Profit_init_list <- ProfitMetric(e1 = e1.init, method = method.in, para1, para2)
  # Profit_mle_list <- ProfitMetric(e1 = e1.mle, method = method.in, para1, para2)
  # Profit_pt_list <- ProfitMetric(e1 = e1.pt, method = method.in, para1, para2)
  Profit_gamma_list <- ProfitMetric(e1 = e1.gamma, method = method.in, para1, para2)
  Profit_tNorm_list <- ProfitMetric(e1 = e1.tNorm, method = method.in, para1, para2)
  
  # Profit_allList$max_price.init[ii] <- Profit_init_list[[1]]
  # Profit_allList$max_price.mle[ii] <- Profit_mle_list[[1]]
  # Profit_allList$max_price.pt[ii] <- Profit_pt_list[[1]]
  Profit_allList$max_price.gamma[ii] <- Profit_gamma_list[[1]]
  Profit_allList$max_price.tNorm[ii] <- Profit_tNorm_list[[1]]
  
  # Profit_allList$Profit.init[ii] <- Profit_init_list[[2]]
  # Profit_allList$Profit.mle[ii] <- Profit_mle_list[[2]]
  # Profit_allList$Profit.pt[ii] <- Profit_pt_list[[2]]
  Profit_allList$Profit.gamma[ii] <- Profit_gamma_list[[2]]
  Profit_allList$Profit.tNorm[ii] <- Profit_tNorm_list[[2]]
  
  # Profit_allList$est_Profit.init[ii] <- Profit_init_list[[3]]
  # Profit_allList$est_Profit.mle[ii] <- Profit_mle_list[[3]]
  # Profit_allList$est_Profit.pt[ii] <- Profit_pt_list[[3]]
  Profit_allList$est_Profit.gamma[ii] <- Profit_gamma_list[[3]]
  Profit_allList$est_Profit.tNorm[ii] <- Profit_tNorm_list[[3]]
  
  print(paste("ii =", ii, "init = ",distanceList$KS.init[ii],"mle = " , distanceList$KS.mle[ii], 
              "pt = ", distanceList$KS.pt[ii], "gamma = ", distanceList$Ks.gamma[ii],
              "tNorm =", distanceList$KS.tNorm[ii]))
}

KS.init.mean <- mean(distanceList$KS.init)
KS.mle.mean <- mean(distanceList$KS.mle)
KS.pt.mean <- mean(distanceList$KS.pt)
KS.gamma.mean <- mean(distanceList$Ks.gamma)
TV.init.mean <- mean(distanceList$TV.init)
TV.mle.mean <- mean(distanceList$TV.mle)
TV.pt.mean <- mean(distanceList$TV.pt)
TV.gamma.mean <- mean(distanceList$TV.gamma)

KS.tNorm.mean <- mean(distanceList$KS.tNorm, na.rm = TRUE)
TV.tNorm.mean <- mean(distanceList$TV.tNorm, na.rm = TRUE)

print(paste("method.in =", method.in))
print(paste("K =",Ka))
print(paste("KS.mle.mean =", KS.mle.mean))
print(paste("KS.init.mean =", KS.init.mean))
print(paste("KS.pt.mean =", KS.pt.mean))
print(paste("KS.gamma.mean =", KS.gamma.mean))
print(paste("TV.mle.mean =", TV.mle.mean))
print(paste("TV.init.mean =", TV.init.mean))
print(paste("TV.pt.mean =", TV.pt.mean))
print(paste("TV.gamma.mean =", TV.gamma.mean))

print(paste("KS.tNorm.mean =", KS.tNorm.mean))
print(paste("TV.tNorm.mean =", TV.tNorm.mean))


save(list="distanceList", file = paste0("SIMU/Table1A_", method.in, "_K=", Ka, ".RData"))

Profit.init.mean <- mean(Profit_allList$Profit.init)
Profit.mle.mean <- mean(Profit_allList$Profit.mle)
Profit.pt.mean <- mean(Profit_allList$Profit.pt)
Profit.gamma.mean <- mean(Profit_allList$Profit.gamma)
Profit.tNorm.mean <- mean(Profit_allList$Profit.tNorm, na.rm = TRUE)

print(paste("Profit.mle.mean =", Profit.mle.mean))
print(paste("Profit.init.mean =", Profit.init.mean))
print(paste("Profit.pt.mean =", Profit.pt.mean))
print(paste("Profit.gamma.mean =", Profit.gamma.mean))
print(paste("Profit.tNorm.mean =", Profit.tNorm.mean))

save(list="Profit_allList", file = paste0("SIMU/Table1A_ProfitRes_", method.in, "_K=", Ka, ".RData"))
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

