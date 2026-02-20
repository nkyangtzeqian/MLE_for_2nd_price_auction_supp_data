## change underlying datageneration mechanisms here

# rm(list = ls())   # Remove everything from the Global environment
# Ka=100
# method.in <- "unif"
# para1 <- 1
# para2 <- 20
# Data_Gen_Ind <- 1 ## 1: Poisson process; 2: Gamma process with shape=2, scale=1; 3: Log-normal process
# source("Simu3.R")

library(tidyverse)
library(GoFKernel) # For the inverse function of a CDF

library(msm) # For the tNorm


# source("0327/EbayFunctions_Ver3.R")
source("Main.R")
source("distributions.R")
source("DataGeneration_more.R")
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

Data_Gen_Raw <- List_Data_Gen_Raw[[Data_Gen_Ind]]

#------------------------------------------------------------------------------------
## Choice of different true valuation distributions F, the corresponding choice of parameters,  
## and reserve prices are listed below:
#------------------------------------------------------------------------------------


distanceList <- list(
  KS.init = rep(NA, replication),
  KS.mle = rep(NA, replication),
  KS.pt = rep(NA, replication),
  Ks.gamma = rep(NA, replication),
  KS.tNorm = rep(NA, replication),
  TV.init = rep(NA, replication),
  TV.mle = rep(NA, replication),
  TV.pt = rep(NA, replication),
  TV.gamma = rep(NA, replication),
  TV.tNorm = rep(NA, replication)
)

Profit_allList <- list(
  max_price.init= rep(NA, replication),
  max_price.mle= rep(NA, replication),
  max_price.pt= rep(NA, replication),
  max_price.gamma= rep(NA, replication),
  max_price.tNorm= rep(NA, replication),
  Profit.init = rep(NA, replication),
  Profit.mle = rep(NA, replication),
  Profit.pt = rep(NA, replication),
  Profit.gamma = rep(NA, replication),
  Profit.tNorm = rep(NA, replication),
  est_Profit.init = rep(NA, replication),
  est_Profit.mle = rep(NA, replication),
  est_Profit.pt = rep(NA, replication),
  est_Profit.gamma = rep(NA, replication),
  est_Profit.tNorm = rep(NA, replication)
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
  ERROR1000=try({Bayes_PT <- MCMC_PT(data_Polya, data$Z_iT_i[,1], x.max=max(x.med.vec))}) ## only happend in K=1000 for Unif, PUnif and Pareto
  if(inherits(ERROR1000, "try-error")){
    print("Error in MCMC_PT.")
  }
  
  
  ## Gamma method
  ERROR1000_2=try({Bayes_gamma <- MCMC_gamma(data_Polya, data$Z_iT_i[,1])})
  if(inherits(ERROR1000_2, "try-error")){
    print("Error in MCMC_gamma.")
  }
  
  ## Truncated Normal method
  ERROR1=try({Bayes_tNorm <- MCMC_tNorm(data_Polya, data$Z_iT_i[,1])}, silent=TRUE)
  if(inherits(ERROR1, "try-error")){
    print("Error in MCMC_tNorm.")
  }
  
  ## Kolmogorov-smirnov distance:
  distanceList$KS.init[ii] <- KS.d(Init.Values$F.y, trueF)
  distanceList$KS.mle[ii] <- KS.d(MLE$F.y, trueF)
  # distanceList$KS.pt[ii] <- KS.d(Bayes_PT$F.y, trueF)
  if(!inherits(ERROR1000, "try-error")){
    distanceList$KS.pt[ii] <- KS.d(Bayes_PT$F.y, trueF)
  }
  if(!inherits(ERROR1000_2, "try-error")){
    distanceList$Ks.gamma[ii] <- KS.d(Bayes_gamma$F.y, trueF)
  }
  if(!inherits(ERROR1, "try-error")){
    distanceList$KS.tNorm[ii] <- KS.d(Bayes_tNorm$F.y, trueF)
  }
  
  
  e1.init <- AbscontDistribution(p = linear_interpolation_cdf(x = Init.Values$F.x, F.x = Init.Values$F.y) )
  e1.mle <- AbscontDistribution(p = linear_interpolation_cdf(x = MLE$F.x, F.x = MLE$F.y) )
  # e1.pt <- AbscontDistribution(p = linear_interpolation_cdf(x = Bayes_PT$F.x, F.x = Bayes_PT$F.y) )
   if(!inherits(ERROR1000, "try-error")){
     e1.pt <- AbscontDistribution(p = linear_interpolation_cdf(x = Bayes_PT$F.x, F.x = Bayes_PT$F.y) )
   }
  if(!inherits(ERROR1000_2, "try-error")){
    e1.gamma <- AbscontDistribution(p = linear_interpolation_cdf(x = Bayes_gamma$F.x, F.x = Bayes_gamma$F.y) )
  }
  if(!inherits(ERROR1, "try-error")){
    ERROR2 <- try({e1.tNorm <- AbscontDistribution(p = linear_interpolation_cdf(x = Bayes_tNorm$F.x, F.x = Bayes_tNorm$F.y) )
    }, silent=TRUE)
    if(class(ERROR2)=="try-error"){ ## strange error might occur in AbscontDistribution function due to randomness
      print("Error in e1.tNorm. Try again")
      e1.tNorm <- AbscontDistribution(p = linear_interpolation_cdf(x = Bayes_tNorm$F.x, F.x = Bayes_tNorm$F.y) )
    }
   }
  
  distanceList$TV.init[ii] <- Setting_Gen_Tab1_TV(e1 = e1.init,data,trueF, method = method.in, para1, para2)[[1]]
  distanceList$TV.mle[ii] <- Setting_Gen_Tab1_TV(e1 = e1.mle,data,trueF, method = method.in, para1, para2)[[1]]
  # distanceList$TV.pt[ii] <- Setting_Gen_Tab1_TV(e1 = e1.pt,data,trueF, method = method.in, para1, para2)[[1]]
  if(!inherits(ERROR1000, "try-error")){
    distanceList$TV.pt[ii] <- Setting_Gen_Tab1_TV(e1 = e1.pt,data,trueF, method = method.in, para1, para2)[[1]]
  }
  if(!inherits(ERROR1000_2, "try-error")){
    distanceList$TV.gamma[ii] <- Setting_Gen_Tab1_TV(e1 = e1.gamma,data,trueF, method = method.in, para1, para2)[[1]]
  }
  if(!inherits(ERROR1, "try-error")){
    distanceList$TV.tNorm[ii] <- Setting_Gen_Tab1_TV(e1 = e1.tNorm,data,trueF, method = method.in, para1, para2)[[1]]
  }
  
  ## Profit Metric:
  Profit_init_list <- ProfitMetric(e1 = e1.init, method = method.in, para1, para2)
  Profit_mle_list <- ProfitMetric(e1 = e1.mle, method = method.in, para1, para2)
  # Profit_pt_list <- ProfitMetric(e1 = e1.pt, method = method.in, para1, para2)
  if(!inherits(ERROR1000, "try-error")){
    Profit_pt_list <- ProfitMetric(e1 = e1.pt, method = method.in, para1, para2)
  }
  if(!inherits(ERROR1000_2, "try-error")){
    Profit_gamma_list <- ProfitMetric(e1 = e1.gamma, method = method.in, para1, para2)
  }
  if(!inherits(ERROR1, "try-error")){
    Profit_tNorm_list <- ProfitMetric(e1 = e1.tNorm, method = method.in, para1, para2)
  }
  
  Profit_allList$max_price.init[ii] <- Profit_init_list[[1]]
  Profit_allList$max_price.mle[ii] <- Profit_mle_list[[1]]
  # Profit_allList$max_price.pt[ii] <- Profit_pt_list[[1]]
  if(!inherits(ERROR1000, "try-error")){
     Profit_allList$max_price.pt[ii] <- Profit_pt_list[[1]]
   }
  if(!inherits(ERROR1000_2, "try-error")){
     Profit_allList$max_price.gamma[ii] <- Profit_gamma_list[[1]]
   }
  if(!inherits(ERROR1, "try-error")){
     Profit_allList$max_price.tNorm[ii] <- Profit_tNorm_list[[1]]
   }
  
  
  Profit_allList$Profit.init[ii] <- Profit_init_list[[2]]
  Profit_allList$Profit.mle[ii] <- Profit_mle_list[[2]]
  # Profit_allList$Profit.pt[ii] <- Profit_pt_list[[2]]
  if(!inherits(ERROR1000, "try-error")){
     Profit_allList$Profit.pt[ii] <- Profit_pt_list[[2]]
   }
  if(!inherits(ERROR1000_2, "try-error")){
     Profit_allList$Profit.gamma[ii] <- Profit_gamma_list[[2]]
   }
  if(!inherits(ERROR1, "try-error")){
     Profit_allList$Profit.tNorm[ii] <- Profit_tNorm_list[[2]]
   }
  
  Profit_allList$est_Profit.init[ii] <- Profit_init_list[[3]]
  Profit_allList$est_Profit.mle[ii] <- Profit_mle_list[[3]]
  # Profit_allList$est_Profit.pt[ii] <- Profit_pt_list[[3]]
  if(!inherits(ERROR1000, "try-error")){
     Profit_allList$est_Profit.pt[ii] <- Profit_pt_list[[3]]
   }
  if(!inherits(ERROR1000_2, "try-error")){
     Profit_allList$est_Profit.gamma[ii] <- Profit_gamma_list[[3]]
   }
  if(!inherits(ERROR1, "try-error")){
     Profit_allList$est_Profit.tNorm[ii] <- Profit_tNorm_list[[3]]
   }
  
  print(paste("ii =", ii, "init = ",distanceList$KS.init[ii],"mle = " , distanceList$KS.mle[ii], 
              "pt = ", distanceList$KS.pt[ii], "gamma = ", distanceList$Ks.gamma[ii]))
}

KS.init.mean <- mean(distanceList$KS.init)
KS.mle.mean <- mean(distanceList$KS.mle)
KS.pt.mean <- mean(distanceList$KS.pt)
KS.gamma.mean <- mean(distanceList$Ks.gamma)
TV.init.mean <- mean(distanceList$TV.init)
TV.mle.mean <- mean(distanceList$TV.mle)
TV.pt.mean <- mean(distanceList$TV.pt)
TV.gamma.mean <- mean(distanceList$TV.gamma)

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


save(list="distanceList", file = paste0("SIMU/Table3_","DG_",Data_Gen_Ind,"_", method.in, "_K=", Ka, ".RData"))

Profit.init.mean <- mean(Profit_allList$Profit.init)
Profit.mle.mean <- mean(Profit_allList$Profit.mle)
Profit.pt.mean <- mean(Profit_allList$Profit.pt)
Profit.gamma.mean <- mean(Profit_allList$Profit.gamma)

print(paste("Profit.mle.mean =", Profit.mle.mean))
print(paste("Profit.init.mean =", Profit.init.mean))
print(paste("Profit.pt.mean =", Profit.pt.mean))
print(paste("Profit.gamma.mean =", Profit.gamma.mean))

save(list="Profit_allList", file = paste0("SIMU/backup_S3/Table3_ProfitRes_","DG_",Data_Gen_Ind,"_",  method.in, "_K=", Ka, ".RData"))
