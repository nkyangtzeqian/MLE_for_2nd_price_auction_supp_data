# rm(list = ls())   # Remove everything from the Global environment
# Ka=100
# method.in <- "unif"
# para1 <- 1
# para2 <- 20
# LAMBDA0_SETTING_INDEX <- 1
# source("Simu1.R")

library(tidyverse)
library(GoFKernel) # For the inverse function of a CDF

library(msm) # for the truncated normal functions

# source("0327/EbayFunctions_Ver3.R")
source("MainL.R")
source("distributionsL.R")
source("DataGenerationL.R")
source("DisCompare.R")
source("SettingGeneration.R")
source("lambda_fun.R")
source("PT_MCMC.R")

source("Para_MCMC.R")
source("ProfitMetric.R")


#------------------------------------------------------------------------------------
## User Input: True lambda, Auction window, Number of independent auctions 'N.auction'
#------------------------------------------------------------------------------------

LAMBDA_SET_LIST=list(
  constant=list(TWO_STAGE=FALSE, h_index=1),
  two_stage=list(TWO_STAGE=TRUE, tau_1=75),
  proportional=list(TWO_STAGE=FALSE, h_index=2)
)

replication <- 100  # replication numbers of each such event i.e. selling an
# identical stuff on ebay with multiple independent auctions.

lambda0_1 <- 1
lambda0_2 <- 5 ## might not be used; 2 is too small to make any difference
lambda0_list <- lambda0_list_init(lambda0_1,lambda0_2,LAMBDA_SET_LIST[[LAMBDA0_SETTING_INDEX]])
rm(lambda0_1,lambda0_2) ## avoid using them
tau_a <- 100

#------------------------------------------------------------------------------------
## Choice of different true valuation distributions F, the corresponding choice of parameters,  
## and reserve prices are listed below:
#------------------------------------------------------------------------------------


distanceList <- list(
  KS.init1 = rep(NA, replication),
  KS.init2 = rep(NA, replication),
  KS.init3 = rep(NA, replication),
  KS.mle1 = rep(NA, replication),
  KS.mle2 = rep(NA, replication),
  KS.mle3 = rep(NA, replication),
  KS.pt = rep(NA, replication),
  KS.gamma = rep(NA, replication),
  KS.tNorm = rep(NA, replication),
  TV.init1 = rep(NA, replication),
  TV.init2 = rep(NA, replication),
  TV.init3 = rep(NA, replication),
  TV.mle1 = rep(NA, replication),
  TV.mle2 = rep(NA, replication),
  TV.mle3 = rep(NA, replication),
  TV.pt = rep(NA, replication),
  TV.gamma = rep(NA, replication),
  TV.tNorm = rep(NA, replication)
)

Profit_allList <- list(
  max_price.init1= rep(NA, replication),
  max_price.init2= rep(NA, replication),
  max_price.init3= rep(NA, replication),
  max_price.mle1= rep(NA, replication),
  max_price.mle2= rep(NA, replication),
  max_price.mle3= rep(NA, replication),
  max_price.pt= rep(NA, replication),
  max_price.gamma= rep(NA, replication),
  max_price.tNorm= rep(NA, replication),
  Profit.init1 = rep(NA, replication),
  Profit.init2 = rep(NA, replication),
  Profit.init3 = rep(NA, replication),
  Profit.mle1 = rep(NA, replication),
  Profit.mle2 = rep(NA, replication),
  Profit.mle3 = rep(NA, replication),
  Profit.pt = rep(NA, replication),
  Profit.gamma = rep(NA, replication),
  Profit.tNorm = rep(NA, replication),
  est_Profit.init1 = rep(NA, replication),
  est_Profit.init2 = rep(NA, replication),
  est_Profit.init3 = rep(NA, replication),
  est_Profit.mle1 = rep(NA, replication),
  est_Profit.mle2 = rep(NA, replication),
  est_Profit.mle3 = rep(NA, replication),
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
  ## Data generation, Initial Estimate, MLEs and so on
  #------------------------------------------------------------------------------------
  
  Raw_data_all_list <- Data_Gen_Raw(Ka, tau_a, lambda0_list,
                                    r_k, method = method.in, para1, para2)
  data <- Data_Gen_Observed(Raw_data_all_list)
  
  Init.Values_1 <- Initialization_lambda_F(data, reserve.price.Cutoff, lambda_set = LAMBDA_SET_LIST[[1]])
  theta.init <- (1-Init.Values_1$F.y)/ c(1,(1-Init.Values_1$F.y)[-length(Init.Values_1$F.y)])
  MLE_1 <- MLE_2ndprice(data, lambda.in = Init.Values_1$lambda, theta.in = theta.init, lambda_set = LAMBDA_SET_LIST[[1]])
  
  Init.Values_2 <- Initialization_lambda_F(data, reserve.price.Cutoff, lambda_set = LAMBDA_SET_LIST[[2]])
  theta.init <- (1-Init.Values_2$F.y)/ c(1,(1-Init.Values_2$F.y)[-length(Init.Values_2$F.y)])
  MLE_2 <- MLE_2ndprice(data, lambda.in = Init.Values_2$lambda, theta.in = theta.init, lambda_set = LAMBDA_SET_LIST[[2]])
  
  Init.Values_3 <- Initialization_lambda_F(data, reserve.price.Cutoff, lambda_set = LAMBDA_SET_LIST[[3]])
  theta.init <- (1-Init.Values_3$F.y)/ c(1,(1-Init.Values_3$F.y)[-length(Init.Values_3$F.y)])
  MLE_3 <- MLE_2ndprice(data, lambda.in = Init.Values_3$lambda, theta.in = theta.init, lambda_set = LAMBDA_SET_LIST[[3]])

  
  ## True F:
  trueF <- Setting_Gen_Tab1_F0(data, method = method.in, para1, para2)[[1]]
  
  
  ## Polya Tree method
  data_Polya <- Data_Gen_Unobserved_PolyaTree(Raw_data_all_list)
  Bayes_PT <- MCMC_PT(data_Polya, data$Z_iT_i[,1], x.max=max(x.med.vec))
  
  ## Gamma method
  Bayes_gamma <- MCMC_gamma(data_Polya, data$Z_iT_i[,1])
  ## Truncated Normal method
  ERROR1=try({Bayes_tNorm <- MCMC_tNorm(data_Polya, data$Z_iT_i[,1])}, silent=TRUE)
  ERROR_TN=F
  if(class(ERROR1)=="try-error"){ ## strange error might occur in MCMC, rbinom
    print("Error in MCMC_tNorm")
    ERROR_TN=T
  }
  
  #------------------------------------------------------------------------------------
  ## record results
  #------------------------------------------------------------------------------------
  
  RES_List <- list(Init.Values_1, Init.Values_2, Init.Values_3, MLE_1, MLE_2, MLE_3, Bayes_PT, Bayes_gamma)
  RES_Name_List <- c("init1", "init2", "init3", "mle1", "mle2", "mle3", "pt", "gamma")
  if(!ERROR_TN){
    RES_List <- c(RES_List, list(Bayes_tNorm))
    RES_Name_List <- c(RES_Name_List, "tNorm")
  } 
  
  for (res in 1:length(RES_List)){
    res_obj <- RES_List[[res]]
    res_name <- RES_Name_List[res]
    
    if(res_name=="tNorm"){ ## might need to try twice for TN
      ERROR2 <- try({e1.obj <- AbscontDistribution(p = linear_interpolation_cdf(x = res_obj$F.x, F.x = res_obj$F.y) )
      }, silent=TRUE)
      if(class(ERROR2)=="try-error"){ ## strange error might occur in AbscontDistribution function occasionally, try again
        print("Error in e1.tNorm. Try again")
        ERROR2 <- try({e1.obj <- AbscontDistribution(p = linear_interpolation_cdf(x = res_obj$F.x, F.x = res_obj$F.y) )
        }, silent=TRUE)
        if(class(ERROR2)=="try-error"){ ## strange error might occur in AbscontDistribution function
          print("Error in e1.tNorm again.") 
          next
        }
      }
    }else{
      e1.obj <- AbscontDistribution(p = linear_interpolation_cdf(x = res_obj$F.x, F.x = res_obj$F.y))
    } 
    
    distanceList[[paste0("KS.", res_name)]][ii] <- KS.d(res_obj$F.y, trueF)
    distanceList[[paste0("TV.", res_name)]][ii] <- Setting_Gen_Tab1_TV(e1 = e1.obj, data, trueF, method = method.in, para1, para2)[[1]]
    
    Profit_list <- ProfitMetric(e1 = e1.obj, method = method.in, para1, para2)
    Profit_allList[[paste0("max_price.", res_name)]][ii] <- Profit_list[[1]]
    Profit_allList[[paste0("Profit.", res_name)]][ii] <- Profit_list[[2]]
    Profit_allList[[paste0("est_Profit.", res_name)]][ii] <- Profit_list[[3]]
  }
  
  print(paste("ii =", ii, "init1 = ",distanceList$KS.init1[ii],"mle1 = " , distanceList$KS.mle1[ii], 
              "init2 = ", distanceList$KS.init2[ii], "mle2 = ", distanceList$KS.mle2[ii],
              "init3 = ", distanceList$KS.init3[ii], "mle3 = ", distanceList$KS.mle3[ii],
              "pt = ", distanceList$KS.pt[ii], "gamma = ", distanceList$KS.gamma[ii],
              "tNorm =", distanceList$KS.tNorm[ii]))
}

KS.init1.mean <- mean(distanceList$KS.init1)
KS.mle1.mean <- mean(distanceList$KS.mle1)
KS.init2.mean <- mean(distanceList$KS.init2)
KS.mle2.mean <- mean(distanceList$KS.mle2)
KS.init3.mean <- mean(distanceList$KS.init3)
KS.mle3.mean <- mean(distanceList$KS.mle3)
KS.pt.mean <- mean(distanceList$KS.pt)
KS.gamma.mean <- mean(distanceList$KS.gamma)
TV.init1.mean <- mean(distanceList$TV.init1)
TV.mle1.mean <- mean(distanceList$TV.mle1)
TV.init2.mean <- mean(distanceList$TV.init2)
TV.mle2.mean <- mean(distanceList$TV.mle2)
TV.init3.mean <- mean(distanceList$TV.init3)
TV.mle3.mean <- mean(distanceList$TV.mle3)
TV.pt.mean <- mean(distanceList$TV.pt)
TV.gamma.mean <- mean(distanceList$TV.gamma)

KS.tNorm.mean <- mean(distanceList$KS.tNorm, na.rm = TRUE)
TV.tNorm.mean <- mean(distanceList$TV.tNorm, na.rm = TRUE)

print(paste("method.in =", method.in))
print(paste("K =",Ka))
print(paste("KS.mle1.mean =", KS.mle1.mean))
print(paste("KS.init1.mean =", KS.init1.mean))
print(paste("KS.mle2.mean =", KS.mle2.mean))
print(paste("KS.init2.mean =", KS.init2.mean))
print(paste("KS.mle3.mean =", KS.mle3.mean))
print(paste("KS.init3.mean =", KS.init3.mean))
print(paste("KS.pt.mean =", KS.pt.mean))
print(paste("KS.gamma.mean =", KS.gamma.mean))
print(paste("KS.tNorm.mean =", KS.tNorm.mean))
print(paste("TV.mle1.mean =", TV.mle1.mean))
print(paste("TV.init1.mean =", TV.init1.mean))
print(paste("TV.mle2.mean =", TV.mle2.mean))
print(paste("TV.init2.mean =", TV.init2.mean))
print(paste("TV.mle3.mean =", TV.mle3.mean))
print(paste("TV.init3.mean =", TV.init3.mean))
print(paste("TV.pt.mean =", TV.pt.mean))
print(paste("TV.gamma.mean =", TV.gamma.mean))
print(paste("TV.tNorm.mean =", TV.tNorm.mean))


save(list="distanceList", file = paste0("SIMU/Table1L_", LAMBDA0_SETTING_INDEX,"_", method.in, "_K=", Ka, ".RData"))

Profit.init1.mean <- mean(Profit_allList$Profit.init1)
Profit.mle1.mean <- mean(Profit_allList$Profit.mle1)
Profit.init2.mean <- mean(Profit_allList$Profit.init2)
Profit.mle2.mean <- mean(Profit_allList$Profit.mle2)
Profit.init3.mean <- mean(Profit_allList$Profit.init3)
Profit.mle3.mean <- mean(Profit_allList$Profit.mle3)
Profit.pt.mean <- mean(Profit_allList$Profit.pt)
Profit.gamma.mean <- mean(Profit_allList$Profit.gamma)
Profit.tNorm.mean <- mean(Profit_allList$Profit.tNorm, na.rm = TRUE)

print(paste("Profit.mle1.mean =", Profit.mle1.mean))
print(paste("Profit.init1.mean =", Profit.init1.mean))
print(paste("Profit.mle2.mean =", Profit.mle2.mean))
print(paste("Profit.init2.mean =", Profit.init2.mean))
print(paste("Profit.mle3.mean =", Profit.mle3.mean))
print(paste("Profit.init3.mean =", Profit.init3.mean))
print(paste("Profit.pt.mean =", Profit.pt.mean))
print(paste("Profit.gamma.mean =", Profit.gamma.mean))
print(paste("Profit.tNorm.mean =", Profit.tNorm.mean))

save(list="Profit_allList", file = paste0("SIMU/Table1L_", LAMBDA0_SETTING_INDEX,"_ProfitRes_", method.in, "_K=", Ka, ".RData"))
