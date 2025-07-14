# rm(list = ls())   # Remove everything from the Global environment
# Ka=100
# NSimu=100
# method.in <- "unif"
# para1 <- 1
# para2 <- 20
# Delta <- 0
# source("Simu2_2_HulC_TruncCoverage.R")



library(tidyverse)
# library(transport) # For calculating Wasserstein distance between two distribution functions.

source("Main.R")
source("distributions.R")
source("DataGeneration.R")
# source("DisCompare.R")
source("SettingGeneration.R")
source("PT_MCMC.R")


source("HulC_main.R")

#------------------------------------------------------------------------------------
## User Input: True lambda, Auction window, Number of independent auctions 'N.auction'
#------------------------------------------------------------------------------------

# repetition_number <- 100  # replication numbers of each such event i.e. selling an 
# identical stuff on ebay with multiple independent auctions.

lambda0 <- 1
tau_a <- 100
# Ka <- 150

#------------------------------------------------------------------------------------
## Choice of different true valuation distributions F, the corresponding choice of parameters,  
## and reserve prices are listed below:
#------------------------------------------------------------------------------------
SETTING_ALL=Setting_Gen_Tab1(Ka, method = method.in, para1, para2)
x.med.vec <- SETTING_ALL$x.med.vec
r_k <- SETTING_ALL$r_k
reserve.price.Cutoff <- SETTING_ALL$reserve.price.Cutoff

## True F:
trueF <- Setting_Gen_Tab1_F0_cdf(method = method.in, para1, para2)[[1]]

trueF_point <- sapply(x.med.vec, trueF)
## For truncated coverage
offbound_x <- which(trueF_point > 0.05 & trueF_point <0.95)
x.med.vec <- x.med.vec[offbound_x]
trueF_point <- trueF_point[offbound_x]


HulC_Coverage_Simu<- function(seed){
  set.seed(seed)
  #------------------------------------------------------------------------------------
  ## Data generation, Initial Estimate, MLE, Polya Tree, HUlC and credible interval
  #------------------------------------------------------------------------------------
  
  Raw_data_all_list <- Data_Gen_Raw(Ka, tau_a, lambda0,
                                    r_k, method = method.in, para1, para2)
  data <- Data_Gen_Observed(Raw_data_all_list)
  
  Init.Values <- Initialization_lambda_F(data, reserve.price.Cutoff)
  theta.init <- (1-Init.Values$F.y)/ c(1,(1-Init.Values$F.y)[-length(Init.Values$F.y)])
  MLE <- MLE_2ndprice(data, lambda.in = Init.Values$lambda, theta.in = theta.init)
  
  Hulc.Conf <- HulC1d_Ebay(raw_data_all_list = Raw_data_all_list, eval.points = x.med.vec, alpha=.1,
                           MLE.Cutoff = reserve.price.Cutoff, diff_bias = 0, Delta=Delta)
  
  ## Polya Tree method
  data_Polya <- Data_Gen_Unobserved_PolyaTree(Raw_data_all_list)
  Bayes_PT <- MCMC_PT(data_Polya, data$Z_iT_i[,1], x.max=max(x.med.vec))
  
  Prediction_PT <- sapply(1:nrow(Bayes_PT$MCMC_Mat), function(ii) {
    cutpt <- c(0, rev(data_Polya$dat[,2]), max(x.med.vec))
    step_est <- approxfun(cutpt, c(0, cumsum(Bayes_PT$MCMC_Mat[ii,])  ))
    return(step_est(x.med.vec))
  })
  
  
  
  Band.CI.Init.lwr = Hulc.Conf$Band.CI.Init$lwr(x.med.vec)
  Band.CI.Init.upr = Hulc.Conf$Band.CI.Init$upr(x.med.vec)
  Band.CI.MLE.lwr = Hulc.Conf$Band.CI.MLE$lwr(x.med.vec)
  Band.CI.MLE.upr = Hulc.Conf$Band.CI.MLE$upr(x.med.vec)
  
  Band.CI.PT.lwr = apply(Prediction_PT, 1, function(x) quantile(x, probs = 0.05))
  Band.CI.PT.upr = apply(Prediction_PT, 1, function(x) quantile(x, probs = 0.95))
  
  ## avoid small numerical errors
  TOL= 1e-3
  Bandwidth.Init = Band.CI.Init.upr - Band.CI.Init.lwr
  Bandwidth.MLE = Band.CI.MLE.upr - Band.CI.MLE.lwr
  Bandwidth.PT = Band.CI.PT.upr - Band.CI.PT.lwr
  Coverage.Init = (trueF_point >= Band.CI.Init.lwr-TOL) & (trueF_point <= Band.CI.Init.upr+TOL)
  Coverage.MLE = (trueF_point >= Band.CI.MLE.lwr-TOL) & (trueF_point <= Band.CI.MLE.upr+TOL)
  Coverage.PT = (trueF_point >= Band.CI.PT.lwr-TOL) & (trueF_point <= Band.CI.PT.upr+TOL)
  
  
  
  
  return(list(
    x.med.vec = x.med.vec,
    Band.CI.Init.lwr = Band.CI.Init.lwr,
    Band.CI.Init.upr = Band.CI.Init.upr,
    Band.CI.MLE.lwr = Band.CI.MLE.lwr,
    Band.CI.MLE.upr = Band.CI.MLE.upr,
    Band.CI.PT.lwr = Band.CI.PT.lwr,
    Band.CI.PT.upr = Band.CI.PT.upr,
    trueF_point = trueF_point,
    Bandwidth.Init = Bandwidth.Init,
    Bandwidth.MLE = Bandwidth.MLE,
    Bandwidth.PT = Bandwidth.PT,
    Coverage.Init = Coverage.Init,
    Coverage.MLE = Coverage.MLE,
    Coverage.PT = Coverage.PT
  ))
  
}

Mat_Bw_Init <- matrix(NA, nrow = NSimu, ncol = length(x.med.vec))
Mat_Bw_MLE <- matrix(NA, nrow = NSimu, ncol = length(x.med.vec))
Mat_Bw_PT <- matrix(NA, nrow = NSimu, ncol = length(x.med.vec))

Mat_Cov_Init <- matrix(NA, nrow = NSimu, ncol = length(x.med.vec))
Mat_Cov_MLE <- matrix(NA, nrow = NSimu, ncol = length(x.med.vec))
Mat_Cov_PT <- matrix(NA, nrow = NSimu, ncol = length(x.med.vec))

for (seed in 1:NSimu) {
  err <- try({res <- HulC_Coverage_Simu(seed);})
  if(inherits(err, "try-error"))
  {
    # print(paste0("Simulation ", seed, " failed."))
    next
  }
  else{
    Mat_Bw_Init[seed, ] <- res$Bandwidth.Init
    Mat_Bw_MLE[seed, ] <- res$Bandwidth.MLE
    Mat_Bw_PT[seed, ] <- res$Bandwidth.PT
    Mat_Cov_Init[seed, ] <- res$Coverage.Init
    Mat_Cov_MLE[seed, ] <- res$Coverage.MLE
    Mat_Cov_PT[seed, ] <- res$Coverage.PT
    print(paste0("Simulation ", seed, " completed."))
  }
}

save(list = ls(),
     file = paste0("SIMU/HulC_TrunCov_", method.in, " With No. Auction = ", Ka, ".RData"))
