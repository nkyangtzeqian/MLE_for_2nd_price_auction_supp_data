#---------------------------------------------------------------------------------------------------
## Do trace plot from 1 repetitions of the training data
rm(list = ls())
source("Main.R")
source("distributions.R")
source("DataGeneration.R")
source("DisCompare.R")
source("PT_MCMC.R")



load(file = "Empirical/Xbox_7day_auctions.RData")

test_X_mk=c()
simu_X_mk.fp=c()
simu_X_mk.init=c()
simu_X_mk.mle=c()
simu_X_mk.pt=c()

NSimu=100
NRep=100

for (replication in 1:NRep) {
  # Training dataset
  # prop_train <- 0.5
  prop_train <- 2/3
  set.seed(replication)
  train.idx <- sample.int(raw_data_list$Ka, ceiling(prop_train * raw_data_list$Ka))
  train.idx <- sort(train.idx)
  
  
  train_data <- Data_Slice_Observed(raw_data_list,train.idx)
  
  reserve.price.Cutoff.train <- ceiling(quantile(train_data$r_k, probs = 0.25))
  
  InitFP.train <- Initialization_lambda_F_fp(data = train_data, reserve.price.Cutoff.train)
  
  Init.Values.train <- Initialization_lambda_F(data = train_data, reserve.price.Cutoff.train)
  theta.init.train <- (1-Init.Values.train$F.y)/ c(1,(1-Init.Values.train$F.y)[-length(Init.Values.train$F.y)])
  MLE.train <- MLE_2ndprice(train_data, lambda.in = Init.Values.train$lambda, theta.in = theta.init.train, tol = 1e-5)
  
  ## Polya Tree method
  data_Polya <- Data_Gen_Unobserved_PolyaTree(Data_Slice_Raw(raw_data_list,train.idx))
  ## Since number of bidders unknown, we have to generate them
  train_N_all <- rpois(train_data$Ka, MLE.train$lambda*train_data$tau_a)
  data_Polya[["dat"]][,1] <- train_N_all
  
  Bayes_PT.train <- MCMC_PT(data_Polya, train_data$Z_iT_i[,1], x.max=500)  
  
  ####################################################################################################
  
  e1.fp.train <- AbscontDistribution(p = linear_interpolation_cdf(x = InitFP.train$F.x,
                                                                F.x = InitFP.train$F.y) )
  
  e1.init.train <- AbscontDistribution(p = linear_interpolation_cdf(x = Init.Values.train$F.x,
                                                                    F.x = Init.Values.train$F.y) )
  
  e1.mle.train <- AbscontDistribution(p = linear_interpolation_cdf(x = MLE.train$F.x,
                                                                   F.x = MLE.train$F.y) )
  
  e1.pt.train <- AbscontDistribution(p = linear_interpolation_cdf(x = Bayes_PT.train$F.x,
                                                                  F.x = Bayes_PT.train$F.y) )
  
  ## Testing dataset
  test.idx <- setdiff(1:raw_data_list$Ka, train.idx)
  
  test_data <- Data_Slice_Observed(raw_data_list,test.idx)
  
  test_X_mk=c(test_X_mk, test_data$X_mk)
  
  ##########################################################################################################
  
  
  for (ii in 1:NSimu) {
    set.seed(ii+replication*NSimu)
    simu.fp.test <- Data_Gen_Observed(Data_Gen_Raw(test_data$Ka, test_data$tau_a, InitFP.train$lambda,
                                                test_data$r_k, method = "customized", customized_CDF = e1.fp.train))
    
    simu.init.test <- Data_Gen_Observed(Data_Gen_Raw(test_data$Ka, test_data$tau_a, Init.Values.train$lambda,
                                                     test_data$r_k, method = "customized", customized_CDF = e1.init.train))
    
    simu.mle.test <- Data_Gen_Observed(Data_Gen_Raw(test_data$Ka, test_data$tau_a, MLE.train$lambda,
                                                    test_data$r_k, method = "customized", customized_CDF = e1.mle.train))
    
    # test_N_all <- sapply(1:length(test.idx), function(ii) nrow(Data_Slice_Raw(raw_data_list,test.idx)$Raw_biddata_list[ii][[1]]))
    test_N_all <- rpois(test_data$Ka, MLE.train$lambda*test_data$tau_a)
    simu.pt.test <- sapply(1:test_data$Ka, function(jj){
      bid_ii=r(e1.pt.train)(test_N_all[jj])
      bid_ii_2nd <- sort(bid_ii, decreasing = T)[2] # 2nd highest bid
      if(max(bid_ii) < test_data$r_k[jj]) # if no bid is above reserve price
        return(NA)
      else if (bid_ii_2nd < test_data$r_k[jj]) # if 2nd highest bid is below reserve price
        return(test_data$r_k[jj])
      else
        return(bid_ii_2nd)
    })
    
    simu_X_mk.fp=c(simu_X_mk.fp, simu.fp.test$X_mk)
    simu_X_mk.init=c(simu_X_mk.init, simu.init.test$X_mk)
    simu_X_mk.mle=c(simu_X_mk.mle, simu.mle.test$X_mk)
    simu_X_mk.pt=c(simu_X_mk.pt, simu.pt.test)
  }
  
  print(paste("Replication", replication, "is done."))

}



# 
# plot(density(test_X_mk,adjust=2), col = "black", main = "Density of simulated final selling prices", xlab = "X_mk",lwd=2)
# lines(density(simu_X_mk.init[!is.na(simu_X_mk.init)],adjust=2), col = "blue", lwd=1)
# lines(density(simu_X_mk.mle[!is.na(simu_X_mk.mle)],adjust=2), col = "red", lwd=1)
# lines(density(simu_X_mk.pt[!is.na(simu_X_mk.pt)],adjust=2), col = "green", lwd=1)



simu_X_mk.fp=matrix(simu_X_mk.fp,ncol = NRep)
simu_X_mk.init=matrix(simu_X_mk.init,ncol = NRep)
simu_X_mk.mle=matrix(simu_X_mk.mle,ncol = NRep)
simu_X_mk.pt=matrix(simu_X_mk.pt,ncol = NRep)


test_X_mk=matrix(test_X_mk,ncol = NRep)

##############################################################################################



MSE_Simu <- function(simuM){
  MSE_rep=c()
  for (ii in 1:NRep) {
    RepM=matrix(simuM[,ii],ncol = NSimu)
    testV=test_X_mk[,ii]
    ## for mean across auctions
    MSEM=apply(RepM, 2, function(x) mean((x-testV)^2 ,na.rm = T)  )
    ## mean across simulations
    MSE_rep=c(MSE_rep, mean(MSEM))
  }
  ## mean across replications
  return(MSE_rep)
}
mean(MSE_Simu(simu_X_mk.fp))
mean(MSE_Simu(simu_X_mk.init))/mean(MSE_Simu(simu_X_mk.fp))
mean(MSE_Simu(simu_X_mk.mle))/mean(MSE_Simu(simu_X_mk.fp))
mean(MSE_Simu(simu_X_mk.pt))/mean(MSE_Simu(simu_X_mk.fp))

save(list=ls(), file = paste0("Empirical/Xbox_simu_", round(prop_train,2),".RData"))
