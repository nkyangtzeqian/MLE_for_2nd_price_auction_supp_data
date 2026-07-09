#---------------------------------------------------------------------------------------------------
## Do trace plot from 1 repetitions of the training data
rm(list = ls())
source("MainL.R")
source("distributionsL.R")
source("DataGenerationL.R")
source("lambda_fun.R")
source("DisCompare.R")
source("PT_MCMC.R")

## to use only the fp conveniently
source("Empirical/Initialization_lambda_F_fp.R")

load(file = "Empirical/Xbox_7day_auctions.RData")

tau_a = raw_data_list$tau_a

LAMBDA_SET_LIST=list(
  constant=list(TWO_STAGE=FALSE, h_index=1),
  two_stage=list(TWO_STAGE=TRUE, tau_1=6),
  proportional=list(TWO_STAGE=FALSE, h_index=2)
)


NSimu=100
NRep=50


# Training dataset
prop_train <- 0.5
# prop_train <- 2/3

NROW <- floor((1-prop_train) * raw_data_list$Ka)*NSimu

res_tab  <- list(test_X_mk=matrix(NA,nrow=NROW,ncol=NRep),
                 simu_X_mk.fp1=matrix(NA,nrow=NROW,ncol=NRep),
                 simu_X_mk.fp2=matrix(NA,nrow=NROW,ncol=NRep),
                 simu_X_mk.fp3=matrix(NA,nrow=NROW,ncol=NRep),
                 simu_X_mk.init1=matrix(NA,nrow=NROW,ncol=NRep),
                 simu_X_mk.init2=matrix(NA,nrow=NROW,ncol=NRep),
                 simu_X_mk.init3=matrix(NA,nrow=NROW,ncol=NRep),
                 simu_X_mk.mle1=matrix(NA,nrow=NROW,ncol=NRep),
                 simu_X_mk.mle2=matrix(NA,nrow=NROW,ncol=NRep),
                 simu_X_mk.mle3=matrix(NA,nrow=NROW,ncol=NRep),
                 simu_X_mk.pt=matrix(NA,nrow=NROW,ncol=NRep)
)
rm(NROW)

for (replication in 1:NRep) {
  set.seed(replication)
  train.idx <- sample.int(raw_data_list$Ka, ceiling(prop_train * raw_data_list$Ka))
  train.idx <- sort(train.idx)
  
  
  train_data <- Data_Slice_Observed(raw_data_list,train.idx)
  
  reserve.price.Cutoff.train <- ceiling(quantile(train_data$r_k, probs = 0.25))
  
  InitFP1.train <- Initialization_lambda_F_fp(data = train_data, reserve.price.Cutoff.train, lambda_set = LAMBDA_SET_LIST[[1]])
  InitFP2.train <- Initialization_lambda_F_fp(data = train_data, reserve.price.Cutoff.train, lambda_set = LAMBDA_SET_LIST[[2]])
  InitFP3.train <- Initialization_lambda_F_fp(data = train_data, reserve.price.Cutoff.train, lambda_set = LAMBDA_SET_LIST[[3]])
  
  Init1.Values.train <- Initialization_lambda_F(data = train_data, reserve.price.Cutoff.train, lambda_set = LAMBDA_SET_LIST[[1]])
  theta.init.train <- (1-Init1.Values.train$F.y)/ c(1,(1-Init1.Values.train$F.y)[-length(Init1.Values.train$F.y)])
  MLE1.train <- MLE_2ndprice_twosteps(train_data, lambda.in = Init1.Values.train$lambda, theta.in = theta.init.train, lambda_set = LAMBDA_SET_LIST[[1]])
  cat(paste("  MLE1.train$lambda:", MLE1.train$lambda))
  
  Init2.Values.train <- Initialization_lambda_F(data = train_data, reserve.price.Cutoff.train, lambda_set = LAMBDA_SET_LIST[[2]])
  theta.init.train <- (1-Init2.Values.train$F.y)/ c(1,(1-Init2.Values.train$F.y)[-length(Init2.Values.train$F.y)])
  MLE2.train <- MLE_2ndprice_twosteps(train_data, lambda.in = Init2.Values.train$lambda, theta.in = theta.init.train, lambda_set = LAMBDA_SET_LIST[[2]])
  cat(paste("  MLE2.train$lambda:", MLE2.train$lambda))
  
  Init3.Values.train <- Initialization_lambda_F(data = train_data, reserve.price.Cutoff.train, lambda_set = LAMBDA_SET_LIST[[3]])
  theta.init.train <- (1-Init3.Values.train$F.y)/ c(1,(1-Init3.Values.train$F.y)[-length(Init3.Values.train$F.y)])
  MLE3.train <- MLE_2ndprice_twosteps(train_data, lambda.in = Init3.Values.train$lambda, theta.in = theta.init.train, lambda_set = LAMBDA_SET_LIST[[3]])
  cat(paste("  MLE3.train$lambda:", MLE3.train$lambda))
  
  ## Polya Tree method
  data_Polya <- Data_Gen_Unobserved_PolyaTree(Data_Slice_Raw(raw_data_list,train.idx))
  ## Since number of bidders unknown, we have to generate them
  train_N_all <- rpois(train_data$Ka, MLE1.train$lambda*train_data$tau_a)
  data_Polya[["dat"]][,1] <- train_N_all
  
  Bayes_PT.train <- MCMC_PT(data_Polya, train_data$Z_iT_i[,1], x.max=500)  
  
  RES_List <- list(fp1=InitFP1.train, fp2=InitFP2.train, fp3=InitFP3.train,
                   init1=Init1.Values.train, init2=Init2.Values.train, init3=Init3.Values.train,
                   mle1=MLE1.train, mle2=MLE2.train, mle3=MLE3.train,
                   pt=Bayes_PT.train)
  RES_Name_List <- c("fp1", "fp2", "fp3",
                     "init1", "init2", "init3",
                     "mle1", "mle2", "mle3",
                     "pt")
  
  ####################################################################################################
  ## Testing dataset
  test.idx <- setdiff(1:raw_data_list$Ka, train.idx)
  
  test_data <- Data_Slice_Observed(raw_data_list,test.idx)

  res_tab$test_X_mk[,replication] <- test_data$X_mk
  
  for(res in 1:length(RES_List)){
    res_name <- RES_Name_List[res]
    res_obj <- RES_List[[res]]
    e1_res.train <- AbscontDistribution(p = linear_interpolation_cdf(x = res_obj$F.x, F.x = res_obj$F.y) )
    
    simu_X_mk.res <- c()
    for (ii in 1:NSimu) {
      set.seed(ii+replication*NSimu)

      
      if(res_name == "pt"){
        # test_N_all <- sapply(1:length(test.idx), function(ii) nrow(Data_Slice_Raw(raw_data_list,test.idx)$Raw_biddata_list[ii][[1]]))
        test_N_all <- rpois(test_data$Ka, MLE1.train$lambda*test_data$tau_a)
        simu_X_mk.res <- sapply(1:test_data$Ka, function(jj){
          bid_ii=r(e1_res.train)(test_N_all[jj])
          bid_ii_2nd <- sort(bid_ii, decreasing = T)[2] # 2nd highest bid
          if(max(bid_ii) < test_data$r_k[jj]) # if no bid is above reserve price
            return(NA)
          else if (bid_ii_2nd < test_data$r_k[jj]) # if 2nd highest bid is below reserve price
            return(test_data$r_k[jj])
          else
            return(bid_ii_2nd)
        })
      }else{
        LAMBDA0_SETTING_INDEX <- as.numeric(substr(res_name, nchar(res_name), nchar(res_name)))
        if(length(res_obj$lambda)== 1) {lambda0_1 <- res_obj$lambda; lambda0_2 <- NULL;}
        if(length(res_obj$lambda)== 2) {lambda0_1 <- res_obj$lambda[1]; lambda0_2 <- res_obj$lambda[2];}
        lambda0_list <- lambda0_list_init(lambda0_1,lambda0_2,LAMBDA_SET_LIST[[LAMBDA0_SETTING_INDEX]])
        simu.res.test <- Data_Gen_Observed(Data_Gen_Raw(test_data$Ka, test_data$tau_a, lambda0_list,
                                                        test_data$r_k, method = "customized", customized_CDF = e1_res.train))
        simu_X_mk.res=c(simu_X_mk.res, simu.res.test$X_mk)
      }
    }
    res_tab[[paste0("simu_X_mk.", res_name)]][,replication] <- simu_X_mk.res
  }
  
  #########################################################################################################
  print(paste("Replication", replication, "is done."))
  
}



# 
# plot(density(test_X_mk,adjust=2), col = "black", main = "Density of simulated final selling prices", xlab = "X_mk",lwd=2)
# lines(density(simu_X_mk.init[!is.na(simu_X_mk.init)],adjust=2), col = "blue", lwd=1)
# lines(density(simu_X_mk.mle[!is.na(simu_X_mk.mle)],adjust=2), col = "red", lwd=1)
# lines(density(simu_X_mk.pt[!is.na(simu_X_mk.pt)],adjust=2), col = "green", lwd=1)



##############################################################################################



MSE_Simu <- function(simuM){
  MSE_rep=c()
  for (ii in 1:NRep) {
    RepM=matrix(simuM[,ii],ncol = NSimu)
    testV=res_tab$test_X_mk[,ii]
    ## for mean across auctions
    MSEM=apply(RepM, 2, function(x) mean((x-testV)^2 ,na.rm = T)  )
    ## mean across simulations
    MSE_rep=c(MSE_rep, mean(MSEM))
  }
  ## mean across replications
  return(MSE_rep)
}
mean(MSE_Simu(res_tab$simu_X_mk.fp1))
mean(MSE_Simu(res_tab$simu_X_mk.fp2))/mean(MSE_Simu(res_tab$simu_X_mk.fp1))
mean(MSE_Simu(res_tab$simu_X_mk.fp3))/mean(MSE_Simu(res_tab$simu_X_mk.fp1))
mean(MSE_Simu(res_tab$simu_X_mk.init1))/mean(MSE_Simu(res_tab$simu_X_mk.fp1))
mean(MSE_Simu(res_tab$simu_X_mk.init2))/mean(MSE_Simu(res_tab$simu_X_mk.fp1))
mean(MSE_Simu(res_tab$simu_X_mk.init3))/mean(MSE_Simu(res_tab$simu_X_mk.fp1))
mean(MSE_Simu(res_tab$simu_X_mk.mle1))/mean(MSE_Simu(res_tab$simu_X_mk.fp1))
mean(MSE_Simu(res_tab$simu_X_mk.mle2))/mean(MSE_Simu(res_tab$simu_X_mk.fp1))
mean(MSE_Simu(res_tab$simu_X_mk.mle3))/mean(MSE_Simu(res_tab$simu_X_mk.fp1))
mean(MSE_Simu(res_tab$simu_X_mk.pt))/mean(MSE_Simu(res_tab$simu_X_mk.fp1))

save(list=ls(), file = paste0("Empirical/Xbox_simuL_", round(prop_train,2),".RData"))
