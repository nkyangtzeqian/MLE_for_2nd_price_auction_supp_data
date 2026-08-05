#---------------------------------------------------------------------------------------------------
## Do trace plot from 1 repetitions of the training data
rm(list = ls())
source("MainL.R")
source("distributionsL.R")
source("DataGenerationL.R")
source("Empirical/lambda_fun_Ebay.R")
source("DisCompare.R")
source("PT_MCMC.R")



load(file = "Empirical/Xbox_7day_auctions.RData")

tau_a <- raw_data_list$tau_a ## 7

tau_1_list <- 1:6

res_est_list <- matrix(NA, nrow = length(tau_1_list), ncol = 3)


train_data <- Data_Gen_Observed(raw_data_list)

reserve.price.Cutoff.train <- ceiling(quantile(train_data$r_k, probs = 0.25))

for (ii in 1:length(tau_1_list)) {
  
  lambda_set <- list(TWO_STAGE=TRUE, tau_1=tau_1_list[ii])

  Init.Values.train <- Initialization_lambda_F(data = train_data, reserve.price.Cutoff.train, lambda_set = lambda_set, MAX_L_EST=10)
  theta.init.train <- (1-Init.Values.train$F.y)/ c(1,(1-Init.Values.train$F.y)[-length(Init.Values.train$F.y)])
  print(Init.Values.train$lambdas)
  
  MLE.train <- MLE_2ndprice_twosteps(train_data, lambda.in = Init.Values.train$lambda, theta.in = theta.init.train, lambda_set = lambda_set)

  res_est_list[ii,1:2] <- MLE.train$lambdas
  res_est_list[ii,3] <- MLE.train$Lik 
  print(MLE.train$Lik)
}

colnames(res_est_list) <- c("lambda_1", "lambda_2", "logLik")
row.names(res_est_list) <- paste0("tau_1=", tau_1_list)

##################################################################################
res_est_list <- rbind(res_est_list, matrix(NA, nrow = 2, ncol = 3))

lambda_set <- list(TWO_STAGE=FALSE,  h_index=1)

Init.Values.train <- Initialization_lambda_F(data = train_data, reserve.price.Cutoff.train, lambda_set = lambda_set, MAX_L_EST=10)
theta.init.train <- (1-Init.Values.train$F.y)/ c(1,(1-Init.Values.train$F.y)[-length(Init.Values.train$F.y)])
print(Init.Values.train$lambdas)

MLE.train <- MLE_2ndprice_twosteps(train_data, lambda.in = Init.Values.train$lambda, theta.in = theta.init.train, lambda_set = lambda_set)
print(MLE.train$Lik)

res_est_list[7,1:2] <- MLE.train$lambdas
res_est_list[7,3] <- MLE.train$Lik 

##################################################################################
lambda_set <- list(TWO_STAGE=FALSE,  h_index=2)

Init.Values.train <- Initialization_lambda_F(data = train_data, reserve.price.Cutoff.train, lambda_set = lambda_set, MAX_L_EST=10)
theta.init.train <- (1-Init.Values.train$F.y)/ c(1,(1-Init.Values.train$F.y)[-length(Init.Values.train$F.y)])
print(Init.Values.train$lambdas)

MLE.train <- MLE_2ndprice_twosteps(train_data, lambda.in = Init.Values.train$lambda, theta.in = theta.init.train, lambda_set = lambda_set)
print(MLE.train$Lik)

res_est_list[8,1:2] <- MLE.train$lambdas
res_est_list[8,3] <- MLE.train$Lik 

row.names(res_est_list)[7:8] <- c("Constant", "Linear")

res_est_list
