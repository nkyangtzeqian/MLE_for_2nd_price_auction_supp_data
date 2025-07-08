rm(list=ls())

## Coordinate wise maximization 

cordwise_step_original <- function( data, lambda, theta){
  
  # Maximizing theta_i's one by one keeping other parameters fixed
  ellPK <- length(theta)
  
  K_sold <- data$K_sold  # This gives the vector of auction indexes for which the item is sold.
  ell <- length(data$Xbar_i)  # Total Number of observed bids for all the auctions.
  
  S_bar_max <- max(data$auxlist$S_bar)
  u_i <- data$auxlist$u_i
  l_i_3p2 <- data$auxlist$l_i_3p2
  Q_i_card <- data$auxlist$Q_i_card
  
  ## Step 3 
  for (ii in 1:S_bar_max){ 
    
    if(ii==1){
      AA <- lambda *  sum(data$Z_iT_i[ii:ellPK,2] * c(1,cumprod( theta[(ii+1):ellPK])) )
    }else if(ii < ellPK){
      AA <- lambda * prod(theta[1:(ii-1)]) * sum(data$Z_iT_i[ii:ellPK,2] * c(1,cumprod( theta[(ii+1):ellPK])) )
    }else if(ii==ellPK){
      AA <- lambda * prod(theta[1:(ii-1)]) * data$Z_iT_i[ellPK,2]
    }
    
    if(is.nan(AA)){
      warning(paste("AA", AA))
      break
    }
    
    BB <- Q_i_card[ii] + ( (1*(ell>0)) * (ell - l_i_3p2[ii]) ) ## for case III, automatically 0
    
    if(ii %in% u_i){ ## Formula 3.9, case I
      denomenator <- 2*AA
      numerator <- (AA+BB+1) - sqrt((AA+BB+1)^2 - (4*AA*BB))
      theta[ii] <- numerator/denomenator
      if ((numerator/denomenator) >= 1){
        print(paste("theta_i =", numerator/denomenator))
        theta[ii] <- 1
      }
    }else{ ## Formula 3.10, case II
      theta[ii] <- min(c(1,BB/AA))
    }
    
  } # END of for(ii in 1:ellPK) loop.
  
  ## Formula 3.11, case III
  if(S_bar_max < ellPK){
    theta[(S_bar_max+1):ellPK] <- 0  
  }
  
  return(list(lambda= lambda, theta = theta))
} # END of cordwise_step_constra function.

MLE_2ndprice_original <- function(data, lambda.in, theta.in, tol){
  
  if(length(data$K_sold)==0){   # if item is not sold in any of the auctions.
    warning("Item is not sold in any of the auctions.")
    theta.in <- rep(0, length(theta.in))
    F.mle <- rep(1,length(theta.in))
    Lik.out <- General_Log_Likelihood_PA(data, lambda.in, theta.in)
    return(list(Lik = Lik.out, theta.mle = theta.in, lambda = lambda.in,
                F.x = data$Z_iT_i[,1], F.y = F.mle, Lik.path = Lik.out))
  }
  
  Lik.in <- -1e100
  
  Lik.out <- General_Log_Likelihood_PA(data, lambda.in, theta.in)
  
  # print(paste("Lik.out initial =", Lik.out))
  
  Lik.path <- NULL
  
  while(Lik.out > Lik.in + tol){
    Lik.in <- Lik.out
    Lik.path <- c(Lik.path, Lik.in )
    
    par.out <- cordwise_step_original(data, lambda.in, theta.in)
    
    # print(par.out$lambda)
    # print(par.out$theta)
    
    theta.in <- par.out$theta
    lambda.in <- lambda_MLE_step(data,theta.in)
    
    # print(theta.in)
    
    Lik.out <- General_Log_Likelihood_PA(data, lambda.in, theta.in)
    
    # print(paste("Lik.out =", Lik.out))
    
    if(is.na(Lik.out) ){
      warning(paste("in", Lik.in, "out", Lik.out))
      Lik.out <- -1e100
    }
  }
  
  G.mle <- cumprod(theta.in)
  F.mle <- 1 - G.mle
  return(list(Lik = Lik.out, theta.mle = theta.in, lambda = lambda.in,
              F.x = data$Z_iT_i[,1], F.y = F.mle, Lik.path = Lik.path))
}


library(tidyverse)
library(GoFKernel) # For the inverse function of a CDF


# source("0327/EbayFunctions_Ver3.R")
source("Main.R")
source("distributions.R")
source("DataGeneration.R")
source("SettingGeneration.R")


########################################################################################################
## Plot
########################################################################################################




lambda0 <- 1
tau_a <- 100

Ka=150
method.in <- "gamma"
para1 <- 10
para2 <- 2


set.seed(123)

## Setting Generation
SETTING_ALL=Setting_Gen_Tab1(Ka, method = method.in, para1, para2)
x.med.vec <- SETTING_ALL$x.med.vec
r_k <- SETTING_ALL$r_k
reserve.price.Cutoff <- SETTING_ALL$reserve.price.Cutoff



## Simulation
Raw_data_all_list <- Data_Gen_Raw(Ka, tau_a, lambda0,
                                  r_k, method = method.in, para1, para2)
data <- Data_Gen_Observed(Raw_data_all_list)

Init.Values <- Initialization_lambda_F(data, reserve.price.Cutoff)
theta.init <- (1-Init.Values$F.y)/ c(1,(1-Init.Values$F.y)[-length(Init.Values$F.y)])
MLE <- MLE_2ndprice(data, lambda.in = Init.Values$lambda, theta.in = theta.init)

MLE_original <- MLE_2ndprice_original(data, lambda.in = Init.Values$lambda, theta.in = theta.init, tol = 1e-3)
MLE_original2 <- MLE_2ndprice_original(data, lambda.in = Init.Values$lambda, theta.in = theta.init, tol = 1e-8)


MLE$Lik
MLE_original$Lik
MLE_original2$Lik

## True F:
trueF <- Setting_Gen_Tab1_F0(data, method = method.in, para1, para2)[[1]]


# # GGplot of F0 vs. F.init vs. F.mle.origin vs F.mle
F.correct.tibble <- tibble(x.cord = c(0,Init.Values$F.x),
                           F.true = c(0,trueF),
                           F.init = c(0,Init.Values$F.y),
                           F.mle_original = c(0,MLE_original$F.y),
                           F.mle.original2 = c(0,MLE_original2$F.y),
                           F.mle = c(0,MLE$F.y)
                         )

( myplot2 <- ggplot(data = F.correct.tibble) +
  geom_line(mapping = aes(x = x.cord, y = F.true), color = "black") +
  geom_line(mapping = aes(x = x.cord, y = F.init), color = "blue") +
  geom_line(mapping = aes(x = x.cord, y = F.mle_original), color = "orange") +
  geom_line(mapping = aes(x = x.cord, y = F.mle.original2), color = "purple") +
  geom_line(mapping = aes(x = x.cord, y = F.mle), color = "red") +
  labs(x = "x", y = "F(x)",
       title = paste("Number of auctions =", Ka, ", True F:", method.in))
)
# ##################################################################################################################

## for high lambda case, we use the two step MLE function to avoid the extremely slow convergence of the MLE function with local optima
rm(list=ls())
lambda0 <- 10
method.in <- "unif"
para1 <- 1
para2 <- 20

## In the more complex case, we have less replications, F_fp for INIT, new order of cord-des and higher tolerance for MLE

library(tidyverse)
library(GoFKernel) # For the inverse function of a CDF


# source("0327/EbayFunctions_Ver3.R")
source("Main.R")
source("distributions.R")
source("DataGeneration.R")
source("DisCompare.R")
source("SettingGeneration.R")
source("PT_MCMC.R")


## when lambda*tau is large, we use this method to avoid numerical issue in F_sp,
Initialization_lambda_F <- Initialization_lambda_F_fp


#------------------------------------------------------------------------------------
## User Input: True lambda, Auction window, Number of independent auctions 'N.auction'
#------------------------------------------------------------------------------------


Ka <- 1000
tau_a <- 100

#------------------------------------------------------------------------------------
## Choice of different true valuation distributions F, the corresponding choice of parameters,
## and reserve prices are listed below:
#------------------------------------------------------------------------------------
set.seed(1)

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

###########################################
## might generate same value even from continuous CDF
TEMP_Z=c((0+min(data$Z_iT_i[,1]))/2,data$Z_iT_i[,1],max(data$Z_iT_i[,1])+1e-6)
for (x in which(TEMP_Z[-length(TEMP_Z)]==TEMP_Z[-1])) {
  TEMP_Z_MIN=max(TEMP_Z[x] - 1e-6,(TEMP_Z[x]+TEMP_Z[x-1])/2)
  TEMP_Z_MAX=min(TEMP_Z[x] + 1e-6,(TEMP_Z[x+1]+TEMP_Z[x+2])/2)
  TEMP_Z[x]=TEMP_Z_MIN
  TEMP_Z[x+1]=TEMP_Z_MAX
  data$Z_iT_i[x-1,1]=TEMP_Z_MIN
  data$Z_iT_i[x,1]=TEMP_Z_MAX
  rm(TEMP_Z_MIN, TEMP_Z_MAX)
}
rm(TEMP_Z)
###########################################

Init.Values <- Initialization_lambda_F(data, reserve.price.Cutoff)
theta.init <- (1-Init.Values$F.y)/ c(1,(1-Init.Values$F.y)[-length(Init.Values$F.y)])

time1=system.time({MLE_original <- MLE_2ndprice(data, lambda.in = Init.Values$lambda, theta.in = theta.init, tol = 8e-4)})

tail(MLE_original$Lik.path)

time2=system.time({MLE <- MLE_2ndprice_twosteps(data, lambda.in = Init.Values$lambda, theta.in = theta.init, tol = 1e-5)})

MLE$Lik.path ## larger than MLE_original$Lik.path


plot(MLE$F.x, MLE$F.y, type = "l", col = "blue", lwd = 2,
     xlab = "x", ylab = "F(x)", main = "MLE of F(x)")
lines(MLE_original$F.x, MLE_original$F.y, col = "orange", lwd = 2)

min(which(diff(MLE$Lik.path)<8e-4))
length(MLE_original$Lik.path)

## same starting point
plot(MLE_original$Lik.path[-1], type = "l", col = "orange", lwd = 2,
     xlab = "Iteration", ylab = "Log-likelihood",
     main = "Log-likelihood path of MLE (orange) and MLE_twosteps (blue)")
lines(MLE$Lik.path[-1], col = "blue", lwd = 2)
# abline(v = min(which(diff(MLE$Lik.path)<8e-4)), col = "red", lwd = 2)
abline(v=min(which(MLE$Lik.path>tail(MLE_original$Lik.path,1))), col = "red", lwd = 1.5)

save(list=ls(), file = "Plots/Example_MLE_2step_for_large_lambda.RData")

#===============================================================================

###################################################################################################################
# 
# #===============================================================================
# ## Illustration of instability of F_{MLE} near 0 using Gamma(10,2) distribution,
# ## when we use old MLE: MLE.old
# #===============================================================================
# 
# legend.vector <- c("True F" = "black", "Initial estimate of F" = "blue", 
#                    "MLE of F" = "orange", "Constrained MLE of F" = "red")
# F.tibble <- tibble("True_F" = pgamma(c(0,data$pooled.data[,1]), shape = para1, rate = para2),
#                    "Init_F" = c(0,Init.Values$F.y),
#                    "MLE_F" = c(0,MLE.old$F.y),
#                    "MLE_F_new" = c(0,MLE$F.y) )
# plot.gamma <- ggplot(data = F.tibble, aes(x = c(0,data$pooled.data[,1]))) +
#   geom_line(aes(y = True_F, color = "True F")) +
#   geom_line(aes(y = Init_F, color = "Initial estimate of F")) +
#   geom_line(aes(y = MLE_F, color = "MLE of F")) +
#   geom_line(aes(y = MLE_F_new, color = "Constrained MLE of F")) +
#   geom_vline(xintercept = min(data$pooled.observed.bids), color = "green") +
#   annotate("text", label = paste("min {all standing prices} = ", round(min(data$pooled.observed.bids),3)), 
#            x = min(data$pooled.observed.bids) + 2, y = -0.03) +
#   labs(x = "x", y = "F(x)", color = "",
#        title = paste("Number of auctions =", N.auction, ", True F is Gamma"))+
#   scale_color_manual(values = legend.vector) +
#   theme(legend.position = "top")
# 
# plot.gamma
# ggsave(filename = paste0(data$true.F,"_plots_",N.auction,".png"), 
#        plot = plot.gamma, height = 6, width = 8)

