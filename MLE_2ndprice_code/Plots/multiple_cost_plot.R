# rm(list = ls())   # Remove everything from the Global environment
# method.in <- "gamma"
# para1 <- 10
# para2 <- 2
# source("Plots/multiple_cost_plot.R") 

Ka <- 100

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

## set cost to be 0.25th quantile of the distribution
if(method.in == "unif"){
  c_list=qunif(seq(0.05,0.95,by=0.05), min=para1, max=para2) # for uniform
}else if( method.in == "piecewise_unif"){
  c_list=sapply(seq(0.05,0.95,by=0.05), 
               function(u) {piecewise.unif.invcdf(u, para1 = para1, para2 = para2)})
}else if (method.in == "pareto"){
  c_list=qpareto(seq(0.05,0.95,by=0.05), m = para1, s = para2) # for pareto
}else if( method.in == "gamma"){
  c_list=qgamma(seq(0.05,0.95,by=0.05), shape = para1, rate=para2) # for gamma
}else if (method.in == "beta"){  
  c_list=qbeta(seq(0.05,0.95,by=0.05), shape1 = para1, shape2 = para2) # for beta
}else{
  stop("The method is not supported. Please use one of the following methods: 
         'unif', 'pareto', 'gamma', 'beta', 'piecewise_unif'.")
}



## copied from ProfitMetric.R
ProfitMetric_test <- function(e1, cost,                              
                         method = c("unif", "pareto", "gamma", "beta",
                                    "piecewise_unif"),
                         para1, para2){
  
  ## codes from George and Hui (2012) to compute Revenue and price
  Compute.Opt.Price <- function(cdf, cost=5, max.price=20){
    result <- list(opt.price=numeric(1), profit=numeric(1))
    tempgrid <- seq(cost, max.price, 0.001)
    val <- sapply(tempgrid, Profit.Func, cdf=cdf, cost=cost)
    result$opt.price <- tempgrid[which.max(val)]
    result$profit <- max(val)
    return(result)
  }
  
  Profit.Func <- function(x, cdf, cost){
    return((1 - cdf(x)) * (x - cost))
  }
  
  
  
  max_cutpt <- cost_def(method, para1, para2)["max_cutpt"]
  # cost <- cost_def(method, para1, para2)["cost"]
  
  e1_cdf <- p(e1)
  
  
  Profit_all <- Compute.Opt.Price(e1_cdf, cost, max.price=max_cutpt)
  max_price <- Profit_all$opt.price
  est_Profit <- Profit_all$profit
  Profit <- Profit.Func(max_price, Setting_Gen_Tab1_F0_cdf(method,para1, para2)[[1]], cost)
  
  return(list(max_price = max_price,
              Profit = Profit,
              est_Profit = est_Profit))
}

## accept multiple costs
ProfitMetric_test_list <- function(e1,                              
                              method = c("unif", "pareto", "gamma", "beta",
                                         "piecewise_unif"),
                              para1, para2){
  max_price_vec <- rep(NA, length(c_list))
  Profit_vec <- rep(NA, length(c_list))
  est_Profit_vec <- rep(NA, length(c_list))
  
  for(i in 1:length(c_list)){
    cost <- c_list[i]
    Profit_all <- ProfitMetric_test(e1, cost, method, para1, para2)
    max_price_vec[i] <- Profit_all$max_price
    Profit_vec[i] <- Profit_all$Profit
    est_Profit_vec[i] <- Profit_all$est_Profit
  }
  return(list(max_price_vec = max_price_vec,
              Profit_vec = Profit_vec,
              est_Profit_vec = est_Profit_vec))
}

Profit_allList <- list(
  max_price.init= matrix(NA, length(c_list), replication),
  max_price.mle=matrix(NA, length(c_list), replication),
  max_price.pt=matrix(NA, length(c_list), replication),
  max_price.gamma=matrix(NA, length(c_list), replication),
  max_price.tNorm=matrix(NA, length(c_list), replication),
  Profit.init =matrix(NA, length(c_list), replication),
  Profit.mle =matrix(NA, length(c_list), replication),
  Profit.pt =matrix(NA, length(c_list), replication),
  Profit.gamma =matrix(NA, length(c_list), replication),
  Profit.tNorm =matrix(NA, length(c_list), replication),
  est_Profit.init =matrix(NA, length(c_list), replication),
  est_Profit.mle =matrix(NA, length(c_list), replication),
  est_Profit.pt =matrix(NA, length(c_list), replication),
  est_Profit.gamma =matrix(NA, length(c_list), replication),
  est_Profit.tNorm =matrix(NA, length(c_list), replication)
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
  
  ## Gamma method
  Bayes_gamma <- MCMC_gamma(data_Polya, data$Z_iT_i[,1])
  ## Truncated Normal method
  ERROR1=try({Bayes_tNorm <- MCMC_tNorm(data_Polya, data$Z_iT_i[,1])})
  if(class(ERROR1)=="try-error"){ ## strange error might occur in AbscontDistribution function due to randomness
    print("Error in MCMC_tNorm")
    next
  }
  
  
  
  e1.init <- AbscontDistribution(p = linear_interpolation_cdf(x = Init.Values$F.x, F.x = Init.Values$F.y) )
  e1.mle <- AbscontDistribution(p = linear_interpolation_cdf(x = MLE$F.x, F.x = MLE$F.y) )
  e1.pt <- AbscontDistribution(p = linear_interpolation_cdf(x = Bayes_PT$F.x, F.x = Bayes_PT$F.y) )
  e1.gamma <- AbscontDistribution(p = linear_interpolation_cdf(x = Bayes_gamma$F.x, F.x = Bayes_gamma$F.y) )
  ERROR2 <- try({e1.tNorm <- AbscontDistribution(p = linear_interpolation_cdf(x = Bayes_tNorm$F.x, F.x = Bayes_tNorm$F.y) )
  })
  if(class(ERROR2)=="try-error"){ ## strange error might occur in AbscontDistribution function due to randomness
    print("Error in e1.tNorm")
    e1.tNorm <- AbscontDistribution(p = linear_interpolation_cdf(x = Bayes_tNorm$F.x, F.x = Bayes_tNorm$F.y) )
  }
  
  print(paste("ii =", ii, "fitted"))
  
  ## Profit Metric:
  Profit_init_list <- ProfitMetric_test_list(e1 = e1.init, method = method.in, para1, para2)
  Profit_mle_list <- ProfitMetric_test_list(e1 = e1.mle, method = method.in, para1, para2)
  Profit_pt_list <- ProfitMetric_test_list(e1 = e1.pt, method = method.in, para1, para2)
  Profit_gamma_list <- ProfitMetric_test_list(e1 = e1.gamma, method = method.in, para1, para2)
  Profit_tNorm_list <- ProfitMetric_test_list(e1 = e1.tNorm, method = method.in, para1, para2)
  
  Profit_allList$max_price.init[ , ii] <- Profit_init_list[[1]]
  Profit_allList$max_price.mle[ , ii] <- Profit_mle_list[[1]]
  Profit_allList$max_price.pt[ , ii] <- Profit_pt_list[[1]]
  Profit_allList$max_price.gamma[ , ii] <- Profit_gamma_list[[1]]
  Profit_allList$max_price.tNorm[ , ii] <- Profit_tNorm_list[[1]]
  
  Profit_allList$Profit.init[ , ii] <- Profit_init_list[[2]]
  Profit_allList$Profit.mle[ , ii] <- Profit_mle_list[[2]]
  Profit_allList$Profit.pt[ , ii] <- Profit_pt_list[[2]]
  Profit_allList$Profit.gamma[ , ii] <- Profit_gamma_list[[2]]
  Profit_allList$Profit.tNorm[ , ii] <- Profit_tNorm_list[[2]]
  
  Profit_allList$est_Profit.init[ , ii] <- Profit_init_list[[3]]
  Profit_allList$est_Profit.mle[ , ii] <- Profit_mle_list[[3]]
  Profit_allList$est_Profit.pt[ , ii] <- Profit_pt_list[[3]]
  Profit_allList$est_Profit.gamma[ , ii] <- Profit_gamma_list[[3]]
  Profit_allList$est_Profit.tNorm[ , ii] <- Profit_tNorm_list[[3]]
  
  print(paste("ii =", ii, "finished"))
}

Profit_mean_list=list(
  max_price.init= rowMeans(Profit_allList$max_price.init, na.rm=TRUE),
  max_price.mle=rowMeans(Profit_allList$max_price.mle, na.rm=TRUE),
  max_price.pt=rowMeans(Profit_allList$max_price.pt, na.rm=TRUE),
  max_price.gamma=rowMeans(Profit_allList$max_price.gamma, na.rm=TRUE),
  max_price.tNorm=rowMeans(Profit_allList$max_price.tNorm, na.rm=TRUE),
  Profit.init =rowMeans(Profit_allList$Profit.init, na.rm=TRUE),
  Profit.mle =rowMeans(Profit_allList$Profit.mle, na.rm=TRUE),
  Profit.pt =rowMeans(Profit_allList$Profit.pt, na.rm=TRUE),
  Profit.gamma =rowMeans(Profit_allList$Profit.gamma, na.rm=TRUE),
  Profit.tNorm =rowMeans(Profit_allList$Profit.tNorm, na.rm=TRUE),
  est_Profit.init =rowMeans(Profit_allList$est_Profit.init, na.rm=TRUE),
  est_Profit.mle =rowMeans(Profit_allList$est_Profit.mle, na.rm=TRUE),
  est_Profit.pt =rowMeans(Profit_allList$est_Profit.pt, na.rm=TRUE),
  est_Profit.gamma =rowMeans(Profit_allList$est_Profit.gamma, na.rm=TRUE),
  est_Profit.tNorm =rowMeans(Profit_allList$est_Profit.tNorm, na.rm=TRUE)
)

## compute default
## copy from Simu1_Plotprofit.R
Profit_default_test <- function(cost,method = c("unif", "pareto", "gamma", "beta",
                                           "piecewise_unif")){
  
  ## codes from George and Hui (2012) to compute Revenue and price
  Compute.Opt.Price <- function(cdf, cost=5, max.price=20){
    result <- list(opt.price=numeric(1), profit=numeric(1))
    tempgrid <- seq(cost, max.price, 0.001)
    val <- sapply(tempgrid, Profit.Func, cdf=cdf, cost=cost)
    result$opt.price <- tempgrid[which.max(val)]
    result$profit <- max(val)
    return(result)
  }
  
  Profit.Func <- function(x, cdf, cost){
    return((1 - cdf(x)) * (x - cost))
  }
  
  source("distributions.R")
  source("SettingGeneration.R")
  para1=Setting_Gen_Tab1(method,K=100,para1=NULL, para2=NULL)$para1
  para2=Setting_Gen_Tab1(method,K=100,para1=NULL, para2=NULL)$para2
  
  source("ProfitMetric.R")
  max_cutpt <- cost_def(method, para1, para2)["max_cutpt"]
  # cost <- cost_def(method, para1, para2)["cost"]
  
  e_cdf <-Setting_Gen_Tab1_F0_cdf(method,para1, para2)[[1]]
  
  
  Profit_all <- Compute.Opt.Price(e_cdf, cost, max.price=max_cutpt)
  max_price <- Profit_all$opt.price
  Profit <- Profit_all$profit
  
  return(list(max_price = max_price,
              Profit = Profit))
}

for(i in 1:length(c_list)){
  cost <- c_list[i]
  Profit_all <- Profit_default_test(cost, method.in)
  Profit_mean_list$max_price.default[i] <- Profit_all$max_price
  Profit_mean_list$Profit.default[i] <- Profit_all$Profit
}


save(list=c("Profit_allList","Profit_mean_list"), file = paste0("Plots/Multiple_costs_trace_",method.in,".RData"))
# #######################################################################################
# ## Plots
# #######################################################################################
# 
# ## max price comparison
# plot(x=seq(0.05,0.95,by=0.05), y=Profit_mean_list$max_price.default, type="l", col="black",
#      xlab="Cost (at quantile of F)", ylab="Optimal Price",
#      main="Optimal Price vs. Cost (Gamma True F)")
# lines(x=seq(0.05,0.95,by=0.05), y=Profit_mean_list$max_price.init, type="l", col="blue")
# lines(x=seq(0.05,0.95,by=0.05), y=Profit_mean_list$max_price.mle, type="l", col="red")
# lines(x=seq(0.05,0.95,by=0.05), y=Profit_mean_list$max_price.pt, type="l", col="green")
# lines(x=seq(0.05,0.95,by=0.05), y=Profit_mean_list$max_price.gamma, type="l", col="purple")
# lines(x=seq(0.05,0.95,by=0.05), y=Profit_mean_list$max_price.tNorm, type="l", col="orange")
# 
# legend("topleft", legend=c("Default","MLE","Initial","Polya Tree","Gamma-F","TN-F"),
#        col=c("black","red","blue","green","purple","orange"), lty=1)
# 
# ## True Profit comparison
# plot(x=seq(0.05,0.95,by=0.05), y=Profit_mean_list$Profit.default, type="l", col="black",
#      xlab="Cost (at quantile of F)", ylab="True Profit",
#      main="True Profit vs. Cost (Gamma True F)")
# lines(x=seq(0.05,0.95,by=0.05), y=Profit_mean_list$Profit.init, type="l", col="blue")
# lines(x=seq(0.05,0.95,by=0.05), y=Profit_mean_list$Profit.mle, type="l", col="red")
# lines(x=seq(0.05,0.95,by=0.05), y=Profit_mean_list$Profit.pt, type="l", col="green")
# lines(x=seq(0.05,0.95,by=0.05), y=Profit_mean_list$Profit.gamma, type="l", col="purple")
# lines(x=seq(0.05,0.95,by=0.05), y=Profit_mean_list$Profit.tNorm, type="l", col="orange")
# 
# legend("topright", legend=c("Default","MLE","Initial","Polya Tree","Gamma-F","TN-F"),
#        col=c("black","red","blue","green","purple","orange"), lty=1)
# 
# ## Deviance comparison
# plot(x=seq(0.05,0.95,by=0.05), y=Profit_mean_list$est_Profit.init-Profit_mean_list$Profit.init, type="l", col="blue",
#      xlab="Cost (at quantile of F)", ylab="Estimated Profit - Profit", ylim=c(-0.5,0.12),
#      main="Estimated Profit - Profit vs. Cost (Gamma True F)")
# lines(x=seq(0.05,0.95,by=0.05), y=Profit_mean_list$est_Profit.mle-Profit_mean_list$Profit.mle, type="l", col="red")
# lines(x=seq(0.05,0.95,by=0.05), y=Profit_mean_list$est_Profit.pt-Profit_mean_list$Profit.pt, type="l", col="green")
# lines(x=seq(0.05,0.95,by=0.05), y=Profit_mean_list$est_Profit.gamma-Profit_mean_list$Profit.gamma, type="l", col="purple")
# lines(x=seq(0.05,0.95,by=0.05), y=Profit_mean_list$est_Profit.tNorm-Profit_mean_list$Profit.tNorm, type="l", col="orange")
# abline(h=0, lty=2)
# 
# legend("bottomright", legend=c("MLE","Initial","Polya Tree","Gamma-F","TN-F"),
#        col=c("red","blue","green","purple","orange"), lty=1)

#######################################################################################
## Plots
#######################################################################################

rm(list=ls())

# library(latex2exp)
library(forcats)



data_list_est_price=data.frame()

# QUANTILE=seq(0.05,0.95,by=0.05)

Namelist=list(
  Setting=rep(c("Uniform","Piecewise-Uniform","Pareto","Gamma","Beta"),times=3),
  Quantile=rep(c(0.1,0.3,0.5),each=5)
)

SimuList=list(
  Setting=rep(c("unif","piecewise_unif","pareto","gamma","beta"),times=3),
  Quantile_P=(Namelist$Quantile)/0.05
)
            
replication=100


for (i in 1:length(SimuList$Setting)) {
  load(paste0("Plots/Multiple_costs_trace_",SimuList$Setting[i],".RData"))
  data_bar_max_price=data.frame(
    Setting=rep(Namelist$Setting[i], replication*5),
    Quantile=rep(Namelist$Quantile[i], replication*5),
    Method=rep(c("Initial","MLE","PT","Gamma-F","TN-F"), each=replication),
    values=c(
      Profit_allList$max_price.init[SimuList$Quantile_P[i],],
      Profit_allList$max_price.mle[SimuList$Quantile_P[i],],
      Profit_allList$max_price.pt[SimuList$Quantile_P[i],],
      Profit_allList$max_price.gamma[SimuList$Quantile_P[i],],
      Profit_allList$max_price.tNorm[SimuList$Quantile_P[i],]
    )
  )
  
  
  data_bar_max_price=cbind(data_bar_max_price,Profit_mean_list$max_price.default[SimuList$Quantile_P[i]])
  colnames(data_bar_max_price)[ncol(data_bar_max_price)]="True_Max_Price"
  
  data_list_est_price=rbind(data_list_est_price,data_bar_max_price)
}

## re_level Method, quantile and Setting
data_list_est_price$Method=factor(data_list_est_price$Method,
                                  levels=c("Initial","MLE","PT","Gamma-F","TN-F"))

data_list_est_price$Setting=factor(data_list_est_price$Setting,
                                   levels=c("Uniform","Piecewise-Uniform","Pareto","Gamma","Beta"))

data_list_est_price$Quantile=paste0("c = ",data_list_est_price$Quantile,"-th quantile")

library(ggplot2)
## need to be modified for different metrics
ggplot_max_price=
  ggplot(data_list_est_price, aes(x=values,y=fct_rev(Method),fill=Method)) +
  facet_wrap(Quantile ~ Setting,ncol=5, scales = "free") +
  theme(strip.text = element_text(size = 8, margin = margin()))+
  geom_boxplot() +
  # scale_fill_manual(values = color6[-c(1,7,11)]) +
  # coord_cartesian(xlim = c(min(data_list[[1]]$rel_values),min(2,max(data_list[[1]]$rel_values,na.rm = T))))+
  geom_vline(aes(xintercept = True_Max_Price), data = data_list_est_price,linetype = "dashed") +
  xlab("Estimated max-profit price") +ylab("Method") +
  ggtitle(paste0("Estimated max-profit price of each method with K=100"))

pdf("Plots/multiple_Max_Profit_Price_1.pdf", width=12, height=4)
ggplot_max_price
dev.off()