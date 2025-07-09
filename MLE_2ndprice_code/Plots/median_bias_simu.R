# rm(list = ls())   # Remove everything from the Global environment
# Ka=20 ## in batches there are fewer auctions
# method.in <- "unif"
# para1 <- 1
# para2 <- 20
# pdf(paste0("Plots/Median_bias_", method.in, " With No. Auction = ", Ka, ".pdf"))
# source("median_bias_simu.R")
# dev.off()

library(tidyverse)
library(GoFKernel) # For the inverse function of a CDF
library(ggplot2)   # PLOT ECDF
library(reshape2)

source("Main.R")
source("distributions.R")
source("DataGeneration.R")
source("DisCompare.R")
source("SettingGeneration.R")
source("PT_MCMC.R")

#------------------------------------------------------------------------------------
## User Input: True lambda, Auction window, Number of independent auctions 'N.auction'
#------------------------------------------------------------------------------------

replication <- 500  # (SOURAV) replication numbers of each such event i.e. selling an
# identical stuff on ebay with multiple independent auctions.

lambda0 <- 1
tau_a <- 100
# Ka=100 ## change with sample size!!

## settings copied from Simulations_with_HulC_Confidence_region.R

SETTING_ALL <- Setting_Gen_Med1(Ka, method = method.in, para1, para2)
y.seq <- SETTING_ALL$y.seq
F.0 <- SETTING_ALL$F.0
r_k <- SETTING_ALL$r_k
reserve.price.Cutoff <- SETTING_ALL$reserve.price.Cutoff


Init_fun_Mat <- matrix(0, nrow = replication, ncol = length(y.seq))
Mle_fun_Mat <- matrix(0, nrow = replication, ncol = length(y.seq))
PT_fun_Mat <- matrix(0, nrow = replication, ncol = length(y.seq))

for(ii in 1:replication){  
  if (ii%%10 ==1) print(paste0("At replication ", ii, " out of ", replication))
  #------------------------------------------------------------------------------------
  ## Data generation, Initial Estimate, MLE, and HUlC
  #------------------------------------------------------------------------------------
  
  Raw_data_all_list <- Data_Gen_Raw(Ka, tau_a, lambda0,
                                    r_k, method = method.in, para1, para2)
  data <- Data_Gen_Observed(Raw_data_all_list)
  
  Init.Values <- Initialization_lambda_F(data, reserve.price.Cutoff)
  theta.init <- (1-Init.Values$F.y)/ c(1,(1-Init.Values$F.y)[-length(Init.Values$F.y)])
  MLE <- MLE_2ndprice(data, lambda.in = Init.Values$lambda, theta.in = theta.init)
  
  ## Polya Tree method
  data_Polya <- Data_Gen_Unobserved_PolyaTree(Raw_data_all_list)
  Bayes_PT <- MCMC_PT(data_Polya, data$Z_iT_i[,1], x.max=max(SETTING_ALL$y.seq))
  
  
  Init_fun <- approxfun(c(0,Init.Values$F.x), y = c(0, Init.Values$F.y), method = "linear", 0, 0.9999999, ties = mean)
  
  Mle_fun <- approxfun(c(0,MLE$F.x), y = c(0, MLE$F.y), method = "linear", 0, 0.9999999, ties = mean)
  
  # PT_fun <- approxfun(c(0, Bayes_PT$F.x), y = c(0, Bayes_PT$F.y), method = "linear", 0, 0.9999999, ties = mean)
  PT_fun <- Bayes_PT$fit_CDF
  
  
  
  # evaluate Initfun and Mlefun at y.seq
  Init_fun_Mat[ii, ] <- Init_fun(y.seq)
  Mle_fun_Mat[ii, ] <- Mle_fun(y.seq)
  PT_fun_Mat[ii, ] <- PT_fun(y.seq)
}
# Init_fun_Mat <- Init_fun_Mat[1:replication,]
# Mle_fun_Mat <- Mle_fun_Mat[1:replication,]
Init_fun_Mat <- t(apply(Init_fun_Mat, 1, function(x)   x-F.0(y.seq)))
Mle_fun_Mat <- t(apply(Mle_fun_Mat, 1, function(x)   x-F.0(y.seq)))
PT_fun_Mat <- t(apply(PT_fun_Mat, 1, function(x)   x-F.0(y.seq)))


med.Init_fun <- apply(Init_fun_Mat, 2, median)
med.Mle_fun <- apply(Mle_fun_Mat, 2, median)
med.PT_fun <- apply(PT_fun_Mat, 2, median)
# Median.bias.MLE <-max(abs(med.Mle_fun - F.0(y.seq)))
# Median.bias.Init <- max(abs(med.Init_fun -F.0(y.seq)))
# Median.bias.PT <- max(abs(med.PT_fun -F.0(y.seq)))

# plot(y.seq,med.Init_fun-F.0(y.seq))
# points(y.seq,med.Mle_fun-F.0(y.seq), col= "red")
# add legend
# plot_data <- data.frame(y.seq, med.Init_fun, med.Mle_fun, med.PT_fun)
# colnames(plot_data) <- c("y", "Initial_estimate", "MLE", "Bayes_PT")
# # melt with id var y
# plot_data_m <- melt(plot_data, id.vars = "y")
# str(plot_data_m)

# pdf(paste0("Plots/Median_bias_", method.in, " With No. Auction = ", Ka, ".pdf"))
# ggplot(plot_data_m) + 
#   geom_line(aes(x = y, y = value, col = variable)) + ylab("Median bias") + ggtitle(paste0("Median bias for MLE, Polya Tree and Init function, for Distribution = ", method.in))
# dev.off()

# (SOURAV) Make the confidence plots wider by the Median.bias.MLE for the MLE and by Median.bias.Init for the Init function. And for each plot, in the caption of the plot add the value of the median bias. 

# Let us know what code you used to compute the confidence bands. That code had an already existing bias correction. We need to remove that and add the Median.bias.MLE as the correction factor. 
# Median_bias_gamma



med.bias.init <-  apply(Init_fun_Mat, 2, function(x)   mean(x>0))
med.bias.mle <-  apply(Mle_fun_Mat, 2, function(x)   mean(x>0))
med.bias.pt <-  apply(PT_fun_Mat, 2, function(x)   mean(x>0))

med_data <- data.frame(y.seq, med.Init_fun, med.Mle_fun, med.PT_fun,
                       med.bias.init, med.bias.mle, med.bias.pt)

med_data=cbind(med_data, med.bias.init, med.bias.mle, med.bias.pt)
colnames(med_data) <- c("y", "diff_med_init", "diff_med_mle", "diff_med_pt",
                         "med_bias_init", "med_bias_mle", "med_bias_pt")

save(med_data, file = paste0("Plots/MdBs/Median_bias_", method.in, " Ka = ", Ka, ".RData"))

#########################################################################################
# plot(med_data$y, abs(med_data$med_bias_init-1/2), type = "l", col = "red",
#      xlab = "y", ylab = "Median bias", main = paste0("Median bias for MLE(G), Polya Tree(B) and Init(R) function, for Distribution = ", method.in))
# lines(med_data$y, abs(med_data$med_bias_mle-1/2), col = "green")
# lines(med_data$y, abs(med_data$med_bias_pt-1/2), col = "blue")

# #########################################################################################
# plot_data <- med_data[,1:4]
# plot_data[,1] <- F.0(y.seq)
# colnames(plot_data) <- c("F_0_y", "Initial_estimate", "MLE", "Bayes_PT")
# # melt with id var y
# plot_data_m <- melt(plot_data, id.vars = "F_0_y")
# # str(plot_data_m)
# 
# 
# myplot_med_diff<- ggplot(plot_data_m) + 
#   geom_line(aes(x = F_0_y, y = value, col = variable)) + ylab("Median bias") + ggtitle(paste0("Difference in Med for MLE, Polya Tree and Init function, for Distribution = ", method.in))
# 
# plot_data2 <- med_data[,c(1,5:7)]
# plot_data2[,1] <- F.0(y.seq)
# plot_data2[,2:4] <- abs(plot_data2[,2:4] - 1/2)
# colnames(plot_data2) <- c("F_0_y", "Initial_estimate", "MLE", "Bayes_PT")
# # melt with id var y
# plot_data_m2 <- melt(plot_data2, id.vars = "F_0_y")
# # str(plot_data_m2)
# myplot_med_bias<- ggplot(plot_data_m2) + 
#   geom_line(aes(x = F_0_y, y = value, col = variable)) + ylab("Median bias") + ggtitle(paste0("Median bias for MLE, Polya Tree and Init function, for Distribution = ", method.in))

#########################################################################################
plot_data <- med_data[,1:4]
colnames(plot_data) <- c("y.seq", "Initial_estimate", "MLE", "Bayes_PT")
# melt with id var y
plot_data_m <- melt(plot_data, id.vars = "y.seq")
# str(plot_data_m)


myplot_med_diff<- ggplot(plot_data_m) + 
  geom_line(aes(x = y.seq, y = value, col = variable)) + ylab("Median bias") + ggtitle(paste0("Difference in Med for MLE, Polya Tree and Init function, for Distribution = ", method.in))

plot_data2 <- med_data[,c(1,5:7)]
plot_data2[,2:4] <- abs(plot_data2[,2:4] - 1/2)
colnames(plot_data2) <- c("y.seq", "Initial_estimate", "MLE", "Bayes_PT")
# melt with id var y
plot_data_m2 <- melt(plot_data2, id.vars = "y.seq")
# str(plot_data_m2)
myplot_med_bias<- ggplot(plot_data_m2) + 
  geom_line(aes(x = y.seq, y = value, col = variable)) + ylab("Median bias") + ggtitle(paste0("Median bias for MLE, Polya Tree and Init function, for Distribution = ", method.in))

