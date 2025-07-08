#---------------------------------------------------------------------------------------------------
## Do trace plot from 1 repetitions of the training data
rm(list = ls())
source("Main.R")
source("distributions.R")
source("DataGeneration.R")
source("DisCompare.R")
source("PT_MCMC.R")



load(file = "Empirical/Xbox_7day_auctions.RData")

data <- Data_Gen_Observed(raw_data_list)
data_Polya <- Data_Gen_Unobserved_PolyaTree(raw_data_list)
## number of observed auctions only for PT
# test_N_all <- sapply(1:length(data$X_mk), function(ii) nrow(raw_data_list$Raw_biddata_list[ii][[1]]))


reserve.price.Cutoff <- ceiling(quantile(data$r_k, probs = 0.25))

true_X_mk=data$X_mk
simu_X_mk.init=c()
simu_X_mk.mle=c()
simu_X_mk.pt=c()


NRep=100

for (replication in 1:NRep) {
  set.seed(replication)
  
  Init.Values <- Initialization_lambda_F(data, reserve.price.Cutoff)
  theta.init <- (1-Init.Values$F.y)/ c(1,(1-Init.Values$F.y)[-length(Init.Values$F.y)])
  MLE <- MLE_2ndprice(data, lambda.in = Init.Values$lambda, theta.in = theta.init, tol = 1e-5)
  
  ## Polya Tree method
  
  ## Since number of bidders unknown, we have to generate them
  train_N_all <- rpois(data$Ka, MLE$lambda*data$tau_a)
  data_Polya[["dat"]][,1] <- train_N_all

  Bayes_PT <- MCMC_PT(data_Polya, data$Z_iT_i[,1], x.max=500)  
  
  ####################################################################################################
  
  e1.init <- AbscontDistribution(p = linear_interpolation_cdf(x = Init.Values$F.x, F.x = Init.Values$F.y) )
  e1.mle <- AbscontDistribution(p = linear_interpolation_cdf(x = MLE$F.x, F.x = MLE$F.y) )
  e1.pt <- AbscontDistribution(p = linear_interpolation_cdf(x = Bayes_PT$F.x, F.x = Bayes_PT$F.y) )
  
  

  
  ##########################################################################################################
  
  

  simu.init <- Data_Gen_Observed(Data_Gen_Raw(data$Ka, data$tau_a, Init.Values$lambda,
                                                   data$r_k, method = "customized", customized_CDF = e1.init))
  
  simu.mle <- Data_Gen_Observed(Data_Gen_Raw(data$Ka, data$tau_a, MLE$lambda,
                                                  data$r_k, method = "customized", customized_CDF = e1.mle))
  
  test_N_all <- rpois(data$Ka, MLE$lambda*data$tau_a)
  simu.pt <- sapply(1:data$Ka, function(jj){
    bid_ii=r(e1.pt)(test_N_all[jj])
    bid_ii_2nd <- sort(bid_ii, decreasing = T)[2] # 2nd highest bid
    if(max(bid_ii) < data$r_k[jj]) # if no bid is above reserve price
      return(NA)
    else if (bid_ii_2nd < data$r_k[jj]) # if 2nd highest bid is below reserve price
      return(data$r_k[jj])
    else
      return(bid_ii_2nd)
  })
  
  simu_X_mk.init=c(simu_X_mk.init, simu.init$X_mk)
  simu_X_mk.mle=c(simu_X_mk.mle, simu.mle$X_mk)
  simu_X_mk.pt=c(simu_X_mk.pt, simu.pt)
  
  print(paste("Replication", replication, "is done."))
  
}

##############################################################################################

simu_X_mk.init=matrix(simu_X_mk.init,ncol = NRep)
simu_X_mk.mle=matrix(simu_X_mk.mle,ncol = NRep)
simu_X_mk.pt=matrix(simu_X_mk.pt,ncol = NRep)

summary(true_X_mk)

min_X_mk_mle <- apply(simu_X_mk.mle, 2, function(x) min(x, na.rm = TRUE))
min_X_mk_init <- apply(simu_X_mk.init, 2, function(x) min(x, na.rm = TRUE))
min_X_mk_pt <- apply(simu_X_mk.pt, 2, function(x) min(x, na.rm = TRUE))

max_X_mk_mle <- apply(simu_X_mk.mle, 2, function(x) max(x, na.rm = TRUE))
max_X_mk_init <- apply(simu_X_mk.init, 2, function(x) max(x, na.rm = TRUE))
max_X_mk_pt <- apply(simu_X_mk.pt, 2, function(x) max(x, na.rm = TRUE))



save(list=ls(), file = paste0("Empirical/Xbox_simu_dist.RData"))

################################################################################################
min_true <- min(true_X_mk)
max_true <- max(true_X_mk)

cov_mle <- sapply(1:NRep, function(i){
  (min(max_X_mk_mle[i], max_true)- max(min_X_mk_mle[i], min_true) )/ (max_true - min_true)
})
cov_init <- sapply(1:NRep, function(i){
  (min(max_X_mk_init[i], max_true)- max(min_X_mk_init[i], min_true) )/ (max_true - min_true)
})
cov_pt <- sapply(1:NRep, function(i){
  (min(max_X_mk_pt[i], max_true)- max(min_X_mk_pt[i], min_true) )/ (max_true - min_true)
})

mean(cov_mle)
mean(cov_init)
mean(cov_pt)

cov_mle2 <- sapply(1:NRep, function(i){
  (max_X_mk_mle[i]- min_X_mk_mle[i])/ (max_true - min_true)
})
cov_init2 <- sapply(1:NRep, function(i){
  (max_X_mk_init[i]- min_X_mk_init[i])/ (max_true - min_true)
})
cov_pt2 <- sapply(1:NRep, function(i){
  (max_X_mk_pt[i]- min_X_mk_pt[i])/ (max_true - min_true)
})

mean(cov_mle2)
mean(cov_init2)
mean(cov_pt2)
# 
# plot1_df <- data.frame(xx=c(cov_mle, cov_init, cov_pt),method=c(rep("mle", NRep), rep("init", NRep), rep("pt", NRep)))
# library(ggplot2)
# ggplot(plot1_df, aes(x=xx,fill=method)) +
#     theme(legend.position='none',text = element_text(size = 10)) +
#     geom_histogram(alpha=0.5,
#                    color="white",
#                    position="identity")+
#   scale_fill_manual(values=c("mle"="red", "init"="blue", "pt"="green")) 
# 
# #####################################################################################################
# mean_X_mk_mle <- apply(simu_X_mk.mle, 2, function(x) mean(x, na.rm = TRUE))
# mean_X_mk_init <- apply(simu_X_mk.init, 2, function(x) mean(x, na.rm = TRUE))
# mean_X_mk_pt <- apply(simu_X_mk.pt, 2, function(x) mean(x, na.rm = TRUE))
# 
# med_X_mk_mle <- apply(simu_X_mk.mle, 2, function(x) median(x, na.rm = TRUE))
# med_X_mk_init <- apply(simu_X_mk.init, 2, function(x) median(x, na.rm = TRUE))
# med_X_mk_pt <- apply(simu_X_mk.pt, 2, function(x) median(x, na.rm = TRUE))
# 
# 
# 
# plot1_df <- data.frame(xx=min_X_mk_mle,xt=min(true_X_mk),method="mle",type="min")
# plot1_df <- rbind(plot1_df,
#                   data.frame(xx=min_X_mk_init,xt=min(true_X_mk),method="init",type="min"))
# plot1_df <- rbind(plot1_df,
#                   data.frame(xx=min_X_mk_pt,xt=min(true_X_mk),method="pt",type="min"))
# 
# plot1_df<- rbind(plot1_df,
#                   data.frame(xx=max_X_mk_mle,xt=max(true_X_mk),method="mle",type="max"))
# plot1_df <- rbind(plot1_df,
#                   data.frame(xx=max_X_mk_init,xt=max(true_X_mk),method="init",type="max"))
# plot1_df <- rbind(plot1_df,
#                   data.frame(xx=max_X_mk_pt,xt=max(true_X_mk),method="pt",type="max"))
# 
# plot1_df <- rbind(plot1_df,
#                   data.frame(xx=mean_X_mk_mle,xt=mean(true_X_mk),method="mle",type="mean"))
# plot1_df <- rbind(plot1_df,
#                   data.frame(xx=mean_X_mk_init,xt=mean(true_X_mk),method="init",type="mean"))
# plot1_df <- rbind(plot1_df,
#                   data.frame(xx=mean_X_mk_pt,xt=mean(true_X_mk),method="pt",type="mean"))
# 
# plot1_df <- rbind(plot1_df,
#                   data.frame(xx=med_X_mk_mle,xt=median(true_X_mk),method="mle",type="med"))
# plot1_df <- rbind(plot1_df,
#                   data.frame(xx=med_X_mk_init,xt=median(true_X_mk),method="init",type="med"))
# plot1_df <- rbind(plot1_df,
#                   data.frame(xx=med_X_mk_pt,xt=median(true_X_mk),method="pt",type="med"))
# 
# get_vline_x <- function(x) {
#   switch(x,
#          "min" = min(true_X_mk),
#          "max" = max(true_X_mk),
#          "mean" = mean(true_X_mk),
#          "med" = median(true_X_mk)
#   )}
# plot1_df$vline_x <- sapply(plot1_df$type, get_vline_x)
# 
# library(ggplot2)
# ggplot(plot1_df, aes(x=xx,fill=method)) +
#   theme(legend.position='none',text = element_text(size = 10)) +
#   geom_histogram(alpha=0.5,
#                  binwidth = 10,
#                  color="white",
#                  position="identity") +
#   facet_grid( method~type,scales="free",labeller = label_parsed)+
#   theme(plot.title = element_text(hjust = 10))+
#   labs(x="value",y="counts") +
#   geom_vline(aes(xintercept = vline_x), color = "black", linetype = "dashed")