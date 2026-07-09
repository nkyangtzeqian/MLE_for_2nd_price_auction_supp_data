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
#####################################################################################

rm(list = ls())   # Remove everything from the Global environment
Ka=100
method.in <- "unif"
para1 <- 1
para2 <- 20
 LAMBDA0_SETTING_INDEX <- 1
source("Simu1L.R")

rm(list = ls())   # Remove everything from the Global environment
Ka=100
method.in <- "piecewise_unif"
para1 <- 2
para2 <- 4
 LAMBDA0_SETTING_INDEX <- 1
source("Simu1L.R")

rm(list = ls())   # Remove everything from the Global environment
Ka=100
method.in <- "pareto"
para1 <- 3
para2 <- 100
 LAMBDA0_SETTING_INDEX <- 1
source("Simu1L.R")

rm(list = ls())   # Remove everything from the Global environment
Ka=100
method.in <- "gamma"
para1 <- 10
para2 <- 2
 LAMBDA0_SETTING_INDEX <- 1
source("Simu1L.R")

rm(list = ls())   # Remove everything from the Global environment
Ka=100
method.in <- "beta"
para1 <- 2
para2 <- 2
LAMBDA0_SETTING_INDEX <- 1
source("Simu1L.R")



rm(list = ls())   # Remove everything from the Global environment
Ka=100
method.in <- "unif"
para1 <- 1
para2 <- 20
LAMBDA0_SETTING_INDEX <- 2
source("Simu1L.R")

rm(list = ls())   # Remove everything from the Global environment
Ka=100
method.in <- "piecewise_unif"
para1 <- 2
para2 <- 4
LAMBDA0_SETTING_INDEX <- 2
source("Simu1L.R")

rm(list = ls())   # Remove everything from the Global environment
Ka=100
method.in <- "pareto"
para1 <- 3
para2 <- 100
LAMBDA0_SETTING_INDEX <- 2
source("Simu1L.R")

rm(list = ls())   # Remove everything from the Global environment
Ka=100
method.in <- "gamma"
para1 <- 10
para2 <- 2
LAMBDA0_SETTING_INDEX <- 2
source("Simu1L.R")

rm(list = ls())   # Remove everything from the Global environment
Ka=100
method.in <- "beta"
para1 <- 2
para2 <- 2
LAMBDA0_SETTING_INDEX <- 2
source("Simu1L.R")



rm(list = ls())   # Remove everything from the Global environment
Ka=100
method.in <- "unif"
para1 <- 1
para2 <- 20
LAMBDA0_SETTING_INDEX <- 3
source("Simu1L.R")

rm(list = ls())   # Remove everything from the Global environment
Ka=100
method.in <- "piecewise_unif"
para1 <- 2
para2 <- 4
LAMBDA0_SETTING_INDEX <- 3
source("Simu1L.R")

rm(list = ls())   # Remove everything from the Global environment
Ka=100
method.in <- "pareto"
para1 <- 3
para2 <- 100
LAMBDA0_SETTING_INDEX <- 3
source("Simu1L.R")

rm(list = ls())   # Remove everything from the Global environment
Ka=100
method.in <- "gamma"
para1 <- 10
para2 <- 2
LAMBDA0_SETTING_INDEX <- 3
source("Simu1L.R")

rm(list = ls())   # Remove everything from the Global environment
Ka=100
method.in <- "beta"
para1 <- 2
para2 <- 2
LAMBDA0_SETTING_INDEX <- 3
source("Simu1L.R")




#####################################################################################
SimuList=list(
  M=rep(c("unif","piecewise_unif","pareto","gamma","beta"),each=3),
  Setting=rep(c("Uniform","P_Unif","Pareto","Gamma","Beta"),each=3),
  l_index=rep(c(1,3,2),times=5),
  lambda=rep(c("Con","Lin","TS"),times=5)
)

Table1=data.frame(Setting=SimuList$Setting,
                  lambda=SimuList$lambda,
                  KS.mle1.mean=NA,
                  KS.mle3.mean=NA,
                  KS.mle2.mean=NA,
                  KS.init.mean=NA,
                  KS.pt.mean=NA,
                  KS.gamma.mean=NA,
                  KS.tNorm.mean=NA,
                  TV.mle1.mean=NA,
                  TV.mle3.mean=NA,
                  TV.mle2.mean=NA,
                  TV.init.mean=NA,
                  TV.pt.mean=NA,
                  TV.gamma.mean=NA,
                  TV.tNorm.mean=NA
)

mean_T=function(x){ ## TN-F always fail
  mean(x[x<1],na.rm = TRUE)
}

for (i in 1:length(SimuList$Setting)) {
  l_index <- SimuList$l_index[i]
  load(paste0("SIMU/Table1L_", l_index,"_", SimuList$M[i], "_K=", 100, ".RData"))
  Table1$KS.mle1.mean[i] <- mean_T(distanceList$KS.mle1)
  Table1$KS.mle2.mean[i] <- mean_T(distanceList$KS.mle2)
  Table1$KS.mle3.mean[i] <- mean_T(distanceList$KS.mle3)
  Table1$KS.init.mean[i] <- mean_T(distanceList[[paste0("KS.init", l_index)]])
  Table1$KS.pt.mean[i] <- mean_T(distanceList$KS.pt)
  Table1$KS.gamma.mean[i] <- mean_T(distanceList$KS.gamma)
  Table1$KS.tNorm.mean[i] <- mean_T(distanceList$KS.tNorm)
  Table1$TV.mle1.mean[i] <- mean_T(distanceList$TV.mle1)
  Table1$TV.mle2.mean[i] <- mean_T(distanceList$TV.mle2)
  Table1$TV.mle3.mean[i] <- mean_T(distanceList$TV.mle3)
  Table1$TV.init.mean[i] <- mean_T(distanceList[[paste0("TV.init", l_index)]])
  Table1$TV.pt.mean[i] <- mean_T(distanceList$TV.pt)
  Table1$TV.gamma.mean[i] <- mean_T(distanceList$TV.gamma)
  Table1$TV.tNorm.mean[i] <- mean_T(distanceList$TV.tNorm)
}

Table1
# cat(print_tex1(Table1[,c(1:2,3:9)]))
# cat(print_tex1(Table1[,c(1:2,10:16)]))

############################################################################################
## Profit 
############################################################################################

rm(list=ls())

# library(latex2exp)
library(forcats)



Profit_default <- function(method = c("unif", "pareto", "gamma", "beta",
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
  
  source("distributionsL.R")
  source("SettingGeneration.R")
  para1=Setting_Gen_Tab1(method,K=100,para1=NULL, para2=NULL)$para1
  para2=Setting_Gen_Tab1(method,K=100,para1=NULL, para2=NULL)$para2
  
  source("ProfitMetric.R")
  max_cutpt <- cost_def(method, para1, para2)["max_cutpt"]
  cost <- cost_def(method, para1, para2)["cost"]
  
  e_cdf <-Setting_Gen_Tab1_F0_cdf(method,para1, para2)[[1]]
  
  
  Profit_all <- Compute.Opt.Price(e_cdf, cost, max.price=max_cutpt)
  max_price <- Profit_all$opt.price
  Profit <- Profit_all$profit
  
  return(list(max_price = max_price,
              Profit = Profit))
}


SimuList=list(
  M=rep(c("unif","piecewise_unif","pareto","gamma","beta"),each=3),
  Setting=rep(c("Uniform","P_Unif","Pareto","Gamma","Beta"),each=3),
  l_index=rep(c(1,3,2),times=5),
  lambda=rep(c("Con","Lin","TS"),times=5)
)

data_list_est_price=data_list_profit=data.frame()


replication=100


for (i in 1:length(SimuList$Setting)) {
  l_index <- SimuList$l_index[i]
  load(paste0("SIMU/Table1L_", l_index,"_ProfitRes_", SimuList$M[i], "_K=", 100, ".RData"))
  data_bar_max_price=data.frame(
    Setting=rep(SimuList$Setting[i], replication*7),
    lambda=rep(SimuList$lambda[i], replication*7),
    Method=rep(c("MLE1","MLE2","MLE3","Initial","PT","Gamma-F","TN-F"), each=replication),
    values=c(
      Profit_allList$max_price.mle1,
      Profit_allList$max_price.mle2,
      Profit_allList$max_price.mle3,
      Profit_allList[[paste0("max_price.init", l_index)]],
      Profit_allList$max_price.pt,
      Profit_allList$max_price.gamma,
      Profit_allList$max_price.tNorm
    )
  )
  data_bar_Profit=data.frame(
    Setting=rep(SimuList$Setting[i], replication*7),
    lambda=rep(SimuList$lambda[i], replication*7),
    Method=rep(c("MLE1","MLE2","MLE3","Initial","PT","Gamma-F","TN-F"), each=replication),
    values=c(
      Profit_allList$Profit.mle1,
      Profit_allList$Profit.mle2,
      Profit_allList$Profit.mle3,
      Profit_allList[[paste0("Profit.init", l_index)]],
      Profit_allList$Profit.pt,
      Profit_allList$Profit.gamma,
      Profit_allList$Profit.tNorm
    )
  )
  data_bar_est_Profit=data.frame(
    Setting=rep(SimuList$Setting[i], replication*7),
    lambda=rep(SimuList$lambda[i], replication*7),
    Method=rep(c("MLE1(est)","MLE2(est)","MLE3(est)","Initial(est)","PT(est)","Gamma-F(est)","TN-F(est)"), each=replication),
    values=c(
      Profit_allList$est_Profit.mle1,
      Profit_allList$est_Profit.mle2,
      Profit_allList$est_Profit.mle3,
      Profit_allList[[paste0("est_Profit.init", l_index)]],
      Profit_allList$est_Profit.pt,
      Profit_allList$est_Profit.gamma,
      Profit_allList$est_Profit.tNorm
    )
  )
  
  Profit_def=Profit_default(method = SimuList$M[i])
  
  data_bar_max_price=cbind(data_bar_max_price,Profit_def$max_price)
  colnames(data_bar_max_price)[ncol(data_bar_max_price)]="True_Max_Price"
  data_bar_Profit=cbind(data_bar_Profit,Profit_def$Profit)
  data_bar_est_Profit=cbind(data_bar_est_Profit,Profit_def$Profit)
  colnames(data_bar_Profit)[ncol(data_bar_Profit)]="True_Max_Profit"
  colnames(data_bar_est_Profit)[ncol(data_bar_est_Profit)]="True_Max_Profit"
  
  data_list_est_price=rbind(data_list_est_price,data_bar_max_price)
  data_list_profit=rbind(data_list_profit,data_bar_Profit,data_bar_est_Profit)
}

## re_level Method and Setting
data_list_est_price$Method=factor(data_list_est_price$Method,
                                  levels=c("MLE1","MLE3","MLE2","Initial","PT","Gamma-F","TN-F"))
data_list_profit$Method=factor(data_list_profit$Method,
                               levels=c("MLE1","MLE1(est)","MLE3","MLE3(est)","MLE2","MLE2(est)","Initial","Initial(est)","PT","PT(est)","Gamma-F","Gamma-F(est)","TN-F","TN-F(est)"))

data_list_est_price$Setting=factor(data_list_est_price$Setting,
                                   levels=c("Uniform","P_Unif","Pareto","Gamma","Beta"))
data_list_profit$Setting=factor(data_list_profit$Setting,
                                levels=c("Uniform","P_Unif","Pareto","Gamma","Beta"))

# library(ggplot2)
# ## need to be modified for different metrics
# ggplot_max_price=
#   ggplot(data_list_est_price, aes(x=values,y=fct_rev(Method),fill=Method)) +
#   facet_wrap(Ka ~ Setting,ncol=5, scales = "free") +
#   theme(strip.text = element_text(size = 8, margin = margin()))+
#   geom_boxplot() +
#   # scale_fill_manual(values = color6[-c(1,7,11)]) +
#   # coord_cartesian(xlim = c(min(data_list[[1]]$rel_values),min(2,max(data_list[[1]]$rel_values,na.rm = T))))+
#   geom_vline(aes(xintercept = True_Max_Price), data = data_list_est_price,linetype = "dashed") +
#   xlab("Estimated max-profit price") +ylab("Method") +
#   ggtitle(paste0("Estimated max-profit price of each method"))
# 
# ggplot_Profit=
#   ggplot(data_list_profit, aes(x=values,y=fct_rev(Method),fill=Method)) +
#   facet_wrap(Ka ~ Setting,ncol=5, scales = "free") +
#   theme(strip.text = element_text(size = 8, margin = margin()))+
#   geom_boxplot() +
#   # scale_fill_manual(values = color6[-c(1,7,11)]) +
#   # coord_cartesian(xlim = c(min(data_list[[1]]$rel_values),min(2,max(data_list[[1]]$rel_values,na.rm = T))))+
#   geom_vline(aes(xintercept = True_Max_Profit), data = data_list_profit,linetype = "dashed") +
#   xlab("Profit from estimated max-profit price") +ylab("Method") +
#   ggtitle(paste0("Profit and Estimated Profit from max-profit price of each method"))
# 
# 
# pdf("Plots/Estimated_Max_Profit_Price_1.pdf", width=16, height=4)
# ggplot_max_price
# dev.off()
# 
# pdf("Plots/Profit_1.pdf", width=16, height=6)
# ggplot_Profit
# dev.off()

###############################################################################################
## Table profit

table_data=matrix(0,nrow=length(SimuList$Setting),ncol=nlevels(data_list_profit$Method)+1)

for (i in 1:length(SimuList$Setting)) {
  for(j in 1:nlevels(data_list_profit$Method)){
    table_data[i,j]=mean(data_list_profit$values[
      data_list_profit$Setting==SimuList$Setting[i] &
        data_list_profit$lambda==SimuList$lambda[i] &
        data_list_profit$Method==levels(data_list_profit$Method)[j]
    ],na.rm=T)
  }
  table_data[i,nlevels(data_list_profit$Method)+1]=data_list_profit$True_Max_Profit[
    min(which(data_list_profit$Setting==SimuList$Setting[i] &
                data_list_profit$lambda==SimuList$lambda[i]))
  ]
}

rownames(table_data)=paste0(SimuList$Setting,"_",SimuList$lambda)
colnames(table_data)=c(levels(data_list_profit$Method),"True_Profit")

Table2=cbind.data.frame(Setting=SimuList$Setting,lambda=SimuList$lambda,
             table_data[,c(2,4,6,8,10,12,14)]-table_data[,c(1,3,5,7,9,11,13)])
# 
# cat(print2_tex(Table2))
