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
  
  source("distributions.R")
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
  Setting=rep(c("unif","piecewise_unif","pareto","gamma","beta"),times=2),
  Ka=rep(c(100,1000),each=5)
)


Namelist=list(
  Setting=rep(c("Uniform","Piecewise-Uniform","Pareto","Gamma","Beta"),times=2),
  Ka=rep(c("K=100","K=1000"),each=5)
)

data_list_est_price=data_list_profit=data.frame()


replication=100


for (i in 1:length(SimuList$Setting)) {
  load(paste0("SIMU/Table1A_ProfitRes_", SimuList$Setting[i], "_K=", SimuList$Ka[i], ".RData"))
  data_bar_max_price=data.frame(
    Setting=rep(Namelist$Setting[i], replication*5),
    Ka=rep(Namelist$Ka[i], replication*5),
    Method=rep(c("Initial","MLE","PT","Gamma-F","TN-F"), each=replication),
    values=c(
      Profit_allList$max_price.init,
      Profit_allList$max_price.mle,
      Profit_allList$max_price.pt,
      Profit_allList$max_price.gamma,
      Profit_allList$max_price.tNorm
    )
  )
  data_bar_Profit=data.frame(
    Setting=rep(Namelist$Setting[i], replication*5),
    Ka=rep(Namelist$Ka[i], replication*5),
    Method=rep(c("Initial","MLE","PT","Gamma-F","TN-F"), each=replication),
    values=c(
      Profit_allList$Profit.init,
      Profit_allList$Profit.mle,
      Profit_allList$Profit.pt,
      Profit_allList$Profit.gamma,
      Profit_allList$Profit.tNorm
    )
  )
  data_bar_est_Profit=data.frame(
    Setting=rep(Namelist$Setting[i], replication*5),
    Ka=rep(Namelist$Ka[i], replication*5),
    Method=rep(c("Initial(est)","MLE(est)","PT(est)","Gamma-F(est)","TN-F(est)"), each=replication),
    values=c(
      Profit_allList$est_Profit.init,
      Profit_allList$est_Profit.mle,
      Profit_allList$est_Profit.pt,
      Profit_allList$est_Profit.gamma,
      Profit_allList$est_Profit.tNorm
    )
  )
  
  Profit_def=Profit_default(method = SimuList$Setting[i])
  
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
                                  levels=c("Initial","MLE","PT","Gamma-F","TN-F"))
data_list_profit$Method=factor(data_list_profit$Method,
                              levels=c("Initial","Initial(est)","MLE","MLE(est)","PT","PT(est)","Gamma-F","Gamma-F(est)","TN-F","TN-F(est)"))

data_list_est_price$Setting=factor(data_list_est_price$Setting,
                                   levels=c("Uniform","Piecewise-Uniform","Pareto","Gamma","Beta"))
data_list_profit$Setting=factor(data_list_profit$Setting,
                               levels=c("Uniform","Piecewise-Uniform","Pareto","Gamma","Beta"))

library(ggplot2)
## need to be modified for different metrics
ggplot_max_price=
  ggplot(data_list_est_price, aes(x=values,y=fct_rev(Method),fill=Method)) +
  facet_wrap(Ka ~ Setting,ncol=5, scales = "free") +
  theme(strip.text = element_text(size = 8, margin = margin()))+
  geom_boxplot() +
  # scale_fill_manual(values = color6[-c(1,7,11)]) +
  # coord_cartesian(xlim = c(min(data_list[[1]]$rel_values),min(2,max(data_list[[1]]$rel_values,na.rm = T))))+
  geom_vline(aes(xintercept = True_Max_Price), data = data_list_est_price,linetype = "dashed") +
  xlab("Estimated max-profit price") +ylab("Method") +
  ggtitle(paste0("Estimated max-profit price of each method"))

ggplot_Profit=
  ggplot(data_list_profit, aes(x=values,y=fct_rev(Method),fill=Method)) +
  facet_wrap(Ka ~ Setting,ncol=5, scales = "free") +
  theme(strip.text = element_text(size = 8, margin = margin()))+
  geom_boxplot() +
  # scale_fill_manual(values = color6[-c(1,7,11)]) +
  # coord_cartesian(xlim = c(min(data_list[[1]]$rel_values),min(2,max(data_list[[1]]$rel_values,na.rm = T))))+
  geom_vline(aes(xintercept = True_Max_Profit), data = data_list_profit,linetype = "dashed") +
  xlab("Profit from estimated max-profit price") +ylab("Method") +
  ggtitle(paste0("Profit and Estimated Profit from max-profit price of each method"))


pdf("Plots/Estimated_Max_Profit_Price_1.pdf", width=16, height=4)
ggplot_max_price
dev.off()

pdf("Plots/Profit_1.pdf", width=16, height=6)
ggplot_Profit
dev.off()

###############################################################################################
## Table profit

table_data=matrix(0,nrow=length(SimuList$Setting),ncol=nlevels(data_list_profit$Method)+1)

for (i in 1:length(SimuList$Setting)) {
  for(j in 1:nlevels(data_list_profit$Method)){
    table_data[i,j]=mean(data_list_profit$values[
        data_list_profit$Setting==Namelist$Setting[i] &
        data_list_profit$Ka==Namelist$Ka[i] &
        data_list_profit$Method==levels(data_list_profit$Method)[j]
    ],na.rm=T)
  }
  table_data[i,nlevels(data_list_profit$Method)+1]=data_list_profit$True_Max_Profit[
    min(which(data_list_profit$Setting==Namelist$Setting[i] &
              data_list_profit$Ka==Namelist$Ka[i]))
  ]
}

rownames(table_data)=paste0(Namelist$Setting,"_",Namelist$Ka)
colnames(table_data)=c(levels(data_list_profit$Method),"True_Profit")

table_data