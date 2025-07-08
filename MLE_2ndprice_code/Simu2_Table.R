library(tidyverse)
library(GoFKernel) # For the inverse function of a CDF


source("Main.R")
source("distributions.R")
source("DataGeneration.R")
# source("DisCompare.R")
source("SettingGeneration.R")
source("PT_MCMC.R")
source("HulC_main.R")

rm(list = ls())   # Remove everything from the Global environment
Ka=100
NSimu=100
Delta=0
method.in <- "unif"
para1 <- 1
para2 <- 20
source("Simu2_HulC_Coverage.R")

rm(list = ls())   # Remove everything from the Global environment
Ka=100
NSimu=100
Delta=0
method.in <- "piecewise_unif"
para1 <- 2
para2 <- 4
source("Simu2_HulC_Coverage.R")

rm(list = ls())   # Remove everything from the Global environment
Ka=100
NSimu=100
Delta=0
method.in <- "pareto"
para1 <- 3
para2 <- 100
source("Simu2_HulC_Coverage.R")

rm(list = ls())   # Remove everything from the Global environment
Ka=100
NSimu=100
Delta=0
method.in <- "gamma"
para1 <- 10
para2 <- 2
source("Simu2_HulC_Coverage.R")

rm(list = ls())   # Remove everything from the Global environment
Ka=100
NSimu=100
Delta=0
method.in <- "beta"
para1 <- 2
para2 <- 2
source("Simu2_HulC_Coverage.R")



rm(list = ls())   # Remove everything from the Global environment
Ka=1000
NSimu=100
Delta=0
method.in <- "unif"
para1 <- 1
para2 <- 20
source("Simu2_HulC_Coverage.R")

rm(list = ls())   # Remove everything from the Global environment
Ka=1000
NSimu=100
Delta=0
method.in <- "piecewise_unif"
para1 <- 2
para2 <- 4
source("Simu2_HulC_Coverage.R")

rm(list = ls())   # Remove everything from the Global environment
Ka=1000
NSimu=100
Delta=0
method.in <- "pareto"
para1 <- 3
para2 <- 100
source("Simu2_HulC_Coverage.R")

rm(list = ls())   # Remove everything from the Global environment
Ka=1000
NSimu=100
Delta=0
method.in <- "gamma"
para1 <- 10
para2 <- 2
source("Simu2_HulC_Coverage.R")

rm(list = ls())   # Remove everything from the Global environment
Ka=1000
NSimu=100
Delta=0
method.in <- "beta"
para1 <- 2
para2 <- 2
source("Simu2_HulC_Coverage.R")


#####################################################################################
SimuList=list(
  Method=rep(c("unif","piecewise_unif","pareto","gamma","beta"),each=2),
  Ka=rep(c(100,1000),times=5)
)

Table2=data.frame(Method=SimuList$Method,
                  Ka=SimuList$Ka,
                  Cov_MLE_HulC=NA,
                  Cov_INIT_HulC=NA,
                  Cov_PT_Cred=NA,
                  Area_MLE_HulC=NA,
                  Area_INIT_HulC=NA,
                  Area_PT_Cred=NA)

function_coverage <- function(Mat_Cov) {
  fun_cov_rep <- apply(Mat_Cov, 1, function(x) {all(x)})
  return(mean(fun_cov_rep,na.rm = TRUE))
}
function_area <- function(Mat_Cov,x.med.vec) {
  diff_x_vec <- diff(c(0,x.med.vec))
  fun_cov_rep <- apply(Mat_Cov, 1, function(x) sum(x*diff_x_vec))
  return(mean(fun_cov_rep, na.rm = TRUE))
}

for (i in 1:length(SimuList$Method)) {
  load(paste0("SIMU/HulC_Cov_", SimuList$Method[i], " With No. Auction = ", SimuList$Ka[i], ".RData"))
  Table2$Cov_MLE_HulC[i] <- function_coverage(Mat_Cov_MLE)
  Table2$Cov_INIT_HulC[i] <- function_coverage(Mat_Cov_Init)
  Table2$Cov_PT_Cred[i] <- function_coverage(Mat_Cov_PT)
  
  Table2$Area_MLE_HulC[i] <- function_area(Mat_Bw_MLE, x.med.vec)
  Table2$Area_INIT_HulC[i] <- function_area(Mat_Bw_Init, x.med.vec)
  Table2$Area_PT_Cred[i] <- function_area(Mat_Bw_PT, x.med.vec)
}
Table2
######################################################################################
Table2=data.frame(Method=SimuList$Method,
                  Ka=SimuList$Ka,
                  Cov_MLE_HulC=NA,
                  Cov_INIT_HulC=NA,
                  Cov_PT_Cred=NA,
                  Area_MLE_HulC=NA,
                  Area_INIT_HulC=NA,
                  Area_PT_Cred=NA,
                  Cov_MLE_HulC_OB=NA,
                  Cov_INIT_HulC_OB=NA,
                  Cov_PT_Cred_OB=NA)

function_coverage_offbound <- function(Mat_Cov) {
  offbound_x <- which(trueF_point > 0.05 & trueF_point <0.95)
  fun_cov_rep <- apply(Mat_Cov, 1, function(x) {all(x[offbound_x])})
  return(mean(fun_cov_rep,na.rm = TRUE))
}

for (i in 1:length(SimuList$Method)) {
  load(paste0("SIMU/HulC_Cov_", SimuList$Method[i], " With No. Auction = ", SimuList$Ka[i], ".RData"))
  Table2$Cov_MLE_HulC[i] <- function_coverage(Mat_Cov_MLE)
  Table2$Cov_INIT_HulC[i] <- function_coverage(Mat_Cov_Init)
  Table2$Cov_PT_Cred[i] <- function_coverage(Mat_Cov_PT)
  
  Table2$Area_MLE_HulC[i] <- function_area(Mat_Bw_MLE, x.med.vec)
  Table2$Area_INIT_HulC[i] <- function_area(Mat_Bw_Init, x.med.vec)
  Table2$Area_PT_Cred[i] <- function_area(Mat_Bw_PT, x.med.vec)
  
  Table2$Cov_MLE_HulC_OB[i] <- function_coverage_offbound(Mat_Cov_MLE)
  Table2$Cov_INIT_HulC_OB[i] <- function_coverage_offbound(Mat_Cov_Init)
  Table2$Cov_PT_Cred_OB[i] <- function_coverage_offbound(Mat_Cov_PT)
}
Table2
