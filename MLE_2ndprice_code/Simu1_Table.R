library(tidyverse)
library(GoFKernel) # For the inverse function of a CDF


# source("0327/EbayFunctions_Ver3.R")
source("Main.R")
source("distributions.R")
source("DataGeneration.R")
source("DisCompare.R")
source("SettingGeneration.R")

source("Para_MCMC.R")
source("ProfitMetric.R")
#####################################################################################

rm(list = ls())   # Remove everything from the Global environment
Ka=100
method.in <- "unif"
para1 <- 1
para2 <- 20
source("Simu1.R")

rm(list = ls())   # Remove everything from the Global environment
Ka=100
method.in <- "piecewise_unif"
para1 <- 2
para2 <- 4
source("Simu1.R")

rm(list = ls())   # Remove everything from the Global environment
Ka=100
method.in <- "pareto"
para1 <- 3
para2 <- 100
source("Simu1.R")

rm(list = ls())   # Remove everything from the Global environment
Ka=100
method.in <- "gamma"
para1 <- 10
para2 <- 2
source("Simu1.R")

rm(list = ls())   # Remove everything from the Global environment
Ka=100
method.in <- "beta"
para1 <- 2
para2 <- 2
source("Simu1.R")




rm(list = ls())   # Remove everything from the Global environment
Ka=1000
method.in <- "unif"
para1 <- 1
para2 <- 20
source("Simu1.R")

rm(list = ls())   # Remove everything from the Global environment
Ka=1000
method.in <- "piecewise_unif"
para1 <- 2
para2 <- 4
source("Simu1.R")

rm(list = ls())   # Remove everything from the Global environment
Ka=1000
method.in <- "pareto"
para1 <- 3
para2 <- 100
source("Simu1.R")

rm(list = ls())   # Remove everything from the Global environment
Ka=1000
method.in <- "gamma"
para1 <- 10
para2 <- 2
source("Simu1.R")

rm(list = ls())   # Remove everything from the Global environment
Ka=1000
method.in <- "beta"
para1 <- 2
para2 <- 2
source("Simu1.R")


#####################################################################################
## add Gamma and truncated Normal appendix

rm(list = ls())   # Remove everything from the Global environment
Ka=100
method.in <- "unif"
para1 <- 1
para2 <- 20
source("simu1_add_ParaF.R")

rm(list = ls())   # Remove everything from the Global environment
Ka=100
method.in <- "piecewise_unif"
para1 <- 2
para2 <- 4
source("simu1_add_ParaF.R")

rm(list = ls())   # Remove everything from the Global environment
Ka=100
method.in <- "pareto"
para1 <- 3
para2 <- 100
source("simu1_add_ParaF.R")

rm(list = ls())   # Remove everything from the Global environment
Ka=100
method.in <- "gamma"
para1 <- 10
para2 <- 2
source("simu1_add_ParaF.R")

rm(list = ls())   # Remove everything from the Global environment
Ka=100
method.in <- "beta"
para1 <- 2
para2 <- 2
source("simu1_add_ParaF.R")




rm(list = ls())   # Remove everything from the Global environment
Ka=1000
method.in <- "unif"
para1 <- 1
para2 <- 20
source("simu1_add_ParaF.R")

rm(list = ls())   # Remove everything from the Global environment
Ka=1000
method.in <- "piecewise_unif"
para1 <- 2
para2 <- 4
source("simu1_add_ParaF.R")

rm(list = ls())   # Remove everything from the Global environment
Ka=1000
method.in <- "pareto"
para1 <- 3
para2 <- 100
source("simu1_add_ParaF.R")

rm(list = ls())   # Remove everything from the Global environment
Ka=1000
method.in <- "gamma"
para1 <- 10
para2 <- 2
source("simu1_add_ParaF.R")

rm(list = ls())   # Remove everything from the Global environment
Ka=1000
method.in <- "beta"
para1 <- 2
para2 <- 2
source("simu1_add_ParaF.R")

#####################################################################################
SimuList=list(
  Method=rep(c("unif","piecewise_unif","pareto","gamma","beta"),each=2),
  Ka=rep(c(100,1000),times=5)
)

Table1=data.frame(Method=SimuList$Method,
                  Ka=SimuList$Ka,
                  KS.mle.mean=NA,
                  KS.init.mean=NA,
                  KS.pt.mean=NA,
                  KS.gamma.mean=NA,
                  KS.tNorm.mean=NA,
                  TV.mle.mean=NA,
                  TV.init.mean=NA,
                  TV.pt.mean=NA,
                  TV.gamma.mean=NA,
                  TV.tNorm.mean=NA
)

mean_T=function(x){ ## TN-F always fail
  mean(x[x<1],na.rm = TRUE)
}

for (i in 1:length(SimuList$Method)) {
  load(paste0("SIMU/Table1A_", SimuList$Method[i], "_K=", SimuList$Ka[i], ".RData"))
  Table1$KS.mle.mean[i] <- mean_T(distanceList$KS.mle)
  Table1$KS.init.mean[i] <- mean_T(distanceList$KS.init)
  Table1$KS.pt.mean[i] <- mean_T(distanceList$KS.pt)
  Table1$KS.gamma.mean[i] <- mean_T(distanceList$Ks.gamma)
  Table1$KS.tNorm.mean[i] <- mean_T(distanceList$KS.tNorm)
  Table1$TV.mle.mean[i] <- mean_T(distanceList$TV.mle)
  Table1$TV.init.mean[i] <- mean_T(distanceList$TV.init)
  Table1$TV.pt.mean[i] <- mean_T(distanceList$TV.pt)
  Table1$TV.gamma.mean[i] <- mean_T(distanceList$TV.gamma)
  Table1$TV.tNorm.mean[i] <- mean_T(distanceList$TV.tNorm)
}

Table1

