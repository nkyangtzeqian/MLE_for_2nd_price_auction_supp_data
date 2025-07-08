library(tidyverse)
library(GoFKernel) # For the inverse function of a CDF


# source("0327/EbayFunctions_Ver3.R")
source("Main.R")
source("distributions.R")
source("DataGeneration.R")
source("DisCompare.R")
source("SettingGeneration.R")


rm(list = ls())   # Remove everything from the Global environment
lambda0=10
method.in <- "unif"
para1 <- 1
para2 <- 20
replication <- 30
source("Simu1_2_morebidder.R")

rm(list = ls())   # Remove everything from the Global environment
lambda0=10
method.in <- "piecewise_unif"
para1 <- 2
para2 <- 4
replication <- 30
source("Simu1_2_morebidder.R")

rm(list = ls())   # Remove everything from the Global environment
lambda0=10
method.in <- "pareto"
para1 <- 3
para2 <- 100
replication <- 30
source("Simu1_2_morebidder.R")

rm(list = ls())   # Remove everything from the Global environment
lambda0=10
method.in <- "gamma"
para1 <- 10
para2 <- 2
replication <- 30
source("Simu1_2_morebidder.R")

rm(list = ls())   # Remove everything from the Global environment
lambda0=10
method.in <- "beta"
para1 <- 2
para2 <- 2
replication <- 30
source("Simu1_2_morebidder.R")



rm(list = ls())   # Remove everything from the Global environment
lambda0=50
method.in <- "unif"
para1 <- 1
para2 <- 20
replication <- 30
source("Simu1_2_morebidder.R")

rm(list = ls())   # Remove everything from the Global environment
lambda0=50
method.in <- "piecewise_unif"
para1 <- 2
para2 <- 4
replication <- 30
source("Simu1_2_morebidder.R")

rm(list = ls())   # Remove everything from the Global environment
lambda0=50
method.in <- "pareto"
para1 <- 3
para2 <- 100
replication <- 30
source("Simu1_2_morebidder.R")

rm(list = ls())   # Remove everything from the Global environment
lambda0=50
method.in <- "gamma"
para1 <- 10
para2 <- 2
replication <- 30
source("Simu1_2_morebidder.R")

rm(list = ls())   # Remove everything from the Global environment
lambda0=50
method.in <- "beta"
para1 <- 2
para2 <- 2
replication <- 30
source("Simu1_2_morebidder.R")
#####################################################################################
SimuList=list(
  Method=rep(c("unif","piecewise_unif","pareto","gamma","beta"),each=2),
  lambda0=rep(c(10,50),times=5)
)



Table12=data.frame(Method=SimuList$Method,
                   lambda0=SimuList$lambda0,
                  KS.mle.mean=NA,
                  KS.init.mean=NA,
                  KS.pt.mean=NA,
                  TV.mle.mean=NA,
                  TV.init.mean=NA,
                  TV.pt.mean=NA)

for (i in 1:length(SimuList$Method)) {
  load(paste0("SIMU/Table1_", SimuList$Method[i], "_Lbd=", SimuList$lambda0[i], ".RData"))
  Table12$KS.mle.mean[i] <- mean(distanceList$KS.mle)
  Table12$KS.init.mean[i] <- mean(distanceList$KS.init)
  Table12$KS.pt.mean[i] <- mean(distanceList$KS.pt)
  Table12$TV.mle.mean[i] <- mean(distanceList$TV.mle)
  Table12$TV.init.mean[i] <- mean(distanceList$TV.init)
  Table12$TV.pt.mean[i] <- mean(distanceList$TV.pt)
}

Table12