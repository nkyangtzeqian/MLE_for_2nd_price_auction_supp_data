library(tidyverse)
# library(transport) # For calculating Wasserstein distance between two distribution functions.

source("Main.R")
source("distributions.R")
source("DataGeneration.R")
# source("DisCompare.R")
source("SettingGeneration.R")
source("PT_MCMC.R")
source("HulC_main.R")

# pdf(paste0("Plots/HulC_Simu_", method.in, " With No. Auction = ", Ka, ".pdf"))
# myplot1000
# dev.off()


rm(list = ls())   # Remove everything from the Global environment
Ka=100
method.in <- "unif"
para1 <- 1
para2 <- 20
Delta <- 0
source("Plots/HulC_plot.R")
myplot100=myplot1000
rm(list = ls()[ls()!="myplot100"])   # Remove everything except for myplot100
Ka=1000
method.in <- "unif"
para1 <- 1
para2 <- 20
Delta <- 0
source("Plots/HulC_plot.R")

library(gridExtra)
pdf(paste0("Plots/HulC_Simu_", method.in, ".pdf"),width=16, height=8)
grid.arrange(myplot100, myplot1000, ncol=2)
dev.off()




rm(list = ls())   # Remove everything from the Global environment
Ka=100
method.in <- "piecewise_unif"
para1 <- 2
para2 <- 4
Delta <- 0
source("Plots/HulC_plot.R")
myplot100=myplot1000
rm(list = ls()[ls()!="myplot100"])   # Remove everything except for myplot100
Ka=1000
method.in <- "piecewise_unif"
para1 <- 2
para2 <- 4
Delta <- 0
source("Plots/HulC_plot.R")

library(gridExtra)
pdf(paste0("Plots/HulC_Simu_", method.in, ".pdf"),width=16, height=8)
grid.arrange(myplot100, myplot1000, ncol=2)
dev.off()



rm(list = ls())   # Remove everything from the Global environment
Ka=100
method.in <- "pareto"
para1 <- 3
para2 <- 100
Delta <- 0
source("Plots/HulC_plot.R")
myplot100=myplot1000
rm(list = ls()[ls()!="myplot100"])   # Remove everything except for myplot100
Ka=1000
method.in <- "pareto"
para1 <- 3
para2 <- 100
Delta <- 0
source("Plots/HulC_plot.R")

library(gridExtra)
pdf(paste0("Plots/HulC_Simu_", method.in, ".pdf"),width=16, height=8)
grid.arrange(myplot100, myplot1000, ncol=2)
dev.off()


rm(list = ls())   # Remove everything from the Global environment
Ka=100
method.in <- "gamma"
para1 <- 10
para2 <- 2
Delta <- 0
source("Plots/HulC_plot.R")
myplot100=myplot1000
rm(list = ls()[ls()!="myplot100"])   # Remove everything except for myplot100
Ka=1000
method.in <- "gamma"
para1 <- 10
para2 <- 2
Delta <- 0
source("Plots/HulC_plot.R")

library(gridExtra)
pdf(paste0("Plots/HulC_Simu_", method.in, ".pdf"),width=16, height=8)
grid.arrange(myplot100, myplot1000, ncol=2)
dev.off()


rm(list = ls())   # Remove everything from the Global environment
Ka=100
method.in <- "beta"
para1 <- 2
para2 <- 2
Delta <- 0
source("Plots/HulC_plot.R")
myplot100=myplot1000
rm(list = ls()[ls()!="myplot100"])   # Remove everything except for myplot100
Ka=1000
method.in <- "beta"
para1 <- 2
para2 <- 2
Delta <- 0
source("Plots/HulC_plot.R")

library(gridExtra)
pdf(paste0("Plots/HulC_Simu_", method.in, ".pdf"),width=16, height=8)
grid.arrange(myplot100, myplot1000, ncol=2)
dev.off()