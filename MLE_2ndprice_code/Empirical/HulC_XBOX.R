rm(list=ls())

library(tidyverse)
# library(transport) # For calculating Wasserstein distance between two distribution functions.

source("Main.R")
source("distributions.R")
source("DataGeneration.R")
# source("DisCompare.R")
source("SettingGeneration.R")
source("PT_MCMC.R")

source("HulC_main.R")

x.med.vec <- seq(0,400,by=20)
Delta=0

load(file = "Empirical/Xbox_7day_auctions.RData")


set.seed(123)


#------------------------------------------------------------------------------------
## On Whole dataset, Initial Estimate, MLE, and HUlC
#------------------------------------------------------------------------------------



data <- Data_Gen_Observed(raw_data_list)

reserve.price.Cutoff <- ceiling(quantile(data$r_k, probs = 0.25))

Init.Values <- Initialization_lambda_F(data, reserve.price.Cutoff)
theta.init <- (1-Init.Values$F.y)/ c(1,(1-Init.Values$F.y)[-length(Init.Values$F.y)])
MLE <- MLE_2ndprice(data, lambda.in = Init.Values$lambda, theta.in = theta.init)

set.seed(123456)
reserve.price.Cutoff <- ceiling(quantile(data$r_k, probs = 0.25))
Hulc.Conf <- HulC1d_Ebay(raw_data_all_list = raw_data_list, eval.points = x.med.vec, alpha=.1,
                         MLE.Cutoff = reserve.price.Cutoff, diff_bias = 0, Delta=Delta)

# set.seed(1)
# 
# reserve.price.Cutoff <- ceiling(quantile(data$r_k, probs = 0.3))
# Hulc.Conf <- HulC1d_Ebay(raw_data_all_list = raw_data_list, eval.points = x.med.vec, alpha=.1,
#                          MLE.Cutoff = reserve.price.Cutoff, diff_bias = 0, Delta=0.1)

#------------------------------------------------------------------------------------
## Location of knots for geom_step() function in ggplot below
#------------------------------------------------------------------------------------

HulC.tibble.long <- tibble(knots.vec = c(0, knots(Hulc.Conf$Band.CI.Init$upr)), 
                           Band.CI.Init.lwr = Hulc.Conf$Band.CI.Init$lwr(knots.vec),
                           Band.CI.Init.upr = Hulc.Conf$Band.CI.Init$upr(knots.vec),
                           Band.CI.MLE.lwr = Hulc.Conf$Band.CI.MLE$lwr(knots.vec),
                           Band.CI.MLE.upr = Hulc.Conf$Band.CI.MLE$upr(knots.vec) ) %>% 
  pivot_longer(colnames(.)[-1], names_to = "type.CI", values_to = "y.CI") %>% 
  mutate(
    final.CI = ifelse(
      grepl(pattern = "Init", type.CI), 
      "HulC CI for Initial estimate",
      "HulC CI for MLE"
    ) 
  )


#------------------------------------------------------------------------------------

#------------------------------------------------------------------------------------

F.tibble.long <- tibble("x.cord" = c(0,data$Z_iT_i[,1]),
                        "Initial estimate of F" = c(0,Init.Values$F.y),
                        "MLE of F" = c(0,MLE$F.y) ) %>% 
  pivot_longer(colnames(.)[-1], names_to = "type.F", values_to = "y.cord")

#-------------------------------------------------------------------------------
## Plot
#-------------------------------------------------------------------------------

( myplot <- ggplot(data = F.tibble.long) +
    geom_line(mapping = aes(x = x.cord, y = y.cord, color = type.F, linetype = type.F)) +
    geom_step(data = filter(HulC.tibble.long, 
                            type.CI == "Band.CI.Init.upr" | type.CI == "Band.CI.MLE.upr"),  
              mapping = aes(x = knots.vec, y = y.CI, color = final.CI, linetype = final.CI),
              direction = "vh") +
    geom_step(data = filter(HulC.tibble.long, 
                            type.CI == "Band.CI.Init.lwr" | type.CI == "Band.CI.MLE.lwr"),  
              mapping = aes(x = knots.vec, y = y.CI, color = final.CI, linetype = final.CI),
              direction = "hv") +
    geom_point(data = HulC.tibble.long[-(1:4), ],  
               mapping = aes(x = knots.vec, y = y.CI, color = final.CI), 
               show.legend = FALSE) +
    scale_color_manual(values = c("blue", "red", "blue", "red")) +
    scale_linetype_manual(values = c(2,2,1,1)) +
    labs(x = "x", y = "F(x)", color = "", linetype = "") +
    theme(legend.position = "top")
)

pdf(paste0("Empirical/HulC_XBOX_all.pdf"), width = 16, height = 8)
myplot
dev.off()