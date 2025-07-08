# rm(list = ls())   # Remove everything from the Global environment
# Ka=100
# method.in <- "unif"
# para1 <- 1
# para2 <- 20
# Delta <- 0
# source("HulC_plot.R")



library(tidyverse)
# library(transport) # For calculating Wasserstein distance between two distribution functions.
# library(gridExtra)

source("Main.R")
source("distributions.R")
source("DataGeneration.R")
# source("DisCompare.R")
source("SettingGeneration.R")
source("PT_MCMC.R")


source("HulC_main.R")

#------------------------------------------------------------------------------------
## User Input: True lambda, Auction window, Number of independent auctions 'N.auction'
#------------------------------------------------------------------------------------

# repetition_number <- 100  # replication numbers of each such event i.e. selling an 
# identical stuff on ebay with multiple independent auctions.

lambda0 <- 1
tau_a <- 100
# Ka <- 150

#------------------------------------------------------------------------------------
## Choice of different true valuation distributions F, the corresponding choice of parameters,  
## and reserve prices are listed below:
#------------------------------------------------------------------------------------

set.seed(123)

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

Hulc.Conf <- HulC1d_Ebay(raw_data_all_list = Raw_data_all_list, eval.points = x.med.vec, alpha=.1,
                         MLE.Cutoff = reserve.price.Cutoff, diff_bias = 0, Delta=Delta)

## Polya Tree method
data_Polya <- Data_Gen_Unobserved_PolyaTree(Raw_data_all_list)
Bayes_PT <- MCMC_PT(data_Polya, data$Z_iT_i[,1], x.max=max(x.med.vec))

Prediction_PT <- sapply(1:nrow(Bayes_PT$MCMC_Mat), function(ii) {
  cutpt <- c(0, rev(data_Polya$dat[,2]), max(x.med.vec))
  step_est <- approxfun(cutpt, c(0, cumsum(Bayes_PT$MCMC_Mat[ii,])  ))
  return(step_est(x.med.vec))
})
Band.CI.PT.lwr = apply(Prediction_PT, 1, function(x) quantile(x, probs = 0.05))
Band.CI.PT.upr = apply(Prediction_PT, 1, function(x) quantile(x, probs = 0.95))

#------------------------------------------------------------------------------------
## Location of knots for geom_step() function in ggplot below
#------------------------------------------------------------------------------------

HulC.tibble.long <- tibble(knots.vec = c(0, knots(Hulc.Conf$Band.CI.Init$upr)), 
                           Band.CI.Init.lwr = Hulc.Conf$Band.CI.Init$lwr(knots.vec),
                           Band.CI.Init.upr = Hulc.Conf$Band.CI.Init$upr(knots.vec),
                           Band.CI.MLE.lwr = Hulc.Conf$Band.CI.MLE$lwr(knots.vec),
                           Band.CI.MLE.upr = Hulc.Conf$Band.CI.MLE$upr(knots.vec),
                           Band.CI.PT.lwr = c(0,Band.CI.PT.lwr),
                           Band.CI.PT.upr = c(0,Band.CI.PT.upr)) %>% 
  pivot_longer(colnames(.)[-1], names_to = "type.CI", values_to = "y.CI") %>% 
  mutate(
    # final.CI = ifelse(
    #   grepl(pattern = "Init", type.CI), 
    #   "HulC CI for Initial estimate",
    #   "HulC CI for MLE"
    # ) 
    final.CI = sapply(type.CI, function(x)
      switch(grepl(pattern = "Init", x)*1 +
               grepl(pattern = "MLE", x)*2 + 
               grepl(pattern = "PT", x)*3,
             "HulC CI for Initial estimate",
             "HulC CI for MLE",
             "HulC CI for Polya Tree"))
  )


#------------------------------------------------------------------------------------
## Choices of True F
#------------------------------------------------------------------------------------

if(method.in == "unif"){
  F.tibble.long <- tibble("x.cord" = c(0,data$Z_iT_i[,1]),
                          "True F" = punif(c(0,data$Z_iT_i[,1]), min = para1, max = para2),
                          "Initial estimate of F" = c(0,Init.Values$F.y),
                          "MLE of F" = c(0,MLE$F.y),
                          "PT estimate of F" = c(0,Bayes_PT$F.y)) %>% 
    pivot_longer(colnames(.)[-1], names_to = "type.F", values_to = "y.cord")
  
}else if( method.in == "piecewise_unif"){
  F0.piecewise.unif_pooled.data <- sapply(c(0,data$Z_iT_i[,1]),
                                          function(y) {piecewise.unif.cdf(y,
                                                                          para1 = para1,
                                                                          para2 = para2
                                          )
                                          }
  )
  
  F.tibble.long <- tibble("x.cord" = c(0,data$Z_iT_i[,1]),
                          "True F" = F0.piecewise.unif_pooled.data,
                          "Initial estimate of F" = c(0,Init.Values$F.y),
                          "MLE of F" = c(0,MLE$F.y),
                          "PT estimate of F" = c(0,Bayes_PT$F.y)) %>% 
    pivot_longer(colnames(.)[-1], names_to = "type.F", values_to = "y.cord")
  
}else if (method.in == "pareto"){
  F.tibble.long <- tibble("x.cord" = c(0,data$Z_iT_i[,1]),
                          "True F" = ppareto(c(0,data$Z_iT_i[,1]), m = para1, s = para2),
                          "Initial estimate of F" = c(0,Init.Values$F.y),
                          "MLE of F" = c(0,MLE$F.y),
                          "PT estimate of F" = c(0,Bayes_PT$F.y)) %>% 
    pivot_longer(colnames(.)[-1], names_to = "type.F", values_to = "y.cord")
  
} else if( method.in == "gamma"){
  F.tibble.long <- tibble("x.cord" = c(0,data$Z_iT_i[,1]),
                          "True F" = pgamma(c(0,data$Z_iT_i[,1]), shape = para1, rate = para2),
                          "Initial estimate of F" = c(0,Init.Values$F.y),
                          "MLE of F" = c(0,MLE$F.y),
                          "PT estimate of F" = c(0,Bayes_PT$F.y)) %>% 
    pivot_longer(colnames(.)[-1], names_to = "type.F", values_to = "y.cord")
  
}else if (method.in == "beta"){  
  F.tibble.long <- tibble("x.cord" = c(0,data$Z_iT_i[,1]),
                          "True F" = pbeta(c(0,data$Z_iT_i[,1]), shape1 = para1, shape2 = para2),
                          "Initial estimate of F" = c(0,Init.Values$F.y),
                          "MLE of F" = c(0,MLE$F.y),
                          "PT estimate of F" = c(0,Bayes_PT$F.y)) %>% 
    pivot_longer(colnames(.)[-1], names_to = "type.F", values_to = "y.cord")
  
}else{
  stop("The method is not supported. Please use one of the following methods: 
       'unif', 'pareto', 'gamma', 'beta', 'piecewise_unif'.")
}

#-------------------------------------------------------------------------------
## Plot
#-------------------------------------------------------------------------------

( myplot1000 <- ggplot(data = F.tibble.long) +
    geom_line(mapping = aes(x = x.cord, y = y.cord, color = type.F, linetype = type.F)) +
    geom_step(data = filter(HulC.tibble.long, 
                            type.CI == "Band.CI.Init.upr" | type.CI == "Band.CI.MLE.upr" | type.CI == "Band.CI.PT.upr"),  
              mapping = aes(x = knots.vec, y = y.CI, color = final.CI, linetype = final.CI),
              direction = "vh") +
    geom_step(data = filter(HulC.tibble.long, 
                            type.CI == "Band.CI.Init.lwr" | type.CI == "Band.CI.MLE.lwr"| type.CI == "Band.CI.PT.lwr"),  
              mapping = aes(x = knots.vec, y = y.CI, color = final.CI, linetype = final.CI),
              direction = "hv") +
    geom_point(data = HulC.tibble.long[-(1:6), ],  
               mapping = aes(x = knots.vec, y = y.CI, color = final.CI), 
               show.legend = FALSE) +
    scale_color_manual(values = c("blue", "red","green4", "blue", "red","green4", "black")) +
    scale_linetype_manual(values = c(2,2,2,1,1,1,1)) +
    labs(x = "x", y = "F(x)", color = "", linetype = "",
         title = paste("Number of auctions =", Ka, ", True F:", method.in)) +
    theme(legend.position = "top")
)


# 
# #-------------------------------------------------------------------------------
# ## plot two images in the same plot
# #-------------------------------------------------------------------------------
# 
# library(grid)
# library(gridExtra)
# plots <- list(myplot100, myplot1000)
# layout <- rbind(c(1,2))
# grid.arrange(grobs = plots, layout_matrix = layout)
# 
# 
# #====================================================================================
# ## Alternative approach to plot Initial estimate, MLE, True F, and HulC confidence regions.
# ## Using Base R (without using ggplot package).
# #====================================================================================
# 
# ## True F: Uniform
# plot(c(0,data$pooled.data[,1]), punif(c(0,data$pooled.data[,1]), min = para1, max = para2),
#      col = "black", xlab = "x", ylab = "F(x)", type = "l",
#      main = paste("No. of auctions =", N.auction, ", True F: Uniform"))
# 
# ## True F: Piecewise Uniform
# F0.piecewise.unif_pooled.data <- sapply(c(0,data$pooled.data[,1]),
#                                         function(y) {piecewise.unif.cdf(y,
#                                                                         para1 = para1,
#                                                                         para2 = para2
#                                         )
#                                         }
# )
# plot(c(0,data$pooled.data[,1]), F0.piecewise.unif_pooled.data,
#      col = "black", xlab = "x", ylab = "F(x)", type = "l",
#      main = paste("No. of auctions =", N.auction))
# 
# ## True F: Pareto
# plot(c(0,data$pooled.data[,1]), ppareto(c(0,data$pooled.data[,1]), m = para1, s = para2),
#      col = "black", xlab = "x", ylab = "F(x)", type = "l",
#      main = paste("No. of auctions =", N.auction))
# 
# ## True F: Gamma
# plot(c(0,data$pooled.data[,1]), pgamma(c(0,data$pooled.data[,1]), shape = para1, rate = para2),
#      col = "black", xlab = "x", ylab = "F(x)", type = "l",
#      main = paste("No. of auctions =", N.auction))
# 
# ## True F: Beta
# plot(c(0,data$pooled.data[,1]), pbeta(c(0,data$pooled.data[,1]), shape1 = para1, shape2 = para2),
#      col = "black", xlab = "x", ylab = "F(x)", type = "l",
#      main = paste("No. of auctions =", N.auction))
# 
# 
# ## Adding the plots of Initial Estimate and MLE:
# lines(c(0,Init.Values$F.x), c(0,Init.Values$F.y), type = "l", col = "blue")
# lines(c(0,MLE$F.x), c(0,MLE$F.y), col = "red")
# 
# 
# ## HulC Confidence band for the whole dataset:
# 
# # points(x.med.vec, as.vector( Hulc.Conf$CI.Init[1,]), col="blue", lwd=1)
# # points(x.med.vec, as.vector( Hulc.Conf$CI.Init[2,]), col="blue", lwd=1)
# # points(x.med.vec, as.vector( Hulc.Conf$CI.MLE[1,]), col="red", lwd=1)
# # points(x.med.vec, as.vector( Hulc.Conf$CI.MLE[2,]), col="red", lwd=1)
# 
# plot(Hulc.Conf$Band.CI.Init$lwr, add=TRUE, col="blue", lwd=1, lty = 2)
# plot(Hulc.Conf$Band.CI.Init$upr, add=TRUE, col="blue", lwd=1, lty = 2)
# plot(Hulc.Conf$Band.CI.MLE$lwr, add=TRUE, col="red", lwd=1, lty = 2)
# plot(Hulc.Conf$Band.CI.MLE$upr, add=TRUE, col="red", lwd=1, lty = 2)
# 
# legend(x = "bottomright", legend = c("True F", "Initial estimate of F", "MLE of F", 
#                                      "HulC CI for Initial estimate", "HulC CI for MLE"),
#        col = c("black", "blue", "red", "blue", "red"), 
#        lty = c(1,1,1,2,2), lwd = 1, cex = 0.5, horiz = F)
# 
