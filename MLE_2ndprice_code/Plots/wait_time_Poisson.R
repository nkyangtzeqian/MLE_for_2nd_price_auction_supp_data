rm(list = ls())   # Remove everything from the Global environment
Ka=100

# method.in <- "unif"
# para1 <- 1
# para2 <- 20


# method.in <- "piecewise_unif"
# para1 <- 2
# para2 <- 4
# source("Main.R")
# 
# 
# method.in <- "pareto"
# para1 <- 3
# para2 <- 100
# 
# 
method.in <- "gamma"
para1 <- 10
para2 <- 2
# 
# 
# method.in <- "beta"
# para1 <- 2
# para2 <- 2


source("distributions.R")
source("DataGeneration.R")
source("DisCompare.R")
source("SettingGeneration.R")


lambda0 <- 1
tau_a <- 1000


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


X=c()
Y=c()

for (i in 1:data$Ka) {
  xx=c(data$data_list_Obs[[i]]$r_k,data$data_list_Obs[[i]]$x_i_p)
  xx=xx[-length(xx)] ## remove final bid price
  yy=data$data_list_Obs[[i]]$t_i_p
  X=c(X,xx)
  Y=c(Y,log(yy))
}

plot(X,Y,xlab="Current selling price",ylab="log(Waiting time)",
     main=paste0("Scatter plot of log(Waiting time) vs. Current selling price"))
     # main=paste0("Scatter plot of waiting time vs. current bidding price for ",method.in," distribution"))

mod=lm(Y~X)
abline(mod,col="red")
# text(1,4, paste0("R2 = ",round(summary(mod)$r.squared,3)), col="red")



## data oberseved and placed bids, modified from DataGeneration.R
Data_Gen_Observed_placed <- function(Raw_data_all_list){
  if(Raw_data_all_list$class != "SecondPriceAuction_Rawdata"){
    stop("This fucntion works only with data in the class 'SecondPriceAuction.Rawdata'.
         See e.g., 'Data_Gen_Raw' where we generate such data.")
  } 
  
  Ka <- Raw_data_all_list$Ka
  tau_a <- Raw_data_all_list$tau_a
  r_k <- Raw_data_all_list$r_k
  
  # # pooled data of reserve prices, observed bids, inter arrival times (including T_0) for all the auctions with M>0
  # pooled.data <- matrix(0, ncol=2, nrow=0) # Pooled data from all the auctions with M>0
  # colnames(pooled.data) <- c("pooled_rec_value", "pooled_wait_time")
  
  Z_iT_i <- matrix(0, ncol=2, nrow=0) # Pooled data from ell + K standing prices from all the auctions (including the reserve prices)
  colnames(Z_iT_i) <- c("pooled_rec_value", "pooled_wait_time")
  Xbar_i <- NULL
  X_1 <- X_mk <- rep(0, Ka) # X_1 is the vector of first jump values / first current selling prices among all the auctions.
  
  T_0 <- rep(0, Ka) # First Time when selling price is strictly larger than reserve price
  K_sold <- rep(0, Ka) # K_sold[which(K_sold!=0)] is the Vector of all the auction indexes where the item has been sold.
  data_list_Obs <- vector("list", Ka) # observed data from all the auctions.
  M_k <- rep(0, Ka)
  
  pp <- c() ## placed bid vector
  pt <- c() ## waiting time vector for placed bid
  
  for (jj in 1: Ka) {
    all_bid_price_j <- Raw_data_all_list$Raw_biddata_list[[jj]][,2] # All the private values of the bidders for the jth process
    # Creating the observed bid vector
    if(sum(r_k[jj]<= all_bid_price_j) == 0){ ## M=O=1
      # xx <- NA
      X_mk[jj] <- NA
      X_1[jj] <- NA # first jump value / first current selling price of the jth auction.
      # pooled.data <- pooled.data   # Since, we don't observe any bid, so pooled.data remains same as the previous step.
      T_0[jj] <- tau_a
      data_list_Obs[[jj]] <- list(r_k = r_k[jj], O_k = 0,
                                  x_i_p = NULL, t_i_p = NULL)
      M_k[jj] <- 0
      Z_iT_i <- rbind(Z_iT_i, cbind(r_k[jj], tau_a))
      
    } else if( sum(r_k[jj]<= all_bid_price_j) == 1){ ## M=0,O=1. SHOULD NOT BE DELTETED: first jump
      
      X_mk[jj] <- r_k[jj]
      X_1[jj] <- NA   # first jump value / first observed bid of the jth auction.
      # pooled.data <- pooled.data   # Since, we don't observe any bid, so pooled.data remains same as the previous step.
      T_0[jj] <- tau_a
      K_sold[jj] <- jj
      data_list_Obs[[jj]] <- list(r_k = r_k[jj], O_k = 1,
                                  x_i_p = NULL, t_i_p = NULL)
      M_k[jj] <- 0
      Z_iT_i <- rbind(Z_iT_i, cbind(r_k[jj], tau_a))
      
    } else if( sum(r_k[jj]<= all_bid_price_j) >= 2){ ## M>=1,O=1
      xx <- c(r_k[jj], rep(0,length(all_bid_price_j))) # Observed (repeated) 2nd price vector of length = no. of bids
      pt_raw <- c() ## time vector for placed bid
      for (ii in 1:length(all_bid_price_j)){
        if(xx[ii] >= all_bid_price_j[ii]){ ## selling price not changed
          xx[ii+1] <- xx[ii]
        } else { 
          temp_2nd_max <- sort(c(r_k[jj],all_bid_price_j[1:ii]), ## from unobserved data
                               decreasing = TRUE)[2]
          xx[ii+1] <- temp_2nd_max
          pp <- c(pp, all_bid_price_j[ii]) ## placed bid price recorded
          pt_raw <- c(pt_raw, Raw_data_all_list$Raw_biddata_list[[jj]][,3][ii]) ##  time for placed bid
        }
      }# END of xx
      pp <- pp[-length(pp)] ## remove the last one, which is not corresponding to next placed bid
      pt <- c(pt,diff(pt_raw)) ## waiting time after each placed bid
      
      inds <- 2:(length(all_bid_price_j)+1) ## reserved not included!
      inds <- inds[diff(xx)>0] # indexes when selling price changes
      sum_ti <- Raw_data_all_list$Raw_biddata_list[[jj]][,3][inds-1] # Vector (T_0, T_0 + T_1, ...., \sum_{i=0}^{i=M-1}T_i)'
      t_i <- diff(c(0,sum_ti, tau_a)) # Vector (T_0, T_1, ...., T_M)'
      # pooled.data <- rbind(pooled.data, 
      #                      cbind(c(r_k[jj], xx[inds]),t_i) )
      Z_iT_i <- rbind(Z_iT_i, 
                      cbind(c(r_k[jj], xx[inds]),t_i) )
      Xbar_i <- c(Xbar_i, xx[inds])  # Pooled observed bids / 2nd price values for all the auctions.
      X_mk[jj] <- tail(xx[inds], n=1)
      X_1[jj] <- xx[inds][1]     # first jump value / first observed bid of the jth auction. For initialization
      T_0[jj] <- t_i[1]
      K_sold[jj] <- jj
      data_list_Obs[[jj]] <- list(r_k = r_k[jj], O_k = 1,
                                  x_i_p = xx[inds],
                                  t_i_p = t_i[-1]) ## DO NOT include T_0, don't use it directly.
      M_k[jj] <- length(xx[inds])
    }# END of M>=1,O=1
  }# END of for(jj in 1: Ka)
  
  # pooled.data <- pooled.data[order(pooled.data[,1]),]
  Z_iT_i <- Z_iT_i[order(Z_iT_i[,1]),,drop=F]
  K_sold <- K_sold[which(K_sold!=0)]
  
  ## help to do the calculation
  # S_bar <- match(X_mk[K_sold], Z_iT_i[ ,1])  # old version, only sold above
  S_bar <- match(X_mk[!is.na(X_mk)], Z_iT_i[ ,1]) # Vector/Set of ranks/positions of the selling prices of all the auctions (where item is sold) in the pooled data set.
  u_i <- sort(match(Xbar_i, Z_iT_i[ ,1]))  # This is vector: (u1, u2,...,uL).
  l_i_3p2 <- sapply(1:(sum(M_k)+Ka), function(ii) {sum(u_i <= ii) })
  Q_i_card <- sapply(1:(sum(M_k)+Ka), function(ii) {sum(S_bar>=ii)})
  
  auxlist <- list( 
    S_bar = S_bar, ## seems that we don't need it
    u_i = u_i,
    l_i_3p2 = l_i_3p2,
    Q_i_card = Q_i_card
  )
  ret <- list(data_list_Obs = data_list_Obs,
              Xbar_i = Xbar_i, #Set of Xi's, not ordered
              X_mk = X_mk,
              T_0 = T_0,
              K_sold = K_sold,
              # pooled.data = pooled.data,
              Z_iT_i = Z_iT_i,
              X_1 = X_1,
              Ka = Raw_data_all_list$Ka,
              tau_a = Raw_data_all_list$tau_a,
              r_k = Raw_data_all_list$r_k,
              F_ture = Raw_data_all_list$True.Method,
              M_k = M_k,
              auxlist = auxlist,
              new_placed=cbind(pp,pt)) ## placed bid data
  
  ret$class = "SecondPriceAuction_Processed"
  return(ret)
}

data2 <- Data_Gen_Observed_placed(Raw_data_all_list)

X=data2$new_placed[,1]
Y=data2$new_placed[,2]
Y=log(Y)

plot(X,Y,xlab="Placed bid price",ylab="log(Waiting time)",
     main=paste0("Scatter plot of log(Waiting time) vs. Placed bid price"))
# main=paste0("Scatter plot of waiting time vs. current bidding price for ",method.in," distribution"))

mod=lm(Y~X)
abline(mod,col="red")