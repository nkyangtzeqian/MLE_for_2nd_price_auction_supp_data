

## Functions to generate the data where $Ka$ is not fixed but the total time
## window is fixed.
## Raw means the data is not observed.

Data_Gen_Raw <- function(Ka, tau_a, lambda0, r_k,
                                   method = c("unif", "pareto", "gamma", "beta",
                                              "piecewise_unif", "customized"),
                                   para1 = NULL, para2 = NULL, customized_CDF = NULL){
  N_all <- rep(0, Ka) # Total number of potential bids for the jth auction, might be smaller than r. Not observed.
  if(is.null(r_k)) r_k <- rep(0, Ka)
  Raw_biddata_list <- vector("list", Ka) # List containing each auction data separately.
  
  for (jj in 1: Ka) { ## for each auctions
    tts <- rexp( (2* tau_a* lambda0), rate = lambda0) ## why not use while loop
    N_all[jj] <- sum(cumsum(tts) <= tau_a) 
    t_otts <- cumsum(tts)[1:N_all[jj]]
    if(method == "unif"){
      yy <- runif(n = N_all[jj], min=para1, max=para2) # All the private values of the bidders for the jth process
    }else if (method == "pareto"){
      yy <- rpareto(n = N_all[jj], m = para1, s = para2) # All the private values of the bidders for the jth process
    } else if( method == "gamma"){
      yy <- rgamma(n = N_all[jj], shape = para1, rate=para2) # All the private values of the bidders for the jth process
    }else if (method == "beta"){  
      yy <- rbeta(n = N_all[jj], shape1 = para1, shape2 = para2) # All the private values of the bidders for the jth process
    }else if (method == "piecewise_unif"){
      uu <- runif(n = N_all[jj])
      yy <- sapply(uu, 
                   function(u) {piecewise.unif.invcdf(u, para1 = para1, para2 = para2)})
    }else if (method == "customized"){
      if(is.null(customized_CDF)){
        stop("For the method = 'customized', please provide the customized CDF as an object of class 'AbscontDistribution'.")
      }else if(class(customized_CDF)[1]!="AbscontDistribution"){
        stop("For the method = 'customized', please provide the customized CDF as an object of class 'AbscontDistribution'.")
      }
      yy <- r(customized_CDF)(N_all[jj])
    }
    
    Raw_biddata_matrix <- matrix(0, nrow = N_all[jj], ncol = 4)
    colnames(Raw_biddata_matrix) <- c("auction_index", "bid_value", "bid_time_cum", "openbid_N_all")
    Raw_biddata_matrix[,1] <- rep(jj, N_all[jj])
    Raw_biddata_matrix[,2] <- yy
    Raw_biddata_matrix[,3] <- t_otts
    Raw_biddata_matrix[,4] <- rep(r_k[jj], N_all[jj])
    
    Raw_biddata_list[[jj]] <- Raw_biddata_matrix
  }# END of for(jj in 1: Ka)
  
  Raw_data_all_list <- list(Raw_biddata_list = Raw_biddata_list,
                        Ka = Ka,
                        tau_a = tau_a,
                        r_k = r_k,
                        True.Method = method) 
  Raw_data_all_list$class = "SecondPriceAuction_Rawdata"
  
  return(Raw_data_all_list)
}


## data oberseved
Data_Gen_Observed <- function(Raw_data_all_list){
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
      for (ii in 1:length(all_bid_price_j)){
        if(xx[ii] >= all_bid_price_j[ii]){ ## selling price not changed
          xx[ii+1] <- xx[ii]
        } else { 
          temp_2nd_max <- sort(c(r_k[jj],all_bid_price_j[1:ii]), ## from unobserved data
                               decreasing = TRUE)[2]
          xx[ii+1] <- temp_2nd_max
        }
      }# END of xx
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
              auxlist = auxlist)
  
  ret$class = "SecondPriceAuction_Processed"
  return(ret)
}

# 
## test
# set.seed(123)
# lambda0 <- 1
# tau_a <- 100
# Ka <- 100
# 
# method.in <- "unif"
# para1 <- 1
# para2 <- 10
# r_k <- runif(Ka, 5, 11)
# 
# Raw_data_all_list <- Data_Gen_Raw(Ka, tau_a, lambda0,
#                                         r_k, method = method.in, para1, para2)
# data <- Data_Gen_Observed(Raw_data_all_list)
# rm(list = ls()[ls() %in% c("Raw_data_all_list")])
# rm(list = ls()[!ls() %in% c("data")])
# 
# 
# sum(data$M_k)+Ka
# length(data$K_sold)
# # which(data$T_0==tau_a)
# sum(is.na(data$X_mk)) ## O=M=0
# sum(is.na(data$X_1)) ## M=0
# 
# source("Main.R")
# source("distributions.R")
# Init.Values <- Initialization_lambda_F(data, 8)
# theta.init <- (1-Init.Values$F.y)/ c(1,(1-Init.Values$F.y)[-length(Init.Values$F.y)])
# MLE <- MLE_2ndprice_constra(data, lambda.in = Init.Values$lambda, theta.in = theta.init, tol = 1e-5)

## slice raw data
Data_Slice_Raw <- function(Raw_data_all_list,index){
  if(Raw_data_all_list$class != "SecondPriceAuction_Rawdata"){
    stop("This fucntion works only with data in the class 'SecondPriceAuction.Rawdata'.
         See e.g., 'Data_Gen_Raw' where we generate such data.")
  } 
  
  ## slicing
  Raw_slice_all_list <- Raw_data_all_list
  Raw_slice_all_list$Raw_biddata_list <- Raw_slice_all_list$Raw_biddata_list[index]
  Raw_slice_all_list$Ka <- length(index)
  Raw_slice_all_list$r_k <- Raw_slice_all_list$r_k[index]

  return(Raw_slice_all_list)
}

## slice and convert to data observed
Data_Slice_Observed <- function(Raw_data_all_list,index){
  Raw_slice_all_list <- Data_Slice_Raw(Raw_data_all_list,index)
  
  res <- Data_Gen_Observed(Raw_slice_all_list)
  return(res)
}
