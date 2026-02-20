source("DataGeneration.R")
## Functions to generate the data where $Ka$ is not fixed but the total time
## window is fixed.
## Raw means the data is not observed.

Data_Gen_Raw_Gamma21 <- function(Ka, tau_a, lambda0, r_k,
                         method = c("unif", "pareto", "gamma", "beta",
                                    "piecewise_unif", "customized"),
                         para1 = NULL, para2 = NULL, customized_CDF = NULL){
  N_all <- rep(0, Ka) # Total number of potential bids for the jth auction, might be smaller than r. Not observed.
  if(is.null(r_k)) r_k <- rep(0, Ka)
  Raw_biddata_list <- vector("list", Ka) # List containing each auction data separately.
  
  for (jj in 1: Ka) { ## for each auctions
    tts <- rgamma( (2* tau_a), shape=2, rate = 1) ## why not use while loop
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

Data_Gen_Raw_LNormal <- function(Ka, tau_a, lambda0, r_k,
                                 method = c("unif", "pareto", "gamma", "beta",
                                            "piecewise_unif", "customized"),
                                 para1 = NULL, para2 = NULL, customized_CDF = NULL){
  N_all <- rep(0, Ka) # Total number of potential bids for the jth auction, might be smaller than r. Not observed.
  if(is.null(r_k)) r_k <- rep(0, Ka)
  Raw_biddata_list <- vector("list", Ka) # List containing each auction data separately.
  
  for (jj in 1: Ka) { ## for each auctions
    tts <- exp(rnorm( (2* tau_a), mean=0, sd = 1)) ## why not use while loop
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

List_Data_Gen_Raw <- list(
  Data_Gen_Raw_Poisson = Data_Gen_Raw,
  Data_Gen_Raw_Gamma21 = Data_Gen_Raw_Gamma21,
  Data_Gen_Raw_LNormal = Data_Gen_Raw_LNormal
)

rm(Data_Gen_Raw)
rm(Data_Gen_Raw_Gamma21)
rm(Data_Gen_Raw_LNormal)