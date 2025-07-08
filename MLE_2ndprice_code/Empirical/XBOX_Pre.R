
## Function to extract organized data from a real (raw) Xbox Auctions data set.

Xbox_data_extract <- function(Raw_unorg_data){
  colnames(Raw_unorg_data) <- c("auctionid", "bid", "bidtime", "bidder",
                                      "bidderrate", "openbid", "price")
  
  # Adding small random noise to all the bids to avoid equality among any two bid values.
  Raw_unorg_data$bid <- Raw_unorg_data$bid + runif(length(Raw_unorg_data$bid), min = 0, max = 0.01)
  # Raw_unorg_data$bid <- Raw_unorg_data$bid + 0.01  # Adding 1 cent to all the bids to avoid matching with the reserve price.
  
  tau_a <- ceiling(max(Raw_unorg_data$bidtime))
  unique.auction.data <- unique(cbind(Raw_unorg_data$auctionid, Raw_unorg_data$openbid))
  Ka <- nrow(unique.auction.data)
  
  N_all<- rep(0, Ka)  # Vector of total number of bids of all the auctions.
  
  r_k <- unique.auction.data[,2]
  if(is.null(r_k)){ r_k <- rep(0, Ka)}
  
  Raw_biddata_list <- vector("list", Ka) # List containing each auction data separately.
  
  for(jj in 1:Ka){
    indexjj <- which(Raw_unorg_data$auctionid == unique.auction.data[jj,1])
    t_otts <- Raw_unorg_data$bidtime[indexjj]
    N_all[jj] <- length(indexjj)   # Total number of bids for the jth auction.
    yy <-  Raw_unorg_data$bid[indexjj]  # All the private values of the bidders for the jth process
    
    Raw_biddata_matrix <- matrix(0, nrow = N_all[jj], ncol = 4)
    colnames(Raw_biddata_matrix) <- c("auction_index", "bid_value", "bid_time_cum", "openbid_N_all")
    Raw_biddata_matrix[,1] <- rep(jj, N_all[jj])
    Raw_biddata_matrix[,2] <- yy
    Raw_biddata_matrix[,3] <- t_otts
    Raw_biddata_matrix[,4] <- rep(r_k[jj], N_all[jj])
    
    Raw_biddata_list[[jj]] <- Raw_biddata_matrix
  }# END of for(jj in 1: Ka)
  
  Raw.data.list <- list(Raw_biddata_list = Raw_biddata_list,
                        Ka = Ka,
                        tau_a = tau_a,
                        r_k = r_k,
                        True.Method = "unknown")
  Raw.data.list$class = "SecondPriceAuction_Rawdata"
  
  return(Raw.data.list)
}

raw_XBOX_7day <- read.csv("Empirical/Xbox_7day_auctions.csv", header = T)
raw_data_list <- Xbox_data_extract(raw_XBOX_7day)

save(raw_data_list, file = "Empirical/Xbox_7day_auctions.RData")