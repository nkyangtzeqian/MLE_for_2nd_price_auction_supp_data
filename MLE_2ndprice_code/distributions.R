## CDF of piecewise uniform distribution

piecewise.unif.cdf <- function(y, para1 = 2, para2 = 4){
  # Considering two uniform distributions: Unif(1, para1] and Unif(para1 + 1, para2] 
  if((para1<=1) || (para2 <= para1+1)) stop("This function works only when para1 > 1 and para2 > para1 + 1, since, by default, we are considering the two uniform distributions to be Unif(1, para1] and Unif(para1 + 1, para2]")
  
  if(y <= 1){
    temp <- 0
  }else if ((1 < y) && (y <= para1)){
    temp <- 0.5 * ((y-1)/(para1-1))
  }else if ((para1 < y) && (y <= (para1 + 1))){
    temp <- 0.5
  }else if (((para1 + 1) < y) && (y <= para2)){
    temp <- 0.5 + 0.5 * ((y - para1 - 1)/ (para2 - para1 - 1) )
  }else{
    temp <- 1
  }
  
  return(temp)
}

## Inverse CDF of piecewise uniform distribution

piecewise.unif.invcdf <- function(u, para1 = 2, para2 = 4){
  if((u<0) || (u >1)){
    stop("This function works only when the first argument u lies within the closed
         interval [0,1]")
  }else if(u == 0.5){
    temp <- para1
  }else{
    temp <- inverse(function(y) {piecewise.unif.cdf(y, para1, para2)},
                    lower = 1, upper = para2)(u)
  }
  return(temp)
}


## Define the general version of 2nd Price likelihood function (considering all the situations) (without any lasso penalty)

General_Log_Likelihood_PA <- function( data, lambda, theta){
  
  # store them in processed.
  # theta is a vector of G(z_i)/G(z_{i-1})
  # The first column of pooled.data is the z vector.
  K_sold <- data$K_sold  # This gives the vector of auction indexes for which the item is sold.
  T_0_sold <- data$T_0[K_sold]
  ell <- length(data$Xbar_i)  # Total Number of observed bids for all the auctions.

  S_bar_max <- max(data$auxlist$S_bar)
  u_i <- data$auxlist$u_i
  l_i_3p2 <- data$auxlist$l_i_3p2
  Q_i_card <- data$auxlist$Q_i_card
  
  # u_i is the vector of ranks of all observed bids in the pooled data Z with the first element being u0 = 0.
  term_without_theta <- ((ell+length(K_sold))*log(lambda)) + (ell*log(2)) + sum(log(T_0_sold))
  # print(paste("term without theta = ", term_without_theta))
  term_with_theta_1st_term <- -(lambda * sum(data$Z_iT_i[ ,2] * cumprod(theta))) + sum( (Q_i_card * log(theta))[1:S_bar_max] )
  # print(paste("term with theta 1st = " ,term_with_theta_1st_term))
  if(ell == 0){
    term_with_theta_2nd_term <- 0
  }else{
    term_with_theta_2nd_term <- sum(log(1 - theta[u_i]) ) + sum( ((ell -  l_i_3p2) * log(theta))[1:S_bar_max] )
  }
  # print(paste("term with theta 2nd = ", term_with_theta_2nd_term))
  # Above if condition takes care of the case when any of the auctions do not have any observed bid i.e. ell = 0 case.
  return(term_without_theta + term_with_theta_1st_term + term_with_theta_2nd_term)
}

## test
# u_i=c(4,6,7,8,9,10,11)
# sapply(1:11, function(kk) {sum(u_i <= kk) })


