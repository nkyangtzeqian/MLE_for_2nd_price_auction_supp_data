### Function to compute Initial values for Lambda and F, with simplified F_SP

## for Lambert W function
# require(pracma)

Initialization_lambda_F_Simp <- function(data, reserve.price.Cutoff){
  ############################################################################################################
  ## Lambda initialization
  ############################################################################################################
  K_r <- which(data$r_k < reserve.price.Cutoff) ## Auctions with reserve price less than reserve.price.Cutoff
  ## K_r should be V(q,e) in the paper
  M_k <- data$M_k[K_r]
  
  inv_g <- function(lambda){
    N_sample=1e4
    N.ran <- rpois(N_sample, lambda * data$tau_a)
    N.ran_thres <- (N.ran>=2)
    k.lam <- 2 * sum( log(N.ran[N.ran_thres] + 0.5) )/ N_sample + 
      2 * (0.5772156649 - 1) * mean(N.ran_thres)
    return(abs(k.lam - mean(M_k)))
  }
  
  lambda_init <- optimize(inv_g, interval = c(0, 10))$minimum
  # print(lambda_init)
  
  ############################################################################################################
  ## F initialization
  ############################################################################################################
  ## Step I
  # G_lambda_r_eta <- function(eta){
  #   if (eta >= 1){
  #     return(1)
  #   } else if (eta <=0){
  #     return(0)
  #   } else{
  #     lambda.tau <- lambda_init * data$tau_a
  #     exp.neg.lambda.tau <- exp(-lambda.tau)
  #     # numerator <- lambda.tau * (1-eta) * ( exp(lambda.tau*eta) - 1) + exp(lambda.tau*eta) - lambda.tau*eta - 1
  #     # numerator <- exp.neg.lambda.tau * numerator
  #     numerator1 <- lambda.tau * (1-eta) * ( exp(-lambda.tau*(1-eta)) - exp.neg.lambda.tau )
  #     numerator2 <- exp(-lambda.tau*(1-eta)) - (lambda.tau*eta*exp.neg.lambda.tau) - exp.neg.lambda.tau
  #     numerator <- numerator1 + numerator2
  #     denominator <- 1 - exp.neg.lambda.tau - lambda.tau * exp.neg.lambda.tau
  #     return(numerator/denominator) 
  #   }
  # }
  trans_Lambert <- function(x){
    ## inverse function of (1+x)exp(-x) where x=lambda*tau*(1-eta)
    ## use the negative branch of Lambert W function, since shoud be decreasing
    return(-1-pracma::lambertWn(-x/exp(1)))
  }
  
  
  # Removing the NA terms and those auction's first jump values whose reserve prices are comparatively higher.
  # We are taking those auctions whose reserve prices are less than an user-specified reserve.price.Cutoff.
  
  X_mk <- data$X_mk[intersect( K_r, which(!is.na(data$X_mk))   )]
  if(length(X_mk)==0){
    stop("Reserve price cutoff is too low for initialization.")
  }
  
  # G_lambda_r_eta_inverse <- inverse(G_lambda_r_eta, lower = 0, upper = 1)
  G_lambda_r_eta_inverse <- function(eta){
    if(length(eta) != 1){
      warning("G_lambda_r_eta_inverse should only be applied to a single value.")
    }
    if (eta >= 1){
      return(1)
    } else if (eta <=0){
      return(0)
    } else{
      return(1-trans_Lambert(eta)/data$tau_a/lambda_init)
    }
  }
  
  
  G_SP <- ecdf(X_mk)
  F_SP <- function(x){
    return(G_lambda_r_eta_inverse(G_SP(x)))
  }
  
  
  ## Step II
  
  # Removing the NA terms and those auction's first jump values whose reserve prices are comparatively higher.
  # We are taking those auctions whose reserve prices are less than an user-specified reserve.price.Cutoff.
  
  X_1 <- data$X_1[intersect( K_r, which(!is.na(data$X_mk))   )]
  G_FP <- ecdf(X_1)
  F_FP <- function(x){
    return(1 - sqrt(1- G_FP(x)))
  }
  
  ## Step III
  
  ## change m1 m2 to p2 p1
  p1 <- max(X_1)
  p2 <- min(X_mk)
  seq_p1p2 <- seq(0, min(p1,p2), by = 0.0001)
  Fc <- sapply(seq_p1p2, F_FP) - F_SP(p2)
  cc <- seq_p1p2[sum(1*(Fc<=0))]
  
  F_0 <- function(x){
    first.term <- F_FP(x) * (x<=cc) + F_SP(x) * (x>p2)
    second.term <- (F_FP(cc) + (  ( (F_SP(p2) - F_FP(cc) ) / (p2 - cc) ) * (x - cc) ) ) * 
      ( ( 1*(x<=p2) ) - ( 1*(x<=cc) ) )
    return(first.term + second.term)
  }
  
  ### deal with numerical issue in inverse function, and avoid log(0) in loglikelihood
  ### actually the Longest Increasing Subsequence algorithm should be applied here. But most of it is increasing
  Fis <- sapply(data$Z_iT_i[ ,1], F_0 )
  x.temp <- data$Z_iT_i[ ,1]
  
  u1 <- max( which(min(data$Xbar_i)== data$Z_iT_i[,1]) )  # index of the min(observed bids) value in the pooled data. u_1
  x.jump <- c(0, x.temp[u1-1], x.temp[diff(c(0,Fis))>0])
  Fis.jump <- c(0, 0, Fis[diff(c(0,Fis))>0])
  
  x.jump[length(x.jump)] <- 1.00001*max(x.temp)
  while (any(diff(c(0,Fis.jump))<0)) {
    warning("F_init is not increasing at some points.")
    x.jump <- c(0, x.temp[u1-1], x.jump[diff(c(0,Fis.jump))>0])
    Fis.jump <- c(0, 0, Fis.jump[diff(c(0,Fis.jump))>0])
  }
  
  Fis.linear <- approx (x.jump, y = Fis.jump, x.temp, method = "linear", yleft=0, yright=.99999)
  
  return(list(lambda = lambda_init, F.x = Fis.linear$x, F.y = Fis.linear$y,
              F_FP= F_FP, F_SP= F_SP  ))
}
