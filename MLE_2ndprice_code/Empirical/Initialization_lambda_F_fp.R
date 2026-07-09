### Function to compute Initial values for Lambda and F

Initialization_lambda_F_fp <- function(data, reserve.price.Cutoff, lambda_set, MAX_L_EST=1000){
  ## MAX_L_EST: avoid too large lambda estimate when necessary, which may cause numerical issue and slow convergence in the optimization.
  ############################################################################################################
  ## Lambda initialization
  ############################################################################################################
  K_r <- which(data$r_k < reserve.price.Cutoff) ## Auctions with reserve price less than reserve.price.Cutoff
  ## K_r should be V(q,e) in the paper
  M_k <- data$M_k[K_r]
  
  ## general inv_g function to find lambda_inits
  inv_g <- function(lt, M){
    N_sample=1e4
    N.ran <- rpois(N_sample, lt)
    N.ran_thres <- (N.ran>=2)
    k.lam <- 2 * sum( log(N.ran[N.ran_thres] + 0.5) )/ N_sample + 
      2 * (0.5772156649 - 1) * mean(N.ran_thres)
    return(abs(k.lam - M))
  }
  
  if(lambda_set$TWO_STAGE){
    tau_1 <- lambda_set$tau_1
    M_k_tau1 <- mean(sapply(K_r, function(k){
      sum( cumsum(c(data$T_0[k],data$data_list_Obs[[k]]$t_i_p))<=tau_1 ) -1 # -1 remove reserve price.
    }))  ## mean of M of K_r auctions with tau_a <= tau_1
    
    lt_init_1 <- optimize(inv_g, interval = c(1e-5, MAX_L_EST*tau_a), M = mean(M_k))$minimum
    lt_init_2 <- optimize(inv_g, interval = c(1e-5, MAX_L_EST*tau_a), M = M_k_tau1)$minimum
    lambda_1_init <- lt_init_2/tau_1
    lambda_2_init <- (lt_init_1 - lt_init_2)/(data$tau_a - tau_1)
    
    lambda_init <- c(lambda_1_init, lambda_2_init)
  }else{
    h_index <- lambda_set$h_index
    Htau <- h_list_gen(h_index)[[2]](data$tau_a)
    lt_init <- optimize(inv_g, interval = c(1e-5, MAX_L_EST*Htau), M = mean(M_k))$minimum
    lambda_init <- lt_init/Htau
  }
  # lambda_init <- optimize(inv_g, interval = c(0, 10))$minimum
  # print(lambda_init)
  
  ############################################################################################################
  ## F initialization
  ############################################################################################################
  ## Step I
  G_lambda_r_eta <- function(eta){
    if (eta >= 1){
      return(1)
    } else if (eta <=0){
      return(0)
    } else{
      lambda.tau <- lambda_list_gen(lambda_init, lambda_set)[[2]](data$tau_a)
      exp.neg.lambda.tau <- exp(-lambda.tau)
      # numerator <- lambda.tau * (1-eta) * ( exp(lambda.tau*eta) - 1) + exp(lambda.tau*eta) - lambda.tau*eta - 1
      # numerator <- exp.neg.lambda.tau * numerator
      numerator1 <- lambda.tau * (1-eta) * ( exp(-lambda.tau*(1-eta)) - exp.neg.lambda.tau )
      numerator2 <- exp(-lambda.tau*(1-eta)) - (lambda.tau*eta*exp.neg.lambda.tau) - exp.neg.lambda.tau
      numerator <- numerator1 + numerator2
      denominator <- 1 - exp.neg.lambda.tau - lambda.tau * exp.neg.lambda.tau
      return(numerator/denominator) 
    }
  }
  
  

  
  
  ## Step II
  
  # Removing the NA terms and those auction's first jump values whose reserve prices are comparatively higher.
  # We are taking those auctions whose reserve prices are less than an user-specified reserve.price.Cutoff.
  
  X_1 <- data$X_1[intersect( K_r, which(!is.na(data$X_mk))   )]
  G_FP <- ecdf(X_1)
  F_FP <- function(x){
    return(1 - sqrt(1- G_FP(x)))
  }
  

  
  F_0 <- F_FP
  
  ### deal with numerical issue in inverse function, and avoid log(0) in loglikelihood
  ### actually the Longest Increasing Subsequence algorithm should be applied here. But most of it is increasing
  Fis <- sapply(data$Z_iT_iA_i[ ,1], F_0 )
  x.temp <- data$Z_iT_iA_i[ ,1]
  
  u1 <- max( which(min(data$Xbar_i)== data$Z_iT_iA_i[,1]) )  # index of the min(observed bids) value in the pooled data. u_1
  x.jump <- c(0, x.temp[u1-1], x.temp[diff(c(0,Fis))>0])
  Fis.jump <- c(0, 0, Fis[diff(c(0,Fis))>0])
  
  x.jump[length(x.jump)] <- 1.00001*max(x.temp)
  while (any(diff(c(0,Fis.jump))<0)) {
    warning("F_init is not increasing at some points.")
    x.jump <- c(0, x.temp[u1-1], x.jump[diff(c(0,Fis.jump))>0])
    Fis.jump <- c(0, 0, Fis.jump[diff(c(0,Fis.jump))>0])
  }
  
  Fis.linear <- approx (x.jump, y = Fis.jump, x.temp, method = "linear", yleft=0, yright=.99999)
  
  return(list(lambdas = lambda_init, F.x = Fis.linear$x, F.y = Fis.linear$y,
              F_FP= F_FP, F_SP= NA  ))
}
