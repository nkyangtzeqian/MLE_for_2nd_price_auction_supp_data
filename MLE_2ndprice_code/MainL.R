library(GoFKernel) # For the inverse function of a CDF


### Function to compute Initial values for Lambda and F

Initialization_lambda_F <- function(data, reserve.price.Cutoff, lambda_set){
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
    
    lt_init_1 <- optimize(inv_g, interval = c(1e-5, 100*tau_a), M = mean(M_k))$minimum
    lt_init_2 <- optimize(inv_g, interval = c(1e-5, 100*tau_a), M = M_k_tau1)$minimum
    lambda_1_init <- lt_init_2/tau_1
    lambda_2_init <- (lt_init_1 - lt_init_2)/(data$tau_a - tau_1)
    
    lambda_init <- c(lambda_1_init, lambda_2_init)
  }else{
    h_index <- lambda_set$h_index
    Htau <- h_list_gen(h_index)[[2]](data$tau_a)
    lt_init <- optimize(inv_g, interval = c(1e-5, 10*Htau), M = mean(M_k))$minimum
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
  
  
  # Removing the NA terms and those auction's first jump values whose reserve prices are comparatively higher.
  # We are taking those auctions whose reserve prices are less than an user-specified reserve.price.Cutoff.
  
  X_mk <- data$X_mk[intersect( K_r, which(!is.na(data$X_mk))   )]
  if(length(X_mk)==0){
    stop("Reserve price cutoff is too low for initialization.")
  }
  
  G_lambda_r_eta_inverse <- inverse(G_lambda_r_eta, lower = 0, upper = 1)
  
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
  p1 <- max(X_1, na.rm = TRUE)
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
              F_FP= F_FP, F_SP= F_SP  ))
}






## Coordinate wise maximization 

cordwise_step <- function( data, lambdas, theta, lambda_set, START = 2){
  
  Lambda=lambda_list_gen(lambdas, lambda_set)[[2]]
  A_m_T <- data$Z_iT_iA_i[ ,3]-data$Z_iT_iA_i[ ,2] 
  
  # Maximizing theta_i's one by one keeping other parameters fixed
  ellPK <- length(theta)
  
  K_sold <- data$K_sold  # This gives the vector of auction indexes for which the item is sold.
  ell <- length(data$Xbar_i)  # Total Number of observed bids for all the auctions.
  
  S_bar_max <- max(data$auxlist$S_bar)
  u_i <- data$auxlist$u_i
  l_i_3p2 <- data$auxlist$l_i_3p2
  Q_i_card <- data$auxlist$Q_i_card
  
  ## F(x) for x<=min(r_k) is not identifiable
  theta[1] <- 1
  
  ## Step 3 
  for (ii in START:S_bar_max){ 
    
    if(ii < ellPK){
      AA <- prod(theta[1:(ii-1)]) * sum(
        (Lambda(data$Z_iT_iA_i[ii:ellPK,3])-Lambda(A_m_T[ii:ellPK])) * c(1,cumprod( theta[(ii+1):ellPK])) 
        )
    }else if(ii==ellPK){
      AA <- prod(theta[1:(ii-1)]) * ( Lambda(data$Z_iT_iA_i[ellPK,3]) - Lambda(A_m_T[ellPK]) )
    }
    
    if(is.nan(AA)){
      warning(paste("AA", AA))
      break
    }
    
    BB <- Q_i_card[ii] + ( (1*(ell>0)) * (ell - l_i_3p2[ii]) ) ## for case III, automatically 0
    
    if(ii %in% u_i){ ## Formula 3.9, case I
      denomenator <- 2*AA
      numerator <- (AA+BB+1) - sqrt((AA+BB+1)^2 - (4*AA*BB))
      theta[ii] <- numerator/denomenator
      if ((numerator/denomenator) >= 1){
        print(paste("theta_i =", numerator/denomenator))
        theta[ii] <- 1
      }
    }else{ ## Formula 3.10, case II
      theta[ii] <- min(c(1,BB/AA))
    }
    
  } # END of for(ii in 1:ellPK) loop.
  
  ## Formula 3.11, case III
  if(S_bar_max < ellPK){
    theta[(S_bar_max+1):ellPK] <- 0  
  }
  
  return(list(lambdas = lambdas, theta = theta))
} # END of cordwise_step function.




## move to lambda_fun.R
# lambda_MLE_step <- function(data,theta){## (alpha-1)/lambda in Gamma pdf
#   alpha.Gamma <- length(data$K_sold)+length(data$Xbar_i) 
#   lambda.Gamma <- sum(data$Z_iT_i[ ,2] * cumprod(theta))
#   if(lambda.Gamma< 1e-5){
#     warning("some thetas are too small")
#   }
#   return(alpha.Gamma/lambda.Gamma)
# }

## MLE finding function, with lambda updated

MLE_2ndprice <- function(data, lambda.in, theta.in, lambda_set, tol=1e-5){
  
  if(length(data$K_sold)==0){   # if item is not sold in any of the auctions.
    warning("Item is not sold in any of the auctions.")
    theta.in <- rep(0, length(theta.in))
    F.mle <- rep(1,length(theta.in))
    Lik.out <- General_Log_Likelihood_PA(data, lambda.in, theta.in)
    return(list(Lik = Lik.out, theta.mle = theta.in, lambdas = lambda.in,
                F.x = data$Z_iT_iA_i[,1], F.y = F.mle, Lik.path = Lik.out))
  }
  
  Lik.in <- -1e100
  
  Lik.out <- General_Log_Likelihood_PA(data, lambda.in, theta.in, lambda_set)
  
  # print(paste("Lik.out initial =", Lik.out))
  
  Lik.path <- Lik.out
  
  while(Lik.out > Lik.in + tol){
    Lik.in <- Lik.out
    
    par.out <- cordwise_step(data, lambda.in, theta.in, lambda_set)
    
    # print(par.out$lambda)
    # print(par.out$theta)
    
    theta.in <- par.out$theta
    lambda.in <- lambda_MLE_step(data, theta.in, lambda.in, lambda_set)
    
    # print(theta.in)
    
    Lik.out <- General_Log_Likelihood_PA(data, lambda.in, theta.in, lambda_set)
    
    # print(paste("Lik.out =", Lik.out))
    
    if(is.na(Lik.out) ){
      warning(paste("in", Lik.in, "out", Lik.out))
      Lik.out <- -1e100
    }
    Lik.path <- c(Lik.path, Lik.out )
  }
  
  G.mle <- cumprod(theta.in)
  F.mle <- 1 - G.mle
  return(list(Lik = Lik.out, theta.mle = theta.in, lambdas = lambda.in,
              F.x = data$Z_iT_iA_i[,1], F.y = F.mle, Lik.path = Lik.path))
}


