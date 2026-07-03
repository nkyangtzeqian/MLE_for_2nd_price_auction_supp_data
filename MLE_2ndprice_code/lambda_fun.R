h_list_gen <- function(h_index){
  force(h_index)
  ## h_list: h and H functions for one stage
  ## h_index=1: constant; h_index=2: proportional ; h_index=3: 
  if(h_index==1){ return(list(function(t){1}, function(t){t})) }
  if(h_index==2){ return(list(
    function(t){1+t/100}, 
    function(t){t + t^2/200}
    )) }
}

H_inv <- function(h_index){
  force(h_index)
  ## h_index=1: constant; h_index=2: proportional ; h_index=3: 
  if(h_index==1){ return(function(t){t}) }
  if(h_index==2){ return(function(t) {100*(-1+sqrt(1+2*t/100))}) }
} 

## definition of lambda and Lambda functions for different h_index and two-stage options
lambda_list_gen <- function(lambdas,lambda_set){
  ## TWO_STAGE = T: ignore the h_index and use the two-step function for lambda
  TWO_STAGE <- lambda_set$TWO_STAGE
  lambda1 <- lambdas[1] ## force lambda
  
  if(TWO_STAGE){
    lambda2 <- lambdas[2]
    tau_1 <- lambda_set$tau_1
    lambda <- function(t){ifelse(t<=tau_1, lambda1, lambda2)}
    Lambda <- function(t){ifelse(t<=tau_1, lambda1*t, lambda1*tau_1 + lambda2*(t-tau_1))}
  }else if(lambda_set$h_index==1){ ## accelerate the reference of function
    lambda <- function(t) sapply(t, function(x) lambda1 * 1)
    Lambda <- function(t) sapply(t, function(x) lambda1 * x)
  }else{
    h_list <- h_list_gen(lambda_set$h_index)
    lambda <- function(t) sapply(t, function(x) lambda1 * h_list[[1]](x))
    Lambda <- function(t) sapply(t, function(x) lambda1 * h_list[[2]](x))
  }
  return(list(lambda, Lambda)) ## list of two functions
}

## only used for initialization of data generation
lambda0_list_init <-  function(lambda1,lambda2=NULL,lambda_set){
  lambdas <- c(lambda1,lambda2)
  lambda0_list_full <- lambda_list_gen(lambdas,lambda_set)
  
  ## TWO_STAGE = T: ignore the h_index and use the two-step function for lambda
  TWO_STAGE <- lambda_set$TWO_STAGE
  
  if(TWO_STAGE){
    tau_1 <- lambda_set$tau_1
    Lambda_inv <- function(t){ifelse(t<=lambda1*tau_1, t/lambda1, tau_1 + (t-lambda1*tau_1)/lambda2)}
  }else{
    h_index <- lambda_set$h_index
    Lambda_inv <- function(xx) sapply(xx, function(x) H_inv(h_index)(x/lambda1))
  }
  lambda0_list_full$Lambda_inv <- Lambda_inv
  return(lambda0_list_full) ## return list of function
}

###########################################################################################



## wrap-up for all cases
lambda_MLE_step <- function(data,theta,lambda_old,lambda_set){
  ## lambda_list_old: a list of lambda and Lambda functions from the previous iteration for two-step H
  TWO_STAGE <- lambda_set$TWO_STAGE
  
  A_m_T <- data$Z_iT_iA_i[ ,3]-data$Z_iT_iA_i[ ,2] 

  # one stage with known h case
  lambda_MLE_step_1 <- function(data,theta,h_index){## (alpha-1)/lambda in Gamma pdf
    if(h_index==1) H <- function(t){t}
    else H <-  h_list_gen(h_index)[[2]]
    
    alpha.Gamma <- length(data$K_sold)+length(data$Xbar_i)
    lambda.Gamma <- sum(cumprod(theta) * ( H(data$Z_iT_iA_i[ ,3]) - H(A_m_T) ) )
    if(lambda.Gamma< 1e-5){
      warning("some thetas are too small")
    }
    return(alpha.Gamma/lambda.Gamma)
  }
  
  # two stage case
  lambda_MLE_step_2 <- function(data,theta,tau_1,lambda_list_old){
    K_sold_tau <- (data$T_0 <= tau_1) & (data$M_k>0)
    K_sold_min_tau <- (data$T_0 > tau_1) & (data$M_k>0)
    T_1 <- (data$Z_iT_iA_i[ ,3] <= tau_1)
    T_2 <- ( (data$Z_iT_iA_i[ ,3] > tau_1) & (A_m_T <= tau_1) )
    T_3 <- (A_m_T > tau_1)
    
    Coe1 <- sum(T_1) + sum(T_2) + sum(K_sold_tau) - data$Ka
    Coe2 <- sum(T_3)
    Coe3 <- sum(data$Z_iT_iA_i[ ,2] * cumprod(theta)*T_1) + 
      sum((tau_1-A_m_T) * cumprod(theta)*T_2)
    Coe4 <- sum(data$Z_iT_iA_i[ ,2] * cumprod(theta)*T_3) +
      sum((data$Z_iT_iA_i[ ,3]-tau_1) * cumprod(theta)*T_2)
    if(length(K_sold_min_tau)==0){ ## all first price changes happens in first stage; quite common
      lambda1 <- Coe1/Coe3
      lambda2 <- Coe2/Coe4
    }else{
      warning("some first price changes happen in the second stage; use numerical optimization to find MLE")
      fr <- function(par) {
        l1 <- par[1]; l2 <- par[2]
        logLik <- Coe1*log(l1)+Coe2*log(l2)-l1*Coe3-l2*Coe4+sum(log(l2*(data$T_0[K_sold_min_tau]-tau_1)+l1*tau_1))
        -logLik
      }
      gr <- function(par) {
        l1 <- par[1]; l2 <- par[2]
        Der <- c(Coe1/l1 - Coe3 + sum(tau_1/(l2*(data$T_0[K_sold_min_tau]-tau_1)+l1*tau_1)),
                 Coe2/l2 - Coe4 + sum((data$T_0[K_sold_min_tau]-tau_1)/(l2*(data$T_0[K_sold_min_tau]-tau_1)+l1*tau_1)))
        -Der
      }
      ## starting point
      l1_old <- lambda_old[1]
      l2_old <- lambda_old[2]
      opim_res <- optim(par=c(l1_old,l2_old), fn=fr,gr=gr, method="L-BFGS-B", lower=c(1e-5, 1e-5))
      
      lambda1 <- opim_res$par[1]
      lambda2 <- opim_res$par[2]
    }
    return(c(lambda1,lambda2))
  }
  
  ##############################################################################
  if(TWO_STAGE){
    tau_1 <- lambda_set$tau_1
    lambdas <- lambda_MLE_step_2(data,theta,tau_1,lambda_list_old)
    lambda1 <- lambdas[1]
    lambda2 <- lambdas[2]
    res <- c(lambda1,lambda2)
  }else{
    res <- lambda_MLE_step_1(data,theta,lambda_set$h_index)
  }
  return(res) 
}
