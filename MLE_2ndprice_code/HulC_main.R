#------------------------------------------------------------------------------------
## HulC1d() uses asymptotic median bias value to construct
## 		convex hull confidence interval for a univariate 
##		parameter. This is Algorithm 1 of the paper.
## data is a data frame.
## estimate is a function that takes a data frame as input 
## 		and returns a one-dimensional estimate.
## alpha is the level.
## Delta is the median bias of the estimate(). 
## randomize is a logical. If TRUE then the number of splits
## 		is randomized. If FALSE, then the larger number of
##		splits is used.
HulC1d_Ebay <- function(raw_data_all_list, eval.points, alpha = 0.1, MLE.Cutoff = 0, Delta=0, diff_bias= NULL){
  if(is.null(diff_bias)) stop("different in bias needs to be input based on the table that computes this quantity.")
  # if(diff_bias==0) warning("diff_bias is zero. Are you sure? This is okay if you are using this code to compute the median bias")
  
  # Init.Values <- Initialization_lambda_F(data = Data_Gen_Observed(raw_data_all_list), 
  #                                           reserve.price.Cutoff = MLE.Cutoff)
  # lambda.initial <- Init.Values$lambda
  
  randomize = TRUE
  # alpha = 0.05; Delta = 0; randomize = TRUE
  nn <- length(raw_data_all_list$Raw_biddata_list)
  # Delta <- 0 ## Median Bias is estimated easily as asymptotic distribution is distribution free.

  corrected.alpha <- alpha/length(eval.points) #Bonferroni Correction
  B1 <- solve_for_B(alpha = corrected.alpha, Delta = Delta, t = 0)
  B <- B1
  if(randomize){
    p1 <- (1/2 + Delta)^B1 + (1/2 - Delta)^B1
    B0 <- B1 - 1
    p0 <- (1/2 + Delta)^B0 + (1/2 - Delta)^B0
    U <- runif(1)
    tau <- (corrected.alpha - p1)/(p0 - p1)
    B <- B0*(U <= tau)+ B1*(U > tau)
  }
  if(B > nn){
    print(paste0("Delta = ", Delta, ", No. of splits = ", B, ", Sample size = ", nn))
    stop("Error: not enough samples for splitting!")
  }
  MLE_est <- matrix(0,nrow=length(eval.points),ncol= B)
  Init_est <- matrix(0,nrow=length(eval.points),ncol= B)
  #Right now we are using the same split but we will later use different splits.
  TMP <- split(sample(nn), sort((1:nn)%%B))
  # idx <- 1
  for(idx in 1:B){
    bat.idx <- sort(TMP[[idx]])
    
    
    bat_data <- Data_Slice_Observed(raw_data_all_list,bat.idx)
    
    # F.local<- MLE.from.raw.data.2ndprice.HulC(bat_data, MLE.Cutoff)
    
    Init.Values <- Initialization_lambda_F(bat_data, MLE.Cutoff)
    theta.init <- (1-Init.Values$F.y)/ c(1,(1-Init.Values$F.y)[-length(Init.Values$F.y)])
    MLE <- MLE_2ndprice(bat_data, lambda.in = Init.Values$lambda, theta.in = theta.init)
    
    MLE_est[,idx] <- approx(x = MLE$F.x, y = MLE$F.y, 
                            xout = eval.points, method = "linear", yleft=0, yright=1)$y
    Init_est[,idx] <- approx(x = Init.Values$F.x, y = Init.Values$F.y, 
                             xout = eval.points, method = "linear", yleft=0, yright=1)$y
    ## Why do I need a linear interpolation?
  }
  CI.Init <- apply(Init_est, 1, range)
  CI.MLE <- apply(MLE_est, 1, range)
  rownames(CI.MLE)<- rownames(CI.Init)<- c("lwr", "upr")
  ## WHY?
  CI.Init[1,] <- 	pmax(CI.Init[1,] - diff_bias,0) # expansion by log(2)/m (split size) to account for binomial dist
  CI.Init[2,] <- 	pmin(CI.Init[2,] + diff_bias, 1) # expansion to account for binomial dist
  CI.MLE[1,] <- 	pmax(CI.MLE[1,] - diff_bias, 0) # expansion by log(2)/m (split size) to account for binomial dist
  CI.MLE[2,] <- 	pmin(CI.MLE[2,] + diff_bias, 1) # expansion to account for binomial dist
  # we can use the monotonicity of $F_0$ to improve provide a confidence band for the whole function, 
  #   suppose $\ell(t_1) \le f_0(t_1) \le u(t_1)$ and $\ell(t_2) \le f_0(t_2) \le u(t_2)$ for two points $t_1 \le t_2\in[0,1]$ in the domain, then using the information $F_0(\cdot)$ is increasing and between 0 and 1, we can conclude that $\overline{\ell}(t) \le f_0(t) \le \overline{u}(t)$ for all $t\in[0, 1]$, where
  #     \begin{equation*}
  #   \overline{\ell}(t) = \begin{cases}
  #   0;&  \mbox{for }t \le t_1,\\
  #   \ell(t_1) &;\mbox{for }t_1 \le t \le t_2,\\
  #   \ell(t_2)&;\mbox{for }t_2 \le t \le 1,\end{cases}
  #   \quad\mbox{and}\quad 
  #   \overline{u}(t) = \begin{cases}u(t_1),&;\mbox{for }t \le t_1,\\
  #   u(t_2),&;\mbox{for }t_1 &l(t); t \le t_2,\\
  #   1, &;\mbox{for }t_2 &l(t); t \le 1.\end{cases}
  #   \label{eq:interval-to-band-monotone}
  #   \end{equation*}
  
  
  extra.info <- list(MLE_est,Init_est) #all the estimators based on the data splits.
  
  MakeBand <- function(temp){
    rval<- NULL
    rval$lwr <- approxfun(eval.points, cummax(temp[1,]), 
                          method = "constant", 
                          yleft = 0, yright = max(temp[1,]), 
                          f = 0, ties = "ordered")
    class(rval$lwr) <- c( "stepfun")
    
    rval$upr <- approxfun(eval.points, rev(cummin(rev(temp[2,]))), 
                          method = "constant", 
                          yleft =  min(temp[2,]), yright = 1, 
                          f = 1, ties = "ordered")
    class(rval$upr) <- c( "stepfun")
    return(rval)
  }
  Band.CI.Init <- MakeBand(CI.Init)
  Band.CI.MLE <- MakeBand(CI.MLE)
  ret <- list(CI.Init = CI.Init, CI.MLE = CI.MLE, B = B, Band.CI.Init = Band.CI.Init, 
              Band.CI.MLE = Band.CI.MLE, extra.info = extra.info)
  return(ret)
}

## For non-negative Delta and t, set
## 	Q(B; Delta, t) = [(1/2 - Delta)^B + (1/2 + Delta)^B]*(1 + t)^{-B+1}
## The following function finds the smallest B for a given t such that
##	Q(B; Delta, t) <= alpha.
solve_for_B <- function(alpha, Delta, t){
  if(Delta == 0.5 && t == 0){
    stop("Delta is 0.5 and t = 0. The estimator lies only on one side of the parameter!")
  }
  B_low <- max(floor(log((2 + 2*t)/alpha, base = 2 + 2*t)), floor(log((1 + t)/alpha, base = (2 + 2*t)/(1 + 2*Delta))))
  B_up <- ceiling(log((2 + 2*t)/alpha, base = (2 + 2*t)/(1 + 2*Delta)))
  Q <- function(B){
    ((1/2 - Delta)^B + (1/2 + Delta)^B)*(1 + t)^(-B + 1)
  }
  for(B in B_low:B_up){
    if(Q(B) <= alpha)
      break
  }
  return(B)
}

# ## For any estimation function estimate() that returns a
# ## univariate estimator, subsamp_median_bias() provides an
# ## estimate of the median bias using subsampling.
# ## The subsample size used is (sample size)^{subsamp_exp}.
# ## The input data is a data frame or a matrix.
# ## nsub is the number of subsamples
# subsamp_median_bias <- function(data, estimate, subsamp_exp = 2/3, nsub = 1000){
#   data <- as.matrix(data)
#   nn <- nrow(data)
#   subsamp_size <- round(nn^subsamp_exp)
#   nsub <- min(nsub, choose(nn, subsamp_size))
#   fulldata_estimate <- estimate(data)
#   Delta <- 0
#   for(b in 1:nsub){
#     TMP <- estimate(data[sample(nn, subsamp_size, replace = FALSE),,drop=FALSE])
#     Delta <- Delta + (TMP - fulldata_estimate <= 0)/nsub
#   }
#   Delta <- abs(Delta - 1/2)
#   return(Delta)
# }
