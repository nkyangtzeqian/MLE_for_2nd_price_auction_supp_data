

## data oberseved
Data_Gen_Unobserved_PolyaTree <- function(Raw_data_all_list){
  if(Raw_data_all_list$class != "SecondPriceAuction_Rawdata"){
    stop("This fucntion works only with data in the class 'SecondPriceAuction.Rawdata'.
         See e.g., 'Data_Gen_Raw' where we generate such data.")
  } 
  
  
  Ka <- Raw_data_all_list$Ka
  dat <- matrix(NA, nrow=Ka, ncol=3)
  
  for (jj in 1: Ka) {
    temp <- Raw_data_all_list$Raw_biddata_list[[jj]][,2]
    ## number of bids
    dat[jj,1] <- length(temp)
    ## 2nd highest
    dat[jj,2] <- sort(temp, decreasing=TRUE)[2]
    ## highest
    dat[jj,3] <- max(temp)
    
  }# END of for(jj in 1: Ka)

  dat <- dat[order(dat[,2], decreasing=TRUE),]
  
  ret <- list(dat=dat)
  
  ret$class = "SecondPriceAuction_PolyaTree"
  
  ret$Ka = Ka
  return(ret)
}


## Wrap up orginal function
MCMC_PT <- function(ret, F.x, x.max){
  if(ret$class != "SecondPriceAuction_PolyaTree"){
    stop("This fucntion works only with data in the class 'SecondPriceAuction_PolyaTree.Rawdata'.
         See e.g., 'Data_Gen_Unobserved_PolyaTree' where we generate such data.")
  } 
  dat <- ret$dat
  Ka <- ret$Ka
  
  # Set "weakly informative" alpha parameters
  alpha.prior <- matrix(0, nrow = Ka, ncol = 2)					
  # x.max <- 20
  alpha.prior[1,] <- c(dat[1,2], x.max - dat[1,2]) / x.max
  for (k in 2:Ka){ 
    alpha.prior[k,] <- (c(dat[k,2], dat[k-1,2] - dat[k,2]) / dat[k-1,2]) * k^2
  }
  alpha.prior <- alpha.prior * exp(-20)
  
  
  result.PT <- EstCDF.PT(dat, alpha.prior, x.max=x.max)
  
  F.y <- sapply(F.x, function(x) result.PT$cdf(x))
  return(list(MCMC_Mat=result.PT$predprob.mat,
              fit_CDF=result.PT$cdf,fit_PDF=result.PT$den,
              F.x = F.x, F.y = F.y))
}

#################################################################################
## Original codes
#################################################################################


# Estimation of CDF Using Polya tree method 
# --------------------------------------------------------------------------------
EstCDF.PT <- function(dat, alpha.prior, x.max, nsim=2000, plotind=FALSE){
  result <- list(cutpt=numeric(1), predprob=numeric(1), predprob.mat=numeric(1), alpha.mat=numeric(1), cdf=function(x){ return(x)}, den=function(x){return(x)});
  result$predprob.mat <- matrix(nrow = nsim, ncol = nrow(dat)+1)
  result$alpha.mat <- matrix(nrow = nsim, ncol = length(c(alpha.prior)))
  
  result$cutpt <- c(0, rev(dat[,2]), x.max)
  Cs <- matrix(0.5, nrow = nrow(alpha.prior), ncol = 2)	# Initialize Cs
  
  # Gibbs sampler for Cs and Z
  for (i in 1:nsim){
    # Draw Z given Cs
    Z <- draw.Z(Cs)
    
    # Draw Cs given Z
    alpha.post <- update.alpha(alpha.prior, dat, Z)
    Cs <- draw.Cs(alpha.post)
    
    # Find P(Bs)
    pred <- rev(pred.dist(Cs))
    
    # Record both P(Bs) and alpha.post of the current draw
    result$predprob.mat[i,] <- pred		# P(Bs)
    result$alpha.mat[i,] <- c(alpha.post)
    
    # if (i %% 100==0) {cat("Iteration ",i," done\n");}
  }
  result$predprob.mat <- result$predprob.mat[-(1:floor(nsim/2)),]
  result$predprob <- colMeans(result$predprob.mat)
  result$cdf <- approxfun(result$cutpt, c(0, cumsum(result$predprob)))
  result$den <- stepfun(result$cutpt[c(-1,-length(result$cutpt))], result$predprob / diff(result$cutpt), f=0)
  
  return(result)
}

# Draw conditional probabilities Cs given the alpha parameters
draw.Cs <- function(alpha){
  Cs <- matrix(nrow = nrow(alpha), ncol = 2)
  Cs[,1] <- rbeta(nrow(alpha), alpha[,1], alpha[,2])
  Cs[,2] <- 1 - Cs[,1]
  return(Cs)
}	

# Draw Z given conditional probabilities Cs
draw.Z <- function(Cs){
  M <- nrow(Cs) ## ADDED BY Z.YANG
  Z.new <- numeric(nrow(Cs))
  baseprob <- numeric(nrow(Cs))
  baseprob[1] <- (1-Cs[1,1])
  for (i in 2:M){
    baseprob[i] <- prod(Cs[1:(i-1),1]) * (1-Cs[i,1])
  }
  for (i in 1:M){
    Z.new[i] <- sample(1:i, 1, prob=baseprob[1:i]/sum(baseprob[1:i]))
  }
  return(Z.new)
}

# Update alpha given N and z
update.alpha <- function(alpha.prior, dat, Z){
  M <- nrow(alpha.prior) ## ADDED BY Z.YANG
  N <- dat[,1] 			# Number of bidders in each auction
  alpha.new <- alpha.prior
  for (i in 1:M){
    alpha.new[i,1] <- alpha.prior[i,1] + sum(N[i:M] - 1) + length(Z[Z>=(i+1)]) - 1 
    alpha.new[i,2] <- alpha.prior[i,2] + length(Z[Z==i]) + 1 
  }
  return(alpha.new)
}

# Compute predictive distribution given Cs
pred.dist <- function(Cs){
  M <- nrow(Cs)
  baseprob <- numeric(M)
  baseprob[1] <- (1-Cs[1,1])
  for (i in 2:M){
    baseprob[i] <- prod(Cs[1:(i-1),1]) * (1-Cs[i,1])
  }
  baseprob <- c(baseprob, 1-sum(baseprob))
  return(baseprob)
}


# Likelihood functions
LL1.PT <- function(cutpt, predprob, n1, y1){
  index <- max((1:length(cutpt))[y1>cutpt])
  fy <- predprob[index] / (cutpt[index+1]-cutpt[index])
  if (index>1){
    Fy <- sum(predprob[1:(index-1)]) + predprob[index] * (y1-cutpt[index]) / (cutpt[index+1] - cutpt[index])
  } else {
    Fy <- predprob[index] * (y1-cutpt[index]) / (cutpt[index+1] - cutpt[index])
  }
  return(log(n1*(n1-1)) + log(1-Fy) + log(fy) + (n1-2)*log(Fy))
}

LL.PT <- function(cutpt, predprob, n, y){
  k <- length(y)
  loglik <- numeric(k)
  
  for (i in 1:k){
    loglik[i] <- LL1.PT(cutpt, predprob, n[i], y[i])
  }
  return(sum(loglik))
}


# # Compute posterior interval for profit
# # ------------------------------------------
# post.int <- function(price, result.PT, cost = 5.2){
#   Niter <- nrow(result.PT$predprob.mat)
#   profit <- matrix(nrow = Niter, ncol = length(price));
#   
#   for (i in 1:Niter){
#     cur.cdf <- approxfun(result.PT$cutpt, c(0, cumsum(result.PT$predprob.mat[i,])))
#     profit[i,] <- Profit.Func(x=price, cdf=cur.cdf, cost=cost)
#   }
#   return(profit)
# }
