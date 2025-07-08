# library(transport) # For calculating Wasserstein distance between two distribution functions.
library(distr)
library(distrEx)   # For calculating Total variation distance between two distributions.


## Function for the linear interpolation cdf used in AbscontDistribution class (in distrEX package):
linear_interpolation_cdf <- function(x, F.x)
{
  xu <- c(0,x)
  F.xu <- c(0,F.x)
  l <- length(xu)
  val <- function(t, lower.tail = T) {
    if(t <= 0){ ## Should be <= correspond to the latter changes. 04/14/25
      F.t <- 0
    }else{
      i <- max(which((xu-t)<0))  # Sourav: It shouldn't have "=" in it. Removed it now. Thanks, Rohit da, for pointing it out.
      if(i==l){
        F.t <- 0.9999999
      }else{
        F.t <- F.xu[i] + (F.xu[i+1] - F.xu[i])*(t - xu[i])/(xu[i+1] - xu[i])
      }
    }
    
    return(F.t)
  }
  
  val.vec <- function(t.vec, lower.tail = T){
    return(sapply(t.vec, FUN = val))
  }
  
  return(val.vec)
}


# Function to calculate KS distance (maximum vertical distance) between two cdf vectors F1 and F2:
KS.d <- function(F1, F2){
  return( max(abs(F1 - F2)) )
}
