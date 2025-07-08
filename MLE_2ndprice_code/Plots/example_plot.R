rm(list=ls())
##New EM
# rep.no<- 50
data <-matrix(0, ncol=4, nrow=0)
para1<- 3
para2<- 100
library(rmutil)

tot.bids<- 100

second.cummax <- function(x, r){
  n <- length(x)
  y <- rep(0,n)
  y[1] <- r
  f <- cummax(x)
  for(i in 2:n){
    if((x[i]< y[i-1]) || (x[i] <r)){
      y[i] <- y[i-1]
    }else if(x[i]>f[i-1]){
      y[i] <- f[i-1]
    }else{
      y[i] <- x[i]
    }
  }
  return(y)
}


tts <- rexp(rate=1, tot.bids)
sum(tts)

yy <- rpareto(tot.bids, m= para1, s=para2)
xx<- second.cummax(yy,2)
inds <- 2:tot.bids
inds<- inds[diff(xx)>0]

inds.bid <- c(1,inds, tot.bids+1)
tts[tot.bids+1] = 1

final.jump.inds <- c(1,inds)

wait.jump.time <- c(diff(cumsum(tts)[inds.bid]))
# wait.jump.time <- diff(cumsum(tts)[inds.bid])
wait.bids <- diff(inds.bid)
temp.d <- cbind(xx[final.jump.inds], wait.jump.time, wait.bids,final.jump.inds)
colnames(temp.d) <- c("rec.value", "wait.time", "wait.bids", "jump.inds")

data <- rbind(data,temp.d)
cbind(yy,xx, tts[-(tot.bids+1)],cumsum(tts[-(tot.bids+1)]))
plot( cumsum(tts[-(tot.bids+1)]), xx, type="l")

time.cum <- cumsum(tts[-(tot.bids+1)])

plot.xx<- xx
plot.yy<- yy
zz<- stepfun(time.cum, c(2,plot.xx), ties = "ordered", right = FALSE)
plot(zz, cex.points= .5, xlim=c(0,sum(tts)), ylim = c(0,max(yy)+1), xlab="Time in minutes", ylab="Asking price in dollars", main="Example bidding data")
for(ii in 1:tot.bids){
  segments(time.cum[ii], 0, time.cum[ii], plot.yy[ii], col="blue")
}
# save.image("bidexample.Rdata")
# 
# 
# load("bidexample.Rdata")



# plot(cumsum(data[,2]),data[,1], xlim=c(0, sum(tts)))

# out.split[[ii]] <- Est.func(temp.d)
# MLE.split[[ii]] <- approxfun(x= out.split[[ii]][,1], y= out.split[[ii]][,2], yleft = 0, yright = 1, ties= mean)

# F.delta.temp <- Del.to.F(temp.d[,3])
#  Delta.split[[ii]] <- approxfun(x= temp.d[,1], y= F.delta.temp, yleft = 0, yright = 1, ties= mean)

# plot(out.split[[ii]], xlim=c(0, 50), ylim=c(0,1))
# lines(1:40, MLE.split[[ii]](1:40))
# lines(1:40, Delta.split[[ii]](1:40), col="red")
# lines(temp.d[,1], ppareto(temp.d[,1], m=para1, s=para2), type="l", col="blue")

############################################################################
## toy examples
############################################################################

rm(list=ls())

second_cummax <- function(x, r){
  if(sum(r <= x) == 0){ ## M=O=1
    # return(NaN)
    return(rep(r, length(x))) # just for plotting purpose
  } else if( sum(r <= x) == 1){ ## M=0,O=1. SHOULD NOT BE DELTETED: first jump
    return(rep(r, length(x)))
  } else if( sum(r<= x) >= 2){ ## M>=1,O=1
    y <- c(r, rep(0,length(x))) # Observed (repeated) 2nd price vector of length = no. of bids
    for (ii in 1:length(x)){
      if(y[ii] >= x[ii]){ ## selling price not changed
        y[ii+1] <- y[ii]
      } else { 
        temp_2nd_max <- sort(c(r,x[1:ii]), ## from unobserved data
                             decreasing = TRUE)[2]
        y[ii+1] <- temp_2nd_max
      }
    }# END of xx
    return(y[-1])
  }
}

# ## lambda=1, tau=10
# tts <- rexp(rate=1, 30)
# U_N <- max(which(cumsum(tts)<10))
# tts <- round(tts[1:U_N],3)
# ## uniform
# xx <- sample(1:25, U_N, replace=TRUE)
# # xx <- sample(1:15, U_N, replace=TRUE)

example_orig=list(auction1=list(r_k=10,xx=c(12,6,9,15,11,19,5,2,23,17),tts=c(1.099,0.264,1.106,0.155,0.184,3.176,0.620,0.268,1.435,0.106)),
                 auction2=list(r_k=5,xx=c(16,1,18,11,16,3,20,2,25,8,14,17,25),tts=c(0.241,0.624,0.819,0.390,1.926,0.723,1.033,0.188,0.994,0.061,1.985,0.197,0.084)),
                 auction3=list(r_k=13,xx=c(8,6,3,3,11,10,6,18,4,3,8,10),tts=c(1.117,0.253,0.462,0.299,0.906,1.836,1.733,0.315,1.281,0.445,0.257,0.163)),
                 auction4=list(r_k=17,xx=c(7,8,4,2,15,12,4,10,3,3,6,15),tts=c(0.090,1.231,0.292,0.827,0.437,0.554,1.604,1.156,1.481,0.357,0.015,1.436)))

pdf("Plots/example_plot.pdf", width=8, height=4)
par(mfrow=c(2,2), mar=c(4,4,2,2), mgp=c(2,0.5,0))

for (i in 1:4) {
  xx <- example_orig[[i]]$xx
  yy <- second_cummax(xx, example_orig[[i]]$r_k)
  time.cum <- cumsum(example_orig[[i]]$tts)
  plot.xx<- xx
  plot.yy<- yy
  zz<- stepfun(time.cum, c(example_orig[[i]]$r_k,plot.yy), ties = "ordered", right = FALSE)
  plot(zz, cex.points= .5, xlim=c(0,10), ylim = c(0,26), xlab="Time in days", ylab="Asking price in dollars", main=paste0("Example bidding data in Acution ",i))
  for(ii in 1:length(time.cum)){
    segments(time.cum[ii], 0, time.cum[ii], plot.xx[ii], col="blue")
  }
}
dev.off()          