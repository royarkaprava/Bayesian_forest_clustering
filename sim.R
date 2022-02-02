if(Sys.info()["sysname"] %in% c("Linux","Darwin") ){
  setwd("~/git/Spectral-clustering/")
}

if(Sys.info()["sysname"] %in% c("Windows") ){
  setwd("\\\\file.phhp.ufl.edu\\home\\ark007\\My Documents\\GitHub\\Spectral-clustering")
}

require("clusterSim")

data<- shapes.two.moon(numObjects = 150)

X<- data$data
label<- data$clusters

plot(X,col=label)

# standardize each col of X
X<- apply(X, 2,  function(x) (x-mean(x))/sd(x))

######################################################################################

source('forestProcess.R')
fit <- clusteringFP(X, Total_itr = 100, burn=50, gamma =1, lambda=1, random_scan_n = 0)


plot(X,col=(fit$C[[1]])[-1])
plot(X,col=(fit$C[[20]])[-1])
plot(X,col=(fit$C[[50]])[-1])


trace_sigma2<- do.call('c',fit$Sig)
acf(trace_sigma2,lag.max = 40)
ts.plot(trace_sigma2)

fit$C