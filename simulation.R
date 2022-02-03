require("clusterSim")

n= 150
data<- shapes.two.moon(numObjects = n/2)

X<- data$data
label<- data$clusters

X<- X+ rnorm(n*2, sd=0.05)

plot(X,col=label)

# standardize each col of X
X<- apply(X, 2,  function(x) (x-mean(x))/sd(x))

source("forestProcessV2.R")

start_time <- Sys.time()
fit <- clusteringFP(X, Total_itr = 200, burn=100)
end_time <- Sys.time()
end_time - start_time

barplot(table(sapply(fit$C,function(x) length(unique(x))))/length(fit$C))

P<- getCoAssignmentProb(fit$C)

library(sommer)
image(P,col=jet.colors(100), zlim=c(0,1))

require("kernlab")

posteriorMode = specc(P,centers=2)

plot(X,col= as.integer(posteriorMode))