require("clusterSim")

n= 400
p<- 2
data<- shapes.two.moon(numObjects = n/2)

y<- data$data
true_membership<- data$clusters

y<- y+ rnorm(n*2, sd=0.1)

plot(y[,1],y[,2],col=true_membership)


source("forestProcess.R")
require("kernlab")
require("mclust")


fc_fit<-forestClust(y,lam = 0.5,n_iter = 1000,burnin = 500)


barplot(table(fc_fit$K)/500)


require("fields")


image.plot(fc_fit$C_mat/500)

point_est_C<- getPointEstC(fc_fit,K=2)

plot(y[,1],y[,2],xlab="",ylab="",col=point_est_C)
