source("forest_class.R")



sep_b =3

n <- 400
p <- 2

# y<- matrix(rnorm(n*p),n,p)
y<- matrix(rt(n*p,df = 5),n,p)


y[201:n,]<- y[201:n,]+ sep_b


true_membership<- c(rep(1,200),rep(3,200))


plot(y[,1],y[,2],xlab="",ylab="",col=true_membership+1)


require("kernlab")
require("mclust")


fc_fit<-forestClust(y,lam = 0.5,n_iter = 1000,burnin = 500)



#uncertainty estimate

C_mat<- getCoAssignmentMat(fc_fit$C)

require("fields")
image.plot(C_mat)

#point estimate
point_est_C<- getPointEstC(C_mat,K=2)
plot(y[,1],y[,2],xlab="",ylab="",col=point_est_C)

# number of clusters
barplot(table(do.call("c",fc_fit$K))/500)

