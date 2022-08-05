source("forestProcess.r")



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



barplot(table(fc_fit$K)/500)


require("fields")


image.plot(fc_fit$C_mat/500)

point_est_C<- getPointEstC(fc_fit,K=2)

plot(y[,1],y[,2],xlab="",ylab="",col=point_est_C)

