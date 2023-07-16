source("forestProcess.r")

load(file="./yaleB10subjects.Rda")

y<- yaleB10subjects$y


true_membership<- as.numeric(yaleB10subjects$label)

n<- nrow(y)
p<- ncol(y)


require("kernlab")
require("mclust")

require("ClusterR")

res<- sapply(c(1:5), function(i){

    km_fit<- kmeans(y,centers = 10,)$cluster

    clusteringAccu(km_fit,true_membership)
    
    })

mean(res)
quantile(res, c(0.025,0.975))


# get sparse representation z_i for each y_i
require("glmnet")

W<- matrix(0,n,n)

for(i in c(1:n)){
    glmnet_fit<- glmnet(t(y[-i,]),y[i,],family = "gaussian",alpha=1,lambda=(1/p),standardize = FALSE,intercept=FALSE)
    W[i,-i]<- as.vector(glmnet_fit$beta*1)
    if(i%%50==0)
    print(i)
    flush.console()
}

Z<- t(apply(W,1, function(x){
    t<- sort(x,decreasing = TRUE)[3]
    x*(x>=t)
    }))

diag(Z)=1

Z2<- (Z%*%t(Z))
Z2<- t(Z2/sqrt(diag(Z2)))/sqrt(diag(Z2))

Z<- Z/sqrt(diag(abs(Z%*%t(Z))))


fixedLogS<- matrix(0,n+1,n+1)
fixedLogS[1:n,1:n]<-  10*Z2 + log(Z2!=0)


fc_fit<- forestClust(Z,useFixedLogS = TRUE, fixedLogS = fixedLogS,n_iter = 500,burnin = 0)

barplot(table(fc_fit$K))

res<- sapply(c(1:20), function(i){
    fc_fit_point_est <- getPointEstC(fc_fit,K=10)
    clusteringAccu(fc_fit_point_est,true_membership)
    })

mean(res)
quantile(res, c(0.025,0.975))

require("fields")

image.plot(fc_fit$C_mat/500, zlim=c(0,1))



D<- as.matrix(dist(y,diag = TRUE,upper = TRUE))

image.plot(D)
