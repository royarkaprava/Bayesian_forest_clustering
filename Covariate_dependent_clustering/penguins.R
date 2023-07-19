source("Main_functions/forest_class.R")



require("palmerpenguins")



table(palmerpenguins::penguins_raw$Species)




data<- palmerpenguins::penguins



true_membership_full<- data$species

yx_full<- cbind(data$bill_depth_mm,data$bill_length_mm, data$body_mass_g,data$flipper_length_mm)


filter = apply(!is.na(yx_full), 1,any)


yx<- yx_full[filter,]
true_membership_label <- as.factor(true_membership_full[filter])



yx<- t((t(yx) - colMeans(yx))/apply(yx,2,sd))



x<- yx[,3:4]
y<- yx[,1:2]



D<- as.matrix(dist(y,upper=TRUE))



filter1<- rowSums(D==0)==1
y<- y[filter1,]
x<- x[filter1,]
yx<- yx[filter1,]
true_membership_label<- true_membership_label[filter1]



n_x<- nrow(x)
p_x<- ncol(x)



colnames(x)<- c("Body.Mass","Flipper.Length")

colnames(y)<- c("Bill.Depth","Bill.Length")



pairs(yx,col=true_membership_label)



true_membership<- as.numeric(true_membership_label)



# filename= paste("~/dropbox/Apps/Overleaf/spanning_tree_clustering/pdf_figs/penguin_data.pdf")

# pdf(filename,4,4)
plot(y[,1] ,y[,2],  col =  as.factor(true_membership),xlab=colnames(y)[1],ylab=colnames(y)[2])
# dev.off()




require("mclust")



mc_fit<- Mclust(y,G = 3)$classification


# pdf(filename,4,4)
plot(y[,1] ,y[,2], col =  as.factor(mc_fit),xlab=colnames(y)[1],ylab=colnames(y)[2])
# dev.off()




n_x - clusteringAccu(mc_fit,true_membership)*n_x 












fc_fit<-forestClust(y, lam = 0.5,n_iter = 1000,burnin = 500)





C_mat<- getCoAssignmentMat(fc_fit$C)

fc_fit_point_est <- getPointEstC(C_mat,K=3)


n_x - clusteringAccu(fc_fit_point_est,true_membership)*n_x

# filename= paste("~/dropbox/Apps/Overleaf/spanning_tree_clustering/pdf_figs/penguin_fc.pdf")

# pdf(filename,4,4)
plot(y[,1] ,y[,2], col =  as.factor(fc_fit_point_est),xlab=colnames(y)[1],ylab=colnames(y)[2])
# dev.off()

barplot(table(do.call("c",fc_fit$K))/500)





length(true_membership)



Sigma<- t(x)%*%x/n_x

SigmaInv<- solve(Sigma)


logPrior<- matrix(0,n_x+1,n_x+1)

for(i in 1:n_x){
  for(j in 1:i){
    d<- x[i,]-x[j,]
    logPrior[i,j]=    -t(d)%*%SigmaInv%*%d/4
    logPrior[j,i] = logPrior[i,j]
  }
  logPrior[i,n_x+1]<-   - t(x[i,])%*% SigmaInv%*% x[i,]/4
  logPrior[n_x+1,i] <- logPrior[i,n_x+1]
}







fc_fit_using_x<-forestClust(y,lam = 0.5
                            ,n_iter = 1000,burnin = 500,logTreePrior = logPrior)

C_mat2 <- getCoAssignmentMat(fc_fit_using_x$C)


fc_fit_point_est2 <- getPointEstC(C_mat2,K=3)


n_x-clusteringAccu(fc_fit_point_est2,true_membership)*n_x

clusteringAccu(fc_fit_point_est2,true_membership)


# filename= paste("~/dropbox/Apps/Overleaf/spanning_tree_clustering/pdf_figs/penguin_fc_w_x.pdf")

# pdf(filename,4,4)
plot(y[,1] ,y[,2], col =  as.factor(fc_fit_point_est2),xlab=colnames(y)[1],ylab=colnames(y)[2])
# dev.off()

barplot(table(do.call("c",fc_fit_using_x$K))/500)





fc_fit_using_full<-forestClust(cbind(y,x),lam = 0.5
                               ,n_iter = 1000,burnin = 500)



# filename= paste("~/dropbox/Apps/Overleaf/spanning_tree_clustering/pdf_figs/penguin_fc_full.pdf")

C_mat3 <- getCoAssignmentMat(fc_fit_using_full$C)
fc_fit_point_est3 <- getPointEstC(C_mat3,K=3)



# pdf(filename,4,4)
plot(y[,1] ,y[,2], col =  as.factor(fc_fit_point_est3),xlab=colnames(y)[1],ylab=colnames(y)[2])

n_x-clusteringAccu(fc_fit_point_est3,true_membership)*n_x

barplot(table(do.call("c",fc_fit_using_full$K))/500)


# dev.off()



require("dirichletprocess")

dp_model <- DirichletProcessMvnormal(y,alphaPriors = c(1/20, 20))
dp_fit <- Fit(dp_model, 1000,progressBar = TRUE)





dp_cmat<- matrix(0,n_x,n_x)

dp_unique_c<- numeric()

for (i in (501:1000)){
  C = dp_fit$labelsChain[[i]]
  dp_cmat = dp_cmat +  outer(C,C,"==")
  dp_unique_c<- c(dp_unique_c, length(unique(C)))
}




barplot(table(dp_unique_c))



# filename= paste("~/dropbox/Apps/Overleaf/spanning_tree_clustering/pdf_figs/penguin_dpmm.pdf")

fc_fit_point_est4 <- getPointEstC( dp_cmat,K=3)


# pdf(filename,4,4)
plot(y[,1] ,y[,2], col =  as.factor(fc_fit_point_est4),xlab=colnames(y)[1],ylab=colnames(y)[2])
# dev.off()


