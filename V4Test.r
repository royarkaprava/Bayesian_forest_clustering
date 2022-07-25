rm(list=ls())

n <- 200
p <- 100

y<- matrix(rt(n*p,df = 5),n,p)

y[101:n,]<- y[101:n,]+ 2

source("forestProcessV4.r")

true_membership<- c(rep(1,100),rep(2,100))



plot(y[,1],y[,2])

fc_fit<-forestClust(y)

barplot(table(fc_fit$K))


require("kernlab")
require("mclust")

image(fc_fit$C_mat)


sc_fit<- as.numeric(specc(as.kernelMatrix(fc_fit$C_mat+0.00001),centers=2))


adjustedRandIndex(sc_fit,true_membership)

# G_D<-graph_from_adjacency_matrix(A_T_mat, mode = "undirected",weighted = TRUE,diag=FALSE)
# G_D$layout=y[,1:2]
# V(G_D)$color = as.numeric(sc_fit)
# 
# G_D<- delete.edges(G_D, E(G_D)[E(G_D)$weight< quantile(E(G_D)$weight,0.8)])
# 
# plot.igraph(G_D,vertex.size=5, vertex.label='')

fit_mclust<- Mclust(y,G = 2,modelNames = "VVV")$classification
adjustedRandIndex(fit_mclust,true_membership)

fit_mclust<- Mclust(y,G = 2,modelNames = "EEI")$classification

adjustedRandIndex(fit_mclust,true_membership)
