# n <- 1000
# 
# tseq <- seq(0, 3, length.out = n)
# Yeu    <- matrix(0, 2, n)
# 
# for(i in 1:n){
#   Yeu[1, i] <- integrate(f=function(x){(cos(x^2))}, 0, tseq[i])$value
#   Yeu[2, i] <- integrate(f=function(x){(sin(x^2))}, 0, tseq[i])$value
# }
# 
# data <- Yeu + matrix(rnorm(2*n, sd = 0.01), 2, n)
# plot(t(data))


setwd("~/git/Spectral-clustering/")
Rcpp::sourceCpp('spectralV1.cpp')

data <- read.csv("./data.txt", header = F)

loc <- as.matrix(data[, 1:2])
label <- data$V3
label[which(label=="Class 1")] <- 1
label[which(label=="Class 2")] <- 2
label <- as.numeric(label)
label

ind1 <- which(label==1)
ind2 <- which(label==2)
idx <- c(sample(ind1,100), sample(ind2, 100))

locR <- loc[idx, ]
label_set = label[idx]
X    <- locR

#plot(X)

# K <- 4
# 
# p     <- ncol(X)
# n     <- nrow(X)
# Xdis  <- as.matrix(exp(-dist(X)/10))
# diag(Xdis) <- 1


######################################################################################
source('spectralV2_withstickbreaking.R')
# source('spectralV2_withQuasiBernoulli.R')

#I am initializing using ComputeMST(), there may be other ways too.

a0 =  1
b0 =  1E-3

fit <- clusteringFP(X,alpha0 =  0.1, a0= a0, b0 = b0, K= 20, Total_itr = 1000, burn=100, gamma=100000000)

ts.plot(fit$lam_ls)
acf(fit$lam_ls,lag.max = 40)

min(dist(X)[dist(X)>0])

rowMeans(fit$clssize_ls) #estimated class sizes
#fit$estiadja #estimated adjacency
rowMeans(fit$stickbrkwts) #Estimated components for stick breaking


barplot(table(apply(fit$clslb_ls,2, function(x)length(unique(x)))))

image(fit$estiadja,col = topo.colors(100, rev=F))

plot(fit$estiadja[1,])
# fit$estiadja

plot(X[,1],X[,2], col = fit$clslb_ls[,100])


