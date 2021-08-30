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

data <- read.csv("\\\\file.phhp.ufl.edu\\home\\ark007\\My Documents\\Spectral clustering\\data.txt", header = F)

loc <- as.matrix(data[, 1:2])
label <- data$V3
label[which(label=="Class 1")] <- 1
label[which(label=="Class 2")] <- 2
label <- as.numeric(label)

ind1 <- which(label==1)
ind2 <- which(label==2)
locR <- loc[c(sample(ind1,100), sample(ind2, 100)), ]
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

#I am initializing using ComputeMST(), there may be other ways too.
fit <- clusteringFP(X, alpha0 = 0.001, a0=0.1, b0 = 0.1, K= 20, Total_itr = 100, burn=50)

##Not run always
rowMeans(fit$clssize_ls) #estimated class sizes
#fit$estiadja #estimated adjacency
rowMeans(fit$stickbrkwts) #Estimated components for stick breaking
