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
locR <- loc[c(sample(ind1,200), sample(ind2, 200)), ]
X    <- locR

plot(X)

# K <- 4
# 
# p     <- ncol(X)
# n     <- nrow(X)
# Xdis  <- as.matrix(exp(-dist(X)/10))
# diag(Xdis) <- 1

library(emstreeR)
library(Rcpp)
library(RcppArmadillo)
library(bmixture)

#clusteringFP <- function(X, K=20, Total_itr = 10000){

Ugamma <- function(gamma){
  fitK <- sum(B[1, ]!=0)
  ret <- -(p*K+2)*log(gamma)-1/gamma 
  
  return(ret)
}

#initialization

lam   <- 1 #variance for other conditional edges
#muprvar <- 0.1 # prior variance for \mu

p     <- ncol(X)
n     <- nrow(X)
# Xdis  <- as.matrix(exp(-dist(X)/10))
# #diag(Xdis) <- 1
# 
# d     <- rowSums(Xdis)
# LaplX <- diag(n) - diag(1/sqrt(d)) %*% Xdis %*% diag(1/sqrt(d))
# Lapei <- eigen(LaplX)$vectors
# Lapei <- Lapei[, n:(n-K+1)]
# sepcl <- kmeans(Lapei, centers = K)

mu    <- c(median(X[, 1]), median(X[, 2]))
Xmu   <- rbind(mu, X)

muprmn  <- mu

d    <- data.frame(Xmu)
out  <- ComputeMST(d)

disX <- array(dist(X)^2)
if(sum(disX==0)){
  disX <- disX[-which(disX==0)] 
}

b0 <- (disX/p) #I will check dist(X)
a0 <- 2

B <- matrix(0, n+1, n)

for(i in n:1){
  B[out[i, 3], i] <- 1
  B[out[i, 4], i] <- -1
}

A <- getA(B)

G <- igraph::graph_from_adjacency_matrix(as.matrix(A[-1,-1]),mode = "undirected")
pathmat <- igraph::shortest.paths(G)

path <- NetworkToolbox::pathlengths(A[-1,-1]) #shortest path length

clsno <- which(B[1, ] != 0)
l <- 1
clslb <- rep(0, n)
for(k in clsno){
  temp <- which(B[-1, k] != 0)
  indcompo <- which(pathmat[temp, ] < Inf)
  
  clslb[indcompo] <- l
  l <- l + 1
}
edges <- out[n:1, 3:4]

Total_itr <- 5000
Ap <- matrix(0, n+1, n+1)
Bp <- matrix(0, n+1, n)
B2inv <- solve(crossprod(B))

itr = 0

B0 =B
A0=A

classMat <- matrix(0, n, K)

classMat[cbind(1:n, clslb[1:n])] <- rep(1, n)

clsmem <- rep(0, K)
for(k in 1:K){
  clsmem[k] <- length(which(clslb==k))
}
alpha = 1/K
wl <- rep(alpha, K)

d <- rep(0, p)

for(i in 1:p){
  d[i] <- max(max(abs(X[,i]))-mean(X[,i]), mean(X[,i])-min(abs(X[,i])))
}

gamma <- 1

gamsd <- 1e-2
argam <- 0

while(itr < Total_itr){
  itr <- itr + 1
  
  #update B
  #beta = rep(0, 401)
  #j=0; k=0 
  # P=1
  # N <- rep(1,200)
  # b <- rep(1, 401)
  
  #update B, here 'b' is useless. I added for debugging.
  BupC(B, B2inv, clslb, clsmem, b, wl, lam, 2^p*gamma^p*prod(d), Xmu)
  A <- getA(B)
  
  ####Update lam
  ind0 <- which(B[1, ]!=0)
  sum1 <- 0
  for(i in 1:p){
    if(length(ind0) < (n-1)){
      sum1 <- sum1 + sum((colSums(B[-1,-ind0]*X[, i]))^2) 
    }
    if(length(ind0) == (n-1)){
      sum1 <- sum1 + sum(((B[-1,-ind0]*X[, i]))^2) 
    }
  }
  ap = a0 + n/2 - length(ind0) / 2
  bp = b0 + sum1/2
  
  lam = 1/rgamma(1, ap, bp)
  
  #update wl
  wl <- rdirichlet(1, rep(alpha, K) + clsmem)
  
  #update gamma using MH
  # temp <- log(sqrt(gamma)) + rnorm(1, sd = gamsd)
  # gammac <- exp(temp^2)
  # #gammac <- max(1.0001, gammac)
  # 
  # D    <- Ugamma(gammac) - Ugamma(gamma)
  # if(is.nan(D)){D=-Inf}
  # if(is.na(D)){D=-Inf}
  # 
  # if(D > log(runif(1))){
  #   argam <- argam + 1
  #   gamma <- gammac
  # }
  
  
  if(itr > 2000){
    Ap <- Ap+A 
    Bp <- Bp + B
  }
  
  # if(itr %% 100 == 0){
  #   ar <- argam/ itr
  #   if(ar<.20){sdgam <- sdgam * (.5)}
  #   if(ar>.40){sdgam <- sdgam * (2)}
  # }
  #update B2inv
  print(itr)
  #print(tau)
}
#}