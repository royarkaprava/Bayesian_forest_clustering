library(emstreeR)
library(Rcpp)
library(RcppArmadillo)
library(bmixture)

#Dirichlet concentration alpha0
#Inverse Gamma shape param a0
#Inverse Gamma scale param b0
#Maximum no. of cluster K

clusteringFP <- function(X, alpha0=0.5, a0=0.1, b0=0.1, K=20, Total_itr = 10000, burn=5000, gamma =1){
  #K=20
  
  Ugamma <- function(gamma){
    fitK <- sum(B[1, ]!=0)
    ret <- -(p*K+2)*log(gamma)-1/gamma 
    
    return(ret)
  }
  
  getA <- function(B){
    A = - B %*% t(B)
    diag(A) <- 0
    
    return(A) 
  }
  #initialization
  
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
  
  
  DisMat2 <- as.matrix(dist(X)^2)
  
  disX <- array(dist(X)^2)
  if(sum(disX==0)){
    disX <- disX[-which(disX==0)] 
  }
  
  
  lam   <- min(disX) #variance for other conditional edges

  
  b0 <- min(disX/p) #I will check dist(X)
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
  clslbp <- matrix(0, n, Total_itr-burn)#clslb
  
  lam_ls<- numeric(length = Total_itr-burn )
  
  for(k in clsno){
    temp <- which(B[-1, k] != 0)
    indcompo <- which(pathmat[temp, ] < Inf)
    
    clslb[indcompo] <- l
    l <- l + 1
  }
  edges <- out[n:1, 3:4]
  
  #Total_itr <- 5000
  Ap <- matrix(0, n+1, n+1)
  Bp <- matrix(0, n+1, n)
  B2inv <- solve(crossprod(B))
  
  itr = 0
  
  B0 =B
  A0=A
  
  classMat <- matrix(0, n, K)
  
  classMat[cbind(1:n, clslb[1:n])] <- rep(1, n)
  
  clsmem <- rep(0, K)
  clsmemp <- matrix(0, K, Total_itr-burn)#clsmem
  for(k in 1:K){
    clsmem[k] <- length(which(clslb==k))
  }
  #alpha0 <- 0.5
  beta  <- rbeta(K-1, 1, alpha0)
  wl <- c(beta, 1)*c(1, exp(cumsum(log(1-beta))))
  
  #wl <- alpha #rdirichlet(1, alpha)#rep(alpha, K)
  
  d <- rep(0, p)
  
  for(i in 1:p){
    d[i] <- (max(X[,i]) - min(X[,i]))/2 #max(max(abs(X[,i]))-mean(X[,i]), mean(X[,i])-min(abs(X[,i])))
  }
  
  # gamma <- 1
  
  gamsd <- 1e-2
  argam <- 0
  
  betap <- matrix(0, K-1, Total_itr-burn)#rep(0, K-1)
  
  pb <- txtProgressBar(min = itr, max = Total_itr, style = 3)
  while(itr < Total_itr){
    itr <- itr + 1
    
    #update B
    #beta = rep(0, 401)
    #j=0; k=0 
    #P=rep(1, 201)
    #N <- rep(1,200)
    b <- rep(1,401)
    b1 <- b
    #update B, here 'b' is useless. I added for debugging.
    BupC(B, B2inv, clslb, clsmem, b1,b, wl, lam, 2^p*gamma^p*prod(d), Xmu)
    A <- getA(B)
    
    ####Update lam
    
    
    # ind0 <- which(B[1, ]!=0)
    # sum1 <- 0
    # for(i in 1:p){
    #   if(length(ind0) < (n-1)){
    #     sum1 <- sum1 + sum((colSums(B[-1,-ind0]*X[, i]))^2) 
    #   }
    #   if(length(ind0) == (n-1)){
    #     sum1 <- sum1 + sum(((B[-1,-ind0]*X[, i]))^2) 
    #   }
    # }
    # ap = a0 + n/2 - length(ind0) / 2
    # bp = b0 + sum1/2
    # 
    # lam = 1/rgamma(1, ap, bp)
    

    A_exclude_root = A[-1,]
    A_exclude_root=  A_exclude_root[,-1]
    A_exclude_root[upper.tri(A_exclude_root)]=0
    ap  = a0 + p * sum(A_exclude_root) /2
    bp  = b0 + sum(DisMat2*A_exclude_root) /2
    
    lam = 1/rgamma(1, ap, bp)
    
    
    
    #update wl
    for(i in 1:(K-1)){
      alpha.k <- 1 + clsmem[i]
      beta.k  <- sum(clsmem[-(1:i)]) + alpha0
      
      beta[i] <- rbeta(1, alpha.k, beta.k)
    }
    
    wl <- c(beta, 1)*c(1, exp(cumsum(log(1-beta))))
    
    #wl <- rdirichlet(1, rep(alpha, K) + clsmem)
    
    #beta of stick-breaking
    
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
    
    
    if(itr > burn){
      Ap <- Ap+A 
      Bp <- Bp + B
      clslbp[, itr-burn]  <- clslb  #Storing class labels from postburn samples
      clsmemp[, itr-burn] <- clsmem #Storing class sizes from postburn samples
      betap[, itr - burn] <- beta   #Storing the K-1 elements of stick-breaking prior after burn-in
      
      lam_ls[itr - burn] <- lam   #Storing the K-1 elements of stick-breaking prior after burn-in
      
    }
    
    # if(itr %% 100 == 0){
    #   ar <- argam/ itr
    #   if(ar<.20){sdgam <- sdgam * (.5)}
    #   if(ar>.40){sdgam <- sdgam * (2)}
    # }
    #update B2inv
    #print(itr)
    #print(tau)
    Sys.sleep(0.1)
    # update progress bar
    setTxtProgressBar(pb, itr)
  }
  close(pb)
  Ap <- Ap / (Total_itr - burn)
  Bp <- Bp / (Total_itr - burn)
  out <- list(estiadja = Ap, estiB = Bp, clslb_ls = clslbp, clssize_ls = clsmemp, stickbrkwts = betap, lam_ls = lam_ls)
  return(out)
}