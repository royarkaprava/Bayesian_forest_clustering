library(emstreeR)
library(Rcpp)
library(RcppArmadillo)
library(bmixture)

#Dirichlet concentration alpha0
#Inverse Gamma shape param a0
#Inverse Gamma scale param b0
#Maximum no. of cluster K


# alpha0=0.5
# a0=0.1
# b0=0.1
# K=20
# Total_itr = 100
# burn=50

clusteringFP <- function(X, p_b=0.1, a0=0.1, b0=0.1, K=20, Total_itr = 10000, burn=5000, gamma =1, random_scan_n = 0){
  #K=20
  
  
  if(random_scan_n==0){
    random_scan_n= nrow(X)
  }
  
  Ugamma <- function(gamma){
    fitK <- sum(B[1, ]!=0)
    ret <- -(p*K+2)*log(gamma)-1/gamma 
    
    return(ret)
  }
  
  
  ###### functions needed for using Quasi-Bernoulli stick-breaking
  # sample from truncated beta supported in (a,b)
  rtbeta <- function(n, alpha, beta, a = 0, b = 1) {
    stopifnot(n > 0 & all(beta > 0) & all(alpha > 0))
    x <- runif(n)
    Fa <- pbeta(a, alpha, beta)
    Fb <- pbeta(b, alpha, beta)
    y <- (1 - x) * Fa + x * Fb
    y[y < 1E-16] = 1E-16
    return(qbeta(y, alpha, beta))
  }
  
  
  updateBeta <- function(n_C, b, alpha_beta = 1, eps = 1E-3, d_beta = 0) {
    
    # n_C <- colSums(C)
    par1_beta = (N - cumsum(n_C)) + alpha_beta + d_beta * c(1:K)
    par2_beta = n_C + 1 - d_beta
    
    beta_1 = rbeta(K, par1_beta, par2_beta)
    if (p_b == 1) {
      beta = beta_1
      if (any(is.na(beta))) {
        # print("beta err")
        beta = updateBeta(n_C, b, alpha_beta, eps, d_beta)
      }
      return(beta)
    }
    beta_2 = rtbeta(K, par1_beta, par2_beta, 0, eps) / eps # truncated Beta distribution
    
    beta = beta_1 * (b == 1) + beta_2 * (b != 1)
    # beta[K]=1E-8
    if (any(is.na(beta))) {
      # print("beta err")
      beta = updateBeta(n_C, b, alpha_beta, eps, d_beta)
    }
    return(beta)
  }
  
  
  updateB <- function(n_C, eps = 1E-3, p_b = 0.9, alpha_beta = 1) {
    
    # n_C <- colSums(C)
    m_C = N - cumsum(n_C)
    
    choice1 = log(p_b)
    choice2 = log(1 - p_b) + pbeta(eps, alpha_beta + m_C, n_C + 1, log.p = T) - alpha_beta * log(eps)
    
    gumbel = -log(-log(runif(K * 2, min = 0, max = 1)))
    
    prob_choice = cbind(choice1, choice2) + gumbel
    b <- colSums((apply(prob_choice, 1, function(x) x == max(x))) * c(1, eps))
    return(b)
  }
  
  updateW <- function(b, beta) {
    v = 1 - b * beta
    w = v * (cumprod(c(1, 1 - v))[1:K])
    return(w)
  }
  
  
  #################
  
  
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
  
  # mu    <- c(median(X[, 1]), median(X[, 2]))
  mu <- apply(X, 2, median)
  Xmu   <- rbind(mu, X)
  
  muprmn  <- mu
  
  d    <- data.frame(Xmu)
  out  <- ComputeMST(d,verbose = FALSE)
  
  DisMat2f <- as.matrix(dist(Xmu)^2)
  DisMat2  <- DisMat2f[-1,-1]
  
  disX <- array(dist(X)^2)
  if(sum(disX==0)){
    disX <- disX[-which(disX==0)] 
  }
  
  sig   <- min(disX) #variance for other conditional edges
  
  # b0 <- min(disX/p) #I will check dist(X)
  # a0 <- 2
  
  B <- matrix(0, n+1, n)
  
  for(i in n:1){
    B[out$from[i], i] <- 1
    B[out$to[i], i] <- -1
  }
  
  A <- getA(B)
  
  G <- igraph::graph_from_adjacency_matrix(as.matrix(A[-1,-1]),mode = "undirected")
  pathmat <- igraph::shortest.paths(G)
  
  path <- NetworkToolbox::pathlengths(A[-1,-1]) #shortest path length
  
  clsno <- which(B[1, ] != 0)
  l <- 1
  clslb <- rep(0, n)
  clslbp <- matrix(0, n, Total_itr-burn)#clslb
  
  sig_ls<- numeric(length = Total_itr-burn )
  
  for(k in clsno){
    temp <- which(B[-1, k] != 0)
    indcompo <- which(pathmat[temp, ] < Inf)
    
    clslb[indcompo] <- l
    l <- l + 1
  }
  edges <- out[n:1, 3:4]
  
  
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
  alpha0 <- 0.5
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
  
  betap <- matrix(0, K, Total_itr-burn)#rep(0, K-1)
  
  pb <- txtProgressBar(min = itr, max = Total_itr, style = 3)
  
  itr=0
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
    
    
    
    BupC(B, B2inv, clslb, clsmem, b1,b, wl, sig, gamma, Xmu, random_scan_n)
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
    
    #A_root = A[1,]
    #A_root=  A_root[,1]
    #A_root[upper.tri(A_exclude_root)]=0
    
    #ap = ap + p*sum(A_root) / 2
    #bp = bp + sum(DisMat2f[1,]*A_root) /(2*kappa)
    
    sig = 1/rgamma(1, ap, bp)
    
    
    
    #update wl
    
    n_C <- clsmem
    N <- n
    p_b = 0.1
    d_beta = 0
    alpha_beta = 1
    eps = 1E-3
    
    b <- updateB(n_C, eps = eps, p_b = p_b, alpha_beta = alpha_beta)
    beta <- updateBeta(n_C, b, alpha_beta = alpha_beta, d_beta = d_beta, eps = eps)
    wl <- updateW(b, beta)
    
    
    # for(i in 1:(K-1)){
    #   alpha.k <- 1 + clsmem[i]
    #   beta.k  <- sum(clsmem[-(1:i)]) + alpha0
    #   
    #   beta[i] <- rbeta(1, alpha.k, beta.k)
    # }
    # 
    # wl <- c(beta, 1)*c(1, exp(cumsum(log(1-beta))))
    
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
      
      sig_ls[itr - burn] <- sig   #Storing the K-1 elements of stick-breaking prior after burn-in
      
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
    # print(n_C)
    
  }
  
  
  close(pb)
  Ap <- Ap / (Total_itr - burn)
  Bp <- Bp / (Total_itr - burn)
  out <- list(estiadja = Ap, estiB = Bp, clslb_ls = clslbp, clssize_ls = clsmemp, stickbrkwts = betap, sig_ls = sig_ls)
  return(out)
}