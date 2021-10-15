#Maximum no. of cluster K


# Find the minimum spanning tree, return the incidence matrix
FindMST <- function(DisMat2){
  n<- nrow(DisMat2)
  B<- matrix(0,n,n-1)
  idx1<- c(1)
  idx2<- c(2:n)
  
  for (l in 1:(n-1)){
    Dsubset <- DisMat2[idx1,idx2]
    idx1mat <- matrix(idx1,length(idx1),length(idx2))
    idx2mat <- matrix(idx2,length(idx1),length(idx2),byrow = T)
    
    picked_i <- idx1mat[Dsubset == min(Dsubset)]
    picked_j <- idx2mat[Dsubset == min(Dsubset)]
    
    B[picked_i,l]= 1
    B[picked_j,l]= -1
    
    idx1<- c(idx1,picked_j)
    idx2<- setdiff(c(1:n),idx1)
  }
  B
}

# function to produce the adjacency matrix from incidence matrix
getA <- function(B){
  A = - B %*% t(B)
  diag(A) <- 0
  return(A) 
}


# function to get submatrix of B2Inv
getSubInverse<- function(B2Inv, k,full_edge_idx){
  sel <- full_edge_idx == k
  sel_not <- full_edge_idx != k
  
  M11 = B2Inv[sel_not, sel_not]
  M12 = B2Inv[sel_not, sel]
  M22 = B2Inv[sel, sel]
  
  return (M11 - M12%*%t(M12)/M22)
}




# log density for f
logden_f <- function(d2, sigma2, p=p){
  - 0.5* d2/sigma2 - 0.5* p*log(sigma2*2*pi) 
}

# log density for r
logden_r <- function(d2, gamma=1, p=p){
  - (1+p)/2* log(1+d2/(gamma^2)) + lgamma((1+p)/2) - (1+p)/2* log(pi) - p *log(gamma)
}


clusteringFP <- function(X, Total_itr = 10000, burn=5000, gamma =1, lambda=1, random_scan_n = 0){
  
  listA<- list()
  listB<- list()
  listC<- list()
  listSig<- list()
  
  # random_scan: randomly pick a subset of edges and update in an iteration  
  if(random_scan_n==0){
    random_scan_n= nrow(X)
  }
  
  
  p     <- ncol(X)
  n     <- nrow(X)
  
  # use mean
  mu <- apply(X, 2, mean)
  
  
  # Compute the MST among y's
  DisMat2 <- as.matrix(dist(X,"euclidean", diag = TRUE, upper = TRUE)^2)
  dist2Xmu<-  colSums((t(X)-mu)**2)
  
  B_without_node0 <- FindMST(DisMat2)
  A_without_node0 <- getA(B_without_node0)
  
  edgeDist2 = (DisMat2*A_without_node0) [lower.tri(DisMat2)]
  edgeDist2 = edgeDist2[edgeDist2>0]
  
  #use the median dist2 on the MST to specify b0_sigma
  b0_sigma = median(edgeDist2)
  a0_sigma = p+1
  
  # initialize sigma2 as b0_sigma/p
  sig   <- b0_sigma/p
  
  # Initialize the B matrix by attaching node 0 to node 1
  # that is, we start with one cluster
  
  B <- cbind( c(1,-1,rep(0,n-1)), rbind(0,B_without_node0))
  A <- getA(B)
  
  
  # initialize the clustering membership
  C <- rep(1,n+1) 
  C[1]<- 0
  
  # S
  
  Sr <- logden_r(dist2Xmu, gamma=1,p=p)
  
  Sf <- logden_f(d2 = DisMat2, sigma2 = sig, p=p)
  S<- rbind(c(0,Sr), cbind(Sr,Sf))
  diag(S)<- 0
  
  ##############
  
  full_edge_idx = c(1:n)
  
  B2Inv <- solve(t(B)%*%B)
  
  itr = 0
  pb <- txtProgressBar(min = itr, max = Total_itr, style = 3)
  
  for (itr in c(1:Total_itr)){
    
    # update every edge
    
    for(k in c(1:n)){
      sel <- full_edge_idx == k
      sel_not <- full_edge_idx != k
      
      B_not_s = B[, sel_not, drop = F]
      
      B_not_s2_inv = getSubInverse(B2Inv, k, full_edge_idx)
      
      beta_k = c(B[,sel] - B_not_s%*%(B_not_s2_inv%*% (t(B_not_s)%*%B[, sel])))
      mid_point = (max(beta_k)+min(beta_k))/2
      
      Ga<- beta_k > mid_point
      # check if Ga includes 0: if not, flip the sign
      if(!Ga[1]){
        Ga <- !Ga
      }
      # the Gb set
      Gb<- !Ga
      
      # the number of clusters within Ga
      
      if(length(C[Ga])>1){
        Kstar = length(unique((C[Ga])[-1]))
      }else{
        Kstar = 1
      }
      
      SGaGb<- S[Ga,Gb]
      
      if(sum(Ga)==1){
        SGaGb<- matrix(SGaGb,nrow=1,ncol=length(SGaGb))
      }
      if(sum(Gb)==1){
        SGaGb<- matrix(SGaGb,ncol=1,nrow=length(SGaGb))
      }
      
      # adjust for the change in Poisson prior if attaching Gb to node 0 (idx 1 in R program here)
      SGaGb[1,] <- SGaGb[1,] + log(lambda) - log(Kstar+1) 
      
      ########
      rgumbel<- -log(-log(runif(length(SGaGb))))
      SGaGbGumbel<- SGaGb+ rgumbel
      
      GaIdx<- c(1:(n+1))[Ga]
      GbIdx<- c(1:(n+1))[Gb]
      
      GaIdx1<- rep(GaIdx,length(GbIdx))
      GbIdx1<- rep(GbIdx,each= length(GaIdx))
      
      picked_i <- GaIdx1[which.max(SGaGbGumbel)]
      picked_j <- GbIdx1[which.max(SGaGbGumbel)]
      
      b_vec<- B[,k]
      B[,k]<- 0
      B[picked_i,k]<- 1
      B[picked_j,k]<- -1
      
      if(sum(abs(B[,k]- b_vec))>0){ # if we change B[,k], update B2Inv and C
        B_s_star = B[, sel]
        b = t(B_not_s)%*%B_s_star  # complexity: O(p^2)
        B_invB_b = B_not_s2_inv%*%b  # complexity: O(p^2)
        M22_star = c(1./(2.0 - t(b)%*%(B_invB_b)))
        M12_star = - B_invB_b * (M22_star)
        
        M11_star = B_not_s2_inv+M12_star%*% t(M12_star)/M22_star
        
        temp = B2Inv[sel_not,]
        temp[, sel_not] = M11_star
        temp[, sel] = M12_star
        B2Inv[sel_not,] = temp
        
        temp = matrix(B2Inv[sel,],nrow = 1)
        temp[,sel_not] = t(M12_star)
        temp[, sel] = M22_star
        B2Inv[sel,] = temp
        
        
        if(picked_i==1){ # add a cluster
          C[Gb]<- min(setdiff(c(0:(n)),C[Ga]))
        }else{
          C[Gb]<- C[picked_i]
        }
      }
    }
    
    A<- getA(B)
    
    
    # update sigma2
    A_without_node0<- A[-1,-1]
    A_without_node0[upper.tri(A_without_node0)]<-0
    
    ap  = a0_sigma + p*sum(A_without_node0)/2
    bp  = b0_sigma + sum(DisMat2*A_without_node0) /2
    sig = 1/rgamma(1, ap, bp)
    
    Sf <- logden_f(d2 = DisMat2, sigma2 = sig, p=p)
    S<- rbind(c(0,Sr), cbind(Sr,Sf))
    diag(S)<- 0
    
    
    
    if(itr > burn){
      listA[[itr-burn]] <- A 
      listB[[itr-burn]] <- B 
      listC[[itr-burn]] <- C 
      listSig[[itr-burn]] <- sig
    }
    
    Sys.sleep(0.1)
    # update progress bar
    setTxtProgressBar(pb, itr)
  }
  close(pb)
  
  
  
  out <- list(A=listA,B=listB,C=listC,Sig= listSig)
  return(out)
}