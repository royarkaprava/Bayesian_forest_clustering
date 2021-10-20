require('torch')
Rcpp::sourceCpp('forest.cpp')

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



# function to reduce label switch

relabelC2<- function(C1,C2){
  ham_dist = sum(C1!=C2)
  
  for( i in 1:length(C1)){
    c1 = C1[i]
    c2 = C2[i]
    if (c1!=c2){
      C2prime = C2
      C2prime[C2==c2] = c1
      C2prime[C2==c1] = c2
      
      prop_dist = sum(C1!=C2prime)
      if(prop_dist<ham_dist){
        C2= C2prime
        ham_dist = prop_dist
      }
      #         print(ham_dist)
    }
  }
  return(C2)
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
    

  #intialize the graph
    
  G <- new(Graph, n+1)
  zero_based_idx<- c(0:n)
    
  for(k in 1:ncol(B)){
    b_edge <- B[,k]
    idx0 <- zero_based_idx[b_edge== -1]
    idx1 <- zero_based_idx[b_edge== 1]
    G$connect(idx0,idx1)    
  }
  
  
  # initialize the clustering membership
  C <- rep(1,n+1) 
  C[1]<- 0
  
  # S
  
  Sr <- logden_r(dist2Xmu, gamma=1,p=p)
  
  Sf <- logden_f(d2 = DisMat2, sigma2 = sig, p=p)
  S<- rbind(c(0,Sr), cbind(Sr,Sf))
  diag(S)<- 0
  
  ##############
  
#   full_edge_idx = c(1:n)
    
  itr = 0
  pb <- txtProgressBar(min = itr, max = Total_itr, style = 3)
  
  for (itr in c(1:Total_itr)){
    
    # create a copy of current C
    
    cur_C <- C
    
    
    # update every edge
    
    for(k in c(1:n)){
#       sel <- full_edge_idx == k
#       sel_not <- full_edge_idx != k
        
      # cut the kth edge
      b_edge <- B[,k]
      idx0 <- zero_based_idx[b_edge== 1]
      idx1 <- zero_based_idx[b_edge== -1]
      G$cut(idx0,idx1)    


    GaIdx <- G$getConnectedNodes(idx1)+1 # add 1 to make the idx 1-based, as used in R
        
    print(length(GaIdx))
        
        
    Ga<- rep(FALSE, n+1)

    Ga[GaIdx]<- TRUE
      
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
        
      G$connect(picked_i-1, picked_j-1)    # connect, reduce idx by 1 to make it 0-based
        
      if(sum(abs(B[,k]- b_vec))>0){ # if we change B[,k], update C
        
        if(picked_i==1){ # add a cluster
          C[Gb]<- min(setdiff(c(1:(n)),C[Ga]))
        }else{
          C[Gb]<- C[picked_i]
        }
      }
    }
    
    # relabeling C to minimize the hamming distance with last C
    
        idxC = min(setdiff(c(1:(n)),C))
        while(idxC< max(C)){
            C[C==max(C)] = idxC
            idxC = min(setdiff(c(1:(n)),C))
        }
    #     C<- relabelC2(cur_C,C)
    
    
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
      listC[[itr-burn]] <- C[-1]
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


getCoAssignmentProb<- function(listC){
  
  coAssign= lapply(listC, function(x)outer(x,x,"==")*1)
  
  coAssignMat = matrix(0,n,n)
  counter = 0
  
  for (l in c(1: length(listC))){
    
    if(length(unique(listC[[l]]))>1)   { 
      coAssignMat = coAssignMat+ coAssign[[l]]
      counter = counter+1
    }
  }
  
  coAssignMat = coAssignMat/counter
  coAssignMat
}

getAssignProb<- function(P,K,learning_rate=0.1){
  
  n<- nrow(P)
  coAssign_tf <- torch_tensor(P)
  
  
  w0 <- torch_randn(c(n,K),requires_grad = TRUE)
  
  W<- torch_exp(w0 - torch_logsumexp(w0,dim = 2,keepdim = TRUE))
  
  
  compLoss<- function(W){
    W2 <- torch_matmul(W, torch_transpose(W,1,2))
    loss <- torch_sum((W2 - coAssign_tf)**2)
    loss
  }
  
  
  optimizer <- optim_adam(w0, lr = learning_rate)
  
  
  for (t in 1:1000) {
    
    ### -------- Forward pass -------- 
    
    ### -------- compute loss -------- 
    W<- torch_exp(w0 - torch_logsumexp(w0,dim = 2,keepdim = TRUE))
    loss <- compLoss(W)
    if (t %% 10 == 0)
      cat("Epoch: ", t, "   Loss: ", loss$item(), "\n")
    
    ### -------- Backpropagation -------- 
    
    # Still need to zero out the gradients before the backward pass, only this time,
    # on the optimizer object
    optimizer$zero_grad()
    
    # gradients are still computed on the loss tensor (no change here)
    loss$backward()
    
    ### -------- Update weights -------- 
    
    # use the optimizer to update model parameters
    optimizer$step()
  }
  
  as_array(W)
  
}
