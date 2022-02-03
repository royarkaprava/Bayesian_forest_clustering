require("igraph")

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
logden_r_elem <- function(d2, gamma=1){
  - log(1+d2/(gamma^2)) - log(pi) - log(gamma)
}


clusteringFP <- function(X, Total_itr = 10000, burn=5000, gamma =1, lambda=1){
  
  listA<- list()
  listB<- list()
  listC<- list()
  listSig<- list()
  

  
  p     <- ncol(X)
  n     <- nrow(X)
  
  # use mean
  mu <- apply(X, 2, mean)
  
  
  # Compute the MST among y's
  DisMat2 <- as.matrix(dist(X,"euclidean", diag = TRUE, upper = TRUE)^2)
  
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
  
    
  G<- graph_from_adjacency_matrix(A,"undirected")
   
  # initialize the clustering membership
  C <- rep(1,n+1) 
  C[1]<- 0
  
  # S
    
  dist2Xmu<-  (t(X)-mu)**2
  
  Sr <- colSums(logden_r_elem(dist2Xmu, gamma=1) )
  
  Sf <- logden_f(d2 = DisMat2, sigma2 = sig, p=p)
  S<- rbind(c(0,Sr), cbind(Sr,Sf))
  diag(S)<- 0
  
  ##############
  
  full_edge_idx = c(1:n)
  
#   B2Inv <- solve(t(B)%*%B)
  
  itr = 0
  pb <- txtProgressBar(min = itr, max = Total_itr, style = 3)
    

  edge_list = ends(G,E(G))

  
  for (itr in c(1:Total_itr)){
    
    # create a copy of current C
    
    cur_C <- C
    
    
    # update every edge
      
    
    for(k in c(1:n)){

        
        e = edge_list[k,]
        G = delete_edges(G, paste(e[1],"|",e[2]))

        Ga = subcomponent(G,v = e[1])
#         Gb = subcomponent(G,v = e[2])
        
        Gb = setdiff(c(1:(n+1)),Ga )

        if(1 %in% Gb){
            tmp <- Gb
            Gb  <- Ga
            Ga  <- tmp
        }

        
      # the number of clusters within Ga
      
      if(length(C[Ga])>1){
        Kstar = length(unique((C[Ga])[-1]))
      }else{
        Kstar = 1
      }
      
      SGaGb<- S[Ga,Gb]
      
      if(length(Ga)==1){
        SGaGb<- matrix(SGaGb,nrow=1,ncol=length(SGaGb))
      }
      if(length(Gb)==1){
        SGaGb<- matrix(SGaGb,ncol=1,nrow=length(SGaGb))
      }
      
      # adjust for the change in Poisson prior if attaching Gb to node 0 (idx 1 in R program here)
      
        
      SGaGb[Ga==1,] <- SGaGb[Ga==1,] + log(lambda) - log(Kstar+1) 
      
      ########
      rgumbel<- -log(-log(runif(length(SGaGb))))
      SGaGbGumbel<- SGaGb+ rgumbel
      
      GaIdx<- c(1:(n+1))[Ga]
      GbIdx<- c(1:(n+1))[Gb]
      
      GaIdx1<- rep(GaIdx,length(GbIdx))
      GbIdx1<- rep(GbIdx,each= length(GaIdx))
      
      max_idx = which.max(SGaGbGumbel)
      picked_i <- GaIdx1[max_idx]
      picked_j <- GbIdx1[max_idx]
      
      b_vec<- B[,k]
      B[,k]<- 0
      B[picked_i,k]<- 1
      B[picked_j,k]<- -1
        
      edge_list[k,1] = picked_i
      edge_list[k,2] = picked_j
        
      G = add_edges(G,edges =c(picked_i,picked_j) )
      
      if(picked_i!=e[1] | picked_j!=e[2]  ){ # if we change B[,k], update B2Inv and C
        
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
    
#     Sys.sleep(0.1)
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

