require(igraph)
require(GIGrvg)

forestClust<- function(y, lam = 0.5,alpha_gamma = 0.5,n_iter = 1000, burnin=500,fastTree= TRUE, updateRootLocation=FALSE,logTreePrior=0){
  n<- nrow(y)
  p<- ncol(y)
  
    
  # get pairwise distance
  D<- as.matrix(dist(y))
  
  mu_sigma_tilde = apply(D, 1, function(x)min(x[x>0]))/sqrt(p)
  
  sigma_tilde<- mu_sigma_tilde
  sigma2<- outer(sigma_tilde, sigma_tilde, "*")
  
  mu_r0 <- colMeans(y)
  mu_r = mu_r0
  
  y_mu_r_2norm<- apply(t(t(y)-mu_r),1,function(x) sqrt(sum(x**2)))
  
  gamma2_r0<-  mean(y_mu_r_2norm**2)/p
  gamma_r = sqrt(gamma2_r0)
  
  dmultinorm<- function(d,sigma2){
    -d**2/2/sigma2 - log(sigma2 * 2*pi)/2*p
  }
  
  dmulticacuhy<- function(d,sigma2){
    -log(1+d**2/sigma2)*(1+p)/2 + lgamma((1+p)/2) -  log(sigma2)/2*p  - (1+p)/2*log(pi)
  }
  
  
  sigma2<- outer(sigma_tilde, sigma_tilde, "*")
  
  logS_data <- dmultinorm(D,sigma2)
  
  logS<- matrix(0,n+1,n+1)
  logS[1:n,1:n]<- logS_data
  
  logS[1:n,n+1]<- log(lam) +  dmulticacuhy(y_mu_r_2norm, gamma_r**2)
  logS[n+1,1:n]<- logS[1:n,n+1]
  
  logS = logS+ logTreePrior
  diag(logS)<- -Inf
                       
    
  
  
  findMST<- function(logS){
    G<- graph_from_adjacency_matrix(-logS, mode = "undirected",weighted = TRUE,diag=FALSE)
    G_mst<- mst(graph = G)
    A_T<- as.matrix(get.adjacency(G_mst))
    A_T
  }
  
  A_T<- findMST(logS)
  
  
  rgumbel<- function(n){
    -log(-log(runif(n)))
  }
  
  
  gumbelMax<- function(logA){
    which.max( logA+ rgumbel(length(logA)))
  }
  
  drawT_approx<- function(logS){
    gumbelMat<- matrix(0,n+1,n+1)
    gumbelMat[lower.tri(gumbelMat,diag = FALSE)]<- rgumbel((n+1)/2*(n))
    gumbelMat<- gumbelMat+ t(gumbelMat)
    A_T<- findMST(logS+gumbelMat)
    A_T
  }
  
  
  drawT_exact<- function(logS){
    
    A_T<- matrix(0,n+1,n+1)
    
    InTree<- list()
    Next <- list()
    
    for (i in 1:(n+1)){
      InTree[[i]]<- FALSE
    }
    
    r<- n+1
    InTree[[r]]<- TRUE
    Next[[r]]<- n+1
    
    for (i in (n+1):1){
      u = i
      while(!InTree[[u]]){
        Next[[u]]<- gumbelMax(logS[u,])
        u <- Next[[u]]
      }
      u = i
      while(!InTree[[u]]){
        InTree[[u]]= TRUE
        u <- Next[[u]]
      }
    }
    
    for (u in 1:n){
      A_T[u, Next[[u]]]=1
      A_T[Next[[u]],u]=1
    } 
    A_T
  }
  
  drawT<- function(logS){
    if(fastTree){
      A_T = drawT_approx(logS)
    }else{
      A_T = drawT_exact(logS)
    }
    return(A_T)
  }
  
  
  updateLogS<- function(sigma_tilde,mu_r,gamma_r){
    sigma2<- outer(sigma_tilde, sigma_tilde, "*")
    logS_data <- dmultinorm(D,sigma2)
    logS<- matrix(0,n+1,n+1)
    logS[1:n,1:n]<- logS_data
    
    y_mu_r_2norm<- apply(t(t(y)-mu_r),1,function(x) sqrt(sum(x**2)))
    logS[1:n,n+1]<- log(lam) +  dmulticacuhy(y_mu_r_2norm, gamma_r**2)
    logS[n+1,1:n]<- logS[1:n,n+1]
    
    diag(logS)<- -Inf
    
    logS+logTreePrior
  }
  
  
  sample_sigma_tilde<- function(A_T){
    
    
    A_forest<- A_T[1:n,1:n] 
    lam_gig<- -rowSums(A_forest)/2*p + (alpha_gamma)
    psi_gig <- 2/(mu_sigma_tilde)
    
    AD2<- A_forest * D**2
    
    for(i in c(1:n)){
      chi_gig_i<- sum(AD2[i,]/sigma_tilde)
      sigma_tilde[i]<- rgig(1, lam_gig[i], chi_gig_i, psi_gig[i] )
    }
    sigma_tilde
  }
  
  
  sample_mu_r_gamma_r<- function(mu_r,gamma_r, A_T){
    y_mu_r_2norm <- apply(t(t(y)-mu_r),1,function(x) sqrt(sum(x**2)))
    
    w_cauchy <- 1/rgamma(n, 1/p + 1/2, rate = (y_mu_r_2norm**2)/2/gamma_r/gamma_r + 1/2)
    
    root_sel <- A_T[n+1,1:n]==1
    
    gamma2_r<- 1/rgamma(1, sum(root_sel)*p/2+ 2,  rate = sum( ((y_mu_r_2norm**2)/2/w_cauchy)[root_sel]) + gamma2_r0)
    gamma_r<- sqrt(gamma2_r)
    
                          
    if(updateRootLocation){
            var_mu_r<- 1/ sum(1/w_cauchy[root_sel]/gamma2_r + 1000)
        if(sum(root_sel)>1){
          mean_mu_r<-  var_mu_r* ( colSums(y[root_sel,]/w_cauchy[root_sel]/gamma2_r) + 1000*mu_r0)
        }else{
          mean_mu_r<-  var_mu_r* ( (y[root_sel,]/w_cauchy[root_sel]/gamma2_r) + 1000*mu_r0)
        }
        mu_r<-  rnorm(p, mean=mean_mu_r, sqrt(var_mu_r))
        }
    return(list(mu_r,gamma_r))
  }
  
  
  param_mu_gamma_r = c(mu_r,gamma_r)
  
  
  trace_K<- numeric()
  trace_A_T<- list()
  
  trace_sigma_tilde<- list()
  
  pb <- txtProgressBar(min = 0,      # Minimum value of the progress bar
                       max = n_iter, # Maximum value of the progress bar
                       style = 3,    # Progress bar style (also available style = 1 and style = 2)
                       width = 50,   # Progress bar width. Defaults to getOption("width")
                       char = "=")   # Character used to create the bar
  
  
  for (step in c(1:n_iter)){
    sigma_tilde<- sample_sigma_tilde(A_T)
    
    res<- sample_mu_r_gamma_r(mu_r,gamma_r,A_T)
    gamma_r<- res[[2]]
    mu_r<- res[[1]]
    
    logS<- updateLogS(sigma_tilde,mu_r,gamma_r)
    
    A_T <- drawT(logS)
    
    setTxtProgressBar(pb, step)
    
    
    if (step> burnin ){
      trace_K<- c(trace_K, sum(A_T[n+1,]))
      trace_A_T[[step-burnin]]<- A_T
      trace_sigma_tilde[[step-burnin]]<- sigma_tilde
    }
  }
  
  # barplot(table(trace_K))
  
  extractC <- function(A_T){
    G<- graph_from_adjacency_matrix(A_T[1:n,1:n],mode = "undirected")
    components(G)$membership
  }
  
  trace_C<- lapply(trace_A_T, extractC)
  
  
  C_mat<- matrix(0,n,n)
  
  for(C in trace_C){
    C_mat = C_mat+ 1*outer(C,C,'==')
  }
  
  
  A_T_mat<- matrix(0,n,n) 
  
  for(A_T in trace_A_T){
    A_T_mat = A_T_mat+ A_T[1:n,1:n]
  }
  A_T_mat = A_T_mat/length(trace_A_T)
  
  
  list(
    "K" = trace_K,
    "A_T" = trace_A_T,
    "sigma_tilde" = trace_sigma_tilde,
    "C" = trace_C,
    "C_mat" = C_mat,
    "logS" = logS
  )
}