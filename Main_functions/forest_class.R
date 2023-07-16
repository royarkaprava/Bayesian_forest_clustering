require("igraph")
require("GIGrvg")
require("gtools")
require("RcppML")

Forest<- setRefClass("Forest",  
                     fields = list(
                       A_T = "matrix", 
                       y="matrix",
                       n="integer",
                       p="integer",
                       D="matrix",
                       mu_sigma_tilde="numeric",
                       logS="matrix",
                       fastTree="logical",
                       updateRootLocation="logical",
                       sigma_tilde = "numeric",
                       gamma_r = "numeric",
                       mu_r = "numeric",
                       alpha_gamma = "numeric",
                       gamma2_r0 = "numeric",
                       lam = "numeric",
                       logTreePrior = "matrix",
                       useFixedLogS = "logical",
                       Z="matrix",
                       Z_d="integer",
                       trace_A_T = "list",
                       trace_C = "list",
                       trace_sigma_tilde = "list",
                       trace_gamma_r = "list",
                       hierachical_prior="logical",
                       eta_sigma = "numeric"
                     ),
                     methods = list(
                       
                       dmultinorm= function(d,sigma2){
                         -d**2/2/sigma2 - log(sigma2 * 2*pi)/2*p
                       },
                       
                       
                       dmulticacuhy=function(d,sigma2){
                         -log(1+d**2/sigma2)*(1+p)/2 + lgamma((1+p)/2) -  log(sigma2)/2*p  - (1+p)/2*log(pi)
                       },
                       
                       
                       findMST= function(logS){
                         G<- graph_from_adjacency_matrix(-logS, mode = "undirected",weighted = TRUE,diag=FALSE)
                         G_mst<- mst(graph = G)
                         A_T_<- as.matrix(get.adjacency(G_mst))
                         A_T_
                       },
                       
                       
                       drawT_approx= function(logS){
                         gumbelMat<- matrix(0,n+1,n+1)
                         gumbelMat[lower.tri(gumbelMat,diag = FALSE)]<- rgumbel((n+1)/2*(n))
                         gumbelMat<- gumbelMat+ t(gumbelMat)
                         A_T_ <- findMST(logS+gumbelMat)
                         A_T_
                       },
                       
                       
                       drawT_exact = function(logS){
                         
                         A_T_<- matrix(0,n+1,n+1)
                         
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
                           A_T_[u, Next[[u]]]=1
                           A_T_[Next[[u]],u]=1
                         } 
                         A_T_
                       },
                       
                       drawT= function(logS){
                         if(fastTree){
                           A_T_ = drawT_approx(logS)
                         }else{
                           A_T_ = drawT_exact(logS)
                         }
                         A_T<<- A_T_
                       },
                       
                       
                       updateLogS = function(sigma_tilde,mu_r,gamma_r){
                         sigma2<- outer(sigma_tilde, sigma_tilde, "*")
                         logS_data <- dmultinorm(D,sigma2)
                         logS_<- matrix(0,n+1,n+1)
                         logS_[1:n,1:n]<- logS_data
                         
                         y_mu_r_2norm<- apply(t(t(y)-mu_r),1,function(x) sqrt(sum(x**2)))
                         logS_[1:n,n+1]<- log(lam) +  dmulticacuhy(y_mu_r_2norm, gamma_r**2)
                         logS_[n+1,1:n]<- logS_[1:n,n+1]
                         
                         diag(logS_)<- -Inf
                         
                         logS<<- logS_+logTreePrior
                       },
                       
                       
                       sample_sigma_tilde = function(A_T_){
                         
                         
                         A_forest<- A_T_[1:n,1:n] 
                         
                         
                         if(hierachical_prior){
                           
                           a_sigma = 100
                           b_sigma = 10
                           xi_sigma = 1
                           
                           beta <- rgamma(1,  n* b_sigma + 1, rate = sum(1/sigma_tilde) + 1/eta_sigma)
                           eta_sigma<<- 1/rgamma(1, a_sigma+ 1,rate = beta+ xi_sigma)
                           
                           lam_gig<- rowSums(A_forest)/2*p + b_sigma
                           
                           AD2<- A_forest * D**2
                           
                           for(i in c(1:n)){
                             chi_gig_i<- sum(AD2[i,]/sigma_tilde)/2 + beta
                             sigma_tilde[i]<<- 1/ rgamma(1,  lam_gig[i],  rate = chi_gig_i)
                           }
                           
                         }
                         else{
                           
                           lam_gig<- -rowSums(A_forest)/2*p + (alpha_gamma)
                           psi_gig <- 2/(mu_sigma_tilde)
                           
                           AD2<- A_forest * D**2
                           
                           for(i in c(1:n)){
                             chi_gig_i<- sum(AD2[i,]/sigma_tilde)
                             
                             
                             sigma_tilde[i]<<- rgig(1, lam_gig[i], chi_gig_i, psi_gig[i] )
                           }
                           
                         }
                       },
                       
                       
                       sample_mu_r_gamma_r = function(mu_r,gamma_r, A_T_){
                         y_mu_r_2norm <- apply(t(t(y)-mu_r),1,function(x) sqrt(sum(x**2)))
                         
                         w_cauchy <- 1/rgamma(n, 1/p + 1/2, rate = (y_mu_r_2norm**2)/2/gamma_r/gamma_r + 1/2)
                         
                         root_sel <- A_T_[n+1,1:n]==1
                         
                         gamma2_r<- 1/rgamma(1, sum(root_sel)*p/2+ 2,  rate = sum( ((y_mu_r_2norm**2)/2/w_cauchy)[root_sel]) + gamma2_r0)
                         gamma_r<<- sqrt(gamma2_r)
                         
                         
                         if(updateRootLocation){
                           var_mu_r<- 1/ sum(1/w_cauchy[root_sel]/gamma2_r + 1000)
                           if(sum(root_sel)>1){
                             mean_mu_r<-  var_mu_r* ( colSums(y[root_sel,]/w_cauchy[root_sel]/gamma2_r) + 1000*mu_r0)
                           }else{
                             mean_mu_r<-  var_mu_r* ( (y[root_sel,]/w_cauchy[root_sel]/gamma2_r) + 1000*mu_r0)
                           }
                           mu_r<<-  rnorm(p, mean=mean_mu_r, sqrt(var_mu_r))
                         }
                       },
                       
                       
                       
                       init = function(y_,use_hierachical_prior=FALSE) {
                         y <<- y_
                         n <<- nrow(y)
                         p <<- ncol(y)
                         D <<- as.matrix(dist(y))
                         mu_sigma_tilde <<- apply(D, 1, function(x)min(x[x>0]))/sqrt(p)
                         eta_sigma<<- min(D[D>0])/sqrt(p)
                         
                         
                         logTreePrior <<- matrix(0,n+1,n+1)
                         useFixedLogS <<- FALSE
                         fastTree <<- TRUE
                         updateRootLocation<<- FALSE
                         hierachical_prior<<- use_hierachical_prior
                         
                         
                         alpha_gamma <<- 0.5
                         lam<<- 0.5
                         
                         sigma_tilde<<- mu_sigma_tilde
                         sigma2<- outer(sigma_tilde, sigma_tilde, "*")
                         
                         mu_r0 <- colMeans(y)
                         mu_r <<- mu_r0
                         
                         
                         y_mu_r_2norm<- apply(t(t(y)-mu_r),1,function(x) sqrt(sum(x**2)))
                         
                         gamma2_r0<<-   mean(y_mu_r_2norm**2)/p
                         gamma_r <<- sqrt(gamma2_r0)
                         
                         sigma2<- outer(sigma_tilde, sigma_tilde, "*")
                         
                         logS_data <- dmultinorm(D,sigma2)
                         
                         logS_<- matrix(0,n+1,n+1)
                         logS_[1:n,1:n]<- logS_data
                         
                         logS_[1:n,n+1]<- log(lam) +  dmulticacuhy(y_mu_r_2norm, gamma_r**2)
                         logS_[n+1,1:n]<- logS_[1:n,n+1]
                         
                         logS_ = logS_+ logTreePrior
                         diag(logS_)<-  -Inf
                         
                         logS<<- logS_
                         
                         if(useFixedLogS){
                           logS<<- fixedLogS
                         }
                         
                         A_T <<- drawT(logS)
                         
                       },
                       
                       MCMC_oneScan= function(){
                         drawT(logS)
                         sample_sigma_tilde(A_T)
                         sample_mu_r_gamma_r(mu_r,gamma_r,A_T)
                         updateLogS(sigma_tilde,mu_r,gamma_r)
                       },
                       
                       MCMC_run_single_graph = function(n_iter= 1000, burnin=500){
                         
                         pb <- txtProgressBar(min = 0,      # Minimum value of the progress bar
                                              max = n_iter, # Maximum value of the progress bar
                                              style = 3,    # Progress bar style (also available style = 1 and style = 2)
                                              width = 50,   # Progress bar width. Defaults to getOption("width")
                                              char = "=")   # Character used to create the bar
                         
                         
                         
                         for(step in 1:n_iter){
                           
                           MCMC_oneScan()
                           if(step>burnin){
                             step1 = step - burnin
                             trace_A_T[[step1]]<<- A_T
                             trace_C[[step1]]<<- extractC(A_T[1:n,1:n])
                             trace_sigma_tilde[[step1]]<<- sigma_tilde
                             trace_gamma_r [[step1]]<<- gamma_r                              
                           }
                           
                           setTxtProgressBar(pb, step)
                           
                           
                         }
                         
                       },
                       
                       randomInitZ=function(Z_d_=2){
                         Z_d<<- as.integer(Z_d_)
                         Z<<- matrix(rnorm((n+1)*Z_d),nrow=n+1)
                         Z[n+1,]<<-0
                       },
                       
                       
                       useZforLogPrior = function(sigma2_z, rho){
                         z<- Z
                         D2_z <- as.matrix(dist(z,upper = TRUE, diag = TRUE))**2
                         
                         z_logdensity<- matrix(0,n+1,n+1)
                         
                         z_logdensity[1:n,1:n] <-  - D2_z[1:n,1:n]/2/rho - Z_d* log(2*pi* rho)/2
                         # z_logdensity[n+1,1:n] <-  - D2_z[n+1,1:n]/2/rho - Z_d* log(2*pi* rho)/2
                         # z_logdensity[1:n,n+1] <-  z_logdensity[n+1,1:n]
                         diag(z_logdensity)<- 0
                         
                         logTreePrior<<- z_logdensity
                       }
                       
                       
                       
                     )
)



rgumbel = function(n){
  -log(-log(runif(n)))
}

gumbelMax = function(logA){
  which.max( logA+ rgumbel(length(logA)))
}

extractC <- function(A_T){
  G<- graph_from_adjacency_matrix(A_T,mode = "undirected")
  components(G)$membership
}


# update eta_star

sample_eta_star<- function(){
  for (l in 1: Z_K_tilde){
    
    Z_sel <- lapply(c(1:S), function(s)forest_objs[[s]]$Z[Z_k[,s]==l,])
    
    Z_sel<- do.call("rbind",Z_sel)
    
    par2 = 1/( nrow(Z_sel)/sigma2_z + 1/10)
    par1 = par2 * (colSums(Z_sel)/sigma2_z)
    eta_star[l,]<<- rnorm(Z_d)* sqrt(par2) + par1
  }}
# update Z

sampleZ<- function(sigma2_z,rho){
  
  
  for(s in 1:S){
    
    a =   forest_objs[[s]]$A_T
    
    a[n+1,] = 0
    a[,n+1] = a[n+1,]
    
    L<- diag(rowSums(a))- a
    L1n<- L[1:n,1:n]
    
    prec<- L1n/rho + diag(1/sigma2_z,n)
    
    chol_prec<- t(chol(prec))
    
    z0<- chol_prec%*% matrix(rnorm(n*Z_d),n) + eta_star[Z_k[,s],]/sigma2_z
    
    forest_objs[[s]]$Z[1:n,]<<- solve(prec,z0)
  }
  
}


# update the assignment Z_k

sampleZ_k<- function(){
  
  for(s in 1:S){
    for(i in 1:n){
      logChoice = -colSums( (t(eta_star) - forest_objs[[s]]$Z[i,])**2)/sigma2_z/2 + log(v[i,])
      Z_k[i,s]<<- gumbelMax(logChoice)
    }
  }
}

# update v:

sampleV<- function(alpha=1){
  count_table = matrix(0, n, Z_K_tilde)
  
  for(l in 1:Z_K_tilde){
    count_table[,l] = rowSums(Z_k==l)
  }
  
  v<<- matrix(rdirichlet(1, c(count_table)+ alpha/Z_K_tilde),nrow = n)
}




# other function


require("clue")

matchAtoB<- function(clusteringA, clusteringB){
  # labels from cluster A will be matched on the labels from cluster B
  minWeightBipartiteMatching <- function(clusteringA, clusteringB) {
    require(clue)
    idsA <- unique(clusteringA)  # distinct cluster ids in a
    idsB <- unique(clusteringB)  # distinct cluster ids in b
    nA <- length(clusteringA)  # number of instances in a
    nB <- length(clusteringB)  # number of instances in b
    if (length(idsA) != length(idsB) || nA != nB) {
      stop("number of cluster or number of instances do not match")
    }
    nC <- length(idsA)
    tupel <- c(1:nA)
    # computing the distance matrix
    assignmentMatrix <- matrix(rep(-1, nC * nC), nrow = nC)
    for (i in 1:nC) {
      tupelClusterI <- tupel[clusteringA == i]
      solRowI <- sapply(1:nC, function(i, clusterIDsB, tupelA_I) {
        nA_I <- length(tupelA_I)  # number of elements in cluster I
        tupelB_I <- tupel[clusterIDsB == i]
        nB_I <- length(tupelB_I)
        nTupelIntersect <- length(intersect(tupelA_I, tupelB_I))
        return((nA_I - nTupelIntersect) + (nB_I - nTupelIntersect))
      }, clusteringB, tupelClusterI)
      assignmentMatrix[i, ] <- solRowI
    }
    # optimization
    result <- solve_LSAP(assignmentMatrix, maximum = FALSE)
    attr(result, "assignmentMatrix") <- assignmentMatrix
    return(result)
  }
  
  matching<- minWeightBipartiteMatching(
    clusteringA,clusteringB
  )
  
  newClusteringA<- clusteringA
  
  tmp <- lapply(1:length(matching), function(i) {
    newClusteringA[which(clusteringA== i)] <<- matching[i]
  })
  
  newClusteringA
}

clusteringAccu<- function(clustering,true_membership){
  n<- length(true_membership)
  sum(matchAtoB(clustering,true_membership)==true_membership)/n
}

getCoAssignmentMat<- function(C_list){
  
  len<- length(C_list)
  n<- length(C_list[[1]])
  mat<- matrix(0,n,n)
  
  for(i in 1:len){
    mat<- mat + outer(C_list[[i]],C_list[[i]],"==")
  }
  mat/len
}


getPointEstC<- function(C_mat, K){
  nmf_fit<-nmf(C_mat, k=K, verbose = F)
  as.numeric(kmeans(nmf_fit$w[,1:K],K)$cluster)
  # as.numeric(specc(as.kernelMatrix( fc_fit$C_mat),centers=K))
}



forestClust<- function(y, lam = 0.5,alpha_gamma = 0.5,n_iter = 1000, burnin=500,fastTree= TRUE, updateRootLocation=FALSE,
                       logTreePrior= 0,
                       useFixedLogS= FALSE,hierachical_prior= FALSE
){
  
  forest<- Forest$new()
  
  forest$init(y)
  forest$useFixedLogS <- useFixedLogS
  forest$fastTree <- fastTree
  forest$updateRootLocation<- updateRootLocation
  forest$hierachical_prior<- hierachical_prior
  forest$lam <- lam
  forest$alpha_gamma<- alpha_gamma
  
  if(length(logTreePrior)==1){
    # do nothing
  }else{
    forest$logTreePrior<- logTreePrior
  }
  
  forest$MCMC_run_single_graph(n_iter,burnin)
  
  
  
  trace_K<- lapply( forest$trace_A_T, function(x){sum(x[forest$n+1,])})
  
  
  res<- list(
    "K" = trace_K,
    "A_T" = forest$trace_A_T,
    "sigma_tilde" = forest$trace_sigma_tilde,
    "C" = forest$trace_C,
    "gamma_r" =  forest$trace_gamma_r
  )
  
  return(res)
}