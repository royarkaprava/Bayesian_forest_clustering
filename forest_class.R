
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
                       Z_d="integer"
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
                         lam_gig<- -rowSums(A_forest)/2*p + (alpha_gamma)
                         psi_gig <- 2/(mu_sigma_tilde)
                         
                         AD2<- A_forest * D**2
                         
                         for(i in c(1:n)){
                           chi_gig_i<- sum(AD2[i,]/sigma_tilde)
                           
                           
                           sigma_tilde[i]<<- rgig(1, lam_gig[i], chi_gig_i, psi_gig[i] )
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
                       
                       
                       
                       init = function(y_) {
                         y <<- y_
                         n <<- nrow(y)
                         p <<- ncol(y)
                         D <<- as.matrix(dist(y))
                         mu_sigma_tilde <<- apply(D, 1, function(x)min(x[x>0]))/sqrt(p)
                         
                         logTreePrior <<- matrix(0,n+1,n+1)
                         useFixedLogS <<- FALSE
                         fastTree <<- TRUE
                         updateRootLocation<<- FALSE
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
                         sample_sigma_tilde(A_T)
                         sample_mu_r_gamma_r(mu_r,gamma_r,A_T)
                         updateLogS(sigma_tilde,mu_r,gamma_r)
                         drawT(logS)
                       },
                       
                       randomInitZ=function(Z_d_=2){
                         Z_d<<- as.integer(Z_d_)
                         Z<<- matrix(rnorm((n+1)*Z_d),nrow=n+1)
                         Z[n+1,]<<-0
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
  G<- graph_from_adjacency_matrix(A_T[1:n,1:n],mode = "undirected")
  components(G)$membership
}

sampleZ<- function(forest, z, rho = 0.01){
  n<-forest$n

  root_z_var<- 1
  
  z_dim = ncol(z)
  
  z_var<- rep(rho,n+1)
  z_var[n+1]<- root_z_var
  
  z[n+1,]<- 0
  
  for(i in 1:n){
      par1<- 1/sum(forest$A_T[i,]/z_var +1)
      par2<- colSums(forest$A_T[i,]*z/z_var )
      z[i,]<- rnorm(z_dim,par1 * par2,sqrt(par1))
  }

  # Laplacian<- function(A){
  #   diag(rowSums(A)) - A
  # }
  # L <- Laplacian(forest$A_T)
  
  # L[n+1,1:n] = L[n+1,1:n]/root_z_var
  # L[,n+1] = L[,n+1]/root_z_var
  # 
  # vari<- solve(L+diag(1,n+1))
  # z<- t(chol(vari))%*%matrix(rnorm(n+1),n+1,1) + mu_mapto_z
  
  D2_z <- as.matrix(dist(z,upper = TRUE, diag = TRUE))**2
  
  z_logdensity<- matrix(0,n+1,n+1)
  
  z_logdensity[1:n,1:n] <-  - D2_z[1:n,1:n]/2/rho - z_dim* log(2*pi* rho)/2
  z_logdensity[n+1,1:n] <-  - D2_z[n+1,1:n]/2/root_z_var - z_dim* log(2*pi* root_z_var)/2
  z_logdensity[1:n,n+1]<- z_logdensity[n+1,1:n]
  diag(z_logdensity)<- 0
  list(z,z_logdensity)
}