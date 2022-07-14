rm(list=ls())

n <- 200
p <- 2

y<- matrix(rnorm(n*2),n,p)
y[101:n,]<- y[101:n,]+5

D<- as.matrix(dist(y))

mu_sigma_tilde = apply(D, 1, function(x)min(x[x>0]))


sigma_tilde<- mu_sigma_tilde
sigma2<- outer(sigma_tilde, sigma_tilde, "*")


mu_r0 <- colMeans(y)

mu_r = mu_r0

y_mu_r_2norm<- apply(t(t(y)-mu_r),1,function(x) sqrt(sum(x**2)))

gamma2_r0<-  mean(y_mu_r_2norm**2)
gamma_r = sqrt(gamma2_r0)

lam = 0.5


sigma2<- outer(sigma_tilde, sigma_tilde, "*")
logS_data <- -D**2/2/sigma2 - log(sigma2 * 2*pi)/2
logS<- matrix(0,n+1,n+1)
logS[1:n,1:n]<- logS_data

logS[1:n,n+1]<- log(lam) +  dt(y_mu_r_2norm /gamma_r, df = 1, log = T) -  log(gamma_r)
logS[n+1,1:n]<- logS[1:n,n+1]

diag(logS)<- -1E8

rgumbel<- function(n){
    -log(-log(runif(n)))
}

gumbelMax<- function(logA){
  which.max( logA+ rgumbel(length(logA)))
}


drawT<- function(logS){
  
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


require(igraph)
require(GIGrvg)



updateLogS<- function(sigma_tilde,mu_r,gamma_r){
  sigma2<- outer(sigma_tilde, sigma_tilde, "*")
  logS_data <- -D**2/2/sigma2 - log(sigma2 * 2*pi)/2
  logS<- matrix(0,n+1,n+1)
  logS[1:n,1:n]<- logS_data
  
  y_mu_r_2norm<- apply(t(t(y)-mu_r),1,function(x) sqrt(sum(x**2)))
  logS[1:n,n+1]<- log(lam) +  dt(y_mu_r_2norm /gamma_r, df = 1, log = T) -  log(gamma_r)
  logS[n+1,1:n]<- logS[1:n,n+1]

  diag(logS)<- -Inf
  
  logS
}


sample_sigma_tilde<- function(A_T){
  
  alpha_gamma =1
  
  A_forest<- A_T[1:n,1:n] 
  lam_gig<- -rowSums(A_forest)/2+1 + (alpha_gamma - 1)
  psi_gig <- 2/(mu_sigma_tilde)
  
  AD2<- A_forest * D**2
  
  for(i in c(1:n)){
    chi_gig_i<- sum(AD2[i,]/sigma_tilde)
    sigma_tilde[i]<- rgig(1, lam_gig[i], chi_gig_i, psi_gig[i] )
  }
  sigma_tilde
}


logCauchyDensity<- function(param,A_T){
  
  mu_r = param[1:p]
  gamma_r = abs(param[p+1])
  
  y_mu_r_2norm <- apply(t(t(y)-mu_r),1,function(x) sqrt(sum(x**2)))
  y_mu_sel <- y_mu_r_2norm[(A_T[n+1,1:n]==1)]
  logDen<- sum(dt(y_mu_sel /gamma_r, df = 1, log = T)- log(gamma_r))
  logPrior<-    dt( sqrt(sum((mu_r- mu_r0)**2)),df = 1,log = T) +
    dgamma(1/gamma_r/gamma_r,shape = 2, rate = gamma2_r0,log = T)
  logDen + logPrior
}


randomWalkMH<- function(param_mu_gamma_r, A_T,eps=0.1){
  cur_loglik = logCauchyDensity(param_mu_gamma_r, A_T)
  cur_param_mu_gamma_r <-  param_mu_gamma_r
  
  param_mu_gamma_r<-   cur_param_mu_gamma_r + rnorm(p+1)*eps
  prop_loglik = logCauchyDensity(param_mu_gamma_r, A_T)
  
  if(log(runif(1))< (prop_loglik - cur_loglik)){
    accept=1
  }else{
    accept=0
    param_mu_gamma_r<- cur_param_mu_gamma_r
  }
  
  return(list(param_mu_gamma_r,accept))
}

param_mu_gamma_r = c(mu_r,gamma_r)

n_iter = 1000

trace_K<- numeric()
trace_A_T<- list()

pb <- txtProgressBar(min = 0,      # Minimum value of the progress bar
                     max = n_iter, # Maximum value of the progress bar
                     style = 3,    # Progress bar style (also available style = 1 and style = 2)
                     width = 50,   # Progress bar width. Defaults to getOption("width")
                     char = "=")   # Character used to create the bar


accept_mh_count = 0
eps_mh = 1

for (step in c(1:n_iter)){
  A_T <-drawT(logS)
  sigma_tilde<- sample_sigma_tilde(A_T)
  # sigma2<- outer(sigma_tilde, sigma_tilde, "*")
  
  
  #MH to update r params:
  res<- randomWalkMH(param_mu_gamma_r,A_T,eps=eps_mh)
  param_mu_gamma_r<- res[[1]]
  accept_mh<- res[[2]]
  mu_r<- param_mu_gamma_r[1:p]
  gamma_r<- abs(param_mu_gamma_r[p+1])
  # print(accept_mh)
  if (step< (n_iter/2)){
    
    #adaptively tune the random walk radius
  
    if(step%%20==0){
      
      eps_mh = eps_mh* exp( accept_mh_count/20 - 0.3)
      
      # print(accept_mh_count/20)
      
      accept_mh_count = 0
    }
    
    accept_mh_count = accept_mh_count+accept_mh
    
  }
  
  
  logS<- updateLogS(sigma_tilde,mu_r,gamma_r)

  setTxtProgressBar(pb, step)
  
  
  if (step> (n_iter/2)){
    trace_K<- c(trace_K, sum(A_T[n+1,]))
    trace_A_T[[step-(n_iter/2)]]<- A_T
  }
}


barplot(table(trace_K))

extractC <- function(A_T){
  G<- graph_from_adjacency_matrix(A_T[1:n,1:n],mode = "undirected")
  components(G)$membership
}

trace_C<- lapply(trace_A_T, extractC)


i= 200
G<- graph_from_adjacency_matrix(trace_A_T[[i]][1:n,1:n],mode = "undirected")
G$layout=y
V(G)$color=trace_C[[i]]
plot.igraph(G,vertex.size=5, vertex.label='')

ts.plot(trace_K)

C_mat<- matrix(0,n,n)

for(C in trace_C){
  C_mat = C_mat+ 1*outer(C,C,'==')
}

image(C_mat)

# plot( (C_mat/length(trace_C))[1,])


# G<- graph_from_adjacency_matrix(A_T[1:n,1:n],mode = "undirected")
# G$layout=y
# V(G)$color=trace_C[[i]]
# plot.igraph(G,vertex.size=5, vertex.label='')

