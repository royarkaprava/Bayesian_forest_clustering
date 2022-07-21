rm(list=ls())

n <- 200
p <- 5

y<- matrix(rt(n*p,df = 10),n,p)

y[101:n,]<- y[101:n,]+2


require(igraph)
require(GIGrvg)

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



lam = 0.5

sigma2<- outer(sigma_tilde, sigma_tilde, "*")

logS_data <- dmultinorm(D,sigma2)

logS<- matrix(0,n+1,n+1)
logS[1:n,1:n]<- logS_data

logS[1:n,n+1]<- log(lam) +  dmulticacuhy(y_mu_r_2norm, gamma_r**2)
logS[n+1,1:n]<- logS[1:n,n+1]

diag(logS)<- -1E8


G<- graph_from_adjacency_matrix(-logS, mode = "undirected",weighted = TRUE,diag=FALSE)
G_mst<- mst(graph = G)
A_T<- as.matrix(get.adjacency(G_mst))

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



updateLogS<- function(sigma_tilde,mu_r,gamma_r){
  sigma2<- outer(sigma_tilde, sigma_tilde, "*")
  logS_data <- dmultinorm(D,sigma2)
  logS<- matrix(0,n+1,n+1)
  logS[1:n,1:n]<- logS_data
  
  y_mu_r_2norm<- apply(t(t(y)-mu_r),1,function(x) sqrt(sum(x**2)))
  logS[1:n,n+1]<- log(lam) +  dmulticacuhy(y_mu_r_2norm, gamma_r**2)
  logS[n+1,1:n]<- logS[1:n,n+1]

  diag(logS)<- -Inf
  
  logS
}


sample_sigma_tilde<- function(A_T){
  
  alpha_gamma = 0.5
  
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
  
  var_mu_r<- 1/ sum(1/w_cauchy[root_sel]/gamma2_r + 1)
  
  if(sum(root_sel)>1){
  mean_mu_r<-  var_mu_r* ( colSums(y[root_sel,]/w_cauchy[root_sel]/gamma2_r) + mu_r0)
  }else{
    mean_mu_r<-  var_mu_r* ( (y[root_sel,]/w_cauchy[root_sel]/gamma2_r) + mu_r0)
  }
  mu_r<-  rnorm(p, mean=mean_mu_r, sqrt(var_mu_r))
  return(list(mu_r,gamma_r))
}


param_mu_gamma_r = c(mu_r,gamma_r)

n_iter = 1000

trace_K<- numeric()
trace_A_T<- list()

trace_sigma_tilde<- list()

pb <- txtProgressBar(min = 0,      # Minimum value of the progress bar
                     max = n_iter, # Maximum value of the progress bar
                     style = 3,    # Progress bar style (also available style = 1 and style = 2)
                     width = 50,   # Progress bar width. Defaults to getOption("width")
                     char = "=")   # Character used to create the bar


accept_mh_count = 0
eps_mh = 1

for (step in c(1:n_iter)){
  sigma_tilde<- sample_sigma_tilde(A_T)
  # sigma2<- outer(sigma_tilde, sigma_tilde, "*")
  
  res<- sample_mu_r_gamma_r(mu_r,gamma_r,A_T)
  mu_r<- res[[1]]
  gamma_r<- res[[2]]
  
  logS<- updateLogS(sigma_tilde,mu_r,gamma_r)

  A_T <- drawT(logS)
  
  
  setTxtProgressBar(pb, step)
  
  
  if (step> (n_iter/2)){
    trace_K<- c(trace_K, sum(A_T[n+1,]))
    trace_A_T[[step-(n_iter/2)]]<- A_T
    trace_sigma_tilde[[step-(n_iter/2)]]<- sigma_tilde
  }
}

barplot(table(trace_K))

extractC <- function(A_T){
  G<- graph_from_adjacency_matrix(A_T[1:n,1:n],mode = "undirected")
  components(G)$membership
}

trace_C<- lapply(trace_A_T, extractC)



ts.plot(trace_K)

C_mat<- matrix(0,n,n)

for(C in trace_C){
    C_mat = C_mat+ 1*outer(C,C,'==')
}

image(C_mat)


A_T_mat<- matrix(0,n,n) 

for(A_T in trace_A_T){
  A_T_mat = A_T_mat+ A_T[1:n,1:n]
}
A_T_mat = A_T_mat/length(trace_A_T)

image(A_T_mat)

require("kernlab")

sc_fit<- specc(as.kernelMatrix(C_mat),centers=2)

G_D<-graph_from_adjacency_matrix(A_T_mat, mode = "undirected",weighted = TRUE,diag=FALSE)
G_D$layout=y[,1:2]
V(G_D)$color = as.numeric(sc_fit)

G_D<- delete.edges(G_D, E(G_D)[E(G_D)$weight< quantile(E(G_D)$weight,0.8)])

plot.igraph(G_D,vertex.size=5, vertex.label='')
