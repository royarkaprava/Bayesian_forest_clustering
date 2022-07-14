n <- 200
p <- 2

y<- matrix(rnorm(n*2),n,p)
y[101:n,]<- y[101:n,]+3

D<- as.matrix(dist(y))

mu_sigma_tilde = apply(D, 1, function(x)min(x[x>0]))


sigma_tilde<- mu_sigma_tilde
sigma2<- outer(sigma_tilde, sigma_tilde, "*")


mu_r <- colMeans(y)
y_mu_r_2norm<- apply(t(t(y)-mu_r),1,function(x) sqrt(sum(x**2)))

gamma_r<-  max(y_mu_r_2norm)

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



updateLogS<- function(sigma_tilde){
  sigma2<- outer(sigma_tilde, sigma_tilde, "*")
  logS_data <- -D**2/2/sigma2 - log(sigma2 * 2*pi)/2
  logS<- matrix(0,n+1,n+1)
  logS[1:n,1:n]<- logS_data
  
  y_mu_r_2norm<- apply(t(t(y)-mu_r),1,function(x) sqrt(sum(x**2)))
  logS[1:n,n+1]<- log(lam) +  dt(y_mu_r_2norm /gamma_r, df = 1, log = T) -  log(gamma_r)
  logS[n+1,1:n]<- logS[1:n,n+1]

  diag(logS)<- -1E8
  
  logS
}


sample_sigma_tilde<- function(A_T){
  A_forest<- A_T[1:n,1:n] 
  lam_gig<- -rowSums(A_forest)/2+1
  psi_gig <- 2/mu_sigma_tilde
  
  AD2<- A_forest * D**2
  
  for(i in c(1:n)){
    chi_gig_i<- sum(AD2[i,]/sigma_tilde)
    sigma_tilde[i]<- rgig(1, lam_gig[i], chi_gig_i, psi_gig[i] )
  }
  sigma_tilde
}



n_iter = 1000

trace_K<- numeric()
trace_A_T<- list()

pb <- txtProgressBar(min = 0,      # Minimum value of the progress bar
                     max = n_iter, # Maximum value of the progress bar
                     style = 3,    # Progress bar style (also available style = 1 and style = 2)
                     width = 50,   # Progress bar width. Defaults to getOption("width")
                     char = "=")   # Character used to create the bar


for (step in c(1:n_iter)){
  A_T <-drawT(logS)
  sigma_tilde<- sample_sigma_tilde(A_T)
  sigma2<- outer(sigma_tilde, sigma_tilde, "*")
  logS<- updateLogS(sigma_tilde)

  setTxtProgressBar(pb, step)
  
  
  if (step> (n_iter/2)){
    trace_K<- c(trace_K, sum(A_T[n+1,]))
    trace_A_T[[step-(n_iter/2)]]<- A_T
  }
}


barplot(table(trace_K))

G<- graph_from_adjacency_matrix(A_T[1:n,1:n],mode = "undirected")
G$layout=y
plot.igraph(G,vertex.size=5)
