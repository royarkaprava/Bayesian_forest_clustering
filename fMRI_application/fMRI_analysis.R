 require("igraph")
require("GIGrvg")
require("gtools")

setwd("//file.phhp.ufl.edu/data/home/ark007/Documents/GitHub/Bayesian_forest_clustering")

load("fMRI_application/time_series_fmri_std.rda")

S <- length(timeseries_fmri_std)
S



load("fMRI_application/fmri_labels.RDa")



labels0<- do.call("c",label_list)



ord_idx<- order(labels0)



labels <- labels0[ord_idx]







source("Main_functions\\forest_class.R")




#create a list of forest objects

forest_objs = list()

y_mat<- numeric()

for(s in 1:S){
  forest_objs[[s]]<- Forest$new()
  forest_objs[[s]]$init(timeseries_fmri_std[[ord_idx[s]]][,1:120])
  y_mat<- rbind(y_mat, (timeseries_fmri_std[[ord_idx[s]]][,1:120]))
  
  forest_objs[[s]]$randomInitZ(2)
}

n<- forest_objs[[1]]$n
Z_d<- forest_objs[[1]]$Z_d
Z_K_tilde<- 3^Z_d

eta_star<- matrix(c(1:3),ncol=1)
for(j in (1:(Z_d-1))){
  eta_star<- rbind(cbind(eta_star,0),cbind(eta_star,1),cbind(eta_star,2))
}


Z_k<- matrix(kmeans(y_mat,centers= Z_K_tilde)$cluster,n,S,byrow = FALSE)
# Z_k<- matrix( as.numeric(sample(c(1:Z_K_tilde),n*S,replace=TRUE)), n)

v<- matrix(1,n,Z_K_tilde)/Z_K_tilde


for(s in 1:S){
  forest_objs[[s]]$Z = eta_star[Z_k[,s],] + rnorm(n*Z_d)*0.1
  forest_objs[[s]]$Z<- rbind(forest_objs[[s]]$Z,0)
}

trace_Z_k<- list()

Cmat_list<- list()
for(s in 1:S){
  Cmat_list[[s]]<- 0
}

n_iter= 100


pb <- txtProgressBar(min = 0,      # Minimum value of the progress bar
                     max = n_iter, # Maximum value of the progress bar
                     style = 3,    # Progress bar style (also available style = 1 and style = 2)
                     width = 50,   # Progress bar width. Defaults to getOption("width")
                     char = "=")   # Character used to create the bar

# rho1= 0.5 #0.5 **2
# rho2= 0.001 #0.1**2

sigma2_z = 0.1**2
rho= 1E-3

for(step in 1:n_iter){
  Z_k_copy<- Z_k
  
  for(s in 1:S){
    forest_objs[[s]]$MCMC_oneScan()
  }
  
  sampleZ(sigma2_z,rho)
  sampleZ_k()
  sampleV()
  
  # sample_eta_star()
  
  for(s in 1:S){
    forest_objs[[s]]$useZforLogPrior(sigma2_z,rho)
  }
  
  
  # print( sum(sum(abs(Z_k_copy-Z_k)>0)))
  
  setTxtProgressBar(pb, step)
  
  trace_Z_k[[step]]<- Z_k
  
  for(s in 1:S){
    mem<- extractC(forest_objs[[s]]$A_T)
    Cmat_list[[s]] = Cmat_list[[s]]+ outer(mem,mem,"==")*1
  }
  
}