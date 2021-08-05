args = commandArgs(trailingOnly=TRUE)

N= as.numeric(args[1])
batchID= as.numeric(args[2])

print(N)

setwd("~/git/quasi-bernoulli-stick-breaking/")
# library(LaplacesDemon)


mu_truth = cbind(c(-4, 0, 5),c(1, 2, 3))

n_truth = round(N* c(0.3, 0.3, 0.4))

y = matrix(0, N,2)



for(i in 1:N){
  y[i,] = mu_truth[sum(i >cumsum(n_truth))+1,] + rnorm(2)
}


y=y[,1]

# plot(y[,1],y[,2])


source("qbsb_gibbs_univariate_gau.R")

res = runQBSBUvGaussian(y, steps= 2E4, burnins=1E4, K=20,p_b = 0.9)
# res = runQBSBUvGaussian(y, steps= 10, burnins=5, K=20,p_b = 0.9)

save(res, file = paste("sim_output/1D_Normal_Sim_N_",N,"_batch_",batchID,".Rda",sep=""))
