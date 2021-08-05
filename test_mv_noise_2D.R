args = commandArgs(trailingOnly = TRUE)

N = as.numeric(args[1])
batchID = as.numeric(args[2])

print(N)

setwd("C://Users//Cheng/OneDrive - University of Florida//Research//Quasi-Bernoulli_Stick-breaking//QBSB-R//simulations")
# library(LaplacesDemon)


mu_truth = cbind(c(-4, 0, 5), c(1, 2, 3))

n_truth = round(N * 0.7 * c(0.3, 0.3, 0.4))
y = matrix(0, N, 2)

for(i in 1:sum(n_truth)){
    y[i, ] = mu_truth[sum(i > cumsum(n_truth)) + 1, ] + rnorm(2)
}
for(i in (sum(n_truth)+1):N) {
    y[i, ] = c(runif(1, -6, 7), runif(1, -1, 5))
}


source("qbsb_gibbs_multivariate_gau.R")

res = runQBSBMvGaussian(y, steps = 2E4, burnins = 1E4, K = 20, p_b = 0.9, eps = 1 / nrow(y) ^ 1.01)
save(y, res, file = paste("sim_output/2D_Noise_Sim_N_", N, "_batch_", batchID, ".Rda", sep = ""))


res = runQBSBMvGaussian(y, steps = 2E4, burnins = 1E4, K = 20, p_b = 1)
save(y, res, file = paste("sim_output/2D_Noise_DP_Sim_N_", N, "_batch_", batchID, ".Rda", sep = ""))


res = runQBSBMvGaussian(y, steps = 2E4, burnins = 1E4, K = 20, p_b = 1, alpha_beta = -0.18, d_beta = 0.35)
save(y, res, file = paste("sim_output/2D_Noise_PY_Sim_N_", N, "_batch_", batchID, ".Rda", sep = ""))