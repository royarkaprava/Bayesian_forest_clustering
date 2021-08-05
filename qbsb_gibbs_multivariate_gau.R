require("LaplacesDemon")

runQBSBMvGaussian <- function(
    y, steps = 2E4, burnins = 1E4, K = 20,
    p_b = 0.9, d_beta = 0, alpha_beta = 1, eps = 1e-3,
    mu_prior_mean = colMeans(y), mu_prior_SigmasInv = solve(cov(y)),
    SigmasInv_prior_nu = ncol(y), SigmasInv_prior_VInv = SigmasInv_prior_nu * cov(y)
    ) {

    stopifnot(eps > 0 & (alpha_beta + d_beta) > 0 & d_beta < 1 & d_beta >= 0 & p_b >= 0 & p_b <= 1)
    stopifnot(p_b == 1 | d_beta == 0)

    N = nrow(y)
    p = ncol(y)
    rtbeta <- function(n, alpha, beta, a = 0, b = 1) {
        stopifnot(n > 0 & all(beta > 0) & all(alpha > 0))
        x <- runif(n)
        Fa <- pbeta(a, alpha, beta)
        Fb <- pbeta(b, alpha, beta)
        y <- (1 - x) * Fa + x * Fb
        y[y < 1E-16] = 1E-16
        return(qbeta(y, alpha, beta))
    }
    
    updateC <- function(mu, SigmasInv, w) {
        
        loglik_normal = matrix(0, N, K)
        
        for(k in c(1:K)) {
            diff = y - rep(mu[k, ], each = N)
            logdet1 = logdet(SigmasInv[, , k])
            loglik_normal[, k] = -diag(diff %*% SigmasInv[, , k] %*% t(diff)) / 2 + logdet1 / 2
        }
        
        loglik = t(t(loglik_normal) + log(c(w))) # R is adding by columns, so I used tranpose's twice 
        gumbel = -log(-log(runif(N * K, min = 0, max = 1)))
        C <- t(apply(loglik + gumbel, 1, function(x) { x == max(x) }))
        return(C)
    }
    
    
    updateMu <- function(C, SigmasInv, mu_prior_mean = rep(0, p), mu_prior_SigmasInv = diag(p)) {
        
        ySum = t(C) %*% y # for each cluster
        n_C = colSums(C)
        
        for(k in c(1:K)) {
            m_part2 = ySum[k, ] %*% SigmasInv[, , k] + mu_prior_mean %*% mu_prior_SigmasInv
            
            m_var = solve(SigmasInv[, , k] * n_C[k] + mu_prior_SigmasInv)
            m_mean = m_var %*% t(m_part2)
            
            mu[k, ] = t(chol(m_var)) %*% rnorm(p) + m_mean
        }
        return(mu)
    }
    
    
    updateSigmasInv <- function(C, mu, SigmasInv_prior_nu = p, SigmasInv_prior_VInv = diag(p)) {
        
        C_label = colSums((t(C) * c(1:K)))
        
        diff = (y - mu[C_label, ])
        
        for(k in c(1:K)) {
            pick = (C_label == k)
            
            if (sum(pick) > 1) {
                par2 = solve(t(diff[pick, ]) %*% (diff[pick, ]) + SigmasInv_prior_VInv)
            }
            
            if (sum(pick) == 1) {
                vec = matrix(diff[pick, ], 1, p)
                par2 = solve(t(vec) %*% (vec)  + SigmasInv_prior_VInv)
            }
            if (sum(pick) == 0) {
                par2 =  solve(SigmasInv_prior_VInv)
            }
            
            par2 = (par2 + t(par2)) / 2
            par1 = sum(C_label == k) + SigmasInv_prior_nu           
            SigmasInv[, , k] = rwishart(par1, par2)
        }        
        return(SigmasInv)
    }
    
    updateBeta <- function(C, b, alpha_beta = 1, eps = 1E-3, d_beta = 0) {

        n_C <- colSums(C)
        par1_beta = (N - cumsum(n_C)) + alpha_beta + d_beta * c(1:K)
        par2_beta = n_C + 1 - d_beta
        
        beta_1 = rbeta(K, par1_beta, par2_beta)
        if (p_b == 1) {
            beta = beta_1
            if (any(is.na(beta))) {
                # print("beta err")
                beta = updateBeta(C, b, alpha_beta, eps, d_beta)
            }
            return(beta)
        }
        beta_2 = rtbeta(K, par1_beta, par2_beta, 0, eps) / eps #truncated Beta distribution
        
        beta = beta_1 * (b == 1) + beta_2 * (b != 1)       
        # beta[K]=1E-8
        if (any(is.na(beta))) {
            # print("beta err")
            beta = updateBeta(C, b, alpha_beta, eps, d_beta)
        }
        return(beta)
    }
    
    
    updateB <- function(C, eps = 1E-3, p_b = 0.9, alpha_beta = 1) {
 
        n_C <- colSums(C)
        m_C = N - cumsum(n_C)
        
        choice1 = log(p_b)
        choice2 = log(1 - p_b) + pbeta(eps, alpha_beta + m_C, n_C + 1, log.p = T) - alpha_beta * log(eps)
        
        gumbel = -log(-log(runif(K * 2, min = 0, max = 1)))
        
        prob_choice = cbind(choice1, choice2) + gumbel
        b <- colSums((apply(prob_choice, 1, function(x) x == max(x))) * c(1, eps))
        return(b)
    }
    
    updateW <- function(b, beta) {
        v = 1 - b * beta
        w = v * (cumprod(c(1, 1 - v))[1:K])
        return(w)
    }
    
    
    partition_prob <- function(n_C, eps = 1E-3, p_b = 0.9, alpha_beta = 1, d_beta = 0) {
        logsumexp <- function(x) log(sum(exp(x - max(x)))) + max(x)
        
        m_C = N - cumsum(n_C)
        part1 = sum(lbeta(m_C + alpha_beta + d_beta * c(1:K), n_C + 1 - d_beta))
        if (p_b == 1) return(part1)

        choice1 = log(p_b)
        choice2 = log(1 - p_b) + pbeta(eps, alpha_beta + m_C, n_C + 1, log.p = T) - alpha_beta * log(eps) 
        part2 = sum(apply(cbind(choice1, choice2), 1, logsumexp))
        return(part1 + part2)
    }
    
    
    MH_order_idx <- function(C, eps = 1E-3, p_b = 0.9, alpha_beta = 1, d_beta = 0) {
        n_C <- colSums(C)
        idx_prop = order(n_C, decreasing = T)
        n_C_prop = n_C[idx_prop]
        
        if (log(runif(1)) < (partition_prob(n_C_prop, eps, p_b, alpha_beta, d_beta) - partition_prob(n_C, eps, p_b, alpha_beta, d_beta))) {
            idx = idx_prop
        }
        else {
            idx = c(1:K)
        }
        return(idx)
    }
       
    # initialize the parameters
    mu = matrix(rnorm(K * p), ncol = p)
    
    I_p = diag(1, p)   
    SigmasInv = array(0, dim = c(p, p, K))    
    for(k in c(1:K)) {
        SigmasInv[, , k] = rwishart(p, I_p)
    }
   
    w = rdirichlet(1, alpha = rep(1, K))
    b = rep(1, K)
    
    # record trace    
    trace_mu = list()
    trace_sigmas_inv = list()
    trace_nC = list()
    trace_label = list()
    trace_b = list()
    
    for (i in 1:steps){
        
    
        C = updateC(mu, SigmasInv, w)
        C_label = colSums((t(C) * c(1:K)))
        n_C <- colSums(C)
        mu = updateMu(C, SigmasInv, mu_prior_mean = mu_prior_mean, mu_prior_SigmasInv = mu_prior_SigmasInv)
        SigmasInv = updateSigmasInv(C, mu, SigmasInv_prior_nu = SigmasInv_prior_nu, SigmasInv_prior_VInv = SigmasInv_prior_VInv)
        
        # MH step to re-order components
        idx = MH_order_idx(C, eps = eps, p_b = p_b, alpha_beta = alpha_beta, d_beta = d_beta)
        mu = mu[idx, ]
        SigmasInv = SigmasInv[, , idx]
        C = C[, idx]
        if (p_b < 1) b <- updateB(C, eps = eps, p_b = p_b, alpha_beta = alpha_beta)
        beta <- updateBeta(C, b, alpha_beta = alpha_beta, d_beta = d_beta, eps = eps)
        w <- updateW(b, beta)
        
        if (i > burnins) {
            idx = i - burnins
            trace_mu[[idx]] = mu
            trace_b[[idx]] = b
            trace_sigmas_inv[[idx]] = SigmasInv
            trace_nC[[idx]] = n_C
            trace_label[[idx]] = C_label
        }
    }
    info = paste("p_b = ", p_b, ", alpha_beta = ", alpha_beta, ", d_beta = ", d_beta, sep = "")
    return(list(info = info, mu = trace_mu, b = trace_b, SigmasInv = trace_sigmas_inv, nC = trace_nC, label = trace_label)) 
}