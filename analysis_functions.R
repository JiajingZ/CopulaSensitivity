## use standarized treatment for analysis
## return matrix B, sigma2_t and A
## Para: Cor: If treatment is Gaussian, C = cov(tr); 
#             If treatment is nongaussian, C is gaussian copula-based correlation matrix
extract_B <- function(C, nP = NA) {
  k <- ncol(C)
  eigen_c <- eigen(C)
  if (is.na(nP)) {nP <- max(which(eigen_c$values > 1))}
  
  if (nP == 1) {
    beta_hat <- eigen_c$vectors[, 1] * sqrt(eigen_c$values[1])
    sigma2_zt_hat <- mean(eigen_c$values[2:k])
    Sigma_Zt_hat <- beta_hat %*% t(beta_hat) + sigma2_zt_hat * diag(k)
    var_u_zt_hat <- as.numeric(diag(nP) - t(beta_hat) %*% solve(Sigma_Zt_hat) %*% beta_hat)
    coef_mu_u_zt_hat <- t(beta_hat) %*% solve(Sigma_Zt_hat)
    return(list(B_hat = beta_hat,
                sigma2_zt_hat = sigma2_zt_hat,
                Sigma_u_zt_hat = var_u_zt_hat,
                coef_mu_u_zt_hat = coef_mu_u_zt_hat))
  } else if (nP > 1) {
    B_hat <- eigen_c$vectors[, 1:nP] %*% diag(sqrt(eigen_c$values[1:nP])) ## columns of B_hat are orthogonal
    sigma2_zt_hat <- mean(eigen_c$values[(nP+1):k])
    Sigma_Zt_hat <- B_hat %*% t(B_hat) + sigma2_zt_hat * diag(k)
    coef_mu_u_zt_hat <- t(B_hat) %*% solve(Sigma_Zt_hat)
    Sigma_u_zt_hat <- diag(nP) - t(B_hat) %*% solve(Sigma_Zt_hat) %*% B_hat
    eigen_Sigma <- eigen(Sigma_u_zt_hat)
    Q_hat <- eigen_Sigma$vectors
    D_hat <- eigen_Sigma$values
    Sigma_inv_half <- Q_hat %*% diag(D_hat^{-1/2}) %*% t(Q_hat)
    A_hat <- Sigma_inv_half %*% coef_mu_u_zt_hat
    return(list(B_hat = B_hat,
                sigma2_zt_hat = sigma2_zt_hat,
                A_hat = A_hat,
                Sigma_u_zt_hat = Sigma_u_zt_hat,
                Sigma_inv_half = Sigma_inv_half,
                coef_mu_u_zt_hat = coef_mu_u_zt_hat))
  } else { cat("Error: unreasonable num. of principals \n") }
}

############################################################
############################################################
# getting bounds of truncated normal for any given t #
tmv_bound <- function(t) {
  nt = length(t)
  for (i in 1:nt) {
    if (t[i] != 1 & t[i] != 0) { stop("Error: t must be binary") }
  }
  a <- rep(NA, nt)
  b <- rep(NA, nt)
  for (i in 1:nt) {
    if (t[i] == 0) {
      a[i] <- -Inf
      b[i] <- 0
    } else if (t[i] == 1) {
      a[i] <- 0
      b[i] <- Inf
    } else { cat("Error: t must be binary") }
  }
  return(list(a=a, b=b))
}

###############################################################
##############################################################
## Y is Gaussian 
## input: t: treatment of interest
#         tau0: regression coefficient of T for Y or Zt for Zy 
#         sigma: standarded deviation of conditional dist. Y|T or Zy|Zt
#         A: 
#         GaussianT: logical, whether treatment is Gaussian or not
#         cov_zt: If GaussianT = F, cov_zt needs to be imputed

est_mean_bound <- function(t, tau0, sigma, A, GaussianT = T, cov_zt) {
  if (GaussianT) {
    biased <- t(tau0) %*% t
    lower <- biased - sigma * sqrt(sum((A %*% t)^2))
    upper <- biased + sigma * sqrt(sum((A %*% t)^2))
    result <- c(biased, lower, upper)
  } else {
    zt_bounds <- tmv_bound(t)
    a <- zt_bounds$a
    b <- zt_bounds$b
    zt_sample <- rtmvnorm(n = 1000, mean = rep(0, ncol(cov_zt)), sigma = cov_zt, 
                          lower = a, upper = b, algorithm = "gibbs")
    biased <- mean(zt_sample %*% tau0)
    Azt <- A%*%t(zt_sample)  ## s by 1000
    norm2_Azt <- sqrt(sum(apply(Azt, 1, mean)^2))
    lower <- biased - sigma * norm2_Azt
    upper <- biased + sigma * norm2_Azt
    result <- c(biased, lower, upper)
  }
  names(result) = c("biased", "lb", "ub")
  return(result)
}
  
  
########################################
#########################################
## Calibrated Bound ##
# bound: returns from "est_mean_bound"， bound matrix with "biased", "lb", "ub" 
# r2: fraction with 2 decimal place, between 0 and 1. partial R^2 for U|T in the outcome model.
est_calibrated_bound <- function(bound, r2) {
  lb.cali <- bound[, 1] - sqrt(r2) * (bound[, 1] - bound[, 2])
  ub.cali <- bound[, 1] + sqrt(r2) * (bound[, 3] - bound[, 1])
  result <- cbind(lb.cali, ub.cali)
  colnames(result) <- c(paste0("lb_", r2), paste0("ub_", r2))
  return(result)
}
  
############################################
############################################
# bound: returns from "est_mean_bound"， bound matrix with "biased", "lb", "ub"
# return: percentage, 100*R^2 when effect = 0
get_zero_effect_r2 <- function(bound) {
  half_width <- bound[, 1] - bound[, 2]
  r2 <- round(100 *(abs(bound[, 1]) / half_width)^2, 2)
  return(r2)
}  


#################################################
################################################
## function: calculate the L2 norm of a vector ##
L2_norm <- function(t) {
  sqrt(sum(t^2))
}
  

## standard vector such that L2-norm = 1 ##
scale_L2 <- function(t) {
  t / sqrt(sum(t^2))
}
  
  
  
  











