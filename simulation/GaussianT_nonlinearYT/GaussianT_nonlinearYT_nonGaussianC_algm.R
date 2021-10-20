library(rstiefel)
source('analysis_functions.R')
library(tidyr)
library(ggplot2)
library(colorspace)
library(tidyverse)
library(colorspace)
library(BART)

# GDP ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
k <- 4
s <- 1
set.seed(123)
B = c(2, 0.5, -0.4, 0.2)
gamma <- 2.8
sigma2_t <- 1
sigma2_y <- 1
tau_l <- c(3, -1, 1, -0.06)
tau_nl <- -4
coef_true <- c(tau_l, tau_nl)
g_yt <- function(t) { # t is n by k matrix
  t[,3] = ifelse(t[,3] > 0, t[,3], 0.7*t[,3])
  t %*% tau_l + t[,1]^2 * tau_nl
}

# model #
set.seed(234)
n = 80000
u = MASS::mvrnorm(n = n, mu = rep(0, s), Sigma = diag(s))
tr = u %*% t(B) + MASS::mvrnorm(n = n, mu = rep(0, k), Sigma = sigma2_t*diag(k))
y =  g_yt(t = tr) + u %*% gamma + rnorm(n, mean = 0, sd = sqrt(sigma2_y))
tr = data.frame(tr); colnames(tr) = c('t1', 't2', 't3')
y = as.numeric(y)

# theoretical values ------------------------------------------------------------------------------------------
coef_mu_u_t <- t(B) %*% solve(B %*% t(B) + sigma2_t * diag(k))
sigma_u_t <- as.numeric(sqrt(1 - t(B) %*% solve(B %*% t(B) + sigma2_t * diag(k)) %*% B))
sigma_y_t <- as.numeric(sqrt( gamma^2*sigma_u_t^2 + sigma2_y ))

t_choice <- diag(k)
t2 = rep(0, k)

# true treatment effect #
effect_true <- c(g_yt(t_choice) - as.numeric(g_yt(t(t2))))
effect_true

# true treatment effect bias #
effect_bias <- c(t_choice %*% t(coef_mu_u_t) %*% gamma)
effect_bias

# true observed treatment effect #
(effect_obs <- effect_true + effect_bias)


# Algorithm
# true copula #
norm_cop <- copula::normalCopula(gamma*sigma_u_t/sigma_y_t)

cal_doMean_algm <- function(t, N = 1e6) {
  message(cat("Case: T = ", t))
  
  ## sample tildet ~ f(t)
  mu_u_tildet <- c(as.matrix(tr[sample(1:n, size = N, replace = TRUE), ]) %*%
                     t(coef_mu_u_t))
  mu_u_t <- c(coef_mu_u_t %*% t)
  ## sample u from f(u|tildet)
  u_sample <- rnorm(N, mean = 0, sd = 1)  
  Uu_sample <- pnorm(u_sample, mean = mu_u_t, sd = sigma_u_t)
  ## sample y from f(y|t)
  mu_y_t <- c(g_yt(t(t)) + t %*% t(coef_mu_u_t) %*% gamma)
  Uy_sample <- copula::cCopula(cbind(Uu_sample, runif(N)), 
                              copula = norm_cop, inverse = TRUE)[,2]
  y_sample <- qnorm(Uy_sample, mean = mu_y_t, sd = sigma_y_t)[Uy_sample!=1]
  
  ## calculate mean(y*f(y|t))
  f_y_t <- dnorm(y_sample, mean = mu_y_t, sd = sigma_y_t)
  sum(y_sample*f_y_t)/sum(f_y_t)
}

cal_doMean_algm <- function(t, N = 1e6) {
  message(cat("Case: T = ", t))
  
  ## sample tildet ~ f(t)
  mu_u_tildet <- c(as.matrix(tr[sample(1:n, size = N, replace = TRUE), ]) %*%
                     t(coef_mu_u_t))
  mu_u_t <- c(coef_mu_u_t %*% t)
  ## sample u from f(u|tildet)
  u_sample <- rnorm(N, mean = mu_u_tildet, sd = sigma_u_t)  
  Uu_sample <- pnorm(u_sample, mean = mu_u_t, sd = sigma_u_t)
  ## sample y from f(y|t)
  mu_y_t <- c(g_yt(t(t)) + t %*% t(coef_mu_u_t) %*% gamma)
  y_sample <- rnorm(N, mean = mu_y_t, sd = sigma_y_t)  
  Uy_sample <- pnorm(y_sample, mean = mu_y_t, sd = sigma_y_t)
  ## calculate mean(y*c)
  C_yu_t <- copula::dCopula(cbind(Uy_sample, Uu_sample), norm_cop)
  sum(y_sample*C_yu_t)/sum(C_yu_t)
}

result_rep <- rep(NA, 10)
for (i in 1:length(result_rep)){
  print(i)
  result_rep[i] <- cal_doMean_algm(c(0,1,0,0))
}
result_rep





# start_time <- Sys.time()
cal_doMean_algm(t = c(0,0,1,0))
# end_time <- Sys.time()
# end_time - start_time
# 
# N_seq <- 10^(1:8)
# doMean_t1_Nseq <- rep(NA, length(N_seq))
# for(i in 1:length(N_seq)){
#   doMean_t1_Nseq[i] <- cal_doMean_algm(t = c(1,0,0,0), N = N_seq[i])
# }
# doMean_t1_Nseq


# Latent Confounder Model ---------------------------------------------------------------------------

fit_ut <- extract_B(cov(tr), nP = 1)
coef_mu_u_t_hat <- fit_ut$coef_mu_u_zt_hat
sigma_u_t_hat <- sqrt(fit_ut$Sigma_u_zt_hat)
u_hat <- as.matrix(tr) %*% t(fit_ut$coef_mu_u_zt_hat) ## n*s

# Observed Outcome model -----------------------------------------
load("simulation/GaussianT_nonlinearYT/bartfit_obs.Rdata")
sigma_y_t_hat <- mean(bartfit_obs$sigma[-(1:100)])
sigma_y_t_hat

cal_ate_obs <- function(t1, bartfit_obs, t2 = rep(0, k)) {
  y1_npost <- predict(bartfit_obs, t(t1))
  y2_npost <- predict(bartfit_obs, t(t2))
  y1_npost - y2_npost
}

ate_obs_npost <- t(apply(t_choice, 1, cal_ate_obs, bartfit_obs = bartfit_obs))

cal_effect_npost_summary <- function(effect_npost) { ## effect_npost is of a single comparison
  mean <- mean(effect_npost)
  q025 <- quantile(effect_npost, prob = 0.025)
  q975 <- quantile(effect_npost, prob = 0.975)
  result <- c(mean, q025, q975)
  names(result) <- c("mean", "q025", "q975")
  return(result)
}

ate_obs_df <- t(apply(ate_obs_npost, 1, cal_effect_npost_summary))
cbind(ate_obs_df, effect_obs)
ate_obs_df[2,1] <- mean(ate_obs_npost[2,1:40])


# Sensitivity Analysis Gaussian Copula---------------------------------------------------------------

ate_cali_mat <- CopSens::gcalibrate(y, tr, t1 = t_choice, t2 = matrix(0, ncol = k, nrow = k),
                                    calitype = "worstcase", mu_y_dt = effect_obs,
                                    sigma_y_t = sigma_y_t,
                                    mu_u_dt = t_choice %*% t(coef_mu_u_t),
                                    cov_u_t = as.matrix(sigma_u_t^2), R2 = 1)$est_df
bound_df <- tibble(x1 = 1:4,
                   y1 = ate_cali_mat[,'R2_1_lwr'],
                   x2 = 1:4,
                   y2 = ate_cali_mat[,'R2_1_upr'])

true_df <- tibble(case = 1:nrow(t_choice),
                  true=c(effect_true),
                  group = rep(8, nrow(t_choice)))
true_df$group <- factor(true_df$group)

plot_nonlinearYT_ate <- tibble(SR1 = ate_cali_mat[,'R2_1_lwr'],
                               SR0 = ate_cali_mat[,'R2_0'],
                               SR_1 = ate_cali_mat[,'R2_1_upr'],
                               case=1:nrow(t_choice)) %>%
  gather(key = "Type", value = "effect", - case) %>%
  ggplot() +
  ungeviz::geom_hpline(aes(x = case, y = effect, col = Type), width = 0.2, size = 1)  +
  geom_point(data = true_df, aes(x = case, y = true, shape = group), size = 2) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_segment(data = bound_df, aes(x=x1,y=y1,xend=x2,yend=y2), size = 0.5) +
  scale_shape_manual(name = "True Effect", values = 8, labels = "") +
  scale_colour_manual(name = "Calibrated",
                      values = c(divergingx_hcl(7,palette = "Zissou 1")[c(7,4,1)]),
                      labels = c(expression(R[paste(tilde(Y), '~', U, '|', T)]^2~'= 1, upper'),
                                 expression(R[paste(tilde(Y), '~', U, '|', T)]^2~'= 0'),
                                 expression(R[paste(tilde(Y), '~', U, '|', T)]^2~'= 1, lower'))) +
  labs(title = bquote(PATE[paste(t[1], ",", t[2])]~"for Gaussian Outcome"),
       y = expression('Causal Effect'), x = 'i') +
  xlim(0.9, 4.4) +
  annotate(geom = "text", x = 1:4 + 0.3, y = c(ate_cali_mat[,'R2_0']), size = 3,
           label = c('0%', 'robust', 'robust', '8.95%'))+
  theme_bw(base_size = 15) +
  theme(plot.title = element_text(hjust = 0.5),
        legend.text.align = 0)
print(plot_nonlinearYT_ate)

# Sensitivity Analysis NonGaussian Copula---------------------------------------------------------------

# Gumbel copula (a.k.a. Gumbel-Hougard copula): an asymmetric Archimedean copula, exhibiting greater dependence in the positive tail than in the negative.
# copula_param of Gumbel: [1, +inf]


cal_doMean_Gumbel <- function(copula_param, t, N = 1e6) {
  
  message(cat("Case: T = ", t))
  
  ## sample tildet ~ f(t)
  mu_u_tildet <- c(as.matrix(tr[sample(1:n, size = N, replace = TRUE), ]) %*%
                     t(coef_mu_u_t))
  mu_u_t <- c(coef_mu_u_t %*% t)
  ## sample u from f(u|tildet)
  u_sample <- rnorm(N, mean = mu_u_tildet, sd = sigma_u_t)  
  Uu_sample <- pnorm(u_sample, mean = mu_u_t, sd = sigma_u_t)
  ## sample y from f(y|t)
  mu_y_t <- c(g_yt(t(t)) + t %*% t(coef_mu_u_t) %*% gamma)
  y_sample <- rnorm(N, mean = mu_y_t, sd = sigma_y_t)  
  Uy_sample <- pnorm(y_sample, mean = mu_y_t, sd = sigma_y_t)
  ## calculate mean(y*c) Gaussian
  C_yu_t_g <- copula::dCopula(cbind(Uy_sample, Uu_sample), norm_cop)
  ## calculate mean(y*c) nonGaussian
  gumbelC <- copula::gumbelCopula(copula_param, dim = 2)
  C_yu_t <- copula::dCopula(cbind(Uy_sample, Uu_sample), gumbelC)
  result <- c(sum(y_sample*C_yu_t_g)/sum(C_yu_t_g), sum(y_sample*C_yu_t)/sum(C_yu_t))
  names(result) <- c("Gaussian", paste0("Gumbel_", copula_param))
  result
}

t_interest <- rbind(t_choice, t2)

cal_doMean_Gumbel_t_interest <- function(copula_param, N = 1e6, n_rep = 10){
  result <- matrix(NA, nrow = nrow(t_interest), ncol = 2)
  colnames(result) <- c("Gaussian", paste0("Gumbel_", copula_param))
  for(case in 1:nrow(t_interest)){
    t <- t_interest[case,]
    result_case_rep <- matrix(NA, nrow = n_rep, ncol = 2)
    for (i in 1:n_rep){
      print(i)
      result_case_rep[i,] <- cal_doMean_Gumbel(copula_param = copula_param, 
                                               t = t)
    }
    result[case,] <- result_case_rep[which.min(abs(result_case_rep[,1]-c(effect_true,0)[case])),]
  }
  result
}

cali_gumbel <- cal_doMean_Gumbel_t_interest(copula_param = 2)

cali_gumbel_df <- tibble(case = 1:nrow(t_choice),
                         cali_gumbel = cali_gumbel[1:4,2] - cali_gumbel[5,2],
                         group = factor(rep(17, nrow(t_choice))))

plot_nonlinearYT_ate_Gumbel <- tibble(SR1 = ate_cali_mat[,'R2_1_lwr'],
                               SR0 = ate_cali_mat[,'R2_0'],
                               SR_1 = ate_cali_mat[,'R2_1_upr'],
                               case=1:nrow(t_choice)) %>%
  gather(key = "Type", value = "effect", - case) %>%
  ggplot() +
  ungeviz::geom_hpline(aes(x = case, y = effect, col = Type), width = 0.2, size = 1)  +
  geom_point(data = true_df, aes(x = case, y = true, shape = group), size = 2) +
  # geom_point(data = cali_gumbel_df, aes(x = case, y = cali_gumbel), size = 2, shape=17) +
  ungeviz::geom_hpline(data = cali_gumbel_df,
                       aes(x = case, y = cali_gumbel, col = group), 
                       width = 0.2, size = 1) + 
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_segment(data = bound_df, aes(x=x1,y=y1,xend=x2,yend=y2), size = 0.5) +
  scale_shape_manual(name = "True Effect", values = 8, labels = "") +
  scale_colour_manual(name = "Calibrated",
                      values = c("#7CAE00", divergingx_hcl(7,palette = "Zissou 1")[c(7,4,1)]),
                      labels = c(expression("Gumbel, "~theta~"= 2"),
                                 expression("Gaussian, "~R[paste(tilde(Y), '~', U, '|', T)]^2~'= 1, upper'),
                                 expression("Gaussian, "~R[paste(tilde(Y), '~', U, '|', T)]^2~'= 0'),
                                 expression("Gaussian, "~R[paste(tilde(Y), '~', U, '|', T)]^2~'= 1, lower'))) +
  labs(title = bquote(PATE[paste(t[1], ",", t[2])]~"for Gaussian Outcome"),
       y = expression('Causal Effect'), x = 'i') +
  xlim(0.9, 4.4) +
  # annotate(geom = "text", x = 1:4 + 0.3, y = c(ate_cali_mat[,'R2_0']), size = 3,
  #          label = c('0%', 'robust', 'robust', '8.95%'))+
  theme_bw(base_size = 15) +
  theme(plot.title = element_text(hjust = 0.5),
        legend.text.align = 0)
print(plot_nonlinearYT_ate_Gumbel)





start_time <- Sys.time()
cal_doMean_Gumbel(copula_param = 2, t = c(1,0,0,0))
end_time <- Sys.time()
end_time - start_time



############### 1 ######################
# cal_doMean_Gumbel <- function(copula_param, case,
#                               M = 1000, N = 100) {
#     mu_y_t <- ate_obs_df[case,1]
#     Uy_sample <- runif(M)
#     y_sample <- qnorm(Uy_sample, mean = mu_y_t, sd = sigma_y_t_hat)
#     mu_u_t_obs <- c(as.matrix(tr) %*% t(coef_mu_u_t_hat)) 
#     mu_u_t <-  c(coef_mu_u_t_hat%*% t_choice[case,])
#     u_sample <- sample(rnorm(n*N, mean = mu_u_t_obs, sd = sigma_u_t_hat), M)
#     Uu_sample <- pnorm(u_sample, mean = mu_u_t, sd = sigma_u_t_hat)
#     gumbelC <- gumbelCopula(copula_param, dim = 2)
#     
#     # W <- rep(NA, M)
#     # for (k in 1:M){
#     #   print(k)
#     #   Cij <- dCopula(u = cbind(rep(Uy_sample[k], n*N), Uu_sample), gumbelC)
#     #   W[k] <- mean(Cij)
#     # }
#     return(mean(y_sample*W))
# }

############### 2 ######################
# cal_doMean_algm <- function(t, N = 100) {
#   message(cat("Case: T = ", t))
#   mu_u_tildet <- c(as.matrix(tr) %*% t(coef_mu_u_t)) ## sample tildet ~ f(t)
#   mu_u_t <- c(coef_mu_u_t %*% t)
#   u_sample <- rnorm(N*n, mean = mu_u_tildet, sd = sigma_u_t)  ## sample u from f(u|tildet)
#   Uu_sample <- pnorm(u_sample, mean = mu_u_t, sd = sigma_u_t)
#   
#   mu_y_t <- c(g_yt(t(t)) + t %*% t(coef_mu_u_t) %*% gamma)
#   y_sample <- rnorm(N*n, mean = mu_y_t, sd = sigma_y_t)  ## sample y from f(y|t)
#   Uy_sample <- pnorm(y_sample, mean = mu_y_t, sd = sigma_y_t)
#   
#   C_yu_t <- copula::dCopula(cbind(Uy_sample, Uu_sample), norm_cop)
#   mean(y_sample*C_yu_t)
# }

############### 3 ######################
# cal_doMean_algm <- function(t, N = 1e6) {
#   message(cat("Case: T = ", t))
#   
#   ## sample tildet ~ f(t)
#   mu_u_tildet <- c(as.matrix(tr[sample(1:n, size = N, replace = TRUE), ]) %*%
#                      t(coef_mu_u_t))
#   mu_u_t <- c(coef_mu_u_t %*% t)
#   ## sample u from f(u|tildet)
#   u_sample <- rnorm(N, mean = mu_u_tildet, sd = sigma_u_t)  
#   Uu_sample <- pnorm(u_sample, mean = mu_u_t, sd = sigma_u_t)
#   ## sample y from f(y|t)
#   mu_y_t <- c(g_yt(t(t)) + t %*% t(coef_mu_u_t) %*% gamma)
#   y_sample <- rnorm(N, mean = mu_y_t, sd = sigma_y_t)  
#   Uy_sample <- pnorm(y_sample, mean = mu_y_t, sd = sigma_y_t)
#   ## calculate mean(y*c)
#   C_yu_t <- copula::dCopula(cbind(Uy_sample, Uu_sample), norm_cop)
#   sum(y_sample*C_yu_t)/sum(C_yu_t)
# }

############### 4 ######################
# cal_doMean_algm <- function(t, N = 1e6) {
#   message(cat("Case: T = ", t))
#   
#   ## sample tildet ~ f(t)
#   mu_u_tildet <- c(as.matrix(tr[sample(1:n, size = N, replace = TRUE), ]) %*%
#                      t(coef_mu_u_t))
#   mu_u_t <- c(coef_mu_u_t %*% t)
#   ## sample u from f(u|tildet)
#   u_sample <- rnorm(N, mean = 0, sd = 1)  
#   Uu_sample <- pnorm(u_sample, mean = mu_u_t, sd = sigma_u_t)
#   ## sample y from f(y|t)
#   mu_y_t <- c(g_yt(t(t)) + t %*% t(coef_mu_u_t) %*% gamma)
#   y_sample <- rnorm(N, mean = mu_y_t, sd = sigma_y_t)  
#   Uy_sample <- pnorm(y_sample, mean = mu_y_t, sd = sigma_y_t)
#   ## calculate mean(y*c)
#   C_yu_t <- copula::dCopula(cbind(Uy_sample, Uu_sample), norm_cop)
#   sum(y_sample*C_yu_t)/sum(C_yu_t)
# }





