library(rstiefel)
source('analysis_functions.R')
library(tidyr)
library(ggplot2)
library(colorspace)
library(tidyverse)
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
y_tilde =  g_yt(t = tr) + u %*% gamma + rnorm(n, mean = 0, sd = sqrt(sigma2_y))
y = ifelse(y_tilde > 0, 1, 0)  ## binary Y
tr = data.frame(tr); colnames(tr) = c('t1', 't2', 't3', 't4')
y = as.numeric(y)

# theoretical values -------------------------------------------------------------
coef_mu_u_t <- t(B) %*% solve(B %*% t(B) + sigma2_t * diag(k))
sigma_u_t <- as.numeric(sqrt(1 - t(B) %*% solve(B %*% t(B) + sigma2_t * diag(k)) %*% B))
sigma_ytilde_t <- as.numeric(sqrt( gamma^2*sigma_u_t^2 + sigma2_y ))
sigma_ytilde_t_do <- as.numeric(sqrt( gamma^2 + sigma2_y ))

t_choice <- diag(k)
t2 = rep(0, k)

# true treatment effect #
# numerically #
cal_y_mean_do_numerical <- function(t) {
  dot_mat <- matrix(rep(t, n), ncol=k, nrow=n, byrow = TRUE)
  ytilde_dot =  g_yt(t = dot_mat) + u %*% gamma + rnorm(n, mean = 0, sd = sqrt(sigma2_y))
  y_dot <- ifelse(ytilde_dot > 0, 1, 0)
  return(mean(y_dot))
}
y_mean_do_numerical <- apply(rbind(t_choice, t2), 1, cal_y_mean_do_numerical)
# analytically #
ytilde_mean_do <- g_yt(rbind(t_choice, t2)) 
y_mean_do <- c(pnorm(ytilde_mean_do/sigma_ytilde_t_do))
y_mean_do
effect_true <- y_mean_do[1:4] - y_mean_do[5]
effect_true

# true treatment effect bias #
ytilde_mean_do_bias <- c(rbind(t_choice, t2) %*% t(coef_mu_u_t) %*% gamma)

# true observed treatment effect #
ytilde_mean_obs <- ytilde_mean_do + ytilde_mean_do_bias
y_mean_obs <- c(pnorm(ytilde_mean_obs/sigma_ytilde_t))
y_mean_obs 
effect_obs <- y_mean_obs[1:4] - y_mean_obs[5]
effect_obs

# Calculate E(y|do(t)) by Algm1 using parameters' value in theory#
nsim = 4000
cal_y_mean_calibrated_algm1 <- function(i) {
  print(i)
  t <- rbind(t_choice, t2)[i,]
  mu_i <- (gamma/sigma_ytilde_t) * (as.matrix(tr) %*% t(coef_mu_u_t) - c(coef_mu_u_t %*% t)) %>% as.numeric()
  ytilde_samples <-rnorm(n = nrow(tr)*nsim, mean = mu_i, sd = 1)
  mu_y_t <- y_prob_obs[i]
  y_samples <- ifelse(pnorm(ytilde_samples) > 1-mu_y_t, 1, 0)
  return(mean(y_samples))
}
mean_sens_algm <- sapply(1:5, cal_y_mean_calibrated_algm1)
mean_sens_algm # 0.3677260 0.3677877 0.6302625 0.4910116 0.4990315
y_prob_true # 0.3683086 0.3683086 0.6316914 0.4919498 0.5000000

