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
y =  g_yt(t = tr) + u %*% gamma + rnorm(n, mean = 0, sd = sqrt(sigma2_y))
tr = data.frame(tr); colnames(tr) = c('t1', 't2', 't3', 't4')
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

# Latent Confounder Model ---------------------------------------------------------------------------

fit_ut <- extract_B(cov(tr), nP = 1)
coef_mu_u_t_hat <- fit_ut$coef_mu_u_zt_hat
sigma_u_t_hat <- sqrt(fit_ut$Sigma_u_zt_hat)
u_hat <- as.matrix(tr) %*% t(fit_ut$coef_mu_u_zt_hat) ## n*s

# Observed Outcome model ---------------------------------------------------------------------------------

bartfit_obs <- BART::wbart(x.train = tr, y.train = y, x.test = rbind(t_choice,t2),
                           sparse = TRUE, rm.const=TRUE, rho = 4,
                           ndpost=1000L, keepevery = 200L) 
save(bartfit_obs, file = "simulation/GaussianT_nonlinearYT/bartfit_obs.Rdata")









