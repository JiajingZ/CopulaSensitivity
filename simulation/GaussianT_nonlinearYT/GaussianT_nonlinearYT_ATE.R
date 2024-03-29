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

# Sensitivity Analysis ---------------------------------------------------------------
ate_cali_mat <- CopSens::gcalibrate(y, tr, t1 = t_choice, t2 = matrix(0, ncol = k, nrow = k),
                                    calitype = "worstcase", mu_y_dt = as.matrix(ate_obs_df[, 'mean']),
                                    sigma_y_t = sigma_y_t_hat,
                                    mu_u_dt = t_choice %*% t(coef_mu_u_t_hat),
                                    cov_u_t = as.matrix(sigma_u_t_hat^2), R2 = 1)$est_df


# robustness value when ATE = 0 #
RV_mean <- CopSens::cal_rv(y, tr, t1 = t_choice, t2 = matrix(0, ncol = k, nrow = k),
                           mu_y_dt = as.matrix(ate_obs_df[, 'mean']), sigma_y_t = sigma_y_t_hat,
                           mu_u_dt = t_choice %*% t(coef_mu_u_t_hat), cov_u_t = as.matrix(sigma_u_t_hat^2))
# cal RV_limit for significant ones #
RV_limit <- CopSens::cal_rv(y, tr, t1 = t_choice, t2 = matrix(0, ncol = k, nrow = k),
                            mu_y_dt = as.matrix(apply(abs(ate_obs_df[, 2:3]), 1, min)), sigma_y_t = sigma_y_t_hat,
                            mu_u_dt = t_choice %*% t(coef_mu_u_t_hat), cov_u_t = as.matrix(sigma_u_t_hat^2))[-1]

# Summarizing Results --------------------------------------------------------------------------------------------
# plot #
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
# ggsave("plot_nonlinearYT_ate.pdf", plot = plot_nonlinearYT_ate,
#        width = 150, height = 100, units = "mm", path = "simulation/GaussianT_nonlinearYT")









