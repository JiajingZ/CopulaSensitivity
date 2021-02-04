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

## theoretical values -------------------------------------------------------------
coef_mu_u_t <- t(B) %*% solve(B %*% t(B) + sigma2_t * diag(k))
sigma_u_t <- as.numeric(sqrt(1 - t(B) %*% solve(B %*% t(B) + sigma2_t * diag(k)) %*% B))
sigma_ytilde_t <- as.numeric(sqrt( gamma^2*sigma_u_t^2 + sigma2_y ))
sigma_ytilde_t_do <- as.numeric(sqrt( gamma^2 + sigma2_y ))

t_choice <- diag(k)
t2 = rep(0, k)

# true Treatment effect #
ytilde_mean_do <- g_yt(rbind(t_choice, t2))
y_mean_do <- c(pnorm(ytilde_mean_do/sigma_ytilde_t_do))
y_mean_do
effect_true <- y_mean_do[1:4]/y_mean_do[5]
effect_true

# true treatment effect bias #
ytilde_mean_do_bias <- c(rbind(t_choice, t2) %*% t(coef_mu_u_t) %*% gamma)

# true observed treatment effect #
ytilde_mean_obs <- ytilde_mean_do + ytilde_mean_do_bias
y_mean_obs <- c(pnorm(ytilde_mean_obs/sigma_ytilde_t))
y_mean_obs
effect_obs <- y_mean_obs[1:4]/y_mean_obs[5]
effect_obs

# Latent Confounder Model --------------------------------------------------------------

fit_ut <- extract_B(cov(tr), nP = 1)
coef_mu_u_t_hat <- fit_ut$coef_mu_u_zt_hat
sigma_u_t_hat <- sqrt(fit_ut$Sigma_u_zt_hat)
u_hat <- as.matrix(tr) %*% t(fit_ut$coef_mu_u_zt_hat) ## n*s


# Observed Outcome model -----------------------------------------
load("simulation/GaussianT_BinaryY_nonlinearYT/bartfit_obs.Rdata")

sigma_ytilde_t <- 1 ## by default

cal_rr_obs <- function(t1, bartfit_obs, t2 = rep(0, k)) {
  prob1_npost <- predict(bartfit_obs, t(t1))$prob.test
  prob2_npost <- predict(bartfit_obs, t(t2))$prob.test
  prob1_npost / prob2_npost
}

rr_obs_npost <- t(apply(t_choice, 1, cal_rr_obs, bartfit_obs = bartfit_obs))

cal_effect_npost_summary <- function(effect_npost) { ## effect_npost is of a single comparison
  mean <- mean(effect_npost)
  q025 <- quantile(effect_npost, prob = 0.025)
  q975 <- quantile(effect_npost, prob = 0.975)
  result <- c(mean, q025, q975)
  names(result) <- c("mean", "q025", "q975")
  return(result)
}

rr_obs_df <- t(apply(rr_obs_npost, 1, cal_effect_npost_summary))
cbind(rr_obs_df, effect_obs)

# Sensitivity Analysis for RR----------------------------------------------------------------------------
gamma_max <- sigma_ytilde_t/sigma_u_t_hat
gamma_seq <- seq(-gamma_max, gamma_max, by  = 0.01)
# prob_cali_results <- CopSens::bcalibrate(y = y, tr = tr,
#                                          t = rbind(t_choice, t2),
#                                          gamma = t(gamma_seq), ## different gamma in columns
#                                          mu_y_t = predict(bartfit_obs, rbind(t_choice, t2))$prob.test.mean,
#                                          mu_u_tr = as.matrix(tr) %*% t(coef_mu_u_t_hat),
#                                          mu_u_t = rbind(t_choice, t2) %*% t(coef_mu_u_t_hat),
#                                          cov_u_t = as.matrix(sigma_u_t_hat^2))
# save(prob_cali_results, file = "simulation/GaussianT_BinaryY_nonlinearYT/prob_cali_results.Rdata")
load("simulation/GaussianT_BinaryY_nonlinearYT/prob_cali_results.Rdata")

rr_cali_mat <- as.matrix(prob_cali_results$est_df[1:4,-1]) /
  matrix(rep(as.numeric(prob_cali_results$est_df[5,-1]), 4), nrow = 4, byrow = TRUE)
rr_cali_bound <- data.frame(lwr = apply(rr_cali_mat, 1, min),
                            upr = apply(rr_cali_mat, 1, max))
R2_seq <- prob_cali_results$R2


# robustness value
cal_RV <- function(rr_est) {
  gamma_null <- gamma_seq[which.min(abs(rr_est - 1))]
  gamma_null^2*sigma_u_t_hat^2
}

RV <- apply(rr_cali_mat[c(1,4),], 1, cal_RV) %>% round(digits = 4)
RV

# R^2 VS RR #
plot_r2_rr1 <- tibble(R2 = sign(gamma_seq)*R2_seq, RR = rr_cali_mat[1,]) %>%
  ggplot(aes(x = R2, y = RR)) +
  geom_line() +
  labs(x =expression(R^2), y = expression(RR[paste(e[1],",","0" )])) +
  theme_bw(base_size = 19)
plot_r2_rr1

plot_r2_rr2 <- tibble(R2 = sign(gamma_seq)*R2_seq, RR = rr_cali_mat[2,]) %>%
  ggplot(aes(x = R2, y = RR)) +
  geom_line() +
  labs(x =expression(R^2), y = expression(RR[paste(e[2],",","0" )])) +
  theme_bw(base_size = 19)
plot_r2_rr2

plot_r2_rr3 <- tibble(R2 = sign(gamma_seq)*R2_seq, RR = rr_cali_mat[3,]) %>%
  ggplot(aes(x = R2, y = RR)) +
  geom_line() +
  labs(x =expression(R^2), y = expression(RR[paste(e[3],",","0" )])) +
  theme_bw(base_size = 19)
plot_r2_rr3

plot_r2_rr4 <- tibble(R2 = sign(gamma_seq)*R2_seq, RR = rr_cali_mat[4,]) %>%
  ggplot(aes(x = R2, y = RR)) +
  geom_line() +
  labs(x =expression(R^2), y = expression(RR[paste(e[4],",","0" )])) +
  theme_bw(base_size = 19)
plot_r2_rr4

# ggsave("plot_r2_rr1.pdf", plot = plot_r2_rr1, width = 110, height = 90,
#        units = "mm", path = "simulation/GaussianT_BinaryY_nonlinearYT")
# ggsave("plot_r2_rr2.pdf", plot = plot_r2_rr2, width = 110, height = 90,
#        units = "mm", path = "simulation/GaussianT_BinaryY_nonlinearYT")
# ggsave("plot_r2_rr3.pdf", plot = plot_r2_rr3, width = 110, height = 90,
#        units = "mm", path = "simulation/GaussianT_BinaryY_nonlinearYT")
# ggsave("plot_r2_rr4.pdf", plot = plot_r2_rr4, width = 110, height = 90,
#        units = "mm", path = "simulation/GaussianT_BinaryY_nonlinearYT")


# Summarising Results --------------------------------------------------------------------------------------------

# plot #
bound_df <- tibble(x1 = 1:4,
                   y1 = rr_cali_bound$lwr,
                   x2 = 1:4,
                   y2 = rr_cali_bound$upr)

true_df <- tibble(case = 1:4,
                  true = effect_true,
                  group = rep(1, 4))
true_df$group = factor(true_df$group)

summarise_df <- tibble(Upper=rr_cali_bound$upr,
                       Naive=rr_obs_df[,'mean'],
                       Lower=rr_cali_bound$lwr,
                       case=1:nrow(t_choice)) %>%
  gather(key = "Type", value = "effect", - case)
summarise_df$Type <- factor(summarise_df$Type,
                            levels = c('Upper', 'Naive', 'Lower'))

plot_nonlinearYT_binaryY_RR <-  ggplot(summarise_df) +
  ungeviz::geom_hpline(aes(x = case, y = effect, col = Type), width = 0.2, size = 1)  +
  geom_hline(yintercept = 1, linetype = "dashed") +
  geom_point(data = true_df, aes(x = case, y = true, shape = group), size = 2) +
  geom_segment(data = bound_df, aes(x=x1,y=y1,xend=x2,yend=y2), size = 0.5) +
  scale_shape_manual(name = "True Effect", values = 8, labels = "") +
  scale_colour_manual(name = "Calibrated",
                      values = c("#F5191C", "#EACB2B", "#3B99B1")) +
  guides(colour = guide_legend(order = 2),
         shape = guide_legend(order = 1)) +
  xlim(0.9,4.4) +
  labs(title = bquote(RR[paste(e[i], ",", 0)]~'for Binary Outcome'),
       y = expression('Causal Effect'), x = 'i') +
  annotate(geom = "text", x = 1:4 + 0.3, y = c(rr_obs_df[,'mean']), size = 3,
           label = c(paste0(RV[1]*100,"%"), 'robust', 'robust', paste0(RV[2]*100,"%")))+
  theme_bw(base_size = 15) +
  theme(plot.title = element_text(hjust = 0.5),
        legend.text.align = 0)
print(plot_nonlinearYT_binaryY_RR)
# ggsave("plot_nonlinearYT_binaryY_RR.pdf", plot = plot_nonlinearYT_binaryY_RR,
#        width = 137, height = 100, units = "mm", path = "simulation/GaussianT_BinaryY_nonlinearYT")












