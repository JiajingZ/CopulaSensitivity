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
gamma_seq <- seq(-gamma_max, gamma_max, by  = 0.001)
cal_rr_calibrated <- function(t1, bartfit_obs, t2 = rep(0, k)) {
  ytilde1_obs <- mean(predict(bartfit_obs, t(t1))$yhat.test)
  ytilde2_obs <- mean(predict(bartfit_obs, t(t2))$yhat.test)
  prob1 <- pnorm((ytilde1_obs - gamma_seq*c(coef_mu_u_t_hat %*% t1)) / sqrt(sigma_ytilde_t^2 + gamma_seq^2*(1-sigma_u_t_hat^2)))
  prob2 <- pnorm( ytilde2_obs / sqrt(sigma_ytilde_t^2 + gamma_seq^2*(1-sigma_u_t_hat^2)))
  prob1 / prob2
}

rr_sens_seq <-  t(apply(t_choice, 1, cal_rr_calibrated, bartfit_obs = bartfit_obs))
rr_sens_bound <- cbind(apply(rr_sens_seq, 1, min),
                        apply(rr_sens_seq, 1, max))
colnames(rr_sens_bound) <- c('lwr', 'upr')
cbind(rr_sens_bound, effect_true)


# robustness value
cal_RV <- function(rr_sens_seq) {
  gamma_null <- gamma_seq[which.min(abs(rr_sens_seq-1))]
  gamma_null^2*sigma_u_t_hat^2
}

RV <- apply(rr_sens_seq[c(1,4),], 1, cal_RV) %>% round(digits = 4)
RV

R2_seq <- gamma_seq^2*sigma_u_t_hat^2 # sigma_{y|t}^2 = gamma^2*sigma_{u|t}^2 + sigma_y^2
## Gamma VS R2 #
plot(gamma_seq, R2_seq, type = "l", xlab = expression(gamma), ylab = expression(R[paste(tilde(Y), '~', U, '|', T)]^2))
# Gamma VS RR #
plot(gamma_seq, rr_sens_seq[1,], type = "l", ylab = "RR", xlab =expression(gamma), main = expression(RR[paste(e[1],",","0" )]))
plot(gamma_seq, rr_sens_seq[2,], type = "l", ylab = "RR", xlab =expression(gamma), main = expression(RR[paste(e[2],",","0" )]))
plot(gamma_seq, rr_sens_seq[3,], type = "l", ylab = "RR", xlab =expression(gamma), main = expression(RR[paste(e[3],",","0" )]))
plot(gamma_seq, rr_sens_seq[4,], type = "l", ylab = "RR", xlab =expression(gamma), main = expression(RR[paste(e[4],",","0" )]))
# R^2 VS RR #
plot_r2_rr1 <- tibble(R2 = sign(gamma_seq)*R2_seq, RR = rr_sens_seq[1,]) %>%
  ggplot(aes(x = R2, y = RR)) + 
  geom_line() + 
  labs(x =expression(R^2), y = expression(RR[paste(e[1],",","0" )])) +
  theme_bw(base_size = 19) 
plot_r2_rr1  

plot_r2_rr2 <- tibble(R2 = sign(gamma_seq)*R2_seq, RR = rr_sens_seq[2,]) %>%
  ggplot(aes(x = R2, y = RR)) + 
  geom_line() + 
  labs(x =expression(R^2), y = expression(RR[paste(e[2],",","0" )])) +
  theme_bw(base_size = 19) 
plot_r2_rr2

plot_r2_rr3 <- tibble(R2 = sign(gamma_seq)*R2_seq, RR = rr_sens_seq[3,]) %>%
  ggplot(aes(x = R2, y = RR)) + 
  geom_line() + 
  labs(x =expression(R^2), y = expression(RR[paste(e[3],",","0" )])) +
  theme_bw(base_size = 19) 
plot_r2_rr3 

plot_r2_rr4 <- tibble(R2 = sign(gamma_seq)*R2_seq, RR = rr_sens_seq[4,]) %>%
  ggplot(aes(x = R2, y = RR)) + 
  geom_line() + 
  labs(x =expression(R^2), y = expression(RR[paste(e[4],",","0" )])) +
  theme_bw(base_size = 19) 
plot_r2_rr4

# ggsave("plot_r2_rr1.pdf", plot = plot_r2_rr1, width = 110, height = 90, units = "mm")
# ggsave("plot_r2_rr2.pdf", plot = plot_r2_rr2, width = 110, height = 90, units = "mm")
# ggsave("plot_r2_rr3.pdf", plot = plot_r2_rr3, width = 110, height = 90, units = "mm")
# ggsave("plot_r2_rr4.pdf", plot = plot_r2_rr4, width = 110, height = 90, units = "mm")

 
# Summarising Results --------------------------------------------------------------------------------------------

# plot #
bound_df <- tibble(x1 = 1:4,
                   y1 = rr_sens_bound[,1],
                   x2 = 1:4,
                   y2 = rr_sens_bound[,2])

true_df <- tibble(case = 1:4,
                  true = effect_true,
                  group = rep(1, 4))
true_df$group = factor(true_df$group)

summarise_df <- tibble(Upper=rr_sens_bound[,2],
                       Naive=rr_obs_df[,'mean'],
                       Lower=rr_sens_bound[,1],
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
           label = c('0.01%', 'robust', 'robust', '15.16%'))+
  theme_bw(base_size = 15) + 
  theme(plot.title = element_text(hjust = 0.5),
        legend.text.align = 0)
print(plot_nonlinearYT_binaryY_RR)
ggsave("plot_nonlinearYT_binaryY_RR.pdf", plot = plot_nonlinearYT_binaryY_RR,
       width = 137, height = 100, units = "mm", path = "simulation/GaussianT_BinaryY_nonlinearYT")














