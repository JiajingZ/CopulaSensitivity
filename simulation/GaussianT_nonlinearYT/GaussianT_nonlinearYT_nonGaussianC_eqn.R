library(rstiefel)
# source('analysis_functions.R')
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


# Sensitivity Analysis NonGaussian Copula---------------------------------------------------------------
# Algorithm
# true copula #
norm_cop <- copula::normalCopula(gamma*sigma_u_t/sigma_y_t)
# Clayton #
readRDS("/Users/jiajing/Desktop/CopulaSensitivity/simulation/GaussianT_nonlinearYT/ate_cali_clayton_tibble.rds") %>%
  gather(key = "Type", value = "effect", - case) -> ate_cali_clayton_narrow
ate_cali_clayton_narrow$Type <- factor(ate_cali_clayton_narrow$Type, 
                                       levels = c("Clayton -0.5","Clayton 2","Clayton 10"),
                                       labels = c("Clayton θ=-0.5","Clayton θ=2","Clayton θ=10"))
# Quadratic #
readRDS("/Users/jiajing/Desktop/CopulaSensitivity/simulation/GaussianT_nonlinearYT/ate_cali_quad_narrow.rds") -> ate_cali_quad_narrow
# Log #
readRDS("/Users/jiajing/Desktop/CopulaSensitivity/simulation/GaussianT_nonlinearYT/ate_cali_log_narrow.rds") -> ate_cali_log_narrow




# cal_doMean_quadratic <- function(t, N = 1e6, a = 1, b = 0) {
#   message(cat("Case: T = ", t))
#   
#   ## sample tildet ~ f(t)
#   mu_u_tildet <- c(as.matrix(tr[sample(1:n, size = N, replace = TRUE), ]) %*%
#                      t(coef_mu_u_t))
#   mu_u_t <- c(coef_mu_u_t %*% t)
#   
#   
#   # u ~ f(u | t)
#   u_sample_cond <- rnorm(N, mean = mu_u_t, sd = sigma_u_t)
#   ## ytilde = f(U) + epsilon, UY relatioship
#   # yu_func <- function(u) {3*sin(10*abs(u)) + rnorm(N, sd=1)}
#   # yu_func <- function(u) {40*abs(u - 0.5) + rnorm(N, sd = 1)}
#   yu_func <- function(u) {a*(u - b)^2}
#   
#   ytilde_sample_cond <- yu_func(u_sample_cond)
#   ytilde_cdf <- ecdf(ytilde_sample_cond)
#   Uy_sample_cond <- ytilde_cdf(ytilde_sample_cond)
#   # u ~ f(u)
#   u_sample_marg <- rnorm(N, mean = mu_u_tildet, sd = sigma_u_t)  
#   ytilde_sample_marg <- yu_func(u_sample_marg)
#   Uy_sample_marg <- ytilde_cdf(ytilde_sample_marg)
#   Uy_sample_marg[Uy_sample_marg > 1-1e-3] <- 1-1e-3
#   Uy_sample_marg[Uy_sample_marg < 1e-3] <- 1e-3
#   
#   ## Uy to Y
#   mu_y_t <- c(g_yt(t(t)) + t %*% t(coef_mu_u_t) %*% gamma)
#   y_sample <- qnorm(Uy_sample_marg, mean = mu_y_t, sd = sigma_y_t)
#   
#   ## calculate mean(y)
#   mean(y_sample)
# }

# plot(pnorm(u_sample_cond, mean = mu_u_t, sd = sigma_u_t), Uy_sample_cond)
# plot(pnorm(u_sample_marg, mean = mu_u_t, sd = sigma_u_t), Uy_sample_marg)


# t_interest <- rbind(t_choice, t2)
# quad_para <- expand_grid(a = c(1, -1), b = c(-2, 0, 2))
# domean_cali_quad_mat <- matrix(NA, nrow = nrow(t_interest), ncol = nrow(quad_para),
#                                dimnames = list(NULL,paste0("quadratic: a = ", quad_para$a, ", b = ", quad_para$b)))
# for(i in 1:nrow(quad_para)){
#   domean_cali_quad_mat[,i] <- apply(t_interest, 1, 
#                                     cal_doMean_quadratic,
#                                     a = quad_para$a[i], 
#                                     b = quad_para$b[i])
# }
# 
# colnames(domean_cali_quad_mat) = paste0("Quadratic: a = ", quad_para$a, ", b = ", quad_para$b)
# ate_cali_quad_narrow <-
#   t(t(domean_cali_quad_mat)[,1:4] - t(domean_cali_quad_mat)[,5]) %>%
#   as_tibble() %>%
#   add_column(case = 1:nrow(t_choice)) %>%
#   gather(key = "Type", value = "effect", - case)
# saveRDS(ate_cali_quad_narrow,
#         file = "/Users/jiajing/Desktop/CopulaSensitivity/simulation/GaussianT_nonlinearYT/ate_cali_quad_narrow.rds")

# cal_doMean_log <- function(t, N = 1e6, a = 1, b = 10) {
#   message(cat("Case: T = ", t))
#   
#   ## sample tildet ~ f(t)
#   mu_u_tildet <- c(as.matrix(tr[sample(1:n, size = N, replace = TRUE), ]) %*%
#                      t(coef_mu_u_t))
#   mu_u_t <- c(coef_mu_u_t %*% t)
#   
#   
#   # u ~ f(u | t)
#   u_sample_cond <- rnorm(N, mean = mu_u_t, sd = sigma_u_t)
#   yu_func <- function(u) {a*log(u + b)}
#   
#   ytilde_sample_cond <- yu_func(u_sample_cond)
#   ytilde_cdf <- ecdf(ytilde_sample_cond)
#   Uy_sample_cond <- ytilde_cdf(ytilde_sample_cond)
#   # u ~ f(u)
#   u_sample_marg <- rnorm(N, mean = mu_u_tildet, sd = sigma_u_t)  
#   ytilde_sample_marg <- yu_func(u_sample_marg)
#   Uy_sample_marg <- ytilde_cdf(ytilde_sample_marg)
#   Uy_sample_marg[Uy_sample_marg > 1-1e-3] <- 1-1e-3
#   Uy_sample_marg[Uy_sample_marg < 1e-3] <- 1e-3
#   
#   ## Uy to Y
#   mu_y_t <- c(g_yt(t(t)) + t %*% t(coef_mu_u_t) %*% gamma)
#   y_sample <- qnorm(Uy_sample_marg, mean = mu_y_t, sd = sigma_y_t)
#   
#   ## calculate mean(y)
#   mean(y_sample)
# }


# domean_cali_log_mat <- cbind(apply(t_interest, 1, cal_doMean_log, a = 1),
#                              apply(t_interest, 1, cal_doMean_log, a = -1)) %>%
#   'colnames<-'(paste0("log: a = ", c(-1, 1)))
# 
# colnames(domean_cali_log_mat) <- paste0("log: a = ", c(-1, 1))
# 
# ate_cali_log_narrow <-
#   t(t(domean_cali_log_mat)[,1:4] - t(domean_cali_log_mat)[,5]) %>%
#   as_tibble() %>%
#   add_column(case = 1:nrow(t_choice)) %>%
#   gather(key = "Type", value = "effect", - case)
# saveRDS(ate_cali_log_narrow,
#         file = "/Users/jiajing/Desktop/CopulaSensitivity/simulation/GaussianT_nonlinearYT/ate_cali_log_narrow.rds")


# plot #
tibble(lwr_g = ate_cali_mat[,'R2_1_lwr'],
       upr_g = ate_cali_mat[,'R2_1_upr'],
       case=1:nrow(t_choice)) %>%
  ggplot() +
  geom_errorbar(aes(ymin = lwr_g, ymax = upr_g, x = case), width = 0.2) + 
  geom_text(data = ate_cali_quad_narrow,
            aes(x = case, y = effect, label = Type), hjust = 0, nudge_x = 0.05) + 
  geom_point(data = ate_cali_quad_narrow,
             aes(x = case, y = effect)) + 
  geom_text(data = ate_cali_log_narrow,
            aes(x = case, y = effect, label = Type), hjust = 0, nudge_x = 0.05) + 
  geom_point(data = ate_cali_log_narrow,
             aes(x = case, y = effect)) + 
  geom_text(data = ate_cali_clayton_narrow,
            aes(x = case, y = effect, label = Type), hjust = 0, nudge_x = 0.05) + 
  geom_point(data = ate_cali_clayton_narrow,
             aes(x = case, y = effect)) +
  # ungeviz::geom_hpline(data = ate_cali_quad_narrow,
  #                      aes(x = case, y = effect, col = Type),
  #                      width = 0.2, size = 1)  +
  # ungeviz::geom_hpline(data = ate_cali_log_narrow,
  #                      aes(x = case, y = effect, col = Type),
  #                      width = 0.2, size = 1)  +
  # ungeviz::geom_hpline(data = ate_cali_clayton_narrow,
  #                      aes(x = case, y = effect, col = Type),
  #                      width = 0.2, size = 1) +
  # scale_colour_manual(name = "Copula Type:",
  #                     values = c(sequential_hcl(11, palette = "Reds 3")[2:7],
  #                                sequential_hcl(7, palette = "Blues")[c(1,3,5)],
  #                                sequential_hcl(7, palette = "Greens")[c(2,4)]),
  #                     labels = c(expression(tilde(y)==-(u+2)^2),
  #                                expression(tilde(y)==-u^2),
  #                                expression(tilde(y)==-(u-2)^2),
  #                                expression(tilde(y)==(u+2)^2),
  #                                expression(tilde(y)== u^2),
  #                                expression(tilde(y)==(u-2)^2),
  #                                bquote("Clayton, "~theta~"= - 0.5"),
  #                                bquote("Clayton, "~theta~"= 2"),
  #                                bquote("Clayton, "~theta~"= 10"),
  #                                expression(tilde(y)==-log(u+10)),
  #                                expression(tilde(y)==log(u+10)))) +
  labs(title = bquote(PATE[paste(t[1], ",", t[2])]~"with Various Copula"),
       y = expression('Causal Effect'), x = 'i') +
  xlim(0.9, 4.4) +
  theme_bw(base_size = 15) +
  theme(plot.title = element_text(hjust = 0.5),
        legend.text.align = 0)






















# draft ---------------------------------------------------------------------------------------------------

ate_cali_nonG_tibble <- tibble(case = 1:nrow(t_choice),
                               cali_nonG = domean_cali_nonGaussian[1:4] - domean_cali_nonGaussian[5],
                               group = factor(rep(17, nrow(t_choice))))

plot_nonlinearYT_ate_nonG <- tibble(SR1 = ate_cali_mat[,'R2_1_lwr'],
                                    SR0 = ate_cali_mat[,'R2_0'],
                                    SR_1 = ate_cali_mat[,'R2_1_upr'],
                                    case=1:nrow(t_choice)) %>%
  gather(key = "Type", value = "effect", - case) %>%
  ggplot() +
  ungeviz::geom_hpline(aes(x = case, y = effect, col = Type), width = 0.2, size = 1)  +
  # geom_point(data = true_df, aes(x = case, y = true, shape = group), size = 2) +
  # geom_point(data = cali_gumbel_df, aes(x = case, y = cali_gumbel), size = 2, shape=17) +
  ungeviz::geom_hpline(data = ate_cali_nonG_tibble,
                       aes(x = case, y = cali_nonG, col = group), 
                       width = 0.2, size = 1) + 
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_segment(data = bound_df, aes(x=x1,y=y1,xend=x2,yend=y2), size = 0.5) +
  scale_shape_manual(name = "True Effect", values = 8, labels = "") +
  scale_colour_manual(name = "Calibrated",
                      values = c("#7CAE00", divergingx_hcl(7,palette = "Zissou 1")[c(7,4,1)]),
                      labels = c(expression("NonGaussian"),
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
print(plot_nonlinearYT_ate_nonG)

geom_errorbar()

# non-monotone
  # quandratic with diff u_cons (-2, 2), +- U^2,perfect dependence
# monotone
  # clayton, 
  # transformation of u, log(u)

