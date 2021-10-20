library(rstiefel)
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
B = c(2, 0.5, -0.4, 0.2)/10
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
  
  # calculate mean(y_sample)
  mean(y_sample)
}

doMean_gaussian_true <- apply(rbind(t_choice, t2), 1, cal_doMean_algm)

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

cal_doMean_Gumbel <- function(t, copula_param, N = 1e6) {
  message(cat("Case: T = ", t))
  
  ## sample tildet ~ f(t)
  mu_u_tildet <- c(as.matrix(tr[sample(1:n, size = N, replace = TRUE), ]) %*%
                     t(coef_mu_u_t))
  mu_u_t <- c(coef_mu_u_t %*% t)
  ## sample u from f(u|tildet)
  u_sample <- rnorm(N, mean = 0, sd = 1)  
  Uu_sample <- pnorm(u_sample, mean = mu_u_t, sd = sigma_u_t)
  Uu_sample[Uu_sample > 1-1e-3] <- 1-1e-3
  Uu_sample[Uu_sample < 1e-3] <- 1e-3
  ## sample y from f(y|t)
  mu_y_t <- c(g_yt(t(t)) + t %*% t(coef_mu_u_t) %*% gamma)
  gumbelC <- copula::gumbelCopula(copula_param, dim = 2)
  Uy_sample <- copula::cCopula(cbind(Uu_sample, runif(N)), 
                               copula = gumbelC, inverse = TRUE)[,2]
  y_sample <- qnorm(Uy_sample, mean = mu_y_t, sd = sigma_y_t)[(Uy_sample!=1)&(Uy_sample!=0)]
  
  # calculate mean(y_sample)
  mean(y_sample)
}

cali_gumbel <- apply(rbind(t_choice, t2), 1, cal_doMean_Gumbel, copula_param = 2)
  
  
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






# copula_param of Clayton: [-1, +inf]\{0}

cal_doMean_Clayton <- function(t, copula_param, N = 1e6) {
  message(cat("Case: T = ", t))
  
  ## sample tildet ~ f(t)
  mu_u_tildet <- c(as.matrix(tr[sample(1:n, size = N, replace = TRUE), ]) %*%
                     t(coef_mu_u_t))
  mu_u_t <- c(coef_mu_u_t %*% t)
  ## sample u from f(u|tildet)
  u_sample <- rnorm(N, mean = 0, sd = 1)  
  Uu_sample <- pnorm(u_sample, mean = mu_u_t, sd = sigma_u_t)
  Uu_sample[Uu_sample > 1-1e-3] <- 1-1e-3
  Uu_sample[Uu_sample < 1e-3] <- 1e-3
  ## sample y from f(y|t)
  mu_y_t <- c(g_yt(t(t)) + t %*% t(coef_mu_u_t) %*% gamma)
  claytonC <- copula::claytonCopula(copula_param, dim = 2)
  Uy_sample <- copula::cCopula(cbind(Uu_sample, runif(N)), 
                               copula = claytonC, inverse = TRUE)[,2]
  y_sample <- qnorm(Uy_sample, mean = mu_y_t, sd = sigma_y_t)
  
  # calculate mean(y_sample)
  mean(y_sample)
}

cali_clayton1 <- apply(rbind(t_choice, t2), 1, cal_doMean_Clayton, copula_param = -0.5)
cali_clayton2 <- apply(rbind(t_choice, t2), 1, cal_doMean_Clayton, copula_param = 2)
cali_clayton3 <- apply(rbind(t_choice, t2), 1, cal_doMean_Clayton, copula_param = 10)


plot_result_clayton <- function(cali_clayton, copula_param){
  cali_clayton_df <- tibble(case = 1:nrow(t_choice),
                            cali_clayton = cali_clayton[1:4] - cali_clayton[5],
                            group = factor(rep(17, nrow(t_choice))))
  plot_nonlinearYT_ate_Clayton <- tibble(SR1 = ate_cali_mat[,'R2_1_lwr'],
                                         SR0 = ate_cali_mat[,'R2_0'],
                                         SR_1 = ate_cali_mat[,'R2_1_upr'],
                                         case=1:nrow(t_choice)) %>%
    gather(key = "Type", value = "effect", - case) %>%
    ggplot() +
    ungeviz::geom_hpline(aes(x = case, y = effect, col = Type), width = 0.2, size = 1)  +
    geom_point(data = true_df, aes(x = case, y = true, shape = group), size = 2) +
    ungeviz::geom_hpline(data = cali_clayton_df,
                         aes(x = case, y = cali_clayton, col = group), 
                         width = 0.2, size = 1) + 
    geom_hline(yintercept = 0, linetype = "dashed") +
    geom_segment(data = bound_df, aes(x=x1,y=y1,xend=x2,yend=y2), size = 0.5) +
    scale_shape_manual(name = "True Effect", values = 8, labels = "") +
    scale_colour_manual(name = "Calibrated",
                        values = c("#7CAE00", divergingx_hcl(7,palette = "Zissou 1")[c(7,4,1)]),
                        labels = c(bquote("Clayton, "~theta~"="~.(copula_param)),
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
  plot_nonlinearYT_ate_Clayton
}

plot_result_clayton(cali_clayton = cali_clayton1, copula_param = -0.5)
plot_result_clayton(cali_clayton = cali_clayton2, copula_param = 2)
plot_result_clayton(cali_clayton = cali_clayton3, copula_param = 10)

ate_cali_clayton_tibble <- cbind(cali_clayton1[1:4]-cali_clayton1[5],
                              cali_clayton2[1:4]-cali_clayton2[5],
                              cali_clayton3[1:4]-cali_clayton3[5]) %>%
  'colnames<-'(paste0("Clayton ", c(-0.5, 2, 10))) %>%
  as_tibble() %>%
  add_column(case = 1:nrow(t_choice))



saveRDS(ate_cali_clayton_tibble,
        file = "/Users/jiajing/Desktop/CopulaSensitivity/simulation/GaussianT_nonlinearYT/ate_cali_clayton_tibble.rds")



























cali_clayton_df <- tibble(case = rep(1:nrow(t_choice), 3),
                         cali_clayton = cali_clayton[1:4] - cali_clayton[5],
                         group = factor(rep(17, nrow(t_choice))))

plot_nonlinearYT_ate_Clayton <- tibble(SR1 = ate_cali_mat[,'R2_1_lwr'],
                                      SR0 = ate_cali_mat[,'R2_0'],
                                      SR_1 = ate_cali_mat[,'R2_1_upr'],
                                      case=1:nrow(t_choice)) %>%
  gather(key = "Type", value = "effect", - case) %>%
  ggplot() +
  ungeviz::geom_hpline(aes(x = case, y = effect, col = Type), width = 0.2, size = 1)  +
  geom_point(data = true_df, aes(x = case, y = true, shape = group), size = 2) +
  ungeviz::geom_hpline(data = cali_clayton_df,
                       aes(x = case, y = cali_clayton, col = group), 
                       width = 0.2, size = 1) + 
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_segment(data = bound_df, aes(x=x1,y=y1,xend=x2,yend=y2), size = 0.5) +
  scale_shape_manual(name = "True Effect", values = 8, labels = "") +
  scale_colour_manual(name = "Calibrated",
                      values = c("#7CAE00", divergingx_hcl(7,palette = "Zissou 1")[c(7,4,1)]),
                      labels = c(expression("Clayton, "~theta~"= 2"),
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
print(plot_nonlinearYT_ate_Clayton)


