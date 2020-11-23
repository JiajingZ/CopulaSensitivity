## ggplot layout ##
library(ggpubr)
#### Gaussian Multi-U #####
library(rstiefel)
library(ggplot2)
library(tidyverse)
library(patchwork)
library(ggridges)
library(colorspace)

source('analysis_functions.R')
# Settings ------------------------------------------------------------
k <- 10
s <- 3
set.seed(123)
B <- rstiefel::rustiefel(10, 3) %*% diag(c(2, -1, 0.5))
tau <- runif(k, min = -1, max = 1)
gamma <- 5 * svd(B)$v[, 1] + 1.5 * svd(B)$v[, 2]
sigma2_t <- 0.75
sigma2_y <- 0.2

#### Data Generation ####
set.seed(234)
n = 10000
u = MASS::mvrnorm(n = n, mu = rep(0, s), Sigma = diag(s))
tr = u%*%t(B) + MASS::mvrnorm(n = n, mu = rep(0, k), Sigma = sigma2_t*diag(k))
y =  tr%*%tau + u%*%gamma + rnorm(n, mean = 0, sd = sqrt(sigma2_y))


#### Analysis-----------------------------------------------------------------
# Extract latent confounder # 
fit_zt <- extract_B(cov(tr))  ## function extract_B used
A_hat <- fit_zt$A_hat

# Regress y on tr #
fit_y_t <- lm(y ~ tr)
summary(fit_y_t)
tau0_hat <- coef(fit_y_t)[-1]
sigma_y_t_hat <- sigma(fit_y_t)

t1_seq <- seq(-2, 3, by = 0.01)
t_test_mat <- cbind(t1_seq, matrix(0.5, nrow = length(t1_seq), ncol = k-1))
bounds <- t(apply(t_test_mat, 1, est_mean_bound,   ## function est_mean_bound used
                tau0 = tau0_hat, sigma = sigma_y_t_hat, A = A_hat, GaussianT = T))
true <- t_test_mat %*% tau
matplot(t1_seq, cbind(true, bounds), type="l", lty = c(2,3,1,1), 
        col = c("red", "blue", "black", "black"),
        ylab = "E(Y|do(t)) - E(Y|do(0))", xlab = "T1",
        main = "Ignorance Region along T1")
legend("topleft", lty = c(2,3,1), col = c("red", "blue", "black"), cex = 0.8,
       legend = c("True", expression(paste(R^2, " = 0")), expression(paste(R^2, "= 1"))))

################### For ||t|| = 1 ################################################
svd_A <- svd(A_hat)
v_A <- svd_A$v
null_A <- rstiefel::NullC(v_A)


dt_max <- v_A[, 1]
dt_min <- null_A[, 1]
theta <- seq(0, pi/2, by = 0.01)
tan_theta <- tan(theta)
delta_t <- dt_min %*% t(tan_theta)
dt_max2min <- (dt_max + delta_t) / matrix(rep(apply(dt_max + delta_t, 2, L2_norm), times = k), nrow = k, byrow = T)

## function for Calculating the Bound of ATE Bias Gieven R^2 ##
cal_ate_bias_bound <- function(dt, r2, sigma, A) {
  quan <- sigma * sqrt(r2) * sqrt(sum((A %*% dt)^2))
  result <- cbind(-quan, quan)
  names(result) <- c(paste0("lb_", r2), paste0("ub_", r2))
  return(result)
}

## Calculating ATE Bias ##
bound_max2min <- t(apply(dt_max2min, 2, cal_ate_bias_bound, 
                         r2 = 1, sigma = sigma_y_t_hat, A = A_hat))
bound_max2min06 <- t(apply(dt_max2min, 2, cal_ate_bias_bound, 
                           r2 = 0.6, sigma = sigma_y_t_hat, A = A_hat))
bound_max2min03 <- t(apply(dt_max2min, 2, cal_ate_bias_bound, 
                           r2 = 0.3, sigma = sigma_y_t_hat, A = A_hat))

#-------------------------plot bound width using ggplot ------------------------------------------------------


bound_df <- tibble(theta = theta,
                   lower = bound_max2min[, "lb_1"],
                   lower06 = bound_max2min06[, "lb_0.6"],
                   lower03 = bound_max2min03[, "lb_0.3"],
                   unconfounded = rep(0, length(theta)),
                   upper03 = bound_max2min03[, "ub_0.3"],
                   upper06 = bound_max2min06[, "ub_0.6"],
                   upper = bound_max2min[, "ub_1"],
                   )
paste0(rep("R2 = ", 7), c(1, 0.6, 0.3, 0, -0.3, -0.6, -1) )
colnames(bound_df)[-1] = paste0(rep("R2=", 7), 
                                c(1, 0.6, 0.3, 0, -0.3, -0.6, -1))


bound_df_narrow <- gather(bound_df, key = R2, value = bound, -theta) 
bound_df_narrow$R2 <- factor(bound_df_narrow$R2, levels = paste0(rep("R2=", 7), c(1, 0.6, 0.3, 0, -0.3, -0.6, -1)),
       labels =  c(1, 0.6, 0.3, 0, -0.3, -0.6, -1))

ate_bias_bound <- bound_df_narrow %>%
  ggplot() + 
  geom_line(aes(x = theta, y = bound, col = R2)) + 
  ##scale_colour_manual(values = c(upper = "black", lower = "black", truth = "red", unconfounded = "blue")) + 
  scale_color_discrete_divergingx(palette = "Zissou 1",
                                  name = expression(R^2~"="),
                                  labels = sapply(100*c(1, 0.6, 0.3, 0, 0.3, 0.6, 1), paste0, '%'))+
                                  # as.character(c(1, 0.6, 0.3, 0, 0.3, 0.6, 1))) + 
  #scale_color_discrete(name = expression(R^2)) + 
  theme_bw(base_size = 16) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  labs(x = expression(theta), y = "Bias", title = bquote("Change of Bias with "~Delta~t)) + 
  scale_x_continuous(breaks = seq(0, pi/2, by = pi/8),
                     labels = c("0", expression(pi/8), expression(pi/4), expression(3*pi/8), expression(pi/2)),
                     limits = c(0, pi/2))
ate_bias_bound
# ggsave("ate_bias_bound.pdf", plot = ate_bias_bound, width = 150, height = 100, units = "mm"s)


#---------------------------- plot distirbution using ggridges ------------------------------------------------

coef_mu_u_t_hat <- fit_zt$coef_mu_u_zt_hat
Sigma_u_t_hat <- fit_zt$Sigma_u_zt_hat
var_u1_t <- Sigma_u_t_hat[1, 1]
useq <- seq(-2, 2, by=0.01)
upop_density_tibble <- tibble(u = useq, 
                              mu = rep(0, length(useq)),
                              t = rep('pop', length(useq)),
                              dens = dnorm(useq, mean=0, sd=1))

mu_u1_t_max <- (coef_mu_u_t_hat %*% dt_max)[1, 1]
density_tibble_max <- as_tibble(rbind(cbind(useq, rep(mu_u1_t_max, length(useq))),
                                      cbind(useq, rep(-mu_u1_t_max, length(useq)))))
colnames(density_tibble_max) <- c("u", "mu")
density_tibble_max$t <- rep(c("A", "B"), each = length(useq))
density_tibble_max <- density_tibble_max %>%
  mutate(dens = dnorm(u, mean=mu, sd=sqrt(var_u1_t))) %>%
  bind_rows(upop_density_tibble)
plot_max <- density_tibble_max %>% 
  ggplot() + 
  ggridges::geom_ridgeline(aes(x=u, height=dens, y=t, fill = t), scale=2, alpha=0.75) + 
  theme_bw(base_size=13) + 
  theme(plot.title = element_text(hjust = 0.5),
        legend.position=c(0.08, 0.975),
        legend.text=element_text(size=11),
        legend.title = element_text(size = 9),
        legend.key.size = unit(0.25, "cm"),
        legend.text.align = 0,
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()) + 
  scale_fill_discrete(labels=c(expression(U[1]~'|'~t[2]),
                               expression(U[1]~'|'~t[1]),
                               expression(U[1])),
                      name = '') +
  guides(fill = guide_legend(reverse = TRUE)) +
  labs(title = expression(Delta~"t" ~ "=" ~ u[1]^B), x=expression(U[1]), 
       y = expression('Disttribution of'~U[1]))
plot_max


t1_min <- c(-1,-1,1,-1,1,1,1,-1,-1,-1)
t2_min <- t1_min - dt_min
mu_u1_t1_min <- (coef_mu_u_t_hat %*% t1_min)[1, 1]
mu_u1_t2_min <- (coef_mu_u_t_hat %*% t2_min)[1, 1]
density_tibble_min <- as_tibble(rbind(cbind(useq, rep(mu_u1_t1_min, length(useq))),
                                      cbind(useq, rep(mu_u1_t2_min, length(useq)))))
colnames(density_tibble_min) <- c("u", "mu")
density_tibble_min$t <- rep(c("A", "B"), each = length(useq))
density_tibble_min <- density_tibble_min %>%
  mutate(dens = dnorm(u, mean=mu, sd=sqrt(var_u1_t))) %>%
  bind_rows(upop_density_tibble)
plot_min <- density_tibble_min %>% 
  ggplot() + 
  ggridges::geom_ridgeline(aes(x=u, height=dens, y=t, fill = t), scale=2, alpha = 0.75) + 
  theme_bw(base_size=13) + 
  theme(plot.title = element_text(hjust = 0.5),
        legend.position=c(0.08, 0.975),
        legend.text=element_text(size=11),
        legend.title = element_text(size = 9),
        legend.key.size = unit(0.3, "cm"),
        legend.text.align = 0,
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()) + 
  scale_fill_discrete(labels=c(expression(U[1]~'|'~t[2]),
                               expression(U[1]~'|'~t[1]),
                               expression(U[1])),
                      name='')+
  guides(fill = guide_legend(reverse = TRUE)) +
  labs(title = expression(Delta ~ "t" ~ "=" ~ n[0]^B), x=expression(U[1]),
       y = expression('Distribution of'~U[1]))
plot_min

plot_max + plot_min
# ggsave("udist.pdf", plot = plot_max + plot_min, width = 150, height = 100, units = "mm")

