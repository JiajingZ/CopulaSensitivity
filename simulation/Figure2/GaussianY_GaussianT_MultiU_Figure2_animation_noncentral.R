## ggplot layout ##
library(ggpubr)
#### Gaussian Multi-U #####
library(rstiefel)
library(ggplot2)
library(tidyverse)
library(patchwork)
library(ggridges)
library(colorspace)
library(gganimate)

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


#### Analysis -----------------------------------------------------------------
# latent confounder model #
fit_zt <- extract_B(cov(tr))  ## function extract_B used
A_hat <- fit_zt$A_hat
coef_mu_u_t_hat <- fit_zt$coef_mu_u_zt_hat
cov_u_t_hat <- fit_zt$Sigma_u_zt_hat


# outcome model #
fit_y_t <- lm(y ~ tr)
summary(fit_y_t)
tau0_hat <- coef(fit_y_t)[-1]
sigma_y_t_hat <- sigma(fit_y_t)

################### For ||t|| = 1 ################################################
svd_A <- svd(A_hat)
v_A <- svd_A$v
null_A <- rstiefel::NullC(v_A)


dt_max <- v_A[, 1]
dt_min <- null_A[, 1]
theta <- seq(0, pi/2, by = 0.01)
tan_theta <- tan(theta)
delta_t <- dt_min %*% t(tan_theta)
dt_max2min <- (dt_max + delta_t) /
  matrix(rep(apply(dt_max + delta_t, 2, L2_norm), times = k), nrow = k, byrow = T)

## Calculating ATE Bias ##
results_cali <- CopSens::gcalibrate(y = 1, tr = tr,
                                    t1 = dt_max2min,
                                    t2 = matrix(0, nrow = nrow(dt_max2min), ncol = k),
                                    calitype = "worstcase",
                                    mu_y_dt = t(dt_max2min) %*% tau0_hat,
                                    sigma_y_t = sigma_y_t_hat,
                                    mu_u_dt = t(coef_mu_u_t_hat %*% dt_max2min),
                                    cov_u_t = cov_u_t_hat,
                                    R2 = c(0.3, 0.6, 1))

#-------------------------plot bound width using ggplot ------------------------------------------------------
bound_df <- tibble(theta = theta,
                   lower = results_cali$est_df[,'R2_1_lwr'] - results_cali$est_df[,'R2_0'],
                   lower06 = results_cali$est_df[,'R2_0.6_lwr'] - results_cali$est_df[,'R2_0'],
                   lower03 = results_cali$est_df[,'R2_0.3_lwr'] - results_cali$est_df[,'R2_0'],
                   unconfounded = rep(0, length(theta)),
                   upper03 = results_cali$est_df[,'R2_0.3_upr'] - results_cali$est_df[,'R2_0'],
                   upper06 = results_cali$est_df[,'R2_0.6_upr'] - results_cali$est_df[,'R2_0'],
                   upper = results_cali$est_df[,'R2_1_upr'] - results_cali$est_df[,'R2_0'])
colnames(bound_df)[-1] = paste0("R2 = ", c(1, 0.6, 0.3, 0, -0.3, -0.6, -1))


bound_df_narrow <- gather(bound_df, key = R2, value = bound, -theta)
bound_df_narrow$R2 <- factor(bound_df_narrow$R2, levels = paste0("R2 = ", c(1, 0.6, 0.3, 0, -0.3, -0.6, -1)),
                             labels =  c(1, 0.6, 0.3, 0, -0.3, -0.6, -1))

ate_bias_bound_animation <- bound_df_narrow %>%
  ggplot() +
  geom_line(aes(x = theta, y = bound, col = R2)) +
  scale_color_discrete_divergingx(palette = "Zissou 1",
                                  name = expression(R^2~"="),
                                  labels = sapply(100*c(1, 0.6, 0.3, 0, 0.3, 0.6, 1), paste0, '%'))+
  labs(x = expression(theta), y = "Bias", title = expression("Change of Bias with "~Delta~t)) +
  theme_bw(base_size = 16) +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_x_continuous(breaks = seq(0, pi/2, by = pi/8),
                     labels = c("0", expression(pi/8), expression(pi/4), expression(3*pi/8), expression(pi/2)),
                     limits = c(0, pi/2)) +
  #---------------------------- add animation ------------------------------#
  geom_vline(aes(xintercept = x_anim), linetype = "dashed",
             data = tibble(x_anim = seq(0, pi/2, by = 0.05))) +
  gganimate::transition_time(time = x_anim)
ate_bias_bound_gif <- gganimate::animate(ate_bias_bound_animation,
                                         rewind=TRUE, nframes=100, fps=20)
ate_bias_bound_gif
gganimate::anim_save("ate_bias_bound.gif", ate_bias_bound_gif,
                     path = "simulation/Figure2")



# ggsave("ate_bias_bound.pdf", plot = ate_bias_bound, width = 150, height = 100, units = "mm"s)


#---------------------------- plot distirbution using ggridges ------------------------------------------------

var_u1_t <- cov_u_t_hat[1, 1]
useq <- seq(-2, 2, by=0.01)
# super-population density #
upop_density_tibble <- as_tibble(expand.grid(useq, rep(0, ncol(dt_max2min)))) %>%
  'colnames<-'(c("u", "mu")) %>%
  mutate(dens = dnorm(u, mean = mu, sd = 1)) %>%
  add_column(t = rep('pop', length(useq)*ncol(dt_max2min)))
# sub-population density #
# t1_max2min <- matrix(rep(c(-1,-1,1,-1,1,1,1,-1,-1,-1), ncol(dt_max2min)), nrow=k)
# t2_max2min <- t1_max2min - dt_max2min
# mu_u1_t1_max2min <- c((coef_mu_u_t_hat %*% t1_max2min)[1, ])
# mu_u1_t2_max2min <- c((coef_mu_u_t_hat %*% t2_max2min)[1, ])


# mu_u1_t_base <- (coef_mu_u_t_hat %*% c(-1,-1,1,-1,1,1,1,-1,-1,-1))[1,]
mu_u1_t_base <- seq(0, (coef_mu_u_t_hat %*% c(-1,-1,1,-1,1,1,1,-1,-1,-1))[1,],
                    length.out = length(mu_u1_dt_half))

mu_u1_dt_half <- c((coef_mu_u_t_hat %*% dt_max2min)[1, ]/2)
density_tibble <- as_tibble(rbind(expand.grid(useq, mu_u1_dt_half + mu_u1_t_base),
                                  expand.grid(useq, -mu_u1_dt_half + mu_u1_t_base))) %>%
  'colnames<-'(c("u", "mu")) %>%
  mutate(dens = dnorm(u, mean = mu, sd = sqrt(var_u1_t))) %>%
  add_column(t = rep(c("A", "B"), each = length(useq)*ncol(dt_max2min))) %>%
  bind_rows(upop_density_tibble) %>%
  add_column(theta = rep(rep(theta, each = length(useq)), 3))

density_animation <-
  density_tibble %>%
  ggplot() +
  ggridges::geom_ridgeline(aes(x=u, height=dens, y = t, fill = t), scale=2, alpha=0.75) +
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
       y = expression('Disttribution of'~U[1])) +
  #---------- add animation ---------------#
  gganimate::transition_time(theta)
density_gif <- animate(density_animation, rewind=TRUE, nframes = 100, fps=20)
density_gif
# gganimate::anim_save("density_gif_noncentral.gif", density_gif, path = "simulation/Figure2")


## Combine the two animation plot plots
ate_bias_bound_mgif <- magick::image_read("simulation/Figure2/ate_bias_bound.gif")
density_mgif <- magick::image_read("simulation/Figure2/density_gif_noncentral.gif")


combined_gif <- magick::image_append(c(ate_bias_bound_mgif[1], density_mgif[1]))
for(i in 2:100){
  combined <- magick::image_append(c(ate_bias_bound_mgif[i], density_mgif[i]))
  combined_gif <- c(combined_gif, combined)
}
combined_gif
# gganimate::anim_save("combined_noncentral.gif", combined_gif, path = "simulation/Figure2")


