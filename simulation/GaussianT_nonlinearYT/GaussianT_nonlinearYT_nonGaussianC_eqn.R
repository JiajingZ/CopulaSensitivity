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
####################################################################


cal_doMean_quadratic <- function(t, N = 1e6, a = 1, b = 0) {
  message(cat("Case: T = ", t))

  ## sample tildet ~ f(t)
  mu_u_tildet <- c(as.matrix(tr[sample(1:n, size = N, replace = TRUE), ]) %*%
                     t(coef_mu_u_t))
  mu_u_t <- c(coef_mu_u_t %*% t)


  # u ~ f(u | t)
  u_sample_cond <- rnorm(N, mean = mu_u_t, sd = sigma_u_t)
  ## ytilde = f(U) + epsilon, UY relatioship
  # yu_func <- function(u) {3*sin(10*abs(u)) + rnorm(N, sd=1)}
  # yu_func <- function(u) {40*abs(u - 0.5) + rnorm(N, sd = 1)}
  yu_func <- function(u) {a*(u - b)^2}

  ytilde_sample_cond <- yu_func(u_sample_cond)
  ytilde_cdf <- ecdf(ytilde_sample_cond)

  ## Should be uniform
  Uy_sample_cond <- ytilde_cdf(ytilde_sample_cond)

                                        ## u ~ f(u)
  u_sample_marg <- rnorm(N, mean = mu_u_tildet, sd = sigma_u_t)
  ytilde_sample_marg <- yu_func(u_sample_marg)
  Uy_sample_marg <- ytilde_cdf(ytilde_sample_marg)
  Uy_sample_marg[Uy_sample_marg > 1-1e-3] <- 1-1e-3
  Uy_sample_marg[Uy_sample_marg < 1e-3] <- 1e-3

  ## Uy to Y
  mu_y_t <- c(g_yt(t(t)) + t %*% t(coef_mu_u_t) %*% gamma)
  y_sample <- qnorm(Uy_sample_marg, mean = mu_y_t, sd = sigma_y_t)

  ## calculate mean(y)
  mean(y_sample)
}

plot(pnorm(u_sample_cond, mean = mu_u_t, sd = sigma_u_t), Uy_sample_cond)
plot(pnorm(u_sample_marg, mean = mu_u_t, sd = sigma_u_t), Uy_sample_marg)


t_interest <- rbind(t_choice, t2)
quad_para <- expand_grid(a = c(1, -1), b = c(-2, 0, 2))
domean_cali_quad_mat <- matrix(NA, nrow = nrow(t_interest), ncol = nrow(quad_para),
                               dimnames = list(NULL, paste0("Quadratic: a = ", quad_para$a, ", b = ", quad_para$b)))
for(i in 1:nrow(quad_para)){
  domean_cali_quad_mat[,i] <- apply(t_interest, 1,
                                    cal_doMean_quadratic,
                                    a = quad_para$a[i],
                                    b = quad_para$b[i])
}

colnames(domean_cali_quad_mat) = paste0("Quadratic: a = ", quad_para$a, ", b = ", quad_para$b)
ate_cali_quad_narrow <-
  t(t(domean_cali_quad_mat)[,1:4] - t(domean_cali_quad_mat)[,5]) %>%
  as_tibble() %>%
  add_column(case = 1:nrow(t_choice)) %>%
  gather(key = "Type", value = "effect", - case)

saveRDS(ate_cali_quad_narrow,
        file = "simulation/GaussianT_nonlinearYT/ate_cali_quad_narrow.rds")

cal_doMean_mix <- function(t, a1, a2, b1, b2, N = 1e6) {
  message(cat("Case: T = ", t))

  ## sample tildet ~ f(t)
  mu_u_tildet <- c(as.matrix(tr[sample(1:n, size = N, replace = TRUE), ]) %*%
                   t(coef_mu_u_t))
  mu_u_t <- c(coef_mu_u_t %*% t)


  ## u ~ f(u | t)
  u_sample_cond <- rnorm(N, mean = mu_u_t, sd = sigma_u_t)
  yu_func <- function(u) {ifelse(runif(n) < 0.5, a1*u+b1, a2*u+b2)}

  ytilde_sample_cond <- yu_func(u_sample_cond)
  ytilde_cdf <- ecdf(ytilde_sample_cond)
  Uy_sample_cond <- ytilde_cdf(ytilde_sample_cond)
  ## u ~ f(u)
  u_sample_marg <- rnorm(N, mean = mu_u_tildet, sd = sigma_u_t)
  ytilde_sample_marg <- yu_func(u_sample_marg)
  Uy_sample_marg <- ytilde_cdf(ytilde_sample_marg)
  Uy_sample_marg[Uy_sample_marg > 1-1e-3] <- 1-1e-3
  Uy_sample_marg[Uy_sample_marg < 1e-3] <- 1e-3

  ## Uy to Y
  mu_y_t <- c(g_yt(t(t)) + t %*% t(coef_mu_u_t) %*% gamma)
  y_sample <- qnorm(Uy_sample_marg, mean = mu_y_t, sd = sigma_y_t)

  ## calculate mean(y)
  mean(y_sample)

}

mix_para <- tibble(a1 = c(1,1,1,-1),
                   a2 = c(-1,-1,1,-1),
                   b1 = c(0,-2,-2,-2),
                   b2 = c(-2,0,2,2))
domean_cali_mix_mat <- matrix(NA, nrow = nrow(t_interest), ncol = nrow(mix_para),
                               dimnames = list(NULL, paste0("Mix: a1=", mix_para$a1, 
                                                            ",b1=", mix_para$b1,
                                                            ",a2=", mix_para$a2,
                                                            ",b2=", mix_para$b2)))
for(i in 1:nrow(mix_para)){
  domean_cali_mix_mat[,i] <- apply(t_interest, 1,
                                   cal_doMean_mix,
                                   a1 = mix_para$a1[i],
                                   a2 = mix_para$a2[i],
                                   b1 = mix_para$b1[i],
                                   b2 = mix_para$b2[i])
}

ate_cali_mix_narrow <-
  t(t(domean_cali_mix_mat)[,1:4] - t(domean_cali_mix_mat)[,5]) %>%
  as_tibble() %>%
  add_column(case = 1:nrow(t_choice)) %>%
  gather(key = "Type", value = "effect", - case)

saveRDS(ate_cali_mix_narrow,
        file = "simulation/GaussianT_nonlinearYT/ate_cali_mix_narrow.rds")



# Sensitivity Analysis NonGaussian Copula---------------------------------------------------------------
# Algorithm
# true copula #
norm_cop <- copula::normalCopula(gamma*sigma_u_t/sigma_y_t)
# Clayton #
readRDS("simulation/GaussianT_nonlinearYT/ate_cali_clayton_tibble.rds") %>%
  gather(key = "Type", value = "effect", - case) -> ate_cali_clayton_narrow
ate_cali_clayton_narrow$Type <- factor(ate_cali_clayton_narrow$Type, 
                                       levels = c("Clayton -0.5","Clayton 2","Clayton 10"),
                                       labels = c("Clayton (theta==-0.5)","Clayton theta==2","Clayton (theta=10)"))
# Quadratic #
readRDS("simulation/GaussianT_nonlinearYT/ate_cali_quad_narrow.rds") -> ate_cali_quad_narrow
# Mix Distributions #
readRDS("simulation/GaussianT_nonlinearYT/ate_cali_mix_narrow.rds") -> ate_cali_mix_narrow


legend_labs <- c(expression((u-2)^2), expression(u^2), expression((u+2)^2),
                 expression(-(u-2)^2), expression(-u^2), expression(-(u+2)^2),
                 expression(Clayton (theta==-0.5)), expression(Clayton (theta==2)), expression(Clayton (theta==10)),
                 expression(ifelse(runif(n) < 0.5, u, -u-2)), 
                 expression(ifelse(runif(n) < 0.5, u-2, -u)),
                 expression(ifelse(runif(n) < 0.5, u-2, u+2)),
                 expression(ifelse(runif(n) < 0.5, -u-2, -u+2)))

## plot #
library(ggrepel)
library(magrittr)
codebook1 <- unique(ate_cali_quad_narrow$Type)
names(codebook1) <- toupper(letters)[1:6]

codebook2 <- as.character(unique(ate_cali_clayton_narrow$Type))
names(codebook2) <- toupper(letters)[7:9]

codebook3 <- as.character(unique(ate_cali_mix_narrow$Type))
names(codebook3) <- toupper(letters)[10:13]


ate_cali_quad_narrow %<>% mutate(Type2 = fct_recode(factor(Type), !!!codebook1)) %>% mutate(Class="Quadratic")
ate_cali_clayton_narrow %<>% mutate(Type2 = fct_recode(factor(Type), !!!codebook2)) %>% mutate(Class = "Clayton")
ate_cali_mix_narrow %<>% mutate(Type2 = fct_recode(factor(Type), !!!codebook3)) %>% mutate(Class = "Mix")

ate_cali_narrow <- bind_rows(ate_cali_quad_narrow, ate_cali_clayton_narrow, ate_cali_mix_narrow)
ate_cali_narrow %<>% mutate(Type3 = paste(Type2, Type, sep=": "))

ate_cali_narrow %<>%
  mutate(Type2=fct_relevel(ate_cali_narrow$Type2, toupper(letters[1:11]))) %>%
  mutate(Type=fct_reorder(factor(ate_cali_narrow$Type), as.numeric(Type2)))

names(codebook1) <- c(expression((u-2)^2), expression(u^2), expression((u+2)^2),
                      expression(-(u-2)^2), expression(-u^2), expression(-(u+2)^2))

ate_cali_narrow <- ate_cali_narrow %>% mutate(Type4 = fct_recode(Type, !!!codebook1))

cols <- colorspace::qualitative_hcl(3)
p <- tibble(lwr_g = ate_cali_mat[,'R2_1_lwr'],
       upr_g = ate_cali_mat[,'R2_1_upr'],
       case=1:nrow(t_choice)) %>%
  ggplot() +
  geom_errorbar(aes(ymin = lwr_g, ymax = upr_g, x = case), width = 0.2) + 
  geom_label_repel(data=ate_cali_narrow, aes(x=case, y=effect, label=Type2, color=Type4),
                   force=10, nudge_x=0.3, parse=TRUE) +
  labs(title = bquote(PATE[paste(t[1], ",", t[2])]~"for alternative Copulas"),
       y = expression('Causal Effect'), x = 'i') +
  xlim(0.9, 4.4) +
  scale_color_manual(values=c(rep(cols[1], 6), rep(cols[2], 3), rep(cols[3], 2)), labels=legend_labs) +
  theme_bw(base_size = 15) +
  theme(plot.title = element_text(hjust = 0.5), legend.text.align = 0)

p <- p + guides(
             col = guide_legend(
               title = "Copula",
               override.aes = list(label = toupper(letters[1:11]))))

p
ggsave(p, file="simulation/GaussianT_nonlinearYT/copula_sim.pdf")










