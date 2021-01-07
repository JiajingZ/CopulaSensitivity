library(tidyverse)
library(glmnet)
library(colorspace)
library(patchwork)
setwd("simulation/Sparse_Effects_Setting")

# import data #
y <- read.csv('y.csv')$X0
tr <- read.csv('tr.csv') %>% as.matrix()
tau <- read.csv('tau.csv')$X0
nontrivial_index <- read.csv('nontrivial_effect_index.csv')$X0 + 1


# fit observed outcome model --------------------------------------------------------------------------------------------------------
## y ~ tr ##
lmfit_y_t <- lm(y ~ tr)

tau_t = coef(lmfit_y_t)[-1]
yhat <- predict(lmfit_y_t)
sigma_y_t <- sigma(lmfit_y_t)
sigma_y_t^2

# Sensitivity Analysis --------------------------------------------------------------------------------------------------
latent_dim <- 3
mu_u_t <- read.csv('LatentDim3/mu_u_t_ise.csv') %>% as.matrix()
cov_u_t <- read.csv('LatentDim3/cov_u_t_ise.csv') %>% as.matrix()
eigen_cov <- eigen(cov_u_t)
cov_halfinv <- eigen_cov$vectors %*% diag(eigen_cov$values^{-1/2}) %*% t(eigen_cov$vectors)
u_t_diff_org_all <- read.csv('LatentDim3/u_t_diff_org_all.csv') %>% as.matrix()
u_t_diff <- u_t_diff_org_all
u_pca <- prcomp(mu_u_t)

# How much confounding components are captured by VAE (mu_u_t_ise)
gamma <- rep(100, 3)
u <- read.csv('u.csv') %>% as.matrix()
y_confound <- u %*% gamma
## gamma'u ~ uhat ##
lmfit_confound <- lm(y_confound ~ mu_u_t)
summary(lmfit_confound) # Multiple R-squared:  0.7989,	Adjusted R-squared:  0.7989


# Optimization ----------------------------------------------------
######## Calibration with Gamma ##########

calibration_opt <- function(objective, init=rep(0, latent_dim), iters=1000) {
  obj_min <- objective(init)
  for (i in 1:iters) {
    gamma0 <- init
    solution <- optim(par = gamma0, fn = objective)
    if (solution$value < obj_min) {
      print("get smaller value!")
      obj_min <- solution$value
      gamma_int_min <- gamma0
      gamma_opted1 <- solution$par
    }
    gamma_opted1 <- solution$par
    gamma0 <- rnorm(latent_dim)*10
  }
  gamma_opted1
}

# Based on gamma #
cal_tau <- function(objective, init=rep(0, latent_dim)) {
  gamma_opted1 <- calibration_opt(objective, init=init)
  return(tau_t - u_t_diff %*% gamma_opted1)
}

# Based on gamma #
cal_tau_gamma <- function(objective, init=rep(0, latent_dim), iters=1000) {
  gamma_opted1 <- calibration_opt(objective, init=init, iters=iters)
  return(list(tau_cali=(tau_t - u_t_diff %*% gamma_opted1), gamma_opt=gamma_opted1))
}

## L infinity norm
tau_cali_max <- cal_tau(function (g) {
  tau_cali <- tau_t - u_t_diff %*% g
  max(abs(tau_cali))
})

## L1 norm
tau_cali_L1 <- cal_tau(function (g) {
  tau_cali <- tau_t - u_t_diff %*% g
  norm(tau_cali, type="1")
})

## L2 norm
tau_cali_L2 <- cal_tau(function (g) {
  tau_cali <- tau_t - u_t_diff %*% g
  norm(tau_cali, type="2")
})


## I had to guess at these.  It will depend on the norm you are using.
## This is L2 which is really similar to L1 in this case
gamma_weights <- c(seq(0.1, 0.9, by=0.005), 10)
tau_mat <- matrix(0, nrow=length(gamma_weights), ncol=length(tau_t))
r2_vec <- numeric(length(gamma_weights))
count <- 1
for(weight in gamma_weights) {

  results <- cal_tau_gamma(function (g) {
    tau_cali <- tau_t - u_t_diff %*% g
    norm(tau_cali, type="2") + weight*norm(g, type="2")
  }, init=10*rnorm(3), iters=10)

  tau_mat[count, ]  <- results$tau_cali
  r2_vec[count]  <- t(results$gamma_opt) %*% cov_u_t %*% results$gamma_opt / sigma_y_t^2
  count <- count + 1
}

plot(r2_vec, apply(tau_mat, 1, function(x) norm(x, type="2")))

tau_ordered <- cbind(tau_mat[, setdiff(1:500, nontrivial_index)], tau_mat[, unique(nontrivial_index)])

tibble(taus=as.numeric(tau_mat), R2 = rep(r2_vec, 500), index=rep(1:500, each=nrow(tau_mat))) %>%
  filter(index >= 400) %>%
  ggplot() + geom_point(aes(x=index, y=taus, col=R2), size=0.2) + theme_bw() +
  scale_color_continuous_sequential(palette="ag_GrnYl", rev=TRUE) + geom_vline(aes(xintercept=455))

plot(r2_vec, apply(tau_mat, 1, function(x) sqrt(mean((x - tau)^2))))

tau_ordered <- cbind(tau_mat[, setdiff(1:500, nontrivial_index)], tau_mat[, unique(nontrivial_index)])

median(abs(tau_ordered[1, 455:500]))
norm(tau_ordered[1, 455:500, drop=FALSE], type='1')
sum(abs(tau_ordered[1, 455:500]))/45


median(abs(tau_ordered[1, 1:455]))
sum(abs(tau_ordered[1, 1:455]))/455
t.test(abs(tau_ordered[1, 455:500]), abs(tau_ordered[1, 1:455]))
