library(tidyverse)
library(glmnet)
library(colorspace)
library(patchwork)
setwd("simulation/Sparse_Effects_Setting")

# import data #
y <- read.csv('y.csv')$X0
tr <- read.csv('tr.csv') %>% as.matrix()
tau <- read.csv('tau.csv')$X0
# nontrivial_index <- read.csv('nontrivial_effect_index.csv')$X0 + 1
nontrivial_index <- which(abs(tau) > 0.1)

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

# histogram of difference in conditional confounder mean -------------------------------------

tibble(u_t_diff = c(u_t_diff[-nontrivial_index, 1], u_t_diff[nontrivial_index, 1]),
       group = c(rep("null", 455), rep("nonnull", 45))) %>%
  ggplot(aes(x = u_t_diff, col = group)) +
  geom_density() +
  labs(x = bquote(mu[paste(u[1], "|", t[1])]~-~mu[paste(u[1], "|", t[2])])) +
  theme_bw()

tibble(u_t_diff = c(u_t_diff[-nontrivial_index, 2], u_t_diff[nontrivial_index, 2]),
       group = c(rep("null", 455), rep("nonnull", 45))) %>%
  ggplot(aes(x = u_t_diff, col = group)) +
  geom_density() +
  labs(x = bquote(mu[paste(u[2], "|", t[1])]~-~mu[paste(u[2], "|", t[2])])) +
  theme_bw()

tibble(u_t_diff = c(u_t_diff[-nontrivial_index, 3], u_t_diff[nontrivial_index, 3]),
       group = c(rep("null", 455), rep("nonnull", 45))) %>%
  ggplot(aes(x = u_t_diff, col = group)) +
  geom_density() +
  labs(x = bquote(mu[paste(u[3], "|", t[1])]~-~mu[paste(u[3], "|", t[2])])) +
  theme_bw()

u_t_diff_projnorm <- apply(u_t_diff %*% cov_halfinv, 1 , function(x) sqrt(sum(x^2)))
tibble(diff = c(u_t_diff_projnorm[-nontrivial_index], u_t_diff_projnorm[nontrivial_index]),
       group = c(rep("null", 455), rep("nonnull", 45))) %>%
  ggplot(aes(x = diff, col = group)) +
  geom_density() +
  labs(x = "L2 norm of projected difference in conditional confounder mean") +
  theme_bw()

# Optimization ----------------------------------------------------
######## Calibration with Gamma ##########
objective <- function(gamma_opt){
  tau_cali <- tau_t - u_t_diff %*% gamma_opt
  norm(tau_cali, type = "1") ## L1 norm
  # mad(tau_cali) ## MAD
}

# optimize #
obj_min <- objective(rep(0, latent_dim))
for (i in 1:1000) {
  gamma0 <- rnorm(latent_dim)*10
  solution <- optim(par = gamma0, fn = objective)
  if (solution$value < obj_min) {
    print("get smaller value!")
    obj_min <- solution$value
    gamma_int_min <- gamma0
    gamma_opted1 <- solution$par
  }
}
objective(rep(0, latent_dim)) # 2525.657
obj_min # 361.9072
gamma_int_min # -3.9380409 -2.7525589  0.1184939
gamma_opted1
# gamma_opted1 <- c(42.783441, 9.247085, -146.155949)
## R^2 ##
(R2 <- c(t(gamma_opted1) %*% cov_u_t %*% gamma_opted1 / sigma_y_t^2) %>% round(digits = 2))



# Calibrating ----------------------------------------------------

# Based on gamma #
cal_tau_calibrated_gamma <- function(gamma = gamma) {
  return(tau_t - u_t_diff %*% gamma)
}

tau_cali <- cal_tau_calibrated_gamma(gamma = gamma_opted1)


## Checking the norm of estimates
norm(as.matrix(tau), type = "1")
norm(as.matrix(tau_t), type = "1")
norm(tau_cali, type = "1")

## Checking whether TEs shrink or not after calibration
sum(tau_t^2/length(tau_t)) # 39.68797
sum(tau_cali^2/length(tau_cali)) # 0.8083098


# RMSE #
# overall #
sqrt(mean((tau_t - tau)^2))
sqrt(mean((tau_cali - tau)^2))
# null #
sqrt(mean(((tau_t - tau)[-nontrivial_index])^2)) # 6.25
sqrt(mean(((tau_cali - tau)[-nontrivial_index])^2)) # 0.82
# non-null #
sqrt(mean(((tau_t - tau)[nontrivial_index])^2)) # 6.58
sqrt(mean(((tau_cali - tau)[nontrivial_index])^2)) # 0.83


## checking shinking direction
diff <- tau  - tau_t
diff_cali <- tau_cali - tau_t
sum(diff*diff_cali > 0) # 472
sum(diff*diff_cali > 0) / ncol(tr)

# how many are expected to decrease
sum(diff < 0)
# among decrese, how many calibrated are aligned
sum((diff < 0) & (diff*diff_cali > 0))
sum((diff < 0) & (diff*diff_cali > 0)) / sum(diff < 0)
# how many are expected to increase
sum(diff > 0)
# among increse, how many calibrated are aligned
sum((diff > 0) & (diff*diff_cali > 0))
sum((diff > 0) & (diff*diff_cali > 0)) / sum(diff > 0)

mean(abs(tau_t)) - mean(abs(tau_cali))
mean(abs(tau_t[nontrivial_index])) - mean(abs(tau_cali[nontrivial_index]))
mean(abs(tau_t[-nontrivial_index])) - mean(abs(tau_cali[-nontrivial_index]))

cbind(tau, tau_t, tau_cali)[nontrivial_index,]
# among nontrivial ones, how many are aligned #
sum((diff*diff_cali > 0)[nontrivial_index])


(nontrivial_count <- sum(abs(tau) > 0.1))
tau_sort <- abs(tau) %>% sort(decreasing = T, index.return = T)
nontrivial_index <- tau_sort$ix[1:nontrivial_count]

############################## ROC #########################################################################
## By P-value ## ------------------------------------------------------------------------------
positive_index <- order(abs(tau), decreasing = T)[1:nontrivial_count]
negative_index <- order(abs(tau), decreasing = T)[(nontrivial_count+1):500]
lmfit_y_t_coef <- summary(lmfit_y_t)$coefficients
pvalue_naive <- lmfit_y_t_coef[-1,4]
sd_tau <- lmfit_y_t_coef[-1,2]
tvalue_cali <- tau_cali /  sd_tau
pvalue_cali <- 2*pt(- abs(tvalue_cali), df = lmfit_y_t$df.residual)

tpr_naive = fpr_naive  = ppv_naive = tpr_cali = fpr_cali = ppv_cali = NULL
for (k in seq(1, 500, by = 1)) {
  #### Naive ####
  positive_index_naive <- order(pvalue_naive)[1:k]
  # TPR  = predicted_positive/true_positive #
  tpr_naive <- c(tpr_naive, mean(positive_index %in% positive_index_naive))
  # FPR = predicted_positive/ true_negative #
  fpr_naive <- c(fpr_naive, mean(negative_index %in% positive_index_naive))
  # PPV = true_positive /  predicted_positive #
  ppv_naive <- c(ppv_naive, mean(positive_index_naive %in% positive_index))
  #### Calibrated ####
  positive_index_cali <- order(pvalue_cali)[1:k]
  # TPR  = predicted_positive/true_positive #
  tpr_cali <- c(tpr_cali, mean(positive_index %in% positive_index_cali))
  # FPR = predicted_positive/ true_negative
  fpr_cali <- c(fpr_cali, mean(negative_index %in% positive_index_cali))
  # PPV = true_positive /  predicted_positive #
  ppv_cali <- c(ppv_cali, mean(positive_index_cali %in% positive_index))
}

plot_roc <- tibble(TPR = c(tpr_naive, tpr_cali),
                   FPR = c(fpr_naive, fpr_cali),
                   Type = rep(c('uncalibrated', 'calibrated'), each = 500)) %>%
  ggplot() +
  geom_line(aes(x = FPR, y = TPR, colour = Type)) +
  scale_colour_manual(name = "",
                      values = divergingx_hcl(5,palette = "Zissou 1")[c(5,1)]) +
  labs(x = "FPR", y = "TPR", title = paste0("ROC (Latent Dim. = ", latent_dim, ')')) +
  theme_bw(base_size = 12) +
  theme(plot.title = element_text(hjust = 0.5),
        legend.title = element_text(size=12),
        legend.text = element_text(size=11))

(auc_naive <- DescTools::AUC(x = fpr_naive, y = tpr_naive)) # 0.5451526
(auc_cali <- DescTools::AUC(x = fpr_cali, y = tpr_cali)) # 0.7302564

plot_roc
# write.csv(fpr_cali, row.names = F, file = 'LatentDim3/fpr_cali_dim3_L1.csv')
# write.csv(tpr_cali, row.names = F, file = 'LatentDim3/tpr_cali_dim3_L1.csv')

#### Scatter Plot ----------------------------------------------------------------

nontrivial_index <- which(abs(tau) > 0.1)
trivial_seleted_index <- sample(1:455, size = 55, replace = FALSE)
scatter_true <- tibble(index = 1:100,
                       effect = c((tau[-nontrivial_index])[trivial_seleted_index],
                                  tau[nontrivial_index]),
                       group = c(rep('null', 55), rep('nonnull', 45))) %>%
  ggplot(aes(x = index, y = effect, colour = group)) +
  geom_point(size = 3) +
  scale_colour_manual(values = c('#C2185B', 'black')) +
  labs(y = expression(tau), x = 'i', title = "True effects:") +
  ylim(-18, 18) +
  theme_bw(base_size = 19) +
  theme(plot.title = element_text(size=19))
# theme(plot.title = element_text(hjust = 0.5))

(rmse_naive <- sqrt(mean((tau_t - tau)^2)) %>% round(digits = 2))

scatter_naive <- tibble(index = 1:100,
                        effect = c((tau_t[-nontrivial_index])[trivial_seleted_index],
                                   tau_t[nontrivial_index]),
                        group = c(rep('null', 55), rep('nonnull', 45))) %>%
  ggplot(aes(x = index, y = effect, colour = group)) +
  geom_point(size = 3) +
  scale_colour_manual(values = c('#C2185B', 'black')) +
  labs(y = expression(tau), x = 'i',
       title = bquote("Estimated effects("~R^2~"= 0 %, RMSE = "~.(rmse_naive)~")")) +
  ylim(-18, 18) +
  theme_bw(base_size = 19) +
  theme(plot.title = element_text(size=19))
# theme(plot.title = element_text(hjust = 0.5))


(rmse_cali <- sqrt(mean((tau_cali - tau)^2)) %>% round(digits = 2))
scatter_cali <- tibble(index = 1:100,
                       effect = c((tau_cali[-nontrivial_index])[trivial_seleted_index],
                                  tau_cali[nontrivial_index]),
                       group = c(rep('null', 55), rep('nonnull', 45))) %>%
  ggplot(aes(x = index, y = effect, colour = group)) +
  scale_colour_manual(values = c('#C2185B', 'black')) +
  geom_point(size = 3) +
  labs(y = expression(tau), x = 'i',
       title = bquote("Estimated effects("~R^2~"= "~.(100*R2)~"%, RMSE = "~.(rmse_cali)~")")) +
  ylim(-18, 18) +
  theme_bw(base_size = 19)  +
  theme(plot.title = element_text(size=19))
# theme(plot.title = element_text(hjust = 0.5))


scatter_coef <- scatter_true/scatter_naive/scatter_cali
# + plot_annotation(title = paste0("Scatter Plots for Coefficients (Latent Dim. = ", latent_dim, ')'))
scatter_coef


ggsave("LatentDim3/scatter_coef_L1.pdf", plot = scatter_coef, width = 300, height = 190, units = "mm")


scatter_est <- tibble(index = rep(1:100, times = 2),
                      effect = c((tau_t[-nontrivial_index])[trivial_seleted_index],
                                 tau_t[nontrivial_index],
                                 (tau_cali[-nontrivial_index])[trivial_seleted_index],
                                 tau_cali[nontrivial_index]),
                      group = rep(c("uncali", "cali"), each = 100)) %>%
  ggplot(aes(x = index, y = effect, colour = group)) +
  geom_point(size = 3, alpha = 0.75) +
  scale_colour_manual(name = "",
                      values = c("#F5191C", "#3B99B1"),
                      labels = c(expression(R^2~" = 74%, RMSE = 0.82"),
                                 expression(R^2~" =  0 %, RMSE = 6.27"))) +
  labs(y = expression(tau), x = 'i') +
  geom_vline(aes(xintercept=55.5), color = "grey") +
  theme_bw(base_size = 19)
scatter_est
ggsave("LatentDim3/scatter_est_L1.pdf", plot = scatter_est, width = 350, height = 140, units = "mm")

############################  Varying R^2 ########################################################
calibration_opt <- function(objective, init=rep(0, latent_dim), iters=1000) {
  obj_min <- objective(init)
  for (i in 1:iters) {
    gamma0 <- init
    solution <- optim(par = gamma0, fn = objective, method = "BFGS")
    if (solution$value < obj_min) {
      print("get smaller value!")
      obj_min <- solution$value
      gamma_opted1 <- solution$par
    }
    gamma_opted1 <- solution$par
    gamma0 <- rnorm(latent_dim)*10
  }
  gamma_opted1
}

cal_tau_gamma <- function(objective, init=rep(0, latent_dim), iters=1000) {
  gamma_opted1 <- calibration_opt(objective, init=init, iters=iters)
  return(list(tau_cali=(tau_t - u_t_diff %*% gamma_opted1), gamma_opt=gamma_opted1))
}

# L1 Norm #
penalty_weights <-  c(seq(0, 8000, by = 500), 20000)
tau_cali_mat <- matrix(0, nrow=length(penalty_weights), ncol=length(tau_t))
R2_vec <- numeric(length(penalty_weights))
count <- 1
for(weight in penalty_weights) {

  results <- cal_tau_gamma(function (g) {
    tau_cali <- tau_t - u_t_diff %*% g
    norm(tau_cali, type = "1") + weight * t(g) %*% cov_u_t %*% g / sigma_y_t^2 ## penalizing R^2
  }, init=10*rnorm(3), iters=100)

  tau_cali_mat[count, ]  <- results$tau_cali
  R2_vec[count]  <- t(results$gamma_opt) %*% cov_u_t %*% results$gamma_opt / sigma_y_t^2
  count <- count + 1
}

plot(R2_vec, apply(tau_cali_mat, 1, function(x) mean(abs(x))), type = "l",
     ylab = expression("m("~tau~")"), xlab = expression(R[paste(Y,'~',U,' | ',T)]^2))
# ylab = expression(paste("||",tau,"||")[1]), xlab = expression(R[paste(Y,'~',U,' | ',T)]^2))
plot(R2_vec, apply(tau_cali_mat[,- nontrivial_index], 1, function(x) mean(abs(x))),
     ylab = expression("m("~tau~")"), xlab = expression(R[paste(Y,'~',U,' | ',T)]^2),
     # ylab = expression(paste("||",tau,"||")[1]), xlab = expression(R[paste(Y,'~',U,' | ',T)]^2),
     ylim = c(0, 5.2))
points(R2_vec, apply(tau_cali_mat[,nontrivial_index], 1, function(x) mean(abs(x))), col = "blue")

plot_r2_norm <- tibble(norm = c(apply(tau_cali_mat[,- nontrivial_index], 1, function(x) mean(abs(x))),
                                apply(tau_cali_mat[,nontrivial_index], 1, function(x) mean(abs(x)))),
                       R2 = rep(R2_vec, times = 2),
                       group = rep(c("null", "nonnull"), each = length(R2_vec))) %>%
  ggplot() +
  geom_smooth(aes(x = R2, y = norm, col = group), se = FALSE) +
  scale_colour_manual(name = "", values = c('#F8766D', '#00BFC4'),
                      guide = guide_legend(reverse=TRUE)) +
  labs(y = expression("m("~tau~")"), x = expression(R[paste(Y,'~',U,' | ',T)]^2)) +
  theme_bw(base_size = 12)
plot_r2_norm
ggsave("LatentDim3/plot_r2_norm.pdf", plot = plot_r2_norm, width = 120, height = 90, units = "mm")

## Scatter plot #
scatter_varyR2 <- tibble(tau = as.numeric(tau_cali_mat[, c(trivial_seleted_index, nontrivial_index)]),
                         R2 = rep(R2_vec, 100),
                         index = rep(1:100, each = nrow(tau_cali_mat))) %>%
  ggplot() + geom_point(aes(x = index, y = tau, col = R2), size = 0.1) +
  scale_color_continuous_sequential(name = expression(R[paste(Y,'~',U,'|',T)]^2~":"),
                                    palette="Sunset", rev=TRUE) +
  # geom_vline(aes(xintercept=55.5), color = "grey", linetype = "dashed") +
  geom_vline(aes(xintercept=55.5), color = "grey") +
  labs(y = expression(tau)) +
  theme_bw()
ggsave("LatentDim3/scatter_varyR2_L1.pdf", plot = scatter_varyR2, width = 350, height = 140, units = "mm")

















