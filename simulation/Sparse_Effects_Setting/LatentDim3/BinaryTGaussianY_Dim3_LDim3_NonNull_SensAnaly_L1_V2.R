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

nontrivial_index <- which(abs(tau) > 0.1)

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
(R2 <- c(t(gamma_opted1) %*% cov_u_t %*% gamma_opted1 / sigma_y_t^2) %>%
    round(digits = 2))



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
sqrt(mean(((tau_t - tau)[nontrivial_index])^2)) # 6.59
sqrt(mean(((tau_cali - tau)[nontrivial_index])^2)) # 0.83


## checking shinking direction
diff <- tau  - tau_t
diff_cali <- tau_cali - tau_t
sum(diff*diff_cali > 0) # 475
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


############################## ROC #########################################################################
## By P-value ## ------------------------------------------------------------------------------
positive_index <- order(abs(tau), decreasing = T)[1:length(nontrivial_index)]
negative_index <- order(abs(tau), decreasing = T)[(length(nontrivial_index)+1):500]
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
scatter_coef
ggsave("LatentDim3/scatter_coef_L1.pdf", plot = scatter_coef, width = 300, height = 190, units = "mm")

############################  Varying R^2 ########################################################
get_opt_gamma <- function(objective, init=rep(0, latent_dim), iters=1000) {
  obj_min <- objective(init)
  for (i in 1:iters) {
    gamma0 <- init
    solution <- optim(par = gamma0, fn = objective, method = "BFGS")
    if (solution$value < obj_min) {
      print("Smaller objective obtained!")
      obj_min <- solution$value
      gamma_opted1 <- solution$par
    }
    gamma_opted1 <- solution$par
    gamma0 <- rnorm(latent_dim)*10
  }
  gamma_opted1
}

cal_tau_gamma <- function(objective, init=rep(0, latent_dim), iters=1000) {
  gamma_opted1 <- get_opt_gamma(objective, init=init, iters=iters)
  return(list(tau_cali=(tau_t - u_t_diff %*% gamma_opted1), gamma_opt=gamma_opted1))
}

# L1 Norm #
penalty_weights <-  c(seq(0, 8000, by = 200), 20000)
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
norm_taucali_nonnull <- apply(tau_cali_mat[,nontrivial_index], 1, function(x) mean(abs(x)))
norm_taucali_null <- apply(tau_cali_mat[, -nontrivial_index], 1, function(x) mean(abs(x)))
norm_ratio <- norm_taucali_nonnull / norm_taucali_null


plot_r2_norm <- tibble(norm = c(norm_taucali_null,
                                norm_taucali_nonnull),
                       R2 = rep(R2_vec, times = 2),
                       group = rep(c("null", "nonnull"), each = length(R2_vec))) %>%
  ggplot() + geom_line(aes(x = R2, y = norm, col = group), size = 1) +
  scale_colour_manual(name = "", values = c('#2C86CA', '#6D6D6D'),
                      guide = guide_legend(reverse=TRUE)) +
  labs(y = expression("E("~"|"~tau~"|"~")"), x = expression(R[paste(Y,'~',U,'|',T)]^2)) +
  theme_bw(base_size = 14)
plot_r2_norm
ggsave("LatentDim3/plot_r2_norm.pdf", plot = plot_r2_norm, width = 120, height = 90, units = "mm")

plot_r2_normratio <- tibble(x = R2_vec, y = norm_ratio) %>%
  ggplot(aes(x = x, y = y)) + geom_line(size = 1) +
  labs(x = expression(R[paste(Y,'~',U,'|',T)]^2),
       y = expression("Ratio of E("~"|"~tau~"|"~") for non-null and null")) +
  theme_bw(base_size = 15) +
  theme(axis.title.x = element_text(size = 16))
plot_r2_normratio
ggsave("LatentDim3/plot_r2_normratio.pdf", plot = plot_r2_normratio,
       width = 115, height = 105, units = "mm")

## Scatter plot #
scatter_varyR2 <- tibble(tau = as.numeric(tau_cali_mat[, c(trivial_seleted_index, nontrivial_index)]),
                         R2 = rep(R2_vec, 100),
                         index = rep(1:100, each = nrow(tau_cali_mat))) %>%
  ggplot() + geom_point(aes(x = index, y = tau, col = R2), size = 1) +
  scale_color_continuous_sequential(name = expression(R[paste(Y,'~',U,'|',T)]^2~":"),
                                    palette="Sunset", rev=TRUE) +
  # geom_vline(aes(xintercept=55.5), color = "grey", linetype = "dashed") +
  geom_vline(aes(xintercept=55.5), color = "grey") +
  labs(y = expression(tau)) +
  theme_bw()
ggsave("LatentDim3/scatter_varyR2_L1.pdf", plot = scatter_varyR2, width = 350, height = 140, units = "mm")


tau_est_selected <- rbind(tau_cali_mat[c(1, 13),], tau_t)
r2_est_selected <- c(round(100*R2_vec[c(1, 13)], digits = 0), 0) %>% as.character()
RMSE_selected <- apply(tau_est_selected, 1, function(x) sqrt(mean((x - tau)^2))) %>% round(digits = 2)

scatter_est_r2 <- tibble(tau = as.numeric(tau_est_selected[, c(trivial_seleted_index, nontrivial_index)]),
                         R2 = factor(rep(r2_est_selected, 100)),
                         index = rep(1:100, each = 3)) %>%
  ggplot() + geom_point(aes(x = index, y = tau, col = R2), size = 2.5, alpha = 0.7) +
  scale_colour_manual(name = "",
                      values = c("#3B99B1", "#FFE93F", "#F5191C"),
                      labels = c(bquote(R^2~"="~.(r2_est_selected[3])~"%, RMSE="~.(RMSE_selected[3])),
                                 bquote(R^2~"="~.(r2_est_selected[2])~"%, RMSE="~.(RMSE_selected[2])),
                                 bquote(R^2~"="~.(r2_est_selected[1])~"%, RMSE="~.(RMSE_selected[1]))))+
  geom_vline(aes(xintercept=55.5), linetype = "dashed") +
  labs(y = expression(tau), x = 'i', title = "Estimates of Treatment Coefficients") +
  ylim(-14.8, 18.3) +
  annotate(geom = "text", x = c(23, 79), y = c(17.5, 17.5), size = 8,
           label = c("null", "non-null")) +
  theme_bw(base_size = 25) +
  theme(plot.title = element_text(hjust = 0.5, size = 23.5),
        legend.position = "bottom")
scatter_est_r2
ggsave("LatentDim3/scatter_est_r2_L1.pdf", plot = scatter_est_r2, width = 300, height = 160, units = "mm")














