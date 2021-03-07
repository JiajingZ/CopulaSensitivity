library(tidyverse)
library(glmnet)
library(colorspace)
library(patchwork)
setwd("simulation/Sparse_Effects_Setting")

# import data #
y <- read.csv('y.csv')$X0
tr <- read.csv('tr.csv') %>% as.matrix()
tau <- read.csv('tau.csv')$X0
nontrivial_index <- read.csv('nontrivial_effect_index.csv')$X0 +1
tau[nontrivial_index]

# fit observed outcome model --------------------------------------------------------------------------------------------------------
## y ~ tr ##
lmfit_y_t <- lm(y ~ tr)
tau_t <- coef(lmfit_y_t)[-1]
yhat <- predict(lmfit_y_t)
sigma_y_t <- sigma(lmfit_y_t)
sigma_y_t^2

# Sensitivity Analysis --------------------------------------------------------------------------------------------------
latent_dim <- 4
mu_u_t <- read.csv('LatentDim4/mu_u_t_ise.csv') %>% as.matrix()
cov_u_t <- read.csv('LatentDim4/cov_u_t_ise.csv') %>% as.matrix()
eigen_cov <- eigen(cov_u_t)
cov_halfinv <- eigen_cov$vectors %*% diag(eigen_cov$values^{-1/2}) %*% t(eigen_cov$vectors)
u_t_diff_org_all <- read.csv('LatentDim4/u_t_diff_org_all.csv') %>% as.matrix()
u_t_diff <- u_t_diff_org_all
u_pca <- prcomp(mu_u_t)

# How much confounding components are captured by VAE (mu_u_t_ise)
gamma <- rep(100, 3)
u <- read.csv('u.csv') %>% as.matrix()
y_confound <- u %*% gamma
## gamma'u ~ uhat ##
lmfit_confound <- lm(y_confound ~ mu_u_t)
summary(lmfit_confound) # Multiple R-squared:  0.7957,	Adjusted R-squared:  0.7957


# Calibrating ----------------------------------------------------
cali_results_R1 <- CopSens::gcalibrate(y, tr, t1 = diag(k),
                                       t2 = matrix(0, ncol = k, nrow = k),
                                       calitype ="multicali",
                                       mu_y_dt = tau_t, sigma_y_t = sigma_y_t,
                                       mu_u_dt = u_t_diff, cov_u_t = cov_u_t,
                                       R2_constr = 1, normtype = "L1")
(R2 <- round(cali_results_R1$R2, digits = 2))
tau_cali <- cali_results_R1$est_df[,2]


## Checking the norm of estimates
norm(as.matrix(tau), type = "1")
norm(as.matrix(tau_t), type = "1")
norm(as.matrix(tau_cali), type = "1")

## Checking whether ATEs shrink or not after calibration
sum(tau_t^2/length(tau_t))
sum(tau_cali^2/length(tau_cali))

# RMSE #
# overall #
sqrt(mean((tau_t - tau)^2))
sqrt(mean((tau_cali - tau)^2))
# null #
sqrt(mean(((tau_t - tau)[-nontrivial_index])^2))
sqrt(mean(((tau_cali - tau)[-nontrivial_index])^2))
# non-null #
sqrt(mean(((tau_t - tau)[nontrivial_index])^2))
sqrt(mean(((tau_cali - tau)[nontrivial_index])^2))

## checking shrinking direction
diff <- tau  - tau_t
diff_cali <- tau_cali - tau_t
sum(diff*diff_cali > 0)
sum(diff*diff_cali > 0) / ncol(tr)

# how many are expected to decrease
sum(diff < 0) # 233
# among decrease, how many calibrated are aligned
sum((diff < 0) & (diff*diff_cali > 0))
sum((diff < 0) & (diff*diff_cali > 0)) / sum(diff < 0)
# how many are expected to increase
sum(diff > 0) # 267
# among increase, how many calibrated are aligned
sum((diff > 0) & (diff*diff_cali > 0))
sum((diff > 0) & (diff*diff_cali > 0)) / sum(diff > 0)

mean(abs(tau_t)) - mean(abs(tau_cali))
mean(abs(tau_t[nontrivial_index])) - mean(abs(tau_cali[nontrivial_index]))
mean(abs(tau_t[-nontrivial_index])) - mean(abs(tau_cali[-nontrivial_index]))

cbind(tau, tau_t, tau_cali)[nontrivial_index,]
# among nontrivial ones, how many are aligned #
sum((diff*diff_cali > 0)[nontrivial_index])

tau_sort <- abs(tau) %>% sort(decreasing = T, index.return = T)
nontrivial_index <- tau_sort$ix[1:45]

############################## ROC #########################################################################
positive_index <- order(abs(tau), decreasing = T)[1:45]
negative_index <- order(abs(tau), decreasing = T)[46:500]

tpr_naive = fpr_naive  = ppv_naive = tpr_cali = fpr_cali = ppv_cali = NULL
for (k in seq(1, 500, by = 1)) {
  #### Naive ####
  positive_index_naive <- order(abs(tau_t), decreasing = TRUE)[1:k]
  # TPR  = predicted_positive/true_positive #
  tpr_naive <- c(tpr_naive, mean(positive_index %in% positive_index_naive))
  # FPR = predicted_positive/ true_negative #
  fpr_naive <- c(fpr_naive, mean(negative_index %in% positive_index_naive))
  # PPV = true_positive /  predicted_positive #
  ppv_naive <- c(ppv_naive, mean(positive_index_naive %in% positive_index))
  #### Calibrated ####
  positive_index_cali <- order(abs(tau_cali), decreasing = TRUE)[1:k]
  # TPR  = predicted_positive/true_positive #
  tpr_cali <- c(tpr_cali, mean(positive_index %in% positive_index_cali))
  # FPR = predicted_positive/ true_negative
  fpr_cali <- c(fpr_cali, mean(negative_index %in% positive_index_cali))
  # PPV = true_positive /  predicted_positive #
  ppv_cali <- c(ppv_cali, mean(positive_index_cali %in% positive_index))
}
# write.csv(fpr_cali, row.names = F, file = 'LatentDim4/fpr_cali_dim4_L1.csv')
# write.csv(tpr_cali, row.names = F, file = 'LatentDim4/tpr_cali_dim4_L1.csv')


