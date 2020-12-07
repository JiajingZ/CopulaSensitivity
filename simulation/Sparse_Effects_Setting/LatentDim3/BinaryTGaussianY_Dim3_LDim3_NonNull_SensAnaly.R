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
objective <- function(gamma_opt){
  return(mad(tau_t - u_t_diff %*% gamma_opt))
  # return(sqrt(mean((tau_t - u_t_diff %*% gamma_opt)^2)))
  # return(var(tau_t - u_t_diff %*% gamma_opt))
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
objective(rep(0, latent_dim)) # 6.252661
obj_min # 0.8590348
gamma_int_min # -20.458805 -13.163856  -4.839688
gamma_opted1 <- c(39.59563, 9.46200, -144.84033)
## R^2 ##
R2 <- t(gamma_opted1) %*% cov_u_t %*% gamma_opted1 / sigma_y_t^2
R2 # 0.72






# Calibrating ----------------------------------------------------

# Based on gamma #
cal_tau_calibrated_gamma <- function(gamma = gamma) {
  return(tau_t - u_t_diff %*% gamma)
} 

tau_cali <- cal_tau_calibrated_gamma(gamma = gamma_opted1)


## Checking the variance of estimates 
mad(tau)
mad(tau_t)
mad(tau_cali) 

## Checking whether TEs shrink or not after calibration
sum(tau_t^2/length(tau_t)) # 39.68797
sum(tau_cali^2/length(tau_cali)) # 0.8269908


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

mad(tau_t) - mad(tau_cali) 
mad(tau_t[nontrivial_index]) - mad(tau_cali[nontrivial_index])
mad(tau_t[-nontrivial_index]) - mad(tau_cali[-nontrivial_index])

cbind(tau, tau_t, tau_cali)[nontrivial_index,]
# among nontrivial ones, how many are aligned #
sum((diff*diff_cali > 0)[nontrivial_index])


nontrivial_count <- sum(abs(tau) > 0.1)
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

(auc_naive <- DescTools::AUC(x = fpr_naive, y = tpr_naive))
(auc_cali <- DescTools::AUC(x = fpr_cali, y = tpr_cali)) 

plot_roc

#### Scatter Plot ----------------------------------------------------------------

nontrivial_index <- which(abs(tau) > 0.1)

scatter_true <- tibble(index = 1:k,
                       effect = c(tau[-nontrivial_index], tau[nontrivial_index]),
                       group = c(rep('null', 455), rep('nonnull', 45))) %>%
  ggplot(aes(x = index, y = effect, colour = group)) + 
  geom_point(alpha = 0.7) + 
  scale_colour_manual(values = c('#F39C12', '#5F6A6A')) + 
  labs(y = expression(tau), x = 'i', title = "True effects:") + 
  ylim(-18, 18) + 
  theme_bw(base_size = 19) +
  theme(plot.title = element_text(size=19))
# theme(plot.title = element_text(hjust = 0.5))

(rmse_naive <- sqrt(mean((tau_t - tau)^2)) %>% round(digits = 2))
scatter_naive <- tibble(index = 1:k,
                        effect = c(tau_t[-nontrivial_index], tau_t[nontrivial_index]),
                        group = c(rep('null', 455), rep('nonnull', 45))) %>%
  ggplot(aes(x = index, y = effect, colour = group)) + 
  geom_point(alpha = 0.7) + 
  scale_colour_manual(values = c('#F39C12', '#5F6A6A')) + 
  labs(y = expression(tau), x = 'i', title=paste0("Uncalibrated effects (RMSE = ", rmse_naive, "):")) + 
  ylim(-18, 18) + 
  theme_bw(base_size = 19) +
  theme(plot.title = element_text(size=19))
# theme(plot.title = element_text(hjust = 0.5))


(rmse_cali <- sqrt(mean((tau_cali - tau)^2)) %>% round(digits = 2))
scatter_cali <- tibble(index = 1:k,
                       effect = c(tau_cali[-nontrivial_index], tau_cali[nontrivial_index]),
                       group = c(rep('null', 455), rep('nonnull', 45))) %>%
  ggplot(aes(x = index, y = effect, colour = group)) + 
  scale_colour_manual(values = c('#F39C12', '#5F6A6A')) + 
  geom_point(alpha = 0.7) + 
  labs(y = expression(tau), x = 'i', title= paste0("Calibrated effects (RMSE = ", rmse_cali, "):")) + 
  ylim(-18, 18) + 
  theme_bw(base_size = 19)  + 
  theme(plot.title = element_text(size=19))
# theme(plot.title = element_text(hjust = 0.5))


scatter_coef <- scatter_true/scatter_naive/scatter_cali
# + plot_annotation(title = paste0("Scatter Plots for Coefficients (Latent Dim. = ", latent_dim, ')'))
scatter_coef


# ggsave("LatentDim3/scatter_coef.pdf", plot = scatter_coef, width = 300, height = 190, units = "mm")





