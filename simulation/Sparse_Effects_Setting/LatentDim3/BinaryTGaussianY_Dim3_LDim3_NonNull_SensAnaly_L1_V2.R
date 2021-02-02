library(tidyverse)
library(glmnet)
library(colorspace)
library(patchwork)
setwd("simulation/Sparse_Effects_Setting")

# import data #
y <- read.csv('y.csv')$X0
tr <- read.csv('tr.csv') %>% as.matrix()
tau <- read.csv('tau.csv')$X0


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


# Calibrating ----------------------------------------------------
cali_results_R1 <- CopSens::gcalibrate(y, tr, t1 = diag(k),
                                   t2 = matrix(0, ncol = k, nrow = k),
                                   calitype ="multicali",
                                   mu_y_dt = tau_t, sigma_y_t = sigma_y_t,
                                   mu_u_dt = u_t_diff, cov_u_t = cov_u_t,
                                   penalty_weight = 0,
                                   n_iter = 1000, normtype = "L1")
(R2 <- round(cali_results_R1$R2, digits = 2))
tau_cali <- cali_results_R1$est_df[,2]


## Checking the norm of estimates
norm(as.matrix(tau), type = "1")
norm(as.matrix(tau_t), type = "1")
norm(as.matrix(tau_cali), type = "1")


## Checking whether ATEs shrink or not after calibration
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

## checking shrinking direction
diff <- tau  - tau_t
diff_cali <- tau_cali - tau_t
sum(diff*diff_cali > 0) # 475
sum(diff*diff_cali > 0) / ncol(tr)

# how many are expected to decrease
sum(diff < 0)
# among decrease, how many calibrated are aligned
sum((diff < 0) & (diff*diff_cali > 0))
sum((diff < 0) & (diff*diff_cali > 0)) / sum(diff < 0)
# how many are expected to increase
sum(diff > 0)
# among increase, how many calibrated are aligned
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
# write.csv(fpr_cali, row.names = F, file = 'LatentDim3/fpr_cali_dim3_L1.csv')
# write.csv(tpr_cali, row.names = F, file = 'LatentDim3/tpr_cali_dim3_L1.csv')



############################  Calibration with various R^2 ########################################################
penalty_weights <-  c(seq(0, 8000, by = 200), 20000)
cali_results <- CopSens::gcalibrate(y, tr, t1 = diag(k),
                                   t2 = matrix(0, ncol = k, nrow = k),
                                   calitype ="multicali",
                                   mu_y_dt = tau_t, sigma_y_t = sigma_y_t,
                                   mu_u_dt = u_t_diff, cov_u_t = cov_u_t,
                                   penalty_weight = penalty_weights,
                                   n_iter = 100, normtype = "L1")
R2_vec <- cali_results$R2
tau_cali_mat <- t(cali_results$est_df[,-1])

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
# ggsave("LatentDim3/plot_r2_norm.pdf", plot = plot_r2_norm, width = 120, height = 90, units = "mm")

plot_r2_normratio <- tibble(x = R2_vec, y = norm_ratio) %>%
  ggplot(aes(x = x, y = y)) + geom_line(size = 1) +
  labs(x = expression(R[paste(Y,'~',U,'|',T)]^2),
       y = expression("Ratio of E("~"|"~tau~"|"~") for non-null and null")) +
  theme_bw(base_size = 15) +
  theme(axis.title.x = element_text(size = 16))
plot_r2_normratio
# ggsave("LatentDim3/plot_r2_normratio.pdf", plot = plot_r2_normratio,
#        width = 125, height = 95, units = "mm")

tau_est_selected <- rbind(tau_cali_mat[c(1, 13),], tau_t)
r2_est_selected <- c(round(100*R2_vec[c(1, 13)], digits = 0), 0) %>% as.character()
RMSE_selected <- apply(tau_est_selected, 1, function(x) sqrt(mean((x - tau)^2))) %>% round(digits = 2)
trivial_seleted_index <- sample(1:455, size = 55, replace = FALSE)

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
# ggsave("LatentDim3/scatter_est_r2_L1_temp.pdf", plot = scatter_est_r2, width = 300, height = 160, units = "mm")














