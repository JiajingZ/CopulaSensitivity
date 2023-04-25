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
                                       R2_constr = 1, normtype = "L1")
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
positive_index <- order(abs(tau), decreasing = T)[1:length(nontrivial_index)]
negative_index <- order(abs(tau), decreasing = T)[(length(nontrivial_index)+1):500]
lmfit_y_t_coef <- summary(lmfit_y_t)$coefficients

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
# write.csv(fpr_cali, row.names = F, file = 'LatentDim3/fpr_cali_dim3_L1.csv')
# write.csv(tpr_cali, row.names = F, file = 'LatentDim3/tpr_cali_dim3_L1.csv')



############################  multivariate calibration(MCC) with various R^2 ########################################################
R2_constr_vec <- seq(1, 0, by = -0.01)

cali_results <- CopSens::gcalibrate(y, tr, t1 = diag(k),
                                    t2 = matrix(0, ncol = k, nrow = k),
                                    calitype ="multicali",
                                    mu_y_dt = tau_t, sigma_y_t = sigma_y_t,
                                    mu_u_dt = u_t_diff, cov_u_t = cov_u_t,
                                    R2_constr = R2_constr_vec, normtype = "L1",
                                    solver="SCS")

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

tau_est_selected <- rbind(tau_cali_mat[c(1, 70),], tau_t)
r2_est_selected <- c(round(100*R2_vec[c(1, 70)], digits = 0), 0) %>% as.character()
RMSE_selected <- apply(tau_est_selected, 1, function(x) sqrt(mean((x - tau)^2))) %>% round(digits = 2)
trivial_index <- (1:k)[!(1:k) %in% nontrivial_index]
set.seed(123)
trivial_seleted_index <- sample(trivial_index, size = 55, replace = FALSE)


## plot the abs(tau -tau.hat), ordered by the magnitude
tau_est_absdiff <- rbind(abs(tau_est_selected[1,] - tau),
                       abs(tau_est_selected[2,] - tau),
                       abs(tau_est_selected[3,] - tau))


## sort abs(tau-tau.hat) by the abs magnitude of tau
# trivial_seleted_index_sorted <- trivial_seleted_index[order(abs(tau[trivial_seleted_index]))]
# nontrivial_index_sorted <- nontrivial_index[order(abs(tau[nontrivial_index]))]

trivial_seleted_index_sorted <- trivial_seleted_index[order(tau_est_absdiff[3,trivial_seleted_index])]
nontrivial_index_sorted <- nontrivial_index[order(tau_est_absdiff[3,nontrivial_index])]

scatter_est_r2 <- tibble(tau = as.numeric(tau_est_absdiff[, c(trivial_seleted_index_sorted, nontrivial_index_sorted)]),
                         R2 = factor(rep(r2_est_selected, 100)),
                         index = rep(1:100, each = 3)) %>%
  ggplot() + geom_point(aes(x = index, y = tau, col = R2), size = 2.5, alpha = 0.7) +
  scale_colour_manual(name = "",
                      values = c("#3B99B1", "#FFE93F", "#F5191C"),
                      labels = c(bquote(R^2~"="~.(r2_est_selected[3])~"%, RMSE="~.(RMSE_selected[3])),
                                 bquote(R^2~"="~.(r2_est_selected[2])~"%, RMSE="~.(RMSE_selected[2])),
                                 bquote(R^2~"="~.(r2_est_selected[1])~"%, RMSE="~.(RMSE_selected[1]))))+
  geom_vline(aes(xintercept=55.5), linetype = "dashed") +
  labs(y = expression("|"~hat(tau)~"-"~tau~"|"), x = 'i', title = "Absolute Difference of the Estimate and True Effect") +
  ylim(-2.5, 18.3) +
  annotate(geom = "text", x = c(23, 79), y = c(17.5, 17.5), size = 8,
           label = c("null", "non-null")) +
  theme_bw(base_size = 25) +
  theme(plot.title = element_text(hjust = 0.5, size = 23.5),
        legend.position = "bottom")
scatter_est_r2
ggsave("LatentDim3/scatter_est_r2_L1_v2.pdf", plot = scatter_est_r2, width = 300, height = 160, units = "mm")




# scatter_est_r2 <- tibble(tau = as.numeric(tau_est_selected[, c(trivial_seleted_index, nontrivial_index)]),
#                          R2 = factor(rep(r2_est_selected, 100)),
#                          index = rep(1:100, each = 3)) %>%
#   ggplot() + geom_point(aes(x = index, y = tau, col = R2), size = 2.5, alpha = 0.7) +
#   scale_colour_manual(name = "",
#                       values = c("#3B99B1", "#FFE93F", "#F5191C"),
#                       labels = c(bquote(R^2~"="~.(r2_est_selected[3])~"%, RMSE="~.(RMSE_selected[3])),
#                                  bquote(R^2~"="~.(r2_est_selected[2])~"%, RMSE="~.(RMSE_selected[2])),
#                                  bquote(R^2~"="~.(r2_est_selected[1])~"%, RMSE="~.(RMSE_selected[1]))))+
#   geom_vline(aes(xintercept=55.5), linetype = "dashed") +
#   labs(y = expression(tau), x = 'i', title = "Estimates of Treatment Coefficients") +
#   ylim(-14.8, 18.3) +
#   annotate(geom = "text", x = c(23, 79), y = c(17.5, 17.5), size = 8,
#            label = c("null", "non-null")) +
#   theme_bw(base_size = 25) +
#   theme(plot.title = element_text(hjust = 0.5, size = 23.5),
#         legend.position = "bottom")
# scatter_est_r2
# ggsave("LatentDim3/scatter_est_r2_L1_v1.pdf", plot = scatter_est_r2, width = 300, height = 160, units = "mm")


###################################### worstcase calibration ############################################################
worstcase_results <- CopSens::gcalibrate(y, tr, t1 = diag(k),
                                         t2 = matrix(0, ncol = k, nrow = k),
                                         calitype ="worstcase",
                                         mu_y_dt = tau_t, sigma_y_t = sigma_y_t,
                                         mu_u_dt = u_t_diff, cov_u_t = cov_u_t,
                                         R2 = c(0.5, 1))
cover_index <- (worstcase_results$est_df[,'R2_1_lwr'] <= tau) & (worstcase_results$est_df[,'R2_1_upr']>= tau)
sum(cover_index)

##### order by true effect ######
## all treatments ##
index <- c(trivial_index[order(tau[trivial_index])],
           nontrivial_index[order(tau[nontrivial_index])])
bound_df <- tibble(x1 = 1:k,
                   y1 = worstcase_results$est_df[index,'R2_1_lwr'],
                   x2 = 1:k,
                   y2 = worstcase_results$est_df[index,'R2_1_upr'])
mean_ig_df <- tibble(case = rep(1:k, times=2),
                     true = rep(tau[index], 2),
                     R2_0 = rep(worstcase_results$est_df[index,'R2_0'], 2),
                     R2_50 = c(worstcase_results$est_df[index, 'R2_0.5_lwr'],
                               worstcase_results$est_df[index, 'R2_0.5_upr']),
                     R2_100 = c(worstcase_results$est_df[index, 'R2_1_lwr'],
                                worstcase_results$est_df[index, 'R2_1_upr'])) %>%
  gather(key = "Type", value = "est", - case)
mean_ig_df$Type <- factor(mean_ig_df$Type,
                          levels = c("true", "R2_0", "R2_50", "R2_100"),
                          labels = c("true", "R2_0", "R2_50", "R2_100"))
plot_nulleffect_worstcase <-
  ggplot() +
  ungeviz::geom_hpline(data = mean_ig_df, aes(x = case, y = est, col = Type), width = 0.8, size = 0.5)  +
  geom_hline(yintercept=0, linetype = "dashed", size = 0.2) +
  geom_vline(aes(xintercept=455.5), linetype = "dashed", size = 0.2) +
  geom_segment(data = bound_df, aes(x=x1,y=y1,xend=x2,yend=y2), size = 0.1) +
  scale_colour_manual(name = "",
                      values = divergingx_hcl(4,palette = "Zissou 1"),
                      labels = c("true effect",
                                 expression(R[paste(tilde(Y), '~', U, '|', T)]^2~'= 0%'),
                                 expression(R[paste(tilde(Y), '~', U, '|', T)]^2~'= 50%'),
                                 expression(R[paste(tilde(Y), '~', U, '|', T)]^2~'= 100%'))) +
  labs(y = expression(tau)) +
  xlim(1,500) +
  theme_bw(base_size = 15) +
  theme(plot.title = element_text(hjust = 0.5),
        # axis.text.x = element_text(size = 13, angle = 75, hjust = 1),
        legend.text.align = 0,
        legend.title = element_text(size=14)) +
  coord_flip()
plot_nulleffect_worstcase
# ggsave("LatentDim3/plot_nulleffect_worstcase.pdf", plot = plot_nulleffect_worstcase,
#        width = 300, height = 450, units = "mm")


## selected treatments ##
uncovered_index <- (1:k)[!cover_index] # 267 307
trivial_index_selected <- c(sample(trivial_index[!(trivial_index %in% uncovered_index)],
                                 size = 53, replace = FALSE),
                            uncovered_index)
index <- c(trivial_index_selected[order(tau[trivial_index_selected])],  ## order by true effect
           nontrivial_index[order(tau[nontrivial_index])])
bound_df <- tibble(x1 = 1:length(index),
                   y1 = worstcase_results$est_df[index,'R2_1_lwr'],
                   x2 = 1:length(index),
                   y2 = worstcase_results$est_df[index,'R2_1_upr'])
mean_ig_df <- tibble(case = rep(1:length(index), times=2),
                     true = rep(tau[index], 2),
                     R2_0 = rep(worstcase_results$est_df[index,'R2_0'], 2),
                     R2_50 = c(worstcase_results$est_df[index, 'R2_0.5_lwr'],
                               worstcase_results$est_df[index, 'R2_0.5_upr']),
                     R2_100 = c(worstcase_results$est_df[index, 'R2_1_lwr'],
                                worstcase_results$est_df[index, 'R2_1_upr'])) %>%
  gather(key = "Type", value = "est", - case)
mean_ig_df$Type <- factor(mean_ig_df$Type,
                          levels = c("true", "R2_0", "R2_50", "R2_100"),
                          labels = c("true", "R2_0", "R2_50", "R2_100"))
rect_x <- c(which(index == 267), which(index == 307))
worstcase_results$est_df[uncovered_index,]
arrow_df <- data.frame(x1 = rect_x,
                       x2 = rect_x,
                       y1 = rep(3,2),
                       y2 = rep(1,2))
plot_nulleffect_worstcase_partial <-
  ggplot() +
  ungeviz::geom_hpline(data = mean_ig_df, aes(x = case, y = est, col = Type), width = 0.9, size = 0.8)  +
  geom_hline(yintercept=0, linetype = "dashed", size = 0.2) +
  geom_vline(aes(xintercept=55.5), linetype = "dashed", size = 0.4) +
  geom_segment(data = bound_df, aes(x=x1,y=y1,xend=x2,yend=y2), size = 0.5) +
  geom_segment(data = arrow_df, aes(x=x1,y=y1,xend=x2,yend=y2), size = 0.8,
               arrow = arrow(length = unit(0.3, "cm")), color = "red") +
  scale_colour_manual(name = "",
                      values = divergingx_hcl(4,palette = "Zissou 1"),
                      labels = c("True Effect",
                                 expression(R[paste(tilde(Y), '~', U, '|', T)]^2~'= 0%'),
                                 expression(R[paste(tilde(Y), '~', U, '|', T)]^2~'= 50%'),
                                 expression(R[paste(tilde(Y), '~', U, '|', T)]^2~'= 100%'))) +
  labs(y = expression(tau)) +
  annotate(geom = "text", x = c(23, 79), y = c(35, 35), size = 6,
           label = c("null", "non-null")) +
  theme_bw(base_size = 15) +
  theme(plot.title = element_text(hjust = 0.5),
        legend.text.align = 0,
        legend.title = element_text(size=14),
        legend.text = element_text(size=14))
plot_nulleffect_worstcase_partial
ggsave("LatentDim3/plot_nulleffect_worstcase_partial_trueorder.pdf",
       plot = plot_nulleffect_worstcase_partial,
       width = 450, height = 180, units = "mm")



