library(CopSens)
# library(patchwork)
library(colorspace)
# library(kableExtra)


# load the data #
y <- micedata[,1]
tr <- micedata[, 2:18]

# treatment model #
nfact <- 1
tr_factanal <- factanal(tr, factors=nfact, scores = "regression")
B_hat <- diag(sqrt(diag(var(tr)))) %*% tr_factanal$loadings
Sigma_t_u_hat <- diag(tr_factanal$uniquenesses * sqrt(diag(var(tr))))
u_hat <- tr_factanal$scores
coef_mu_u_t_hat <- t(B_hat) %*% solve(B_hat %*% t(B_hat) + Sigma_t_u_hat)
cov_u_t_hat <- diag(nfact) - t(B_hat) %*% solve(B_hat %*% t(B_hat) + Sigma_t_u_hat) %*% B_hat

# outcome model #
lmfit_y_t <- lm(y ~ ., data = micedata[,1:18])
beta_t <- coef(lmfit_y_t)[-1]
names(beta_t) <- colnames(tr)
sigma_y_t_hat <- sigma(lmfit_y_t)


# Sensitivity Analysis -----------------------------------------------------------------------------------------------------

k <- ncol(tr)
t1 <- diag(k)
t2 <- matrix(0, ncol = k, nrow = k)
u_t_diff <- (t1 - t2) %*% t(coef_mu_u_t_hat)

# worst-case calibration #
R2 <- c(0.15, 0.5, 1)
worstcase_results <- gcalibrate(y, tr, t1 = t1, t2 = t2, calitype = "worstcase",
                                mu_y_dt = as.matrix(beta_t), sigma_y_t =  sigma_y_t_hat,
                                mu_u_dt = u_t_diff, cov_u_t = cov_u_t_hat, R2 = R2)
rownames(worstcase_results$est_df) <- names(beta_t)
names(worstcase_results$rv) <- names(beta_t)

## multivariate calibration ##
# with L1 norm #
multcali_results_L1 <- gcalibrate(y, tr, t1 = t1, t2 = t2, calitype = "multicali",
                                  mu_y_dt = as.matrix(beta_t), sigma_y_t =  sigma_y_t_hat,
                                  mu_u_dt = u_t_diff, cov_u_t = cov_u_t_hat, normtype = "L1")
# with L2 norm #
multcali_results_L2 <- gcalibrate(y, tr, t1 = t1, t2 = t2, calitype = "multicali", 
                                  mu_y_dt = as.matrix(beta_t), sigma_y_t =  sigma_y_t_hat,
                                  mu_u_dt = u_t_diff, cov_u_t = cov_u_t_hat, normtype = "L2")



order_name <- rownames(multcali_results_L2$est_df)[order(multcali_results_L2$est_df[,2])]
summary_df <- data.frame(uncali = round(multcali_results_L2$est_df[order_name, 1], 3),
                         multicali_L1 = round(multcali_results_L1$est_df[order_name, 2], 3),
                         multicali_L2 = round(multcali_results_L2$est_df[order_name, 2], 3),
                         miao_nulltr =  mice_est_nulltr[order_name,]$esti,
                         miao_nulltr_sig = mice_est_nulltr[order_name,]$signif,
                         worstcase_lwr = worstcase_results$est_df[order_name, 'R2_1_lwr'],
                         worstcase_upr = worstcase_results$est_df[order_name, 'R2_1_upr'])
rownames(summary_df) <- order_name
summary_df_narrow <- data.frame(summary_df[,c(1:4)], case = 1:nrow(summary_df)) %>%
  gather(key = "Type", value = "effect", - case)
summary_df_narrow$Type <- factor(summary_df_narrow$Type, 
                                 levels = c("uncali", "miao_nulltr", "multicali_L1", "multicali_L2"),
                                 labels = c("uncali", "miao_nulltr", "multicali_L1", "multicali_L2"))
plot_L1L2Null <- summary_df_narrow %>%
  ggplot() +
  geom_errorbar(aes(ymin = lwr, ymax = upr, x = case), width = 0.5, size = 0.6,
                data = tibble(lwr = summary_df$worstcase_lwr, 
                              upr = summary_df$worstcase_upr,
                              case = 1:nrow(summary_df))) + 
  ungeviz::geom_hpline(aes(x = case, y = effect, col = Type), width = 0.4, size = 1)  +
  scale_colour_manual(name = "",
                      values = c("#F5191C", "#3B99B1", "#7CBA96", "#FFC300"),
                      # divergingx_hcl(5, palette = "Zissou 1")[c(1, 2, 3, 5)],
                      labels = c("naive",
                                 "Miao",
                                 bquote("MCC (L1),"~R^2~"="~
                                          .(round(multcali_results_L1$R2*100,0))~"%"),
                                 bquote("MCC (L2),"~R^2~"="~
                                          .(round(multcali_results_L2$R2*100,0))~"%"))) +
  scale_x_continuous(breaks = 1:k, labels = order_name,
                     limits = c(0.5,k + 0.5)) +
  labs(y = "Causal Effect", x = "") +
  theme_bw(base_size = 14) +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(size = 13, angle = 75, hjust = 1),
        legend.text.align = 0,
        legend.title = element_text(size=10))
print(plot_L1L2Null)


ggsave("mice_L1L2Null.pdf", plot = plot_L1L2Null,
       width = 270, height=130, units = "mm", 
       path = "~/Desktop/CopulaSensitivity/mice_analysis")




















