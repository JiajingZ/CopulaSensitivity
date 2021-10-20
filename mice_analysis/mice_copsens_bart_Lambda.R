library(CopSens)
library(BART)
library(patchwork)
library(colorspace)

## load the data #
y <- micedata[,1]
tr <- micedata[, 2:38] ## 2-18 are nonnull tr
# treatment contrast of interest #
percentiles  <- seq(0.05, 0.95, by=0.05)
tr_test <- matrix(0, nrow=0, ncol=ncol(tr))
for(i in 1:ncol(tr)) {
  gene_mat <- matrix(apply(tr, 2, median), nrow=length(percentiles), ncol=ncol(tr), byrow=TRUE)  ## median, repeated in rows
  gene_mat[, i]  <- as.numeric(quantile(tr[,i], probs=percentiles)) ## typo
  tr_test <- rbind(tr_test, gene_mat)
}
colnames(tr_test)  <-  colnames(tr)
rownames(tr_test) <- rep(colnames(tr), each=length(percentiles))



# treatment model #
# ut_cv <- pcaMethods::kEstimate(tr, method = "ppca", allVariables = TRUE)
nfact <- 3 ## by scree plot of the singular values of tr
ut_factanal <- factanal(tr, factors = nfact)
B_hat <- diag(sqrt(diag(var(tr)))) %*% ut_factanal$loadings
Lambda_hat <- diag(ut_factanal$uniquenesses * sqrt(diag(var(tr))))
cov_u_t_hat <- diag(nfact) -
  t(B_hat) %*% solve(B_hat %*% t(B_hat) + Lambda_hat) %*% B_hat
coef_mu_u_t_hat <- t(B_hat) %*% solve(B_hat %*% t(B_hat) + Lambda_hat)


# outcome model #
## Evaluate at quantiles of the treatment
bartfit <- BART::gbart(y.train=y, x.train=tr, x.test=tr_test)

samples_tibble <- as_tibble(expand.grid(sample=1:1000, quantile=percentiles, gene=colnames(tr))) ## each row of tr_test has 1000 post draws
samples_tibble$yhat <- as.numeric(bartfit$yhat.test)
# naive estimate mu_y_dt #
samples_tibble %>%
  group_by(gene) %>%
  mutate(yhat = yhat - yhat[quantile==0.5]) %>%
  group_by(gene, quantile) %>%
  summarize(lower=quantile(yhat, 0.025), upper=quantile(yhat, 0.975), mid=median(yhat)) ->
  grouped_ests
# visualize the naive estimate mu_y_dt #
grouped_ests %>%
  ggplot() +
  geom_ribbon(aes(x=quantile, ymin=lower, ymax=upper), alpha=0.5, col="gray") +
  geom_line(aes(x=quantile, y=mid)) +
  ylim(-8, 5) + 
  facet_wrap(~ gene) + theme_bw() + geom_hline(yintercept=0, linetype="dashed")

gene_sig <- grouped_ests %>%
  group_by(gene) %>%
  mutate(significant = any(lower*upper > 0)) %>% 
  select(gene, significant) %>%
  distinct() %>%
  filter(significant==TRUE) %>%
  select(gene)
gene_sig 

# Sensitivity Analysis -----------------------------------------------------------------------------------

t1 <- tr_test
t2 <- matrix(apply(tr, 2, median), nrow=length(percentiles)*ncol(tr),
             ncol = ncol(tr), byrow=TRUE)
R2 <- c(0.15, 0.5, 1)
worstcase_results <- CopSens::gcalibrate(y = y, tr = tr,
                                t1 = t1, t2 = t2,
                                calitype = "worstcase",
                                mu_y_dt = grouped_ests$mid,
                                sigma_y_t =  bartfit$sigma.mean,
                                mu_u_dt = as.matrix(t1 - t2) %*% t(coef_mu_u_t_hat), 
                                cov_u_t = cov_u_t_hat,
                                R2 = R2)
worstcase_lwr <- CopSens::gcalibrate(y = y, tr = tr,
                                     t1 = t1, t2 = t2,
                                     calitype = "worstcase",
                                     mu_y_dt = grouped_ests$lower,
                                     sigma_y_t =  bartfit$sigma.mean,
                                     mu_u_dt = as.matrix(t1 - t2) %*% t(coef_mu_u_t_hat), 
                                     cov_u_t = cov_u_t_hat,
                                     R2 = R2)
worstcase_upr <- CopSens::gcalibrate(y = y, tr = tr,
                                     t1 = t1, t2 = t2,
                                     calitype = "worstcase",
                                     mu_y_dt = grouped_ests$upper,
                                     sigma_y_t =  bartfit$sigma.mean,
                                     mu_u_dt = as.matrix(t1 - t2) %*% t(coef_mu_u_t_hat), 
                                     cov_u_t = cov_u_t_hat,
                                     R2 = R2)

calibrated_narrow <- as_tibble(expand.grid(quantile=percentiles, gene=colnames(tr))[rep(1:nrow(grouped_ests), 4),]) %>%
  add_column(lwr = c(worstcase_lwr$est_df[,'R2_0'],
                     as.matrix(worstcase_lwr$est_df)[,c(2,4,6)]),
             upr = c(worstcase_upr$est_df[,'R2_0'],
                     as.matrix(worstcase_upr$est_df)[,c(3,5,7)]),
             R2 = as.character(rep(c(0,R2), each = nrow(grouped_ests))))


plot_gene <- function(target_gene){
  calibrated_narrow %>%
    filter(gene == target_gene) %>%
    ggplot() +
    geom_ribbon(aes(x=quantile, ymin=lwr, ymax=upr, colour = R2), alpha=0.1) +
    geom_hline(yintercept=0, linetype="dashed") +
    scale_colour_manual(name = expression(R[paste(tilde(Y),'~',U,'|',T)]^2~":"),
                        values = divergingx_hcl(4,palette = "Zissou 1"),
                        labels = sapply(c(0,R2*100), paste0, "%")) +
    labs(title = target_gene, y = "causal effect") +
    theme_bw(base_size = 15) +
    theme(plot.title = element_text(hjust = 0.5))
}

plot_gene <- function(target_gene){
  calibrated_narrow %>%
    filter(gene == target_gene) %>%
    ggplot() +
    geom_ribbon(aes(x=quantile, ymin=lwr, ymax=upr, colour = R2), alpha=0.1) +
    geom_hline(yintercept=0, linetype="dashed") +
    scale_colour_manual(name = expression(R[paste(tilde(Y),'~',U,'|',T)]^2~":"),
                        values = divergingx_hcl(4,palette = "Zissou 1"),
                        labels = sapply(c(0,R2*100), paste0, "%")) +
    labs(title = target_gene, y = "causal effect") +
    theme_bw(base_size = 15) +
    theme(plot.title = element_text(hjust = 0.5))
}

plot_gene("Igfbp2")
ggsave("mice_igfbp2_bart.pdf", plot = plot_gene("Igfbp2"),
       width = 145, height = 100, units = "mm", 
       path = "/Users/jiajing/Documents/Documents /Research-Papers/miao2020/code-miao2020identifying/Results")


# Robustness Value #
RV_null_limit <- CopSens::cal_rv(y = y, tr = tr,
                                 t1 = t1, t2 = t2,
                                 mu_y_dt = apply(abs(grouped_ests[,-c(1,2)]), 1, min),
                                 sigma_y_t =  bartfit$sigma.mean,
                                 mu_u_dt = as.matrix(t1 - t2) %*% t(coef_mu_u_t_hat), 
                                 cov_u_t = cov_u_t_hat)
RV_null_limit[is.nan(RV_null_limit)]  <- 0


target_gene = "Igfbp2"
sig_index <- (grouped_ests$lower*grouped_ests$upper > 0)[grouped_ests$gene==target_gene]
(RV_null_limit[names(RV_null_limit) == target_gene])[sig_index]
# credible intervals of naive estimates include zero, not robust to estimation uncertainty
# credible intervals of naive estimates exclude zero, all robust to confoundedness.

# MCC ------------------------------------------------------------------------------
naive_est_tibble <- samples_tibble %>%
  group_by(gene) %>%
  mutate(yhat = yhat - yhat[quantile==0.5])

# L1 #
mcc_L1_tibble <- as_tibble(expand.grid(sample=1:1000, quantile=percentiles, gene=colnames(tr))) %>%
  add_column(effect_cali = rep(NA, 1000*length(percentiles)*ncol(tr)))
mcc_L1_R2 <- rep(NA, 1000)
for(draw_i in 1:1000){
  print(draw_i)
  mcc_results_L1 <- gcalibrate(y, tr, t1 = t1, t2 = t2, calitype = "multicali",
                               mu_y_dt = (naive_est_tibble %>% filter(sample==draw_i))$yhat , 
                               sigma_y_t =  bartfit$sigma[100+draw_i],
                               mu_u_dt = as.matrix(t1 - t2) %*% t(coef_mu_u_t_hat), 
                               cov_u_t = cov_u_t_hat, normtype = "L1")
  mcc_L1_tibble$effect_cali[mcc_L1_tibble$sample==draw_i] <- mcc_results_L1$est_df[,2]
  mcc_L1_R2[draw_i] <- mcc_results_L1$R2
}

round(mcc_L1_R2, 2)

mcc_L1_tibble %>%
  group_by(gene, quantile) %>%
  summarize(lower=quantile(effect_cali, 0.025), 
            upper=quantile(effect_cali, 0.975), 
            mid=median(effect_cali)) %>%
  ggplot() +
  geom_ribbon(aes(x=quantile, ymin=lower, ymax=upper), alpha=0.5, col="gray") +
  geom_line(aes(x=quantile, y=mid)) +
  ylim(-8, 5) + 
  facet_wrap(~ gene) + theme_bw() + geom_hline(yintercept=0, linetype="dashed")

# L2 #
mcc_L2_tibble <- as_tibble(expand.grid(sample=1:1000, quantile=percentiles, gene=colnames(tr))) %>%
  add_column(effect_cali = rep(NA, 1000*length(percentiles)*ncol(tr)))
mcc_L2_R2 <- rep(NA, 1000)
for(draw_i in 1:1000){
  print(draw_i)
  mcc_results_L2 <- gcalibrate(y, tr, t1 = t1, t2 = t2, calitype = "multicali",
                               mu_y_dt = (naive_est_tibble %>% filter(sample==draw_i))$yhat , 
                               sigma_y_t =  bartfit$sigma[100+draw_i],
                               mu_u_dt = as.matrix(t1 - t2) %*% t(coef_mu_u_t_hat), 
                               cov_u_t = cov_u_t_hat, normtype = "L2")
  mcc_L2_tibble$effect_cali[mcc_L2_tibble$sample==draw_i] <- mcc_results_L2$est_df[,2]
  mcc_L2_R2[draw_i] <- mcc_results_L2$R2
}
mcc_L2_R2
mcc_L2_tibble %>%
  group_by(gene, quantile) %>%
  summarize(lower=quantile(effect_cali, 0.025), 
            upper=quantile(effect_cali, 0.975), 
            mid=median(effect_cali)) %>%
  ggplot() +
  geom_ribbon(aes(x=quantile, ymin=lower, ymax=upper), alpha=0.5, col="gray") +
  geom_line(aes(x=quantile, y=mid)) +
  ylim(-8, 5) + 
  facet_wrap(~ gene) + theme_bw() + geom_hline(yintercept=0, linetype="dashed")
mcc_L2_R2 %>% round(2)






















# draft -------------------------------------------------------------------------------------------------


worstcase_results$rv[is.nan(worstcase_results$rv)]  <- 0
plot_gene <- function(gene){
  gene_index <- which(colnames(tr) == gene)
  result_index <- ((gene_index-1)*length(percentiles) + 1):(gene_index*length(percentiles))
  est <- list(est_df = worstcase_results$est_df[result_index,],
              rv = worstcase_results$rv[result_index],
              R2 = worstcase_results$R2)
  CopSens::plot_estimates(est, order="asis",
                          labels = percentiles) +
    geom_hline(yintercept=0, linetype="dashed") + ggtitle(gene)
}

for(i in (1:ncol(tr))){
  plt <- plot_gene(colnames(tr)[i])
  print(plt)
}

gene_sig
plot_gene("Igfbp2")
plot_gene("Apoa4")

# (plot_gene("Igfbp2") + plot_gene("Lamc1") + plot_gene("Sirpa")) / 
#   (plot_gene("Gstm2") + plot_gene("Ccnl2") + plot_gene("Glcci1")) /
#   (plot_gene("Vwf") + plot_gene("Irx3") + plot_gene("Apoa4"))
# 
# (plot_gene("Socs2") + plot_gene("Avpr1a") + plot_gene("Abca8a")) / 
#   (plot_gene("Gpld1") + plot_gene("Fam105a") + plot_gene("Dscam")) /
#   (plot_gene("Slc22a3") + plot_gene("N04Rik") + plot_gene("Emr4"))

