library(CopSens)
library(patchwork)
library(colorspace)
library(kableExtra)

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
data_lm <- cbind(y, rep(1, nrow(micedata)), tr) %>%
  as.data.frame()
colnames(data_lm)[-1] <- paste0("X", 1:(ncol(tr)+1))
lmfit_y_t <- rstanarm::stan_glm(y ~.-1, data = data_lm, seed = 111)

coef_mu_y_t_draws <- as.matrix(lmfit_y_t)[3001:4000,c(-1,-39)]
colnames(coef_mu_y_t_draws) <- colnames(tr)
sigma_y_t_hat_draws <- as.matrix(lmfit_y_t)[3001:4000,39]
tr_contrasts <- tr_test - matrix(apply(tr, 2, median),
                                 nrow=length(percentiles)*ncol(tr),
                                 ncol = ncol(tr), byrow=TRUE)
# x_test <- cbind(rep(0,nrow(tr_test)), tr_contrasts)
# colnames(x_test) <- paste0("X", 1:(ncol(tr)+1))

samples_tibble <- as_tibble(expand.grid(sample=1:1000, quantile=percentiles, gene=colnames(tr))) ## each row of tr_test has 1000 post draws
samples_tibble$yhat <- as.numeric(coef_mu_y_t_draws %*% t(tr_contrasts))
  
# yhat_test <- predict(lmfit_y_t, newdata = data.frame(x_test),
#                      interval = 'confidence')
# naive_grouped <- as_tibble(expand.grid(quantile=percentiles, gene=colnames(tr))) %>%
#   add_column(data.frame(yhat_test))

samples_tibble %>%
  # group_by(gene) %>%
  # mutate(yhat = yhat - yhat[quantile==0.5]) %>%
  group_by(gene, quantile) %>%
  summarize(lwr=quantile(yhat, 0.025), upr=quantile(yhat, 0.975), mid=median(yhat)) ->
  naive_grouped
# visualize the naive estimate mu_y_dt #
ribbon_naive <- naive_grouped %>%
  ggplot() +
  geom_ribbon(aes(x=quantile, ymin=lwr, ymax=upr), alpha=0.5, col="gray") +
  geom_line(aes(x=quantile, y=mid)) +
  facet_wrap(~ gene) + theme_bw() + geom_hline(yintercept=0, linetype="dashed")


gene_sig <- naive_grouped %>%
  group_by(gene) %>%
  mutate(significant = any(lwr*upr > 0)) %>%
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
                                         mu_y_dt = naive_grouped$mid,
                                         sigma_y_t =  mean(sigma_y_t_hat_draws),
                                         mu_u_dt = as.matrix(t1 - t2) %*% t(coef_mu_u_t_hat), 
                                         cov_u_t = cov_u_t_hat,
                                         R2 = R2)
worstcase_lwr <- CopSens::gcalibrate(y = y, tr = tr,
                                     t1 = t1, t2 = t2,
                                     calitype = "worstcase",
                                     mu_y_dt = naive_grouped$lwr,
                                     sigma_y_t =  mean(sigma_y_t_hat_draws),
                                     mu_u_dt = as.matrix(t1 - t2) %*% t(coef_mu_u_t_hat), 
                                     cov_u_t = cov_u_t_hat,
                                     R2 = R2)
worstcase_upr <- CopSens::gcalibrate(y = y, tr = tr,
                                     t1 = t1, t2 = t2,
                                     calitype = "worstcase",
                                     mu_y_dt = naive_grouped$upr,
                                     sigma_y_t =  mean(sigma_y_t_hat_draws),
                                     mu_u_dt = as.matrix(t1 - t2) %*% t(coef_mu_u_t_hat), 
                                     cov_u_t = cov_u_t_hat,
                                     R2 = R2)

calibrated_narrow <- as_tibble(expand.grid(quantile=percentiles, gene=colnames(tr))[rep(1:nrow(naive_grouped), 4),]) %>%
  add_column(lwr = c(worstcase_lwr$est_df[,'R2_0'],
                     as.matrix(worstcase_lwr$est_df)[,c(2,4,6)]),
             upr = c(worstcase_upr$est_df[,'R2_0'],
                     as.matrix(worstcase_upr$est_df)[,c(3,5,7)]),
             R2 = as.character(rep(c(0,R2), each = nrow(naive_grouped))))

plot_gene <- function(target_gene){
  calibrated_narrow %>%
    filter(gene == target_gene) %>%
    ggplot() +
    geom_ribbon(aes(x=quantile, ymin=lwr, ymax=upr, colour = R2), alpha=0.1) +
    geom_hline(yintercept=0, linetype="dashed") +
    # scale_fill_discrete_sequential(name = expression(R[paste(tilde(Y),'~',U,'|',T)]^2~":"),
    #                                 labels = sapply(c(0,R2*100), paste0, "%"),
    #                                 palette = "Oranges", rev = TRUE) +
    # scale_colour_discrete_sequential(name = expression(R[paste(tilde(Y),'~',U,'|',T)]^2~":"),
    #                                   labels = sapply(c(0,R2*100), paste0, "%"),
    #                                   palette = "Oranges", rev = TRUE) +
    scale_colour_manual(name = expression(R[paste(tilde(Y),'~',U,'|',T)]^2~":"),
                        values = divergingx_hcl(4,palette = "Zissou 1"),
                        labels = sapply(c(0,R2*100), paste0, "%")) +
    # scale_fill_manual(name = expression(R[paste(tilde(Y),'~',U,'|',T)]^2~":"),
    #                     values = divergingx_hcl(4,palette = "Zissou 1"),
    #                     labels = sapply(c(0,R2*100), paste0, "%")) +
    labs(title = target_gene, y = "causal effect") +
    theme_bw(base_size = 15) +
    theme(plot.title = element_text(hjust = 0.5))
}

plot_gene("Igfbp2")


# plot_gene("Igfbp2") + plot_gene("Ccnl2") + plot_gene("Vwf") +
#   plot_gene("Irx3") + plot_gene("Abca8a") + plot_gene("Fam105a") +
#   plot_gene("Slc22a3") + plot_layout(nrow = 3, byrow = TRUE)
# 
# plot_gene("2010002N04Rik") + plot_gene("Serpina6") + plot_gene("Ear11") +
#   plot_gene("Ndrg1") + plot_gene("Gapdh") + plot_gene("Kdm4a") +
#   plot_gene("Mest") + plot_layout(nrow = 3, byrow = TRUE)


# Robustness Value #
RV_null_limit_all <- CopSens::cal_rv(y = y, tr = tr,
                                 t1 = t1, t2 = t2,
                                 mu_y_dt = apply(abs(naive_grouped[,-c(1,2)]), 1, min),
                                 sigma_y_t =  mean(sigma_y_t_hat_draws),
                                 mu_u_dt = as.matrix(t1 - t2) %*% t(coef_mu_u_t_hat), 
                                 cov_u_t = cov_u_t_hat)
RV_null_limit_all[is.nan(RV_null_limit_all)]  <- 0
RV_null_limit <- RV_null_limit_all[19*(1:ncol(tr))]
RV_null_mean_all <- CopSens::cal_rv(y = y, tr = tr,
                                     t1 = t1, t2 = t2,
                                     mu_y_dt = naive_grouped$mid,
                                     sigma_y_t =  mean(sigma_y_t_hat_draws),
                                     mu_u_dt = as.matrix(t1 - t2) %*% t(coef_mu_u_t_hat), 
                                     cov_u_t = cov_u_t_hat)
RV_null_mean_all[is.nan(RV_null_mean_all)]  <- 0
RV_null_mean <- RV_null_mean_all[19*(1:ncol(tr))]
RV_null_mean[gene_sig$gene]

# plot_gene_rv <- function(target_gene){
#   tibble(quantile = percentiles,
#          rv = RV_null_limit[names(RV_null_limit) == target_gene])
# }
## y_mu_dt, proprtional to dt, proportional to mu_u_dt, rv same for each gene

rv_siggene <- (RV_null_limit_all[19*(1:ncol(tr))])[gene_sig$gene]
# significant genes that are robust to confoundedness
gene_robust <- which(is.na(rv_siggene)) %>% names()
gene_robust
# significant, not robust to confoundedness
gene_sig$gene[!(as.character(gene_sig$gene) %in% gene_robust)] %>% as.character()

# plot_gene("Vwf") + plot_gene("2010002N04Rik") + plot_gene("Gapdh")
# 
# plot_gene("Abca8a") + plot_gene("Slc22a3") + plot_gene("Serpina6") + 
#   plot_gene("Ear11") + plot_gene("Kdm4a") +plot_gene("Mest")
# 
# plot_gene("Irx3")

######### Robustness Value Table ##########
coef_mu_y_t_summary <- tibble(mean = apply(coef_mu_y_t_draws, 2, mean),
                              lwr = apply(coef_mu_y_t_draws, 2, quantile, probs = 0.025),
                              upr = apply(coef_mu_y_t_draws, 2, quantile, probs = 0.975)) %>%
  round(digits = 2)
dt <- tibble(Name = gene_sig$gene,
                 naive_mean  = (coef_mu_y_t_summary$mean)[gene_sig$gene],
                 naive_limit = (ifelse(abs(coef_mu_y_t_summary$lwr) < abs(coef_mu_y_t_summary$upr),
                                      coef_mu_y_t_summary$lwr,
                                      coef_mu_y_t_summary$upr))[gene_sig$gene],
                 rv_limit = 100*RV_null_limit[gene_sig$gene]) %>%
  arrange(desc(naive_mean)) %>%
  mutate_at(c("rv_limit"), funs(replace_na(., 'robust'))) %>%
  column_to_rownames(var = "Name") %>%
  'colnames<-'(c('tau_{mean}', 'tau_{limit}', 'RV_limit(%)')) 
dt   
kable(dt, 'latex', booktabs = T, linesep = "", align = "c",
      caption = "Robustness Value for Significant Genes") %>%
  kable_styling(position = "center")

# MCC unit ------------------------------------------------------------------------------
mcc_L1_R2_ <- readRDS("~/Documents/Documents /Research-Papers/miao2020/code-miao2020identifying/Results/mcc_L1_R2_unit.rds")
mcc_L1_draws <- readRDS("~/Documents/Documents /Research-Papers/miao2020/code-miao2020identifying/Results/mcc_L1_tibble_unit.rds")
mcc_L2_R2 <- readRDS("~/Documents/Documents /Research-Papers/miao2020/code-miao2020identifying/Results/mcc_L2_R2_unit.rds")
mcc_L2_draws <- readRDS("~/Documents/Documents /Research-Papers/miao2020/code-miao2020identifying/Results/mcc_L2_tibble_unit.rds")

# mcc_L1_tibble <- as_tibble(expand.grid(sample=1:1000, quantile=percentiles, gene=colnames(tr))) %>%
#   add_column(effect_cali = rep(NA, 1000*length(percentiles)*ncol(tr)))
# mcc_L1_R2 <- rep(NA, 1000)
# for(draw_i in 1:1000){
#   print(draw_i)
#   mcc_results_L1 <- gcalibrate(y, tr, t1 = t1, t2 = t2, calitype = "multicali",
#                                mu_y_dt = (samples_tibble %>% filter(sample==draw_i))$yhat , 
#                                sigma_y_t =  sigma_y_t_hat_draws[draw_i],
#                                mu_u_dt = as.matrix(t1 - t2) %*% t(coef_mu_u_t_hat), 
#                                cov_u_t = cov_u_t_hat, normtype = "L1")
#   mcc_L1_tibble$effect_cali[mcc_L1_tibble$sample==draw_i] <- mcc_results_L1$est_df[,2]
#   mcc_L1_R2[draw_i] <- mcc_results_L1$R2
# }

# round(mcc_L1_R2, 2)

# L1 #
# mcc_L1_draws <- matrix(NA, ncol = ncol(coef_mu_y_t_draws),
#                        nrow = nrow(coef_mu_y_t_draws),
#                        dimnames = list(NULL, colnames(coef_mu_y_t_draws)))
# mcc_L1_R2 <- rep(NA, 1000)
# for(draw_i in 1:1000){
#   print(draw_i)
#   mcc_results_L1 <- gcalibrate(y, tr, t1 = t1, t2 = t2, calitype = "multicali",
#                                mu_y_dt =  coef_mu_y_t_draws[draw_i,],
#                                sigma_y_t =  sigma_y_t_hat_draws[draw_i],
#                                mu_u_dt = diag(ncol(tr)) %*% t(coef_mu_u_t_hat),
#                                cov_u_t = cov_u_t_hat, normtype = "L1")
#   mcc_L1_draws[draw_i,] <- mcc_results_L1$est_df[,2]
#   mcc_L1_R2[draw_i] <- mcc_results_L1$R2
# }

mcc_L1Null_tibble <- tibble(gene = rownames(mice_est_nulltr),
                            lwr = apply(mcc_L1_draws[,rownames(mice_est_nulltr)],2,quantile,probs=0.025),
                            mid = apply(mcc_L1_draws[,rownames(mice_est_nulltr)],2,median),
                            upr = apply(mcc_L1_draws[,rownames(mice_est_nulltr)],2,quantile,probs=0.975),
                            miao_nulltr = mice_est_nulltr$esti) %>%
  arrange(mid) %>%
  add_column(case = 1:nrow(mice_est_nulltr))


mcc_L1Null_tibble %>%
  select(-c(gene,mid)) %>%
  round(digits = 2) %>%
  gather(key = "Type", value = "effect", - case) %>%
  ggplot() +
  ungeviz::geom_hpline(aes(x = case, y = effect, col = Type), width = 0.5, size = 1)  +
  geom_segment(aes(x=x1,y=y1,xend=x2,yend=y2), size = 0.5, 
               data = tibble(x1 = mcc_L1Null_tibble$case,
                             y1 = mcc_L1Null_tibble$lwr,
                             x2 = mcc_L1Null_tibble$case,
                             y2 = mcc_L1Null_tibble$upr)) + 
  scale_colour_manual(name = "",
                      values = c("#7CBA96", "#3B99B1", "#7CBA96"),
                      labels = c("MCC_L1_upr","miao_nulltr", "MCC_L1_lwr")) + 
  scale_x_continuous(breaks = 1:nrow(mice_est_nulltr), labels = mcc_L1Null_tibble$gene,
                     limits = c(0.5,nrow(mice_est_nulltr) + 0.5)) +
  labs(y = "Causal Effect", x = "") +
  theme_bw(base_size = 14) +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(size = 13, angle = 75, hjust = 1),
        legend.text.align = 0,
        legend.title = element_text(size=10))

# L2 #
# mcc_L2_draws <- matrix(NA, ncol = ncol(coef_mu_y_t_draws),
#                        nrow = nrow(coef_mu_y_t_draws),
#                        dimnames = list(NULL, colnames(coef_mu_y_t_draws)))
# mcc_L2_R2 <- rep(NA, 1000)
# for(draw_i in 1:1000){
#   print(draw_i)
#   mcc_results_L2 <- gcalibrate(y, tr, t1 = t1, t2 = t2, calitype = "multicali",
#                                mu_y_dt =  coef_mu_y_t_draws[draw_i,],
#                                sigma_y_t =  sigma_y_t_hat_draws[draw_i],
#                                mu_u_dt = diag(ncol(tr)) %*% t(coef_mu_u_t_hat),
#                                cov_u_t = cov_u_t_hat, normtype = "L2")
#   mcc_L2_draws[draw_i,] <- mcc_results_L2$est_df[,2]
#   mcc_L2_R2[draw_i] <- mcc_results_L2$R2
# }

mcc_L2Null_tibble <- tibble(gene = rownames(mice_est_nulltr),
                            lwr = apply(mcc_L2_draws[,rownames(mice_est_nulltr)],2,quantile,probs=0.025),
                            mid = apply(mcc_L2_draws[,rownames(mice_est_nulltr)],2,median),
                            upr = apply(mcc_L2_draws[,rownames(mice_est_nulltr)],2,quantile,probs=0.975),
                            miao_nulltr = mice_est_nulltr$esti) %>%
  arrange(mid) %>%
  add_column(case = 1:nrow(mice_est_nulltr))


mcc_L2Null_tibble %>%
  select(-c(gene,mid)) %>%
  round(digits = 2) %>%
  gather(key = "Type", value = "effect", - case) %>%
  ggplot() +
  ungeviz::geom_hpline(aes(x = case, y = effect, col = Type), width = 0.5, size = 1)  +
  geom_segment(aes(x=x1,y=y1,xend=x2,yend=y2), size = 0.5, 
               data = tibble(x1 = mcc_L2Null_tibble$case,
                             y1 = mcc_L2Null_tibble$lwr,
                             x2 = mcc_L2Null_tibble$case,
                             y2 = mcc_L2Null_tibble$upr)) + 
  scale_colour_manual(name = "",
                      values = c("#7CBA96", "#3B99B1", "#7CBA96"),
                      labels = c("MCC_L2_upr","miao_nulltr", "MCC_L2_lwr")) + 
  scale_x_continuous(breaks = 1:nrow(mice_est_nulltr), labels = mcc_L2Null_tibble$gene,
                     limits = c(0.5,nrow(mice_est_nulltr) + 0.5)) +
  labs(y = "Causal Effect", x = "") +
  theme_bw(base_size = 14) +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(size = 13, angle = 75, hjust = 1),
        legend.text.align = 0,
        legend.title = element_text(size=10))



# saveRDS(mcc_L1_R2,
#         file = "~/Documents/Documents /Research-Papers/miao2020/code-miao2020identifying/Results/mcc_L1_R2_unit.rds")
# saveRDS(mcc_L1_draws,
#         file = "~/Documents/Documents /Research-Papers/miao2020/code-miao2020identifying/Results/mcc_L1_tibble_unit.rds")
# saveRDS(mcc_L2_R2,
#         file = "~/Documents/Documents /Research-Papers/miao2020/code-miao2020identifying/Results/mcc_L2_R2_unit.rds")
# saveRDS(mcc_L2_draws,
#         file = "~/Documents/Documents /Research-Papers/miao2020/code-miao2020identifying/Results/mcc_L2_tibble_unit.rds")





# MCC quantile ------------------------------------------------------------------------------

mcc_L1_R2<- readRDS("~/Documents/Documents /Research-Papers/miao2020/code-miao2020identifying/Results/mcc_L1_R2_quantile.rds")
mcc_L1_tibble <- readRDS("~/Documents/Documents /Research-Papers/miao2020/code-miao2020identifying/Results/mcc_L1_tibble_quantile.rds")
mcc_L2_R2 <- readRDS("~/Documents/Documents /Research-Papers/miao2020/code-miao2020identifying/Results/mcc_L2_R2_quantile.rds")
mcc_L2_tibble <- readRDS("~/Documents/Documents /Research-Papers/miao2020/code-miao2020identifying/Results/mcc_L2_tibble_quantile.rds")

############ L1 ###############
# mcc_L1_tibble <- as_tibble(expand.grid(sample=1:1000, quantile=percentiles, gene=colnames(tr))) %>%
#   add_column(effect_cali = rep(NA, 1000*length(percentiles)*ncol(tr)))
# mcc_L1_R2 <- rep(NA, 1000)
# for(draw_i in 1:1000){
#   print(draw_i)
#   mcc_results_L1 <- gcalibrate(y, tr, t1 = t1, t2 = t2, calitype = "multicali",
#                                mu_y_dt = (samples_tibble %>% filter(sample==draw_i))$yhat , 
#                                sigma_y_t =  sigma_y_t_hat_draws[draw_i],
#                                mu_u_dt = as.matrix(t1 - t2) %*% t(coef_mu_u_t_hat), 
#                                cov_u_t = cov_u_t_hat, normtype = "L1")
#   mcc_L1_tibble$effect_cali[mcc_L1_tibble$sample==draw_i] <- mcc_results_L1$est_df[,2]
#   mcc_L1_R2[draw_i] <- mcc_results_L1$R2
# }

round(mcc_L1_R2, 2)

ribbon_mccl1 <- mcc_L1_tibble %>%
  group_by(gene, quantile) %>%
  summarize(lower=quantile(effect_cali, 0.025), 
            upper=quantile(effect_cali, 0.975), 
            mid=median(effect_cali)) %>%
  ggplot() +
  geom_ribbon(aes(x=quantile, ymin=lower, ymax=upper), alpha=0.5, col="gray") +
  geom_line(aes(x=quantile, y=mid)) +
  ylim(-8, 5) + 
  facet_wrap(~ gene) + theme_bw() + geom_hline(yintercept=0, linetype="dashed")


gene_sig_mccl1 <- mcc_L1_tibble %>%
  group_by(gene, quantile) %>%
  summarize(lwr=quantile(effect_cali, 0.025), 
            upr=quantile(effect_cali, 0.975), 
            mid=median(effect_cali)) %>%
  group_by(gene) %>%
  mutate(significant = any(lwr*upr > 0)) %>%
  select(gene, significant) %>%
  distinct() %>%
  filter(significant==TRUE) %>%
  select(gene)
# insignificant genes become significant
gene_sig_mccl1$gene[!(gene_sig_mccl1$gene %in% gene_sig$gene)]
# significant genes become insignificant
gene_sig$gene[!(gene_sig$gene %in% gene_sig_mccl1$gene)]

plot_contrast_mcc <- function(gene_interested, mcc_tibble, 
                              ylim_lrw = -4, ylim_upr = 3){
  cali_plot <- mcc_tibble %>%
    filter(gene==gene_interested) %>%
    group_by(quantile) %>%
    summarize(lower=quantile(effect_cali, 0.025), 
              upper=quantile(effect_cali, 0.975), 
              mid=median(effect_cali)) %>%
    ggplot() +
    geom_ribbon(aes(x=quantile, ymin=lower, ymax=upper), alpha=0.5, col="gray") +
    geom_line(aes(x=quantile, y=mid)) +
    ylim(ylim_lrw, ylim_upr) + labs(y="causal effect", title = "MCC Calibrated") + 
    theme_bw() + geom_hline(yintercept=0, linetype="dashed") + 
    theme(plot.title = element_text(hjust = 0.5))
  naive_plot <- naive_grouped %>%
    filter(gene==gene_interested) %>%
    ggplot() +
    geom_ribbon(aes(x=quantile, ymin=lwr, ymax=upr), alpha=0.5, col="gray") +
    geom_line(aes(x=quantile, y=mid)) +
    ylim(ylim_lrw, ylim_upr) + labs(y="causal effect", title = "Naive") +
    theme_bw() + geom_hline(yintercept=0, linetype="dashed") + 
    theme(plot.title = element_text(hjust = 0.5))
  naive_plot + cali_plot + 
    plot_annotation(title = paste0(gene_interested, ":"),
                    theme = theme(plot.title = element_text(size = 18)))
}

(plot_contrast_mcc(gene_interested = "Gpld1",
                  mcc_tibble = mcc_L1_tibble) -> plot_mccl1_Gpld1) %>% print()
(plot_contrast_mcc(gene_interested = "Igfbp2",
                   mcc_tibble = mcc_L1_tibble,
                   ylim_lrw = -5,
                   ylim_upr =  4) -> plot_mccl1_Igfbp2) %>% print()
(plot_contrast_mcc(gene_interested = "Mest",
                   mcc_tibble = mcc_L1_tibble) -> plot_mccl1_mest) %>% print()
(plot_contrast_mcc(gene_interested = "Kdm4a",
                   mcc_tibble = mcc_L1_tibble) -> plot_mccl1_kdm4a) %>% print()


############ L2 ##################
# mcc_L2_tibble <- as_tibble(expand.grid(sample=1:1000, quantile=percentiles, gene=colnames(tr))) %>%
#   add_column(effect_cali = rep(NA, 1000*length(percentiles)*ncol(tr)))
# mcc_L2_R2 <- rep(NA, 1000)
# for(draw_i in 1:1000){
#   print(draw_i)
#   mcc_results_L2 <- gcalibrate(y, tr, t1 = t1, t2 = t2, calitype = "multicali",
#                                mu_y_dt = (samples_tibble %>% filter(sample==draw_i))$yhat , 
#                                sigma_y_t =  sigma_y_t_hat_draws[draw_i],
#                                mu_u_dt = as.matrix(t1 - t2) %*% t(coef_mu_u_t_hat), 
#                                cov_u_t = cov_u_t_hat, normtype = "L2")
#   mcc_L2_tibble$effect_cali[mcc_L2_tibble$sample==draw_i] <- mcc_results_L2$est_df[,2]
#   mcc_L2_R2[draw_i] <- mcc_results_L2$R2
# }

mcc_L2_R2 %>% round(2)

ribbon_mccl2 <- mcc_L2_tibble %>%
  group_by(gene, quantile) %>%
  summarize(lower=quantile(effect_cali, 0.025), 
            upper=quantile(effect_cali, 0.975), 
            mid=median(effect_cali)) %>%
  ggplot() +
  geom_ribbon(aes(x=quantile, ymin=lower, ymax=upper), alpha=0.5, col="gray") +
  geom_line(aes(x=quantile, y=mid)) +
  ylim(-8, 5) + 
  facet_wrap(~ gene) + theme_bw() + geom_hline(yintercept=0, linetype="dashed")

gene_sig_mccl2 <- mcc_L2_tibble %>%
  group_by(gene, quantile) %>%
  summarize(lwr=quantile(effect_cali, 0.025), 
            upr=quantile(effect_cali, 0.975), 
            mid=median(effect_cali)) %>%
  group_by(gene) %>%
  mutate(significant = any(lwr*upr > 0)) %>%
  select(gene, significant) %>%
  distinct() %>%
  filter(significant==TRUE) %>%
  select(gene)
# insignificant genes become significant
gene_sig_mccl2$gene[!(gene_sig_mccl2$gene %in% gene_sig$gene)]
# significant genes become insignificant
gene_sig$gene[!(gene_sig$gene %in% gene_sig_mccl2$gene)]


# plot_contrast_mcc <- function(gene_interested, mcc_tibble){
#   cali_plot <- mcc_tibble %>%
#     filter(gene==gene_interested) %>%
#     group_by(quantile) %>%
#     summarize(lower=quantile(effect_cali, 0.025), 
#               upper=quantile(effect_cali, 0.975), 
#               mid=median(effect_cali)) %>%
#     ggplot() +
#     geom_ribbon(aes(x=quantile, ymin=lower, ymax=upper), alpha=0.5, col="gray") +
#     geom_line(aes(x=quantile, y=mid)) +
#     ylim(-4, 3) + labs(y="causal effect", title = "MCC Calibrated") + 
#     theme_bw() + geom_hline(yintercept=0, linetype="dashed") + 
#     theme(plot.title = element_text(hjust = 0.5))
#   naive_plot <- naive_grouped %>%
#     filter(gene==gene_interested) %>%
#     ggplot() +
#     geom_ribbon(aes(x=quantile, ymin=lwr, ymax=upr), alpha=0.5, col="gray") +
#     geom_line(aes(x=quantile, y=mid)) +
#     ylim(-4, 3) + labs(y="causal effect", title = "Naive") +
#     theme_bw() + geom_hline(yintercept=0, linetype="dashed") + 
#     theme(plot.title = element_text(hjust = 0.5))
#   naive_plot + cali_plot + 
#     plot_annotation(title = paste0(gene_interested, ":"),
#                     theme = theme(plot.title = element_text(size = 18)))
# }

(plot_contrast_mcc(gene_interested = "Mest",
                   mcc_tibble = mcc_L2_tibble) -> plot_mccl2_mest) %>% print()
(plot_contrast_mcc(gene_interested = "Kdm4a",
                   mcc_tibble = mcc_L2_tibble) -> plot_mccl2_kdm4a) %>% print()

ggsave("plot_mccl2_mest.pdf", plot = plot_mccl2_mest, width = 170, height = 87, units = "mm",
       path = "~/Documents/Documents /Research-Papers/miao2020/code-miao2020identifying/Results")
ggsave("plot_mccl2_kdm4a.pdf", plot = plot_mccl2_kdm4a, width = 170, height = 87, units = "mm",
       path = "~/Documents/Documents /Research-Papers/miao2020/code-miao2020identifying/Results")


# saveRDS(mcc_L1_R2,
#         file = "~/Documents/Documents /Research-Papers/miao2020/code-miao2020identifying/Results/mcc_L1_R2_quantile.rds")
# saveRDS(mcc_L1_tibble,
#         file = "~/Documents/Documents /Research-Papers/miao2020/code-miao2020identifying/Results/mcc_L1_tibble_quantile.rds")
# saveRDS(mcc_L2_R2,
#         file = "~/Documents/Documents /Research-Papers/miao2020/code-miao2020identifying/Results/mcc_L2_R2_quantile.rds")
# saveRDS(mcc_L2_tibble,
#         file = "~/Documents/Documents /Research-Papers/miao2020/code-miao2020identifying/Results/mcc_L2_tibble_quantile.rds")


