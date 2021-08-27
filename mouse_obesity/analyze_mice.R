library(CopSens)
library(BART)

## load the data #
ytrain <- micedata[,1]
xtrain <- micedata[, 2:31]

apply(xtrain, 2, function(x) quantile(x,


# treatment model #
nfact <- 3
tr_factanal <- factanal(xtrain, factors=nfact, scores = "regression")
B_hat <- diag(sqrt(diag(var(xtrain)))) %*% tr_factanal$loadings
Sigma_t_u_hat <- diag(tr_factanal$uniquenesses * sqrt(diag(var(xtrain))))
u_hat <- tr_factanal$scores
coef_mu_u_t_hat <- t(B_hat) %*% solve(B_hat %*% t(B_hat) + Sigma_t_u_hat)
cov_u_t_hat <- diag(nfact) - t(B_hat) %*% solve(B_hat %*% t(B_hat) + Sigma_t_u_hat) %*% B_hat




## outcome model
## Evaluate at quantiles of the treatment
percentiles  <- seq(0.05, 0.95, by=0.05)
xtest <- matrix(0, nrow=0, ncol=ncol(xtrain))
for(i in 1:ncol(xtrain)) {
  gene_mat <- matrix(apply(xtrain, 2, median), nrow=length(percentiles), ncol=ncol(xtrain), byrow=TRUE)

  gene_mat[, i]  <- as.numeric(quantile(xtrain[i,], probs=percentiles))
  xtest <- rbind(xtest, gene_mat)
}
colnames(xtest)  <-  colnames(xtrain)
rownames(xtest) <- rep(colnames(xtrain), each=length(percentiles))
bart_result <- gbart(y.train=ytrain, x.train=xtrain, x.test=xtest)

samples_tibble <- as_tibble(expand.grid(sample=1:1000, quantile=percentiles, gene=colnames(xtrain)))
samples_tibble$yhat <- as.numeric(bart_result$yhat.test)

samples_tibble %>%
  group_by(gene) %>%
  mutate(yhat = yhat - yhat[quantile==0.5]) %>%
  group_by(gene, quantile) %>%
  summarize(lower=quantile(yhat, 0.025), upper=quantile(yhat, 0.975), mid=median(yhat)) %>%
  ggplot() +
  geom_ribbon(aes(x=quantile, ymin=lower, ymax=upper), alpha=0.5, col="gray") +
  geom_line(aes(x=quantile, y=mid)) +
  facet_wrap(~ gene) + theme_bw() + geom_hline(yintercept=0, linetype="dashed")




samples_tibble %>%
  group_by(gene) %>%
  mutate(yhat = yhat - yhat[quantile==0.5]) %>%
  group_by(gene, quantile) %>%
  summarize(lower=quantile(yhat, 0.025), upper=quantile(yhat, 0.975), mid=median(yhat)) ->
  grouped_ests

## target_gene <- "Igfbp2"
## target_gene <- "Fam105a"
target_gene <- "Apoa4"
k <- ncol(tr)
t1 <- xtest[rownames(xtest) == target_gene, ]
t2 <- t1
t2[, target_gene] <- median(t2[, target_gene])
u_t_diff <- (t1 - t2) %*% t(coef_mu_u_t_hat)
mu_y_dt <- grouped_ests %>% filter(gene==target_gene) %>% pull(lower)
# worst-case calibration #
R2 <- c(0.15, 0.5, 1)
worstcase_results <- gcalibrate(y=ytrain, tr=xtrain,
                                t1=t1, t2=t2,
                                calitype = "worstcase",
                                mu_y_dt = grouped_ests$mid[1:length(percentiles)], sigma_y_t =  bart_result$sigma.mean,
                                mu_u_dt = u_t_diff, cov_u_t = cov_u_t_hat, R2 = R2)
worstcase_results$rv[is.nan(worstcase_results$rv)]  <- 0
rownames(worstcase_results$est_df) <- c(0.025, 0.25, 0.5, 0.75, 0.975)
names(worstcase_results$rv) <- c(0.025, 0.25, 0.5, 0.75, 0.975)
plt <- plot_estimates(worstcase_results, labels=c(0.025, 0.25, 0.5, 0.75, 0.975), order="asis",
                      axis.text.x = element_text(size = 10, angle = 75, hjust = 1))
plt + geom_hline(yintercept=0, linetype="dashed") + ggtitle(target_gene)
