library(patchwork)
library(tidyverse)
library(knitr)
library(kableExtra)
library(colorspace)
library(stringr)
library(patchwork)



movie <- read.csv("movie_analysis/movie.csv")
y <- movie %>% select(log_revenue) %>% as.matrix
t <- movie %>% select(Zoe.Saldana:Ethan.Hawke) %>% as.matrix
colnames(movie) <- make.names(colnames(movie))
n_movie <- read.csv("movie_analysis/n_movie.csv",
                    row.names = 1)

load('movie_analysis/vaedata+sensanaly/lmfit_y_t.Rdata')

# Observed Outcome Model ---------------------------------------------------------------------------------------------

# observed outcome model #
tau_t <- coef(lmfit_y_t)[-1]
names(tau_t) <- rownames(n_movie)
sigma_y_t <- sigma(lmfit_y_t)
yhat <- predict(lmfit_y_t)

# Compute Cor(yhat, log_budget) #
cor.test(yhat, movie$log_budget, method = "spearman")

## Sensitivity Analysis ------------------------------------------------------------------------

# setwd("~/Movie_Application/movie_sensanaly_dim20")
mu_u_t = read.csv("movie_analysis/vaedata+sensanaly/mu_u_t_ise.csv") %>% as.matrix()
cov_u_t = read.csv("movie_analysis/vaedata+sensanaly/cov_u_t_ise.csv") %>% as.matrix()
eigen_cov = eigen(cov_u_t)
cov_halfinv = eigen_cov$vectors %*% diag(eigen_cov$values^{-1/2}) %*% t(eigen_cov$vectors)
u_t_diff <- read.csv('movie_analysis/vaedata+sensanaly/u_t_diff_ise_sigactors.csv') %>% as.matrix()
sig_actors_index <- read.csv('movie_analysis/vaedata+sensanaly/sig_actors_index.csv')$x
sig_actors <- rownames(n_movie)[sig_actors_index]
rownames(u_t_diff) <- sig_actors
u_pca = prcomp(mu_u_t)


## Multivariate Calibration Criteria  ----------------------------------------------------
# R2 = 1 #
tau_t_sig <- tau_t[sig_actors_index]
tau_benchmark <- 0

results_multicali_R1 <- CopSens::gcalibrate(y = y, tr = t, t1 = diag(length(tau_t_sig)),
                                            t2 = matrix(0, ncol = length(tau_t_sig), ncol = length(tau_t_sig)),
                                            calitype = "multicali",
                                            mu_y_dt = as.matrix(tau_t_sig), sigma_y_t = sigma_y_t,
                                            mu_u_dt = u_t_diff, cov_u_t = cov_u_t)
gamma_opt <- as.numeric(results_multicali_R1$gamma)
results_multicali_R1$R2


# Compute Cor(gamma'uhat, log_budget) #
cor.test(mu_u_t %*% gamma_opt, movie$log_budget) # cor = 0.2008618, p-value < 2.2e-16


# Varying R^2 #
penalty_weights <-  seq(0, 10, by = 0.1)
results_multicali_vR2 <- CopSens::gcalibrate(y = y, tr = t, t1 = diag(length(tau_t_sig)),
                                             t2 = matrix(0, ncol = length(tau_t_sig), ncol = length(tau_t_sig)),
                                             calitype = "multicali",
                                             mu_y_dt = as.matrix(tau_t_sig), sigma_y_t = sigma_y_t,
                                             mu_u_dt = u_t_diff, cov_u_t = cov_u_t,
                                             penalty_weight = penalty_weights)

plot(results_multicali_vR2$R2, apply(results_multicali_vR2$est_df[,-1], 2, norm, type = "2"), type = "l",
     ylab = expression(paste("||",tau,"||")[2]), xlab = expression(R[paste(Y,'~',U,'|',T)]^2))


## Analyze multivariate calibration results ##
# how many has estimates' magnitude shrinking #
sum(abs(results_multicali_R1$est_df[,'R2_0.74']) < abs(tau_t_sig))

## checking shrinking direction
diff <- tau_benchmark  - tau_t_sig
diff_cali <- results_multicali_R1$est_df[,'R2_0.74'] - tau_t_sig
# overall #
sum(diff*diff_cali > 0)
sum(diff*diff_cali > 0) / length(diff)
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

# ntau
tau_cali_mat <- t(results_multicali_vR2$est[,-1])
ntau_obs_df <- tau_t_sig * (n_movie$x)[sig_actors_index]
ntau_cali74 <- tau_cali_mat[1,] * (n_movie$x)[sig_actors_index]
ntau_cali30 <- tau_cali_mat[16,] * (n_movie$x)[sig_actors_index]
ntau_cali10 <- tau_cali_mat[76,] * (n_movie$x)[sig_actors_index]

plot_summary_con <- function(cast) {
  k = length(cast)
  cast = names(sort(ntau_obs_df[cast]))
  bound_df <- tibble(x1 = 2*(1:k),
                     y1 = ntau_obs_df[cast],
                     x2 = 2*(1:k),
                     y2 = ntau_cali74[cast])
  mean_ig_df <- tibble(uncali = ntau_obs_df[cast],
                       cali74 = ntau_cali74[cast],
                       cali30 = ntau_cali30[cast],
                       cali10 = ntau_cali10[cast],
                       case = 2*(1:k)) %>%
    gather(key = "Type", value = "est", - case)
  mean_ig_df$Type = factor(mean_ig_df$Type, levels = c("uncali", "cali10", "cali30","cali74"),
                           labels = c("uncali", "cali10", "cali30","cali74"))
  plot = ggplot() +
    ungeviz::geom_hpline(data = mean_ig_df, aes(x = case, y = est, col = Type), width = 0.8, size = 1.5)  +
    # geom_segment(data = bound_df, aes(x=x1,y=y1,xend=x2,yend=y2),
    #              arrow = arrow(length = unit(0.2, "cm"))) +
    geom_hline(yintercept=tau_benchmark, linetype = "dashed") +
    scale_colour_manual(name = expression(R[paste(Y,'~',U,'|',T)]^2~":"),
                        values = divergingx_hcl(7,palette = "Zissou 1")[c(1,2,4,7)],
                          # divergingx_hcl(5, palette = "Zissou 1")[c(1, 2, 3, 5)],
                        labels = c("0 %",
                                   paste0(round(results_multicali_vR2$R2[76]*100, digits = 0), "%"),
                                   paste0(round(results_multicali_vR2$R2[16]*100, digits = 0), "%"),
                                   paste0(round(results_multicali_vR2$R2[1]*100, digits = 0), "%"))) +
    scale_x_continuous("Actor j", breaks = 2*(1:length(cast)), labels = gsub('.', ' ', x = cast, fixed =TRUE)) +
    labs(y = expression(eta[j])) +
    annotate(geom="text", x=10.3, y=55, label="L2 minimization of effects",size = 5.6) +
    theme_bw(base_size = 15) +
    theme(plot.title = element_text(hjust = 0.5),
          axis.text.x = element_text(size = 13, angle = 75, hjust = 1),
          legend.text.align = 0,
          legend.title = element_text(size=14))
  return(plot)
}

plot_movie_multicali_r2 <- plot_summary_con(cast = sig_actors)
plot_movie_multicali_r2

# ggsave("plot_movie_multicali_r2.pdf", plot_movie_multicali_r2, width = 300, height=180, units = "mm",
#        path = "movie_analysis/Figures")


## Calibration, worst case --------------------------------------------------------------------------------------------------

R2_seq <- seq(0.2,0.5, by=0.1)
results_worstcase <- CopSens::gcalibrate(y = y, tr = t, t1 = diag(length(tau_t_sig)),
                                         t2 = matrix(0, ncol = length(tau_t_sig), ncol = length(tau_t_sig)),
                                         calitype = "worstcase",
                                         mu_y_dt = as.matrix(tau_t_sig), sigma_y_t = sigma_y_t,
                                         mu_u_dt = u_t_diff, cov_u_t = cov_u_t,
                                         R2 = R2_seq)

ntau_cali_worstcase <- results_worstcase$est_df *
                       matrix(rep((n_movie$x)[sig_actors_index],
                                  times = ncol(results_worstcase$est_df)),
                              nrow = length(tau_t_sig))


# plot #
plot_summary <- function(cast) {
  k = length(cast)
  cast = names(sort(ntau_obs_df[cast]))
  bound_df <- tibble(x1 = 2*(1:k),
                     y1 = ntau_cali_worstcase[cast,'R2_0.5_lwr'],
                     x2 = 2*(1:k),
                     y2 = ntau_cali_worstcase[cast,'R2_0.5_upr'])
  mean_ig_df <- tibble(case = rep(2*(1:k), times=2),
                       R2_0 = rep(ntau_cali_worstcase[cast,'R2_0'], 2),
                       R2_20 = c(ntau_cali_worstcase[cast, 'R2_0.2_lwr'],
                                 ntau_cali_worstcase[cast, 'R2_0.2_upr']),
                       R2_30 = c(ntau_cali_worstcase[cast, 'R2_0.3_lwr'],
                                 ntau_cali_worstcase[cast, 'R2_0.3_upr']),
                       R2_40 = c(ntau_cali_worstcase[cast, 'R2_0.4_lwr'],
                                 ntau_cali_worstcase[cast, 'R2_0.4_upr']),
                       R2_50 = c(ntau_cali_worstcase[cast, 'R2_0.5_lwr'],
                                 ntau_cali_worstcase[cast, 'R2_0.5_upr'])) %>%
    gather(key = "Type", value = "est", - case)
  plot = ggplot() +
    ungeviz::geom_hpline(data = mean_ig_df, aes(x = case, y = est, col = Type), width = 0.6, size = 1.5)  +
    geom_hline(yintercept=tau_benchmark, linetype = "dashed") +
    geom_segment(data = bound_df, aes(x=x1,y=y1,xend=x2,yend=y2), size = 0.5) +
    scale_colour_manual(name = expression(R[paste(tilde(Y),'~',U,'|',T)]^2~":"),
                        values = divergingx_hcl(7,palette = "Zissou 1")[c(1,3,4,5,6)],
                        labels = sapply(c(0,R2_seq*100), paste0, "%")) +
    scale_x_continuous("", breaks = 2*(1:length(cast)), labels = gsub('.', ' ', x = cast, fixed =TRUE)) +
    labs(y = expression(eta[j])) +
    annotate(geom="text", x=12.8, y=220, label="worst-case ignorance regions", size = 5.6) +
    theme_bw(base_size = 15) +
    theme(plot.title = element_text(hjust = 0.5),
          # axis.text.x = element_text(size = 13, angle = 75, hjust = 1),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          legend.text.align = 0,
          legend.title = element_text(size=14))
  return(plot)
}


plot_movie_worstcase <- plot_summary(sig_actors)
plot_movie_worstcase
# ggsave("plot_movie_worstcase.pdf", plot_movie_worstcase, width = 300, height=180, units = "mm",
#        path = "movie_analysis/Figures")

plot_movie_worstmulti <- plot_movie_worstcase / plot_movie_multicali_r2 & theme(plot.margin = margin(0.1, 0.1, 0.1, 0.1, "cm"))
plot_movie_worstmulti
ggsave("plot_movie_worstmulti.pdf", plot = plot_movie_worstmulti, 
       width = 300, height=200, units = "mm", path = "movie_analysis/Figures")

# Robustness Value Table ----------------------------------------------------------------------------------------------------

RV_mean <- CopSens::cal_rv(y = y, tr = t, t1 = diag(length(tau_t_sig)),
                            t2 = matrix(0, ncol = length(tau_t_sig), ncol = length(tau_t_sig)),
                            mu_y_dt = as.matrix(tau_t_sig), sigma_y_t = sigma_y_t,
                            mu_u_dt = u_t_diff, cov_u_t = cov_u_t) %>% as.numeric()*100

tau_sig_confint <- confint(lmfit_y_t)[-1,][sig_actors_index,] %>% as.data.frame()
rownames(tau_sig_confint) <- sig_actors
RV_null_limit <- CopSens::cal_rv(y = y, tr = t, t1 = diag(length(tau_t_sig)),
                                 t2 = matrix(0, ncol = length(tau_t_sig), ncol = length(tau_t_sig)),
                                 mu_y_dt = apply(abs(tau_sig_confint), 1, min), sigma_y_t = sigma_y_t,
                                 mu_u_dt = u_t_diff, cov_u_t = cov_u_t) %>% as.numeric()*100

dt <- data.frame(Name = gsub('.', ' ', x = sig_actors, fixed =TRUE),
                 Effect = as.numeric(round(ntau_obs_df[sig_actors],2)),
                 RV_mean = RV_mean,
                 RV_null_limit = RV_null_limit) %>%
  arrange(desc(Effect)) %>%
  mutate_at(c("RV_null_limit"), funs(replace_na(., 'NS'))) %>%
  column_to_rownames(var = "Name") %>%
  'colnames<-'(c('Effect', 'RV_mean(%)', 'RV_limit(%)'))
dt
kable(dt, 'latex', booktabs = T, linesep = "", align = "c",
      caption = "Robustness Value for Significant Actors") %>%
  kable_styling(position = "center")




