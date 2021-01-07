library(tidyverse)
library(knitr) 
library(kableExtra)
library(colorspace)
library(stringr)


movie <- read.csv("movie_analysis/movie.csv")
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

# sig_actors_index <- which(summary(lmfit_y_t)$coefficients[-1,4] < 0.05)
# sig_actors <- names(tau_t)[sig_actors_index]

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


# Multivariate Calibration Criteria  ----------------------------------------------------
tau_t_sig <- tau_t[sig_actors_index]
tau_benchmark <- 0

objective <- function(gamma_opt){
  return(sqrt(sum((tau_t_sig - u_t_diff %*% gamma_opt)^2))) ## L2(tau_cali)
}

gamma0 <- u_pca$rotation[,1]  # initial gamma
objective(gamma0) # initial objective
# optimize #
solution <- optim(par = gamma0, fn = objective, method = "BFGS")
gamma_opt <- solution$par
gamma_opt
solution$value # objective with calibration
objective(rep(0,ncol(mu_u_t))) # objective without calibration 
# Check R^2 #
R2 <- c(t(gamma_opt) %*% cov_u_t %*% gamma_opt /  sigma_y_t^2)
R2 # 0.7376014

# Compute Cor(gamma'uhat, log_budget) #
cor.test(mu_u_t %*% gamma_opt, movie$log_budget) # cor = 0.2008618, p-value < 2.2e-16

# Varying R^2 #
calibration_opt <- function(objective, init=rep(0, ncol(u_t_diff)), iters=1000) {
  obj_min <- objective(init)
  for (i in 1:iters) {
    gamma0 <- init
    solution <- optim(par = gamma0, fn = objective, method = "BFGS")
    if (solution$value < obj_min) {
      print("got smaller value!")
      obj_min <- solution$value
      gamma_opted1 <- solution$par
    }
    gamma_opted1 <- solution$par
    gamma0 <- ncol(u_t_diff)
  }
  gamma_opted1
}

cal_tau_gamma <- function(objective, init=rep(0, ncol(u_t_diff)), iters=1000) {
  gamma_opted1 <- calibration_opt(objective, init=init, iters=iters)
  return(list(tau_cali=(tau_t_sig - u_t_diff %*% gamma_opted1), gamma_opt=gamma_opted1))
}

penalty_weights <-  seq(0, 10, by = 0.1)
tau_cali_mat <- matrix(NA, nrow=length(penalty_weights), ncol=length(tau_t_sig)) 
colnames(tau_cali_mat) <- sig_actors
gamma_opt_mat <- matrix(NA, nrow=length(penalty_weights), ncol=ncol(u_t_diff) )
R2_vec <- numeric(length(penalty_weights))
count <- 1
for(weight in penalty_weights) {
  results <- cal_tau_gamma(function (g) {
    tau_cali <- tau_t_sig - u_t_diff %*% g
    norm(tau_cali, type = "2") + weight * t(g) %*% cov_u_t %*% g / sigma_y_t^2 ## penalizing R^2
  }, init = u_pca$rotation[,1] , iters=1)
  
  tau_cali_mat[count, ]  <- results$tau_cali
  gamma_opt_mat[count, ] <- results$gamma_opt
  R2_vec[count]  <- (t(results$gamma_opt) %*% cov_u_t %*% results$gamma_opt) / sigma_y_t^2
  count <- count + 1
}

plot(R2_vec, apply(tau_cali_mat, 1, norm, type = "2"), type = "l", 
     ylab = expression(paste("||",tau,"||")[2]), xlab = expression(R[paste(Y,'~',U,'|',T)]^2))
############################### Calibration #############################-#############################
## Calibration with mutlivariate calibration criteria ------------------------------------------------------
cal_tau_calibrated_con <- function(cast, gamma = gamma) {
  return(tau_t[cast] - u_t_diff[cast,] %*% gamma)
}

tau_sens_p_con <- sapply(sig_actors, cal_tau_calibrated_con, gamma = gamma_opt)

# how many has estimates' magnitude shrinking #
sum(abs(tau_sens_p_con) < abs(tau_t_sig))

## checking shinking direction
diff <- tau_benchmark  - tau_t_sig
diff_cali <- tau_sens_p_con - tau_t_sig
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
                        values = c("#3B99B1", "#7CBA96", "#FFC300", "#F5191C"),
                          # divergingx_hcl(5, palette = "Zissou 1")[c(1, 2, 3, 5)],
                        labels = c("0 %",
                                   paste0(round(R2_vec[76]*100, digits = 0), "%"),
                                   paste0(round(R2_vec[16]*100, digits = 0), "%"),
                                   paste0(round(R2_vec[1]*100, digits = 0), "%"))) + 
    scale_x_continuous("Actor j", breaks = 2*(1:length(cast)), labels = gsub('.', ' ', x = cast, fixed =TRUE)) + 
    labs(y = expression(eta[j])) +
    theme_bw(base_size = 15) + 
    theme(plot.title = element_text(hjust = 0.5),
          axis.text.x = element_text(size = 13, angle = 75, hjust = 1),
          legend.text.align = 0,
          legend.title = element_text(size=10))
  return(plot)
}

plot_movie_multicali_r2 <- plot_summary_con(cast = sig_actors)
plot_movie_multicali_r2

# ggsave("plot_movie_multicali_r2.pdf", plot_movie_multicali_r2, width = 300, height=180, units = "mm")
ggsave("plot_movie_multicali_r2.pdf", plot_movie_multicali_r2, width = 300, height=180, units = "mm",
       path = "movie_analysis/Figures")


## Calibration, worst case --------------------------------------------------------------------------------------------------
cal_tau_calibrated <- function(cast, S = 1, R2 = 0) {
  bias <- S*sqrt(R2)*sigma_y_t*sqrt(sum((cov_halfinv %*% u_t_diff[cast,])^2))
  result <- as.numeric(tau_t[cast] - bias)
  return(result)
}

R2_seq <- seq(0.2,0.5, by=0.1)

cal_tau_sens <- function(R2) {
  ntau_sens_p <- sapply(sig_actors, cal_tau_calibrated, S = 1, R2 = R2) * (n_movie$x)[sig_actors_index]
  ntau_sens_n <- sapply(sig_actors, cal_tau_calibrated, S = -1, R2 = R2) * (n_movie$x)[sig_actors_index]
  result <- cbind(ntau_sens_p, ntau_sens_n)
  colnames(result) <- c('p', 'n')
  return(result)
}

ntau_obs_df <- tau_t_sig * (n_movie$x)[sig_actors_index]
ntau_sens_20 <- cal_ntau_sens(0.2)
ntau_sens_30 <- cal_ntau_sens(0.3)
ntau_sens_40 <- cal_ntau_sens(0.4)
ntau_sens_50 <- cal_ntau_sens(0.5)


# plot #
plot_summary <- function(cast) {
  k = length(cast)
  cast = names(sort(ntau_obs_df[cast]))
  bound_df <- tibble(x1 = 2*(1:k),
                     y1 = ntau_sens_50[cast,'p'],
                     x2 = 2*(1:k),
                     y2 = ntau_sens_50[cast,'n'])
  mean_ig_df <- tibble(case = rep(2*(1:k), times=2),
                       R2_0 = c(ntau_obs_df[cast], ntau_obs_df[cast]),
                       R2_20 = c(ntau_sens_20[cast, 'p'], ntau_sens_20[cast, 'n']),
                       R2_30 = c(ntau_sens_30[cast, 'p'], ntau_sens_30[cast, 'n']),
                       R2_40 = c(ntau_sens_40[cast, 'p'], ntau_sens_40[cast, 'n']),
                       R2_50 = c(ntau_sens_50[cast, 'p'], ntau_sens_50[cast, 'n'])) %>% 
    gather(key = "Type", value = "est", - case)
  plot = ggplot() + 
    ungeviz::geom_hpline(data = mean_ig_df, aes(x = case, y = est, col = Type), width = 0.6, size = 1.5)  +
    geom_segment(data = bound_df, aes(x=x1,y=y1,xend=x2,yend=y2), size = 0.5) +
    scale_colour_manual(name = expression(R[paste(tilde(Y),'~',U,'|',T)]^2~":"),
                        values = divergingx_hcl(5,palette = "Zissou 1"),
                        labels = sapply(c(0,R2_seq*100), paste0, "%")) +
    scale_x_continuous("Actor j", breaks = 2*(1:length(cast)), labels = gsub('.', ' ', x = cast, fixed =TRUE)) + 
    labs(y = expression(eta[j])) +
    theme_bw(base_size = 15) + 
    theme(plot.title = element_text(hjust = 0.5),
          axis.text.x = element_text(size = 13, angle = 75, hjust = 1),
          legend.text.align = 0,
          legend.title = element_text(size=10))
  return(plot)
}



plot_movie_worstcase <- plot_summary(sig_actors)
plot_movie_worstcase
ggsave("plot_movie_worstcase.pdf", plot_movie_worstcase, width = 300, height=180, units = "mm",
       path = "movie_analysis/Figures")



# Robustness Value Table ----------------------------------------------------------------------------------------------------
ig_halfwidth <- sapply(sig_actors, cal_tau_calibrated, S = 1, R2 = 1) - tau_t_sig

# cal_RV function #
# estimate: = mean for RV of observed mean 
#            = cbind(lwr, upr) for RV of limit that is closer to the null
# type = 'mean' or 'limit'
cal_RV <- function(estimate, type = c('mean', 'limit')) {
  if (type=='mean'){
    df = abs(estimate)
  } else if (type=='limit') {
    estimate <- estimate %>%
      rownames_to_column() %>%
      filter( (lwr>0 & upr>0) | (lwr<0 & upr<0) ) %>%
      column_to_rownames(var = 'rowname')
    df = apply(abs(estimate), 1, min)
  } else {
    stop("Please choose type to be 'mean' or 'limit' ")
  }
  (df / ig_halfwidth)^2
}


RV_mean <- cal_RV(estimate = tau_t_sig, type='mean')

tau_sig_confint <- confint(lmfit_y_t)[-1,][sig_actors_index,] %>% as.data.frame()
rownames(tau_sig_confint) <- sig_actors
colnames(tau_sig_confint) <- c('lwr', 'upr')
RV_null_limit <- cal_RV(tau_sig_confint, type = 'limit')

dt <- data.frame(Name = gsub('.', ' ', x = sig_actors, fixed =TRUE),
                 Effect = round(ntau_obs_df[sig_actors],2),
                 RV_mean = round(100*RV_mean ,2),
                 RV_null_limit = round(100*RV_null_limit,2)) %>%
  arrange(desc(Effect)) %>%
  mutate_at(c("RV_null_limit"), funs(replace_na(., 'NS'))) %>%
  column_to_rownames(var = "Name") %>%
  'colnames<-'(c('Effect', 'RV_mean(%)', 'RV_limit(%)')) 

dt

kable(dt, 'latex', booktabs = T, linesep = "", align = "c",
      caption = "Robustness Value for Significant Actors") %>%
  kable_styling(position = "center")




