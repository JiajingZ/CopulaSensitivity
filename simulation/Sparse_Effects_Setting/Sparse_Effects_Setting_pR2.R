library(tidyverse)
library(glmnet)
library(colorspace)
library(patchwork)
setwd("simulation/Sparse_Effects_Setting")

# import data #
y <- read.csv('y.csv')$X0
tr <- read.csv('tr.csv') %>% as.matrix()
# k <- ncol(tr)
# tau <- read.csv('tau.csv')$X0
# nontrivial_index <- read.csv('nontrivial_effect_index.csv')$X0 +1

## Partial R^2 based on model: y ~ t, R^2_{Y ~ T_j | T_{-j}} ------------------------------------------------
## y ~ tr ##
lmfit_y_t <- lm(y ~ tr)
# tau_t <- coef(lmfit_y_t)[-1]
yhat_full <- predict(lmfit_y_t)
ss_res_full <- sum((y - yhat_full)^2)
# sigma_y_t <- sigma(lmfit_y_t)
# sigma_y_t^2

## Partial R^2 based on model: y ~ t, R^2_{Y ~ T_j | T_{-j}} ----------------------------------------------------
# lmfit_castfull <- lm(y ~ t)
# yhat_castfull <- predict(lmfit_castfull)
# ss_res_castfull <- sum((y - yhat_castfull)^2)

cal_partial_R2 <- function(var) {
  cat(var, "\n")
  tr_sub <- tr %>% as_tibble() %>% select(-one_of(var)) %>% as.matrix
  fit_sub <- lm(y ~ tr_sub)
  yhat_sub <- predict(fit_sub)
  ss_res_sub <- sum((y - yhat_sub)^2)
  (ss_res_sub - ss_res_full)/ss_res_sub
}

pR2_tj <- sapply(colnames(tr), cal_partial_R2)
pR2_tj %>% sort(decreasing = TRUE)

# write.csv(pR2_tj, row.names = T, file = 'pR2_tj.csv')
# pR2_tj <- read.csv("pR2_tj.csv", row.names = 1)
