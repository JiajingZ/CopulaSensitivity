library(tidyverse)

# import data #
movie <- read_csv("movie_analysis/movie.csv")
colnames(movie) <- make.names(colnames(movie))
y <- movie %>% select(log_revenue) %>% as.matrix
t <- movie %>% select(Zoe.Saldana:Ethan.Hawke) %>% as.matrix
tb <- movie %>% select(c(log_budget, Zoe.Saldana:Ethan.Hawke)) %>% as.matrix
budget <- movie %>% select(log_budget) %>% as.matrix

# fitting lm -----------------------------------------------------------------------------------------------
lmfit_y_t <- lm(y ~ t)
summary(lmfit_y_t)
# save(lmfit_y_t,file = 'movie_analysis/outcome_model/lmfit_y_t.Rdata')
beta_t_lm <- coef(lmfit_y_t)[-1]

mean(beta_t_lm^2) # 0.3251733
mean(beta_t^2) # 0.08598833
var(beta_t_lm) # 0.2616276
var(beta_t) # 0.0683734 