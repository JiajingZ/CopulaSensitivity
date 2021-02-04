library(tidyverse)


movie <- read_csv("movie_analysis/movie.csv")
colnames(movie) <- make.names(colnames(movie))
y <- movie %>% select(log_revenue) %>% as.matrix
x <- movie %>% select(-c(log_revenue, id, TV.Movie)) %>% as.matrix
t <- movie %>% select(Zoe.Saldana:Ethan.Hawke) %>% as.matrix

## Partial R^2 based on model: y ~ t + x --------------------------------------------------------------------
lmfit_full <- lm(y ~ x)
summary(lmfit_full)
yhat_full <- predict(lmfit_full)
ss_res_full <- sum((y - yhat_full)^2)

cal_partial_R2 <- function(var) {
  x_sub <- x %>% as_tibble() %>% select(-one_of(var)) %>% as.matrix
  fit_sub <- lm(y ~ x_sub)
  yhat_sub <- predict(fit_sub)
  ss_res_sub <- sum((y - yhat_sub)^2)
  (ss_res_sub - ss_res_full)/ss_res_sub
}

pR2_singles <- sapply(c("log_budget", "runtime", "release_month", "release_year"), cal_partial_R2)
pR2_singles

# genre #
genre_all <- colnames(x)[5:23]
pR2_genre <- cal_partial_R2(genre_all)
pR2_genre

# cast #
cast_all <- colnames(x)[24:350]
pR2_cast <- cal_partial_R2(cast_all)
pR2_cast

## Barplot of Partial R^2 for Observed Covariates ##
pR2_df <- tibble(pR2 = c(pR2_singles, pR2_genre, pR2_cast),
                 factor = c(names(pR2_singles), "genre", "cast")) %>%
  arrange(pR2)
pR2_df$factor <- factor(pR2_df$factor, levels = pR2_df$factor)

# Barplot for all factors ##
pR2_plot_all <- pR2_df %>%
  ggplot() +
  geom_bar(aes(x = factor, y = pR2), stat = "identity") +
  # geom_hline(yintercept = 0.015, linetype = "dashed", color = "red") +
  theme_bw(base_size = 20) +
  labs(y = expression("Partial "~R^2), x = "",
       title = expression("Partial "~R^2~" for Observed Factors"))+
  theme(plot.title = element_text(hjust = 0.5))
print(pR2_plot_all)


# Barplot for factors with cast excluded ##
pR2_plot_full <- pR2_df %>% filter(factor != "cast") %>%
  ggplot() +
  geom_bar(aes(x = factor, y = pR2), stat = "identity") +
  # geom_hline(yintercept = 0.015, linetype = "dashed", color = "red") +
  theme_bw(base_size = 24) +
  scale_y_continuous(breaks = seq(0, 0.3, by=0.1),
                     labels = sapply(seq(0, 0.3, by=0.1)*100, paste0, "%")) +
  labs(y = expression("Partial "~R^2), x = "",
       title = expression("Partial "~R^2~" for Observed Factors"))+
  theme(plot.title = element_text(hjust = 0.5))
print(pR2_plot_full)
ggsave("movie_pR2_full_lm.pdf", plot = pR2_plot_full, width = 300, height = 150, units = "mm",
       path = "movie_analysis/Figures")


# Barplot for factors with cast, log_budget excluded ##
pR2_plot <- pR2_df %>% filter(factor != "log_budget" & factor != "cast") %>%
  ggplot() +
  geom_bar(aes(x = factor, y = pR2), stat = "identity") +
  # geom_hline(yintercept = 0.015, linetype = "dashed", color = "red") +
  theme_bw(base_size = 20) +
  scale_y_continuous(breaks = seq(0, 0.03, by=0.005),
                     labels = sapply(c('0.0','0.5','1.0','1.5','2.0','2.5','3.0'), paste0, "%")) +
  labs(y = expression("Partial "~R^2), x = "",
       title = expression("Partial "~R^2~" for Observed Factors"))+
  theme(plot.title = element_text(hjust = 0.5))
print(pR2_plot)
# ggsave("movie_pR2_lm.pdf", plot = pR2_plot, width = 180, height = 120, units = "mm",
#        path = "/movie_analysis/Figures")

## Partial R^2 based on model: y ~ t, R^2_{Y ~ T_j | T_{-j}} -----------------------------------------------------------------
lmfit_castfull <- lm(y ~ t)
yhat_castfull <- predict(lmfit_castfull)
ss_res_castfull <- sum((y - yhat_castfull)^2)

cal_partial_R2_cast <- function(var) {
  cat(var, "\n")
  t_sub <- t %>% as_tibble() %>% select(-one_of(var)) %>% as.matrix
  fit_sub <- lm(y ~ t_sub)
  yhat_sub <- predict(fit_sub)
  ss_res_sub <- sum((y - yhat_sub)^2)
  (ss_res_sub - ss_res_castfull)/ss_res_sub
}

pR2_tj <- sapply(colnames(t), cal_partial_R2_cast)
pR2_tj %>% sort(decreasing = TRUE)

