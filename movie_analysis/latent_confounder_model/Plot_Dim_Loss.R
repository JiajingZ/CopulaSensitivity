library(tidyverse)
setwd("/Users/jiajing/Desktop/copula-sensitivity/tmdb-movie-metadata/LatentConfounderAnalysis/Cast_VAE_Python/NumLatentVarCheck/LabelSmooth")
dim_loss_df <- read.csv('dim_loss_df.csv')
dim(dim_loss_df)
dim <- c(seq(1,7,by = 2),seq(10, 40, by=5))
loss <- rowMeans(dim_loss_df)

tibble(dim = c(seq(1,7,by = 2),seq(10, 40, by=5)),
       loss = rowMeans(dim_loss_df)) %>%
  ggplot(aes(x = dim, y = loss)) + 
  geom_line() + 
  labs(x = '# latent Dimension', y = "test loss")


plot_dim_loss <- tibble(dim = c(seq(1,7,by = 2),seq(10, 40, by=5)),
       loss = rowMeans(dim_loss_df)) %>%
  ggplot(aes(x = dim, y = loss)) + 
  geom_smooth(se = FALSE, method = "loess", color = "black") + 
  labs(x = '# latent dimension', y = "validation loss",
       title = "Latent Dimension Selection") + 
  theme_bw(base_size = 15) +
  theme(plot.title = element_text(hjust = 0.5))

plot_dim_loss

ggsave("plot_dim_loss.pdf", plot = plot_dim_loss,
       width = 127, height = 100, units = "mm")
