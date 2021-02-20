library(tidyverse)

fpr_naive <- read.csv('simulation/Sparse_Effects_Setting/fpr_naive.csv')$x
tpr_naive <- read.csv('simulation/Sparse_Effects_Setting/tpr_naive.csv')$x
fpr_cali_dim2 <- read.csv("simulation/Sparse_Effects_Setting/LatentDim2/fpr_cali_dim2_L1.csv")$x
tpr_cali_dim2 <- read.csv("simulation/Sparse_Effects_Setting/LatentDim2/tpr_cali_dim2_L1.csv")$x
fpr_cali_dim3 <- read.csv("simulation/Sparse_Effects_Setting/LatentDim3/fpr_cali_dim3_L1.csv")$x
tpr_cali_dim3 <- read.csv("simulation/Sparse_Effects_Setting/LatentDim3/tpr_cali_dim3_L1.csv")$x
fpr_cali_dim4 <- read.csv("simulation/Sparse_Effects_Setting/LatentDim4/fpr_cali_dim4_L1.csv")$x
tpr_cali_dim4 <- read.csv("simulation/Sparse_Effects_Setting/LatentDim4/tpr_cali_dim4_L1.csv")$x

plot_roc <- tibble(TPR = c(tpr_naive, tpr_cali_dim2, tpr_cali_dim3, tpr_cali_dim4),
                   FPR = c(fpr_naive, fpr_cali_dim2, fpr_cali_dim3, fpr_cali_dim4),
                   Type = rep(c('uncali', 'cali_dim2', 'cali_dim3', 'cali_dim4'), each = 500)) %>%
  ggplot() + 
  geom_line(aes(x = FPR, y = TPR, colour = Type)) +
  scale_colour_manual(name = "",
                      values = colorspace::divergingx_hcl(7,palette = "Fall")[c(2,7,6,1)],
                      # values = colorspace::divergingx_hcl(4,palette = "PuOr")[c(2,4,3,1)],
                      labels = c(expression(bolditalic(R^2)~" = 1 (latent dim. = 2)"),
                                 expression(bolditalic(R^2)~" = 1 (latent dim. = 3)"),
                                 expression(bolditalic(R^2)~" = 1 (latent dim. = 4)"),
                                 expression(bolditalic(R^2)~" = 0"))) +
  labs(x = "FPR", y = "TPR", title = 'ROC Curves') + 
  theme_bw(base_size = 12) + 
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.y= element_text(size = 12),
        axis.text.x= element_text(size = 12),
        legend.title = element_text(size=12),
        legend.text = element_text(size=12),
        legend.position = "bottom",
        legend.direction = "vertical",
        legend.text.align = 0)
plot_roc

ggsave("simulation/Sparse_Effects_Setting/plot_roc_curve_L1.pdf", 
       plot = plot_roc, width = 87, height = 115, units = "mm")


DescTools::AUC(x = fpr_cali_dim2, y = tpr_cali_dim2) # 0.61
DescTools::AUC(x = fpr_cali_dim3, y = tpr_cali_dim3) # 0.73
DescTools::AUC(x = fpr_cali_dim4, y = tpr_cali_dim4) # 0.64









