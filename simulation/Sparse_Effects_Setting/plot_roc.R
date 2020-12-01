fpr_naive <- read.csv('simulation/Sparse_Effects_Setting/fpr_naive.csv')$x
tpr_naive <- read.csv('simulation/Sparse_Effects_Setting/tpr_naive.csv')$x

fpr_cali_dim2 <- read.csv("simulation/Sparse_Effects_Setting/LatentDim2/fpr_cali_dim2.csv")$x
tpr_cali_dim2 <- read.csv("simulation/Sparse_Effects_Setting/LatentDim2/tpr_cali_dim2.csv")$x
fpr_cali_dim3 <- read.csv("simulation/Sparse_Effects_Setting/LatentDim3/fpr_cali_dim3.csv")$x
tpr_cali_dim3 <- read.csv("simulation/Sparse_Effects_Setting/LatentDim3/tpr_cali_dim3.csv")$x
fpr_cali_dim4 <- read.csv("simulation/Sparse_Effects_Setting/LatentDim4/fpr_cali_dim4.csv")$x
tpr_cali_dim4 <- read.csv("simulation/Sparse_Effects_Setting/LatentDim4/tpr_cali_dim4.csv")$x


plot_roc <- tibble(TPR = c(tpr_naive, tpr_cali_dim2, tpr_cali_dim3, tpr_cali_dim4),
                   FPR = c(fpr_naive, fpr_cali_dim2, fpr_cali_dim3, fpr_cali_dim4),
                   Type = rep(c('uncali', 'cali_dim2', 'cali_dim3', 'cali_dim4'), each = 500)) %>%
  ggplot() + 
  geom_line(aes(x = FPR, y = TPR, colour = Type)) +
  scale_colour_manual(name = "",
                      values = divergingx_hcl(4,palette = "Zissou 1")[c(2,4,3,1)],
                      labels = c('calibrated (latent dim. = 2)',
                                 'calibrated (latent dim. = 3)',
                                 'calibrated (latent dim. = 4)',
                                 'uncalibrated')) +
  labs(x = "FPR", y = "TPR", title = 'ROC Curves') + 
  theme_bw(base_size = 12) + 
  theme(plot.title = element_text(hjust = 0.5),
        legend.title = element_text(size=12),
        legend.text = element_text(size=11),
        legend.position = "bottom",
        legend.direction = "vertical")
plot_roc


ggsave("plot_roc_curve.pdf", plot = plot_roc, width = 90, height = 120, units = "mm")


DescTools::AUC(x = fpr_cali_dim2, y = tpr_cali_dim2)
DescTools::AUC(x = fpr_cali_dim4, y = tpr_cali_dim4)












