library(ggplot2)
library(ggridges)

###############################################################################################


R_indirect_est <- 50
xmin <- 0
xmax <- 3
gdp_mu <- 1 
confidence_level <- 0.95 
n <- 100
nSIM <- 10000
B <- 200
gdp_mu_list <- c(0.1,0.3,1,3,10)
seed_num <- 0
estimate_lists <- c()
for (gdp_mu in gdp_mu_list) {
  txt_name <- paste("estimates/paramboot_normal_meansd", "-clamp_", xmin, "_", xmax, 
                    "-gdp_mu=", gdp_mu, 
                    "-N=", n, "-nSIM=", nSIM, 
                    "-R=", R_indirect_est, 
                    "-B=", B, "-seed=", seed_num,
                    ".adaptive_indirect.estimates.csv", sep='')
  estimate_list <- read.csv(txt_name, sep=",")
  estimate_lists <- rbind(estimate_lists, estimate_list)
}
estimate_lists <- rbind(estimate_lists, estimate_lists)
estimate_lists <- data.frame(estimate_lists)
total_len <- length(estimate_lists[,1])
estimate_lists$estimate <- c(estimate_lists$mean[1:(total_len/2)],
                             estimate_lists$std[(total_len/2+1):total_len])
estimate_lists$parameter <- rep(c('estimated mean', 'estimated standard deviation'), each=total_len/2)
estimate_lists$gdp <- factor(estimate_lists$gdp)

ggplot(estimate_lists) +
  geom_density_ridges(aes(x = estimate, y = gdp, fill = gdp), alpha=0.9, quantile_lines = TRUE, quantiles = c(0.5)) + #, scale=1. , panel_scaling=FALSE
  facet_grid(. ~ parameter) +
  geom_text(label = "*", x = 1, y = 0.4, size = 8, color = "#C28E0E") +
  theme_bw() +
  labs(x = "Value", y = "GDP (epsilon)", color = "Method") +
  theme(legend.position = "none") 
ggsave(paste("bias_comparison_gdp.pdf", sep=""), width=8, height=3)




R_indirect_est <- 50
gdp_mu <- 1
xmax_list <- c(0.1, 0.5, 1, 3, 5)
seed_num <- 0
estimate_lists <- c()
for (xmax in xmax_list) {
  txt_name <- paste("estimates/paramboot_normal_meansd", "-clamp_", xmin, "_", xmax, 
                    "-gdp_mu=", gdp_mu, 
                    "-N=", n, "-nSIM=", nSIM, 
                    "-R=", R_indirect_est, 
                    "-B=", B, "-seed=", seed_num,
                    ".adaptive_indirect.estimates.csv", sep='')
  estimate_list <- read.csv(txt_name, sep=",")
  estimate_lists <- rbind(estimate_lists, estimate_list)
}
estimate_lists <- rbind(estimate_lists, estimate_lists)
estimate_lists <- data.frame(estimate_lists)
total_len <- length(estimate_lists[,1])
estimate_lists$estimate <- c(estimate_lists$mean[1:(total_len/2)],
                             estimate_lists$std[(total_len/2+1):total_len])
estimate_lists$parameter <- rep(c('estimated mean', 'estimated standard deviation'), each=total_len/2)
estimate_lists$clip <- factor(estimate_lists$clip)

ggplot(estimate_lists) +
  geom_density_ridges(aes(x = estimate, y = clip, fill = clip), alpha=0.9, quantile_lines = TRUE, quantiles = c(0.5)) + #, scale=1. , panel_scaling=FALSE
  facet_grid(. ~ parameter) +
  geom_text(label = "*", x = 1, y = 0.4, size = 8, color = "#C28E0E") +
  theme_bw() +
  labs(x = "Value", y = "Clamping (upper bound)", color = "Method") +
  theme(legend.position = "none") 
ggsave(paste("bias_comparison_clamp.pdf", sep=""), width=8, height=3)




xmax <- 3
R_indirect_est_list <- c(10,20,50,100,200)
seed_num <- 0
estimate_lists <- c()
for (R_indirect_est in R_indirect_est_list) {
  txt_name <- paste("estimates/paramboot_normal_meansd", "-clamp_", xmin, "_", xmax, 
                    "-gdp_mu=", gdp_mu, 
                    "-N=", n, "-nSIM=", nSIM, 
                    "-R=", R_indirect_est, 
                    "-B=", B, "-seed=", seed_num,
                    ".adaptive_indirect.estimates.csv", sep='')
  estimate_list <- read.csv(txt_name, sep=",")
  print(c(R_indirect_est, sd(estimate_list[,1]), sd(estimate_list[,2])))
  estimate_lists <- rbind(estimate_lists, estimate_list)
}
estimate_lists <- rbind(estimate_lists, estimate_lists)
estimate_lists <- data.frame(estimate_lists)
total_len <- length(estimate_lists[,1])
estimate_lists$estimate <- c(estimate_lists$mean[1:(total_len/2)],
                             estimate_lists$std[(total_len/2+1):total_len])
estimate_lists$parameter <- rep(c('estimated mean', 'estimated standard deviation'), each=total_len/2)
estimate_lists$R <- factor(estimate_lists$R)

ggplot(estimate_lists) +
  geom_density_ridges(aes(x = estimate, y = R, fill = R), alpha=0.9, quantile_lines = TRUE, quantiles = c(0.5)) + #, scale=1. , panel_scaling=FALSE
  facet_grid(. ~ parameter) +
  geom_text(label = "*", x = 1, y = 0.4, size = 8, color = "#C28E0E") +
  theme_bw() +
  labs(x = "Value", y = "Indirect samples (R)", color = "Method") +
  theme(legend.position = "none") 
ggsave(paste("bias_comparison_indirectR.pdf", sep=""), width=8, height=3)



##################################################################################################################


R_indirect_est <- 50
gdp_mu <- 1
xmax_list <- c(0.1, 0.5, 1, 3, 5)
seed_num <- 0
estimate_lists <- c()
for (xmax in xmax_list) {
  txt_name <- paste("estimates/paramboot_normal_meansd", "-clamp_", xmin, "_", xmax, 
                    "-gdp_mu=", gdp_mu, 
                    "-N=", n, "-nSIM=", nSIM, 
                    "-R=", R_indirect_est, 
                    "-B=", B, "-seed=", seed_num,
                    ".IND.estimates.csv", sep='')
  estimate_list <- read.csv(txt_name, sep=",")
  estimate_lists <- rbind(estimate_lists, estimate_list)
}
estimate_lists <- rbind(estimate_lists, estimate_lists)
estimate_lists <- data.frame(estimate_lists)
total_len <- length(estimate_lists[,1])
estimate_lists$estimate <- c(estimate_lists$mean[1:(total_len/2)],
                             estimate_lists$std[(total_len/2+1):total_len])
estimate_lists$parameter <- rep(c('estimated mean', 'estimated standard deviation'), each=total_len/2)
estimate_lists$clip <- factor(estimate_lists$clip)

ggplot(estimate_lists) +
  geom_density_ridges(aes(x = estimate, y = clip, fill = clip), alpha=0.9, quantile_lines = TRUE, quantiles = c(0.5)) + #, scale=1. , panel_scaling=FALSE
  facet_grid(. ~ parameter) +
  geom_text(label = "*", x = 1, y = 0.4, size = 8, color = "#C28E0E") +
  theme_bw() +
  labs(x = "Value", y = "Clamping (upper bound)", color = "Method") +
  theme(legend.position = "none") 
ggsave(paste("bias_comparison_clamp_IND.pdf", sep=""), width=8, height=3)

