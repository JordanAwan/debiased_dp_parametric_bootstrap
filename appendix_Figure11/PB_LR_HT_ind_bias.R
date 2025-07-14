Delta = 2
ep = 1
mu = .5
tau = 1
# X ~ N(mu, tau^2)
beta = c(-.5,1) # first entry is intercept, second is slope.
sa = .5
# Y ~ N(X * beta[2] + beta[1], sa^2)
seed_num <- 0
theta = c(beta[2],beta[1],mu,tau,sa)


alpha=.05
R=200 # number of paramboot samples
reps=10000 # num of simulation for computing coverage

clamp_list <- c(0.5, 0.8, 1.0, 1.5, 2.0, 5.0, 10.0)
n_list <- c(100,200,300,400,500,1000)
n_list <- c(100,200,300,500,1000,2000)
n_list <- c(2000,1000,500,300,200,100)

ncol <- length(clamp_list)
nrow <- length(n_list)
result_table1 <- matrix(nrow = nrow, ncol = ncol)
rownames(result_table1) <- n_list
colnames(result_table1) <- clamp_list

result_table2 <- result_table1
result_table3 <- result_table1
result_table0 <- result_table1


estimate_lists <- c()
for (j in c(1:ncol)){
  for (i in c(1:nrow)) {
    n <- n_list[i]
    Delta <- clamp_list[j]
    txt_name <- paste("./results/PB_LR_HT/estden/PB_LR_HT_ada_o4", "-clamp_", Delta, "-n_", n,
                      "-gdp_ep=", ep, "-alpha=", alpha,
                      "-X_mu=", mu, "-tau=", tau,
                      "-beta=(", beta[1], ", b1)",
                      "-sa=", sa,
                      "-R=", R, "-reps=", reps, "_checkclamp_estden.csv", sep='')
    est_results <- read.csv(txt_name)
    result_table0[i,j] <- mean(est_results$V1 == 0)
    result_table1[i,j] <- mean(est_results$V1 == 1)
    result_table2[i,j] <- mean(est_results$V1 == 2)
    result_table3[i,j] <- mean(est_results$V1 == 3)
    
    estimate_list <- c()
    estimate_list_part <- data.frame(value=est_results$V3)
    estimate_list_part$param <- c('beta0')
    estimate_list <- rbind(estimate_list, estimate_list_part)
    estimate_list_part <- data.frame(value=est_results$V4)
    estimate_list_part$param <- c('E[X]')
    estimate_list <- rbind(estimate_list, estimate_list_part)
    estimate_list_part <- data.frame(value=est_results$V5)
    estimate_list_part$param <- c('SD(X)')
    estimate_list <- rbind(estimate_list, estimate_list_part)
    estimate_list_part <- data.frame(value=est_results$V6)
    estimate_list_part$param <- c('SD(e)')
    estimate_list <- rbind(estimate_list, estimate_list_part)
    estimate_list$n <- n
    estimate_list$Delta <- paste("Delta=",Delta, sep="")
    estimate_lists <- rbind(estimate_lists, estimate_list)
  }
}


library(ggplot2)
library(ggridges)
estimate_lists_new <- estimate_lists
estimate_lists_new$param <- factor(estimate_lists$param, levels=c('beta0', 'E[X]', 'SD(X)', 'SD(e)'))
estimate_lists_new$n <- factor(estimate_lists$n, levels=n_list)
estimate_lists_new$Delta <- factor(estimate_lists$Delta, levels=c("Delta=0.5","Delta=0.8","Delta=1","Delta=1.5","Delta=2","Delta=5","Delta=10"))
alpha_val <- 1
num_sample_size <- length(unique(estimate_lists_new$n))
override.alpha <- rep(alpha_val, num_sample_size)
override.linetype <- c(1:num_sample_size)
ggplot(estimate_lists_new) +
  # stat_density(aes(x = value, linetype = n, color = n), geom="line", position="identity") +
  # scale_linetype_manual(values=c('solid','dashed','dotted','dotdash','longdash', 'twodash')[1:num_sample_size])+
  # # xlim(c(-1.5, 2.5)) +
  geom_density_ridges(aes(x = value, y = n, fill = n), alpha=0.5, quantile_lines = TRUE, quantiles = c(0.5)) + #, scale=1. , panel_scaling=FALSE
  facet_grid(Delta ~ param, scales="free") +
  theme_bw(base_size = 25) +
  labs(x = "Estimate Value", y = "n", color = "") +
  theme(legend.position = "none") 
  # theme(legend.position = "top", panel.spacing = unit(0.7, "lines")) +
  # # scale_color_discrete(labels=c('adaptive indirect', 'naïve', 'Ferrando et al. (2022)')) +
  # guides(color=guide_legend(title="sample size n", nrow=1, 
  #                           override.aes = list(alpha = override.alpha, linetype = override.linetype)), 
  #        linetype="none", alpha="none") 

ggsave(paste("bias_comparison_HT_LR_null.pdf", sep=""), width=17, height=22)




result_table0 <- t(result_table0)
result_table1 <- t(result_table1)
result_table2 <- t(result_table2)
result_table3 <- t(result_table3)
result_table0 # reject
result_table1 # accept as T is NaN
result_table2 # accept as too many NaN are in synthetic_T
result_table3 # accept as T is not large enough

