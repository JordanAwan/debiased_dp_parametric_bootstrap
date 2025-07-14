library(ggplot2)
library(ggridges)

estimate_lists <- read.csv("results/paramboot_normal_comparison_estimates.csv", sep=",")
estimate_lists2 <- read.csv("results/paramboot_normal_comparison_adaptive_indirect_estimates.csv", sep=",")
estimate_lists2$std <- sqrt(estimate_lists2$std)

estimate_lists <- rbind(estimate_lists, estimate_lists2)
estimate_lists$method <- factor(estimate_lists$method, levels=c("adaptive indirect", "simplified t", 
                                                                "naive", "Efron's BC", "automatic percentile"))
estimate_lists_new <- data.frame(estimate=c(estimate_lists$mean, estimate_lists$std), 
                                 method=c(estimate_lists$method, estimate_lists$method),
                                 parameter=c(rep(c('estimated mean', 'estimated standard deviation'), each=length(estimate_lists$std))))

estimate_lists_new2 <- estimate_lists_new[estimate_lists_new$method != "Efron's BC",]
estimate_lists_new <- estimate_lists_new2[estimate_lists_new2$method != "automatic percentile",]

method1 <- estimate_lists_new$method
method1 <- factor(method1, level=c("adaptive indirect", "simplified t", "naive", "naïve"))
method1 <- replace(method1, method1=="naive", 'naïve')
estimate_lists_new$method <- method1

ggplot(estimate_lists_new) +
  geom_density_ridges(aes(x = estimate, y = method, fill = method), alpha=0.9, quantile_lines = TRUE, quantiles = c(0.5)) + #, scale=1. , panel_scaling=FALSE
  facet_grid(. ~ parameter, scales="free") +
  geom_text(label = "*", x = 1, y = 0.4, size = 8, color = "#C28E0E") +
  theme_bw() +
  labs(x = "Value", y = "Density", color = "Method") +
  theme(legend.position = "none") 
ggsave(paste("results/bias_comparison.pdf", sep=""), width=7, height=1.5)
