# experiment settings ####
nSIM <- 1000 # num of simulation for computing coverage
B_paramboot <- 200 # number of paramboot samples
R_indirect_est <- 50

n <- 100  # sample size
true_normal_mean <- 1
true_normal_sd <- 1
gdp_mu <- 1
confidence_level_list <- c(0.95)
seed_num <- 0

xmin <- 0
xmax <- 3

result_digits <- 6

# automatic install of packages if they are not installed already ####
list.of.packages <- c(
  "foreach",
  "doSNOW",
  "ddalpha"
)

new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]

if(length(new.packages) > 0){
  install.packages(new.packages, dep=TRUE)
}

#loading packages
for(package.i in list.of.packages){
  suppressPackageStartupMessages(
    library(
      package.i,
      character.only = TRUE
    )
  )
}

n.cores <- min(124, parallel::detectCores() - 1)
cl <- makeSOCKcluster(n.cores)
registerDoSNOW(cl)

# show progress bar
pb <- txtProgressBar(max = nSIM, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress=progress)

####################################
check_coverage <- function(true_normal_mean, true_normal_sd, 
                           gdp_mu, confidence_level_list, xmin, xmax, n, B_paramboot,
                           R_indirect_est, seed_num=0) {
  
  upper_clamp <- xmax
  lower_clamp <- xmin
  
  set.seed(seed_num)
  
  # privacy noise
  sensitivity <- (upper_clamp - lower_clamp)/n
  sd_of_noise_mean <- sensitivity / gdp_mu
  sensitivity_var <- (upper_clamp - lower_clamp)^2/n
  sd_of_noise_var <- sensitivity_var / gdp_mu
  
  # the clamped statistic (private statistic before adding noise)
  clean_clamp_meanvar <- function(x) {
    clamp_x <- pmax(lower_clamp, pmin(upper_clamp, x))
    return(c(mean(clamp_x), var(clamp_x)))
  }
  
  # compute the private statistic 
  sdp_vec <- function(data_randomness, privacy_noises, sa, mu) {
    # data_randomness should be (R_indirect_est x n)
    # privacy_noises should be (R_indirect_est x 2)
    n <- dim(data_randomness)[2]
    data <- sa * data_randomness + mu
    return(t(apply(data, 1, clean_clamp_meanvar)) + privacy_noises)
  }
  
  # objective function to minimize
  score <- function(optim_par, data_randomness, privacy_noises, dp_statistic) {
    mu <- optim_par[1]
    sa <- optim_par[2]
    synth <- sdp_vec(data_randomness, privacy_noises, sa, mu)
    D_synth <- depth.Mahalanobis(dp_statistic, synth)
    return(-D_synth) 
  }
  
  # indirect inference (for getting unbiased estimate)
  solve_meanstd_from_clamp <- function(clamped_meanvar) {
    initialized_value <- clamped_meanvar
    initialized_value[2] <- sqrt(max(1e-12, initialized_value[2]))
    
    data_randomness <- matrix(rnorm(n * R_indirect_est), ncol = n, nrow = R_indirect_est)
    privacy_noises <- matrix(rnorm(2 * R_indirect_est), ncol = 2, nrow = R_indirect_est)
    privacy_noises <- t(t(privacy_noises) * c(sd_of_noise_mean, sd_of_noise_var))
    
    opt = optim(
      par = initialized_value,
      fn = score,
      lower = c(-2, 1e-6),
      upper = c(10, 10),
      method = "L-BFGS-B",
      data_randomness = data_randomness,
      privacy_noises = privacy_noises,
      dp_statistic = clamped_meanvar,
    )
    return(opt$par)
  }
  
  start_time <- Sys.time()
  set.seed(42 + seed_num)
  
  # precomputed nSIM statistics for nSIM times of simulation
  clean_means_vars <- foreach(
    i = 1:nSIM,
    .combine = 'rbind',
    .options.snow=opts
  ) %dopar% {
    set.seed(i + seed_num)
    result <- clean_clamp_meanvar(rnorm(n = n, mean=true_normal_mean, 
                                        sd=true_normal_sd))
  }
  noisy_means_vars <- 
    clean_means_vars + t(t(matrix(rnorm(2 * nSIM), ncol = 2, nrow = nSIM)) * c(sd_of_noise_mean, sd_of_noise_var))
  
  ##### compute adaptive inference estimate & CI #####
  method_result <- foreach(
    i = 1:nSIM,
    .combine = 'cbind',
    .packages = c("ddalpha"),
    .options.snow=opts
  ) %dopar% {
    set.seed(i + seed_num)
    # solve unbiased estimate from precomputed one of nSIM statistics (satisfying DP)
    ADI_meansd <- solve_meanstd_from_clamp(noisy_means_vars[i,])
    
    ##### parametric bootstrap #####
    clean_means_vars_new <- sapply(seq_len(B_paramboot), 
                                   function(x) clean_clamp_meanvar(rnorm(n = n, mean=ADI_meansd[1], sd=ADI_meansd[2])))
    
    privacy_noises <- matrix(rnorm(2 * B_paramboot), ncol = 2, nrow = B_paramboot)
    privacy_noises <- t(t(privacy_noises) * c(sd_of_noise_mean, sd_of_noise_var))
    
    noisy_means_vars_new <- t(clean_means_vars_new) + privacy_noises
    
    # solve PB unbiased estimates from PB statistics (computed following the same DP procedure)
    ADI_meansd_pb <- sapply(seq_len(B_paramboot),
                                    function(i) solve_meanstd_from_clamp(noisy_means_vars_new[i,]))
    ADI_meansd_pb <- t(ADI_meansd_pb)
    
    # build CI from PB unbiased estimates
    current_result <- c()
    for (confidence_level in confidence_level_list) {
      CI_mean_ends <- quantile(2*ADI_meansd[1] - ADI_meansd_pb[,1],
                               probs=c((1-confidence_level)/2, (1+confidence_level)/2),
                               names=FALSE)
      CI_std_ends <- quantile(pmax(0, 2*ADI_meansd[2] - ADI_meansd_pb[,2]),
                              probs=c((1-confidence_level)/2, (1+confidence_level)/2),
                              names=FALSE)
      current_result <- rbind(current_result, 
                              CI_mean_ends[1], CI_mean_ends[2], 
                              CI_std_ends[1], CI_std_ends[2])
    }
    current_result <- rbind(current_result, ADI_meansd[1], ADI_meansd[2])
    current_result
  }
  all_ADI_meansd <- t(method_result[c(5,6),]) # estimate
  estimate_lists <- data.frame(cbind(all_ADI_meansd,
                                    rep('adaptive indirect', each=nSIM)))
  colnames(estimate_lists) <- c('mean', 'std', 'method')
  write.table(estimate_lists, "results/paramboot_normal_comparison_adaptive_indirect_estimates.csv", sep=",", row.names=TRUE, col.names=TRUE)
  
  method_result <- method_result[c(1:4),] # CI 
  CI_infos <- c()
  for (i in c(1:length(confidence_level_list))) {
    covered <- (method_result[4*i-3,] <= true_normal_mean) & 
      (method_result[4*i-2,] >= true_normal_mean)
    coverage_mean <- mean(covered)
    std_coverage_mean <- round(sd(covered) / sqrt(nSIM), digits = result_digits)
    
    covered <- (method_result[4*i-1,] <= true_normal_sd^2) & 
      (method_result[4*i,] >= true_normal_sd^2)
    coverage_sd <- mean(covered)
    std_coverage_sd <- round(sd(covered) / sqrt(nSIM), digits = result_digits)
    
    width <- c(method_result[4*i-2,] - method_result[4*i-3,])
    mean_width_mean <- round(mean(width), digits = result_digits)
    std_width_mean <- round(sd(width) / sqrt(nSIM), digits = result_digits)
    
    width <- c(method_result[4*i,] - method_result[4*i-1,])
    mean_width_sd <- round(mean(width), digits = result_digits)
    std_width_sd <- round(sd(width) / sqrt(nSIM), digits = result_digits)
    CI_infos <- c(CI_infos, coverage_mean, std_coverage_mean, coverage_sd, std_coverage_sd, mean_width_mean, std_width_mean, mean_width_sd, std_width_sd)
  }
  method_stat <- data.frame(matrix(CI_infos, nrow=1))
  method_stat$time <- round(Sys.time() - start_time, 2)
  print(method_stat)
  return(method_stat)
}

method_stat <- check_coverage(true_normal_mean, true_normal_sd,
                              gdp_mu, confidence_level_list, xmin, xmax, n,
                              B_paramboot, R_indirect_est=R_indirect_est, 
                              seed_num=seed_num)
print(method_stat)
txt_name <- paste("results/paramboot_normal_meansd", "-clamp_", xmin, "_", xmax, 
                  "-gdp_mu=", gdp_mu, "-conf=", confidence_level_list[1],
                  "-N=", n, "-nSIM=", nSIM, 
                  "-R=", R_indirect_est, 
                  "-B=", B_paramboot, "-seed=", seed_num,
                  ".adaptive_indirect.csv", sep='')
write.table(method_stat, txt_name, sep=",", row.names=TRUE, col.names=TRUE)



stopCluster(cl)

