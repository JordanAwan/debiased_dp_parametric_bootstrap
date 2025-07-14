
dir.create(file.path('.', 'estimates'), showWarnings = FALSE, recursive = TRUE)

# experiment settings ####
nSIM <- 10000 # num of simulation for computing coverage
B_paramboot <- 200 # number of paramboot samples
R_indirect_est <- 50
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



n <- 100  # sample size
true_normal_mean <- 1
true_normal_sd <- 1
gdp_mu <- 1
confidence_level_list <- c(0.95)
seed_num <- 0

####################################
check_coverage <- function(true_normal_mean, true_normal_sd, 
                           gdp_mu, xmin, xmax, n, B_paramboot,
                           R_indirect_est, seed_num=0) {
  
  upper_clamp <- xmax
  lower_clamp <- xmin
  
  set.seed(seed_num)
  
  # privacy noise
  sensitivity <- (upper_clamp - lower_clamp)/n
  sd_of_noise_mean <- sensitivity / gdp_mu
  sensitivity_var <- (upper_clamp - lower_clamp)^2/n
  sd_of_noise_var <- sensitivity_var / gdp_mu
  
  clean_clamp_meanvar <- function(x) {
    clamp_x <- pmax(lower_clamp, pmin(upper_clamp, x))
    return(c(mean(clamp_x), var(clamp_x)))
  }
  
  sdp_vec <- function(data_randomness, privacy_noises, sa, mu) {
    # data_randomness should be (R_indirect_est x n)
    # privacy_noises should be (R_indirect_est x 2)
    n <- dim(data_randomness)[2]
    data <- sa * data_randomness + mu
    return(t(apply(data, 1, clean_clamp_meanvar)) + privacy_noises)
  }
  
  score <- function(optim_par, data_randomness, privacy_noises, dp_statistic) {
    mu <- optim_par[1]
    sa <- optim_par[2]
    synth <- sdp_vec(data_randomness, privacy_noises, sa, mu)
    synth <- rbind(synth, dp_statistic) # both are (noisy_clamp_mean, noisy_clamp_var)
    D_synth <- depth.Mahalanobis(synth, synth)
    return(-D_synth[R_indirect_est + 1]) #(rank-(R+1)/2)^2)
  }
  
  solve_meanvar_from_clamp <- function(clamped_meanvar) {
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
  # solve_meanvar_from_clamp(noisy_means_vars[1,])
  
  start_time <- Sys.time()
  set.seed(42 + seed_num)
  
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
  
  ##### simulations #####
  method_result <- foreach(
    i = 1:nSIM,
    .combine = 'cbind',
    .packages = c("ddalpha"),
    .options.snow=opts
  ) %dopar% {
    set.seed(i + seed_num)
    ADI_meansd <- solve_meanvar_from_clamp(noisy_means_vars[i,])    
  }
  all_ADI_meansd <- t(method_result[c(1,2),])
  estimate_lists <- data.frame(cbind(all_ADI_meansd,
                                    rep(xmax, each=nSIM),
                                    rep(gdp_mu, each=nSIM),
                                    rep(R_indirect_est, each=nSIM)))
  colnames(estimate_lists) <- c('mean', 'std', 'clip', 'gdp', 'R')
  return(estimate_lists)
}

true_normal_mean <- 1 # true_normal_mean_list[1]
true_normal_sd <- 1 # true_normal_sd_list[1]
gdp_mu <- 1 # gdp_mu_list[1]
R_indirect_est <- 50

xmax <- 3
gdp_mu_list <- c(0.1,0.3,1,3,10)
seed_num <- 0
for (gdp_mu in gdp_mu_list) {
  estimate_lists <- check_coverage(true_normal_mean, true_normal_sd,
                    gdp_mu, xmin, xmax, n,
                    B_paramboot, R_indirect_est=R_indirect_est, 
                    seed_num=seed_num)
  txt_name <- paste("estimates/paramboot_normal_meansd", "-clamp_", xmin, "_", xmax, 
                    "-gdp_mu=", gdp_mu, 
                    "-N=", n, "-nSIM=", nSIM, 
                    "-R=", R_indirect_est, 
                    "-B=", B_paramboot, "-seed=", seed_num,
                    ".adaptive_indirect.estimates.csv", sep='')
  write.table(estimate_lists, txt_name, sep=",", row.names=FALSE, col.names=TRUE)
}
gdp_mu <- 1


xmax_list <- c(0.1, 0.5, 1, 3, 5)
seed_num <- 0
for (xmax in xmax_list) {
  estimate_lists <- check_coverage(true_normal_mean, true_normal_sd,
                    gdp_mu, xmin, xmax, n,
                    B_paramboot, R_indirect_est=R_indirect_est, 
                    seed_num=seed_num)
  txt_name <- paste("estimates/paramboot_normal_meansd", "-clamp_", xmin, "_", xmax, 
                    "-gdp_mu=", gdp_mu, 
                    "-N=", n, "-nSIM=", nSIM, 
                    "-R=", R_indirect_est, 
                    "-B=", B_paramboot, "-seed=", seed_num,
                    ".adaptive_indirect.estimates.csv", sep='')
  write.table(estimate_lists, txt_name, sep=",", row.names=FALSE, col.names=TRUE)
}
xmax <- 3


R_indirect_est_list <- c(10,20,100,200)
seed_num <- 0
for (R_indirect_est in R_indirect_est_list) {
  estimate_lists <- check_coverage(true_normal_mean, true_normal_sd,
                    gdp_mu, xmin, xmax, n,
                    B_paramboot, R_indirect_est=R_indirect_est, 
                    seed_num=seed_num)
  txt_name <- paste("estimates/paramboot_normal_meansd", "-clamp_", xmin, "_", xmax, 
                    "-gdp_mu=", gdp_mu, 
                    "-N=", n, "-nSIM=", nSIM, 
                    "-R=", R_indirect_est, 
                    "-B=", B_paramboot, "-seed=", seed_num,
                    ".adaptive_indirect.estimates.csv", sep='')
  write.table(estimate_lists, txt_name, sep=",", row.names=FALSE, col.names=TRUE)
}
R_indirect_est <- 50

stopCluster(cl)

