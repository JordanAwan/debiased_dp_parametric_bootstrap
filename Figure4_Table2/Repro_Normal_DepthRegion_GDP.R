# Set CRAN mirror
options(repos = c(CRAN = "https://cloud.r-project.org"))

# --- Global Simulation Settings ---
nSIM <- 1000           # Number of simulation trials (Start with 50 to test timing)
value_r <- 50         # Grid density (50x50 = 2500 cells per trial)
n <- 100              # Sample size
R_synthetic <- 200    # Number of synthetic samples per depth calculation
ep <- 1             # Privacy parameter (epsilon)
alpha <- 0.05         # Significance level (for 95% Confidence Region)
tol <- 10^-8

population_mu <- 1    # True population mean
population_sigma <- 1 # True population SD
upper_clamp <- 3
lower_clamp <- 0

###################################################
# Package Installation and Loading
list.of.packages <- c("foreach", "doSNOW", "ddalpha", "parallelly", "rmutil")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages) > 0) install.packages(new.packages, dep=TRUE)

for(package.i in list.of.packages) library(package.i, character.only = TRUE)

# Parallel Computing Setup
# Use 90% of cores instead of 100%
n.cores <- floor(parallelly::availableCores() * 0.4)
cl <- makeSOCKcluster(n.cores)

# n.cores <- parallelly::availableCores(constraints = "connections")
# cl <- makeSOCKcluster(n.cores)
registerDoSNOW(cl)

# Progress Bar for the SIMULATION trials
pb <- txtProgressBar(max = nSIM, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress=progress)

#################################################
# Helper Functions (Data Generation & Scoring)

sdp_vec <- function(data_randomness, privacy_noises, sa, mu) {
  n_inner <- dim(data_randomness)[2]
  data <- sa * data_randomness + mu
  data_clamp <- pmax(pmin(data, upper_clamp), lower_clamp)
  s1 <- apply(data_clamp, 1, mean) + (upper_clamp - lower_clamp) / (n_inner * ep) * privacy_noises[, 1]
  s2 <- apply(data_clamp, 1, var) + (upper_clamp - lower_clamp)^2 / (n_inner * ep) * privacy_noises[, 2]
  return(cbind(s1, s2))
}

score_mu_sa <- function(optim_par, data_randomness, privacy_noises, dp_statistic, depth_type){
  mu <- optim_par[1]; sa <- optim_par[2]
  synth <- sdp_vec(data_randomness, privacy_noises, sa, mu)
  synth <- rbind(synth, dp_statistic)
  
  D_synth <- switch(depth_type,
    'mahalanobis' = depth.Mahalanobis(synth, synth),
    'halfspace'   = depth.halfspace(synth, synth),
    'simplicial'  = depth.simplicial(synth, synth),
    'spatial'     = depth.spatial(synth, synth)
  )
  
  r = rank(D_synth, ties.method="max")[R_synthetic+1]
  s = r + D_synth[R_synthetic+1]
  return(-s)
}

score_sa_mu <- function(optim_par, data_randomness, privacy_noises, dp_statistic, depth_type){
  # Swaps indices for the bisection search on the second parameter
  return(score_mu_sa(c(optim_par[2], optim_par[1]), data_randomness, privacy_noises, dp_statistic, depth_type))
}

accept <- function(optim_par, data_randomness, privacy_noises, dp_statistic, 
                   search_lower, search_upper, nuisance_lower, nuisance_upper, depth_type, score_func) {
  
  proposed_result <- score_func(optim_par, data_randomness, privacy_noises, dp_statistic, depth_type)
  if ((-proposed_result) >= (floor((alpha)*(R_synthetic+1))+1)) return(optim_par)

  opt <- optim(par = optim_par, fn = score_func, method = "L-BFGS-B",
               lower = c(search_lower, nuisance_lower), upper = c(search_upper, nuisance_upper),
               data_randomness = data_randomness, privacy_noises = privacy_noises,
               dp_statistic = dp_statistic, depth_type = depth_type)
  
  if ((-opt$value) >= (floor((alpha)*(R_synthetic+1))+1)) return(opt$par)
  return(c(NA, NA))
}

getConfidenceInterval <- function(optim_par, dp_statistic, data_randomness, privacy_noises, 
                                  search_lower, search_upper, nuisance_lower, nuisance_upper, depth_type, score_func) {
  optim_par <- accept(optim_par, data_randomness, privacy_noises, dp_statistic, 
                      search_lower, search_upper, nuisance_lower, nuisance_upper, depth_type, score_func)
  if (is.na(optim_par[1])) return(c(NA, NA))
  
  # Bisection Logic
  t_low <- search_lower; t_up <- search_upper; t_mid_val <- optim_par[1]
  
  # Left Boundary
  l_low <- t_low; l_up <- t_mid_val - tol
  while(l_up - l_low > 0.1){ # Larger tolerance for faster grid bounding
    l_mid <- (l_low + l_up) / 2
    if(!is.na(accept(c(l_mid, optim_par[2]), data_randomness, privacy_noises, dp_statistic, 
                     l_low, l_mid, nuisance_lower, nuisance_upper, depth_type, score_func)[1])) l_up <- l_mid - tol else l_low <- l_mid
  }
  
  # Right Boundary
  r_low <- t_mid_val + tol; r_up <- t_up
  while(r_up - r_low > 0.1){
    r_mid <- (r_low + r_up) / 2
    if(!is.na(accept(c(r_mid, optim_par[2]), data_randomness, privacy_noises, dp_statistic, 
                     r_mid, r_up, nuisance_lower, nuisance_upper, depth_type, score_func)[1])) r_low <- r_mid + tol else r_up <- r_mid
  }
  return(c(l_low, r_up))
}

################################################################
# --- Main Simulation Loop ---

depth_type <- 'mahalanobis' # Change to run others
cat("Starting Simulation for:", depth_type, "\n")

sim_results <- foreach(s = 1:nSIM, .combine = 'rbind', .packages = c('ddalpha'), .options.snow = opts) %dopar% {
  
  set.seed(s + 123)
  # 1. Generate Observed Private Data
  raw_data <- rnorm(n, mean = population_mu, sd = population_sigma)
  data_c <- pmax(pmin(raw_data, upper_clamp), lower_clamp)
  obs_s1 <- mean(data_c) + (upper_clamp - lower_clamp) / (n * ep) * rnorm(1)
  obs_s2 <- var(data_c) + (upper_clamp - lower_clamp)^2 / (n * ep) * rnorm(1)
  obs_stats <- c(obs_s1, obs_s2)
  
  # 2. Setup Repro Internal Randomness
  data_rand <- matrix(rnorm(n * R_synthetic), ncol = n, nrow = R_synthetic)
  priv_noise <- matrix(rnorm(R_synthetic * 2), ncol = 2, nrow = R_synthetic)
  
  # 3. Find Search Bounds
  d1_rng <- getConfidenceInterval(c(1,1), obs_stats, data_rand, priv_noise, -2, 5, 0.1, 5, depth_type, score_mu_sa)
  d2_rng <- getConfidenceInterval(c(1,1), obs_stats, data_rand, priv_noise, 0.1, 5, -2, 5, depth_type, score_sa_mu)
  
  if(any(is.na(d1_rng)) || any(is.na(d2_rng))) return(c(0, 0))
  
  # 4. Grid Search for Area and Coverage
  mu_vals <- seq(d1_rng[1], d1_rng[2], length.out = value_r)
  sa_vals <- seq(d2_rng[1], d2_rng[2], length.out = value_r)
  step_mu <- mu_vals[2] - mu_vals[1]
  step_sa <- sa_vals[2] - sa_vals[1]
  
  total_area <- 0
  covered <- 0
  
  for(m_i in 1:(value_r-1)){
    for(s_i in 1:(value_r-1)){
      mid_pt <- c((mu_vals[m_i] + mu_vals[m_i+1])/2, (sa_vals[s_i] + sa_vals[s_i+1])/2)
      # Check acceptance
      if(!is.na(accept(mid_pt, data_rand, priv_noise, obs_stats, mu_vals[m_i], mu_vals[m_i+1], 
                       sa_vals[s_i], sa_vals[s_i+1], depth_type, score_mu_sa)[1])){
        total_area <- total_area + (step_mu * step_sa)
        # Check if truth is in this cell
        if(population_mu >= mu_vals[m_i] && population_mu <= mu_vals[m_i+1] &&
           population_sigma >= sa_vals[s_i] && population_sigma <= sa_vals[s_i+1]) covered <- 1
      }
    }
  }
  c(covered, total_area)
}

# --- Final Output ---
stopCluster(cl)
final_coverage <- mean(sim_results[,1])
final_area <- mean(sim_results[,2])
final_coverage_se <- sd(sim_results[,1]) / sqrt(nSIM)
final_area_se <- sd(sim_results[,2]) / sqrt(nSIM)

cat("\n--- Final Simulation Results ---\n")
cat("Depth Type:      ", depth_type, "\n")
cat("Joint Coverage:  ", round(final_coverage, 4), "(SE:", round(final_coverage_se, 4), ")\n")
cat("Average Area:    ", round(final_area, 4), "(SE:", round(final_area_se, 4), ")\n")
cat("Total Time:      ", round(difftime(Sys.time(), start_time, units="mins"), 2), "mins\n")