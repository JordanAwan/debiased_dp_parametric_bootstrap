nsim_sub <- 200
nSIM_range <- 1000
nSIM_start <- 1
target_value <- 2
# bias_correct <- 1
R_synthetic <- 200

# Initialize empty matrices to store coverage and width values
coverage_matrix <- matrix(0, nrow=5, ncol=5)
width_matrix <- matrix(0, nrow=5, ncol=5)

# Define epsilon and N values in order they'll appear in loop
epsilon_vals <- c(0.1, 0.3, 1.0, 3.0, 10)
# N_vals <- c(100, 500, 1000, 2000, 5000, 10000, 20000, 50000)
N_vals <- c(100, 200, 500, 1000, 2000)

for(epsilon in epsilon_vals){
    for(n in N_vals){
        
    # results/small_problem_compare_CI_results-lapldp_ep=0.1-n=100-R_synthetic=20-nsim_sub=200-nSIM=100_1bias_correct=1.rdata
        path <- paste("results/logistic_naive_compare_CI_results-gdp_ep=", epsilon, "-n=", n, 
                  "-R_synthetic=", R_synthetic, "-nsim_sub=", nsim_sub, "-nSIM=", nSIM_range, "_", nSIM_start, 
                  ".rdata", sep='')
                  
        load(path)
        contains_target <- apply(CI_results, 1, function(row) {
          target_value >= row[1] && target_value <= row[2]
        })
        coverage <- mean(contains_target) 
        cat("coverage for epsilon: ", epsilon, "and N: ", n, "is", coverage, "\n")
        
        width <- CI_results[,2] - CI_results[,1]
        width_mean <- mean(width)
        cat("width for epsilon: ", epsilon, "and N: ", n, "is", width_mean, "\n")

        # Store coverage and width values in respective matrices
        epsilon_idx <- which(epsilon_vals == epsilon)
        N_idx <- which(N_vals == n)   

        coverage_matrix[epsilon_idx, N_idx] <- coverage
        width_matrix[epsilon_idx, N_idx] <- width_mean
    }
}

coverage_matrix
width_matrix

# Save coverage and width matrices
file_path <- paste("results/coverage_width_matrices_objpert_nsim_sub=", nsim_sub, "-nSIM_range=1000", "_", nSIM_start, "-R_synthetic=", R_synthetic, ".rdata", sep='')
# results/coverage_width_matrices_lapldp_nsim_sub=200-nSIM_range=100_1-R_synthetic=20.rdata
save(coverage_matrix, width_matrix, file=file_path)
