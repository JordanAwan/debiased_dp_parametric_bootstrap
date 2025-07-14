### parametric bootstrap
nSIM = 1000
B_paramboot = 200

n = 100
true_normal_mean = 1
true_normal_sd = 1
ep = 1
confidence_level = 0.95
alpha = 1 - confidence_level

seed_num = 1000

xmin = 0
xmax = 3

# automatic install of packages if they are not installed already ####
list.of.packages <- c(
  "foreach",
  "doSNOW"
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
opts <- list()
opts <- list(progress=progress)

###################################
set.seed(seed_num)

clamp = function(x){
  return(pmin(pmax(x, xmin), xmax))
}
generate_X_func = function(u, m, s){
  # u ~ N(0,1), m is mean, s is std
  return(s * u + m)
}

results <- foreach(
    r = 1:nSIM,
    .combine = 'rbind',
    .options.snow=opts
  ) %dopar% {
  set.seed(r + seed_num)
  X_0 = rnorm(n, true_normal_mean, true_normal_sd)
  
  # (naive) private estimate of mean and std
  x_clamp = clamp(X_0)
  z_m = mean(x_clamp) + rnorm(n=1,m=0,s=(b-a)/ep) / n
  z_s = var(x_clamp) + rnorm(n=1,m=0,s=(b-a)^2/ep) / n 
  z_s_sqrt=sqrt(max(z_s, 0))
  estimate_list <- c(z_m, z_s_sqrt)
  
  # from the (naive) private estimate, we conduct parametric bootstrap
  ## fix seeds for data_generation (u), privacy_noise_generation (w1, w2)
  u = matrix(rnorm(B_paramboot * n, m=0, s=1), nrow=B_paramboot, ncol=n)
  w1 = rnorm(n=B_paramboot, m=0, s=(b-a)/ep)
  w2 = rnorm(n=B_paramboot, m=0, s=(b-a)^2/ep)
  m_star = m_star2 = m_star3 = m_star5 = m_starPB = rep(0,B_paramboot)
  s_star = s_star2 = s_star3 = s_star5 = s_starPB = rep(0,B_paramboot)

  for(i in 1:R){
    data = generate_X_func(u[i,], m=z_m, s=z_s_sqrt)
    data_clamp = clamp(data)
    z_star1 = mean(data_clamp) + w1[i] / n
    z_star2 = var(data_clamp) + w2[i] / n

    # naive PB estimates
    m_star[i] = z_star1
    s_star[i] = sqrt(max(z_star2,0))
    
    # simplified-t PB estimates
    m_star2[i]= 2*z_m - z_star1
    s_star2[i] = 2*z_s_sqrt - s_star[i]
  }
  estimate_list2 <- c(mean(m_star2), mean(s_star2))
  
  # Ferrando et al. estimates
  m_star3 = m_star - (mean(m_star) - z_m)
  s_star3 = s_star - (mean(s_star) - z_s_sqrt)
  estimate_list3 <- c(z_m - (mean(m_star) - z_m), 
                            z_s_sqrt - (mean(s_star) - z_s_sqrt))
  
  # CI for mean
  # naive percentile CI
  interval = quantile(m_star,c(alpha/2,1-alpha/2))
  if(interval[1]<=true_normal_mean & interval[2]>=true_normal_mean){
    coverage = 1
  }
  else{coverage = 0}
  width = interval[2]-interval[1]
  
  # simplified-t CI
  interval2 = quantile(m_star2,c(alpha/2,1-alpha/2))
  if(interval2[1]<=true_normal_mean & interval2[2]>=true_normal_mean){
    coverage2 = 1
  }
  else{coverage2 = 0}
  width2 = interval2[2]-interval2[1]
  
  # Ferrando et al. CI
  interval3 = quantile(m_star3,c(alpha/2,1-alpha/2))
  if(interval3[1]<=true_normal_mean & interval3[2]>=true_normal_mean){
    coverage3 = 1
  }
  else{coverage3 = 0}
  width3 = interval3[2]-interval3[1]
  
  # CI for std
  # naive percentile CI
  intervalS = quantile(s_star,c(alpha/2,1-alpha/2))
  if(intervalS[1]<=true_normal_sd & intervalS[2]>=true_normal_sd){
    coverageS = 1
  }
  else{coverageS = 0}
  widthS = intervalS[2]-intervalS[1]
  
  # simplified-t CI
  intervalS2 = quantile(s_star2,c(alpha/2,1-alpha/2))
  if(intervalS2[1]<=true_normal_sd & intervalS2[2]>=true_normal_sd){
    coverageS2 = 1
  }
  else{coverageS2 = 0}
  widthS2 = intervalS2[2]-intervalS2[1]
  
  # Ferrando et al. CI
  intervalS3 = quantile(s_star3,c(alpha/2,1-alpha/2))
  if(intervalS3[1]<=true_normal_sd & intervalS3[2]>=true_normal_sd){
    coverageS3 = 1
  }
  else{coverageS3 = 0}
  widthS3 = intervalS3[2]-intervalS3[1]
  
  # Efron's BC estimates 
  # for mean
  hat_z_0 = qnorm(sum(m_star < z_m) / B_paramboot)
  z_alpha_lower = 2 * hat_z_0 - qnorm(alpha/2)
  z_alpha_upper = 2 * hat_z_0 - qnorm(1-alpha/2)
  m_star_sorted = sort(m_star)
  theta_upper = m_star_sorted[as.integer(R * pnorm(z_alpha_lower))+1]
  theta_lower = m_star_sorted[as.integer(R * pnorm(z_alpha_upper))+1]
  interval4 = c(theta_lower, theta_upper)
  corrected_m <- m_star_sorted[max(1, as.integer(B_paramboot * pnorm(2 * hat_z_0)))]
  
  # for std
  hat_z_0 = qnorm(sum(s_star < z_s_sqrt) / B_paramboot)
  z_alpha_lower = 2 * hat_z_0 - qnorm(alpha/2)
  z_alpha_upper = 2 * hat_z_0 - qnorm(1-alpha/2)
  s_star_sorted = sort(s_star)
  theta_upper = s_star_sorted[max(1, as.integer(B_paramboot * pnorm(z_alpha_lower)))]
  theta_lower = s_star_sorted[max(1, as.integer(B_paramboot * pnorm(z_alpha_upper)))]
  intervalS4 = c(theta_lower, theta_upper)
  corrected_s <- s_star_sorted[max(1, as.integer(B_paramboot * pnorm(2 * hat_z_0)))]
  
  estimate_list4 <- c(corrected_m, corrected_s)
  
  # Efron's BC CI
  # for mean
  if(interval4[1]<=true_normal_mean & interval4[2]>=true_normal_mean){
    coverage4 = 1
  }
  else{coverage4 = 0}
  width4 = interval4[2]-interval4[1] 
  
  # for std
  if(intervalS4[1]<=true_normal_sd & intervalS4[2]>=true_normal_sd){
    coverageS4 = 1
  }
  else{coverageS4 = 0}
  widthS4 = intervalS4[2]-intervalS4[1]
  
  # automatic percentile 
  
  ## for mean: lower bound for CI
  mu_0 = quantile(m_star, alpha/2)
  for(i in 1:B_paramboot){
    data = generate_X_func(u[i,], m=mu_0, s=z_s_sqrt)
    data_clamp = clamp(data)
    z_star12 = mean(data_clamp) + w1[i] / n
    z_star22 = var(data_clamp) + w2[i] / n
    
    m_star5[i] = z_star12
    s_star5[i] = sqrt(max(z_star22,0))
  }
  mu_0_prime = quantile(m_star5, 1-alpha/2)
  for(i in 1:B_paramboot){
    data = generate_X_func(u[i,], m=mu_0_prime, s=z_s_sqrt)
    data_clamp = clamp(data)
    z_star12 = mean(data_clamp) + w1[i] / n
    z_star22 = var(data_clamp) + w2[i] / n
    
    m_star5[i] = z_star12
    s_star5[i] = sqrt(max(z_star22, 0))
  }
  m_star5_sorted = sort(m_star5)
  new_quantile = mean(m_star5_sorted <= mu_0)
  theta_lower = quantile(m_star, new_quantile)
  
  ## for mean: upper bound for CI
  mu_0 = quantile(m_star, 1-alpha/2)
  for(i in 1:B_paramboot){
    data = generate_X_func(u[i,], m=mu_0, s=z_s_sqrt)
    data_clamp = clamp(data)
    z_star12 = mean(data_clamp) + w1[i] / n
    z_star22 = var(data_clamp) + w2[i] / n
    
    m_star5[i] = z_star12
    s_star5[i] = sqrt(max(z_star22, 0))
  }
  mu_0_prime = quantile(m_star5, alpha/2)
  for(i in 1:B_paramboot){
    data = generate_X_func(u[i,], m=mu_0_prime, s=z_s_sqrt)
    data_clamp = clamp(data)
    z_star12 = mean(data_clamp) + w1[i] / n
    z_star22 = var(data_clamp) + w2[i] / n
    
    m_star5[i] = z_star12
    s_star5[i] = sqrt(max(z_star22, 0))
  }
  m_star5_sorted = sort(m_star5)
  new_quantile = mean(m_star5_sorted <= mu_0)
  theta_upper = quantile(m_star, new_quantile)
  
  # automatic percentile CI for mean
  interval5 = c(theta_lower, theta_upper)
  if(interval5[1]<=true_normal_mean & interval5[2]>=true_normal_mean){
    coverage5 = 1
  }
  else{coverage5 = 0}
  width5 = interval5[2]-interval5[1]
  
  # automatic percentile type of estimate for mean (0% CI)
  mu_0 = quantile(m_star, 0.5)
  w1 = rnorm(n=B_paramboot, m=0, s=(b-a)/ep)
  w2 = rnorm(n=B_paramboot, m=0, s=(b-a)^2/ep)
  for(i in 1:B_paramboot){
    data = generate_X_func(u[i,],m=mu_0,s=z_s_sqrt)
    data_clamp = clamp(data)
    z_star12 = mean(data_clamp) + w1[i] / n
    z_star22 = var(data_clamp) + w2[i] / n
    
    m_star5[i] = z_star12
    s_star5[i] = sqrt(max(z_star22, 0))
  }
  mu_0_prime = quantile(m_star5, 0.5)
  for(i in 1:B_paramboot){
    data = generate_X_func(u[i,], m=mu_0_prime, s=z_s_sqrt)
    data_clamp = clamp(data)
    z_star12 = mean(data_clamp) + w1[i] / n
    z_star22 = var(data_clamp) + w2[i] / n
    
    m_star5[i] = z_star12
    s_star5[i] = sqrt(max(z_star22,0))
  }
  m_star5_sorted = sort(m_star5)
  new_quantile = mean(m_star5_sorted <= mu_0)
  corrected_m = quantile(m_star, new_quantile)
  
  
  ## for std: lower bound for CI
  sa_0 = quantile(s_star, alpha/2)
  for(i in 1:R){
    data = generate_X_func(u[i,],m=z_m,s=sa_0)
    data_clamp = clamp(data)
    z_star12 = mean(data_clamp) + w1[i] / n
    z_star22 = var(data_clamp) + w2[i] / n
    
    m_star5[i] = z_star12
    s_star5[i] = sqrt(max(z_star22,0))
  }
  sa_0_prime = quantile(s_star5, 1-alpha/2)
  for(i in 1:R){
    data = generate_X_func(u[i,],m=z_m,s=sa_0_prime)
    data_clamp = clamp(data)
    z_star12 = mean(data_clamp) + w1[i] / n
    z_star22 = var(data_clamp) + w2[i] / n
    
    m_star5[i] = z_star12
    s_star5[i] = sqrt(max(z_star22,0))
  }
  s_star5_sorted = sort(s_star5)
  new_quantile = mean(s_star5_sorted <= sa_0)
  theta_lower = quantile(s_star, new_quantile)
  
  ## for std: upper bound for CI
  sa_0 = quantile(s_star, 1-alpha/2)
  for(i in 1:R){
    data = generate_X_func(u[i,],m=z_m,s=sa_0)
    data_clamp = clamp(data)
    z_star12 = mean(data_clamp) + w1[i] / n
    z_star22 = var(data_clamp) + w2[i] / n
    
    m_star5[i] = z_star12
    s_star5[i] = sqrt(max(z_star22,0))
  }
  sa_0_prime = quantile(s_star5, alpha/2)
  for(i in 1:R){
    data = generate_X_func(u[i,],m=z_m,s=sa_0_prime)
    data_clamp = clamp(data)
    z_star12 = mean(data_clamp) + w1[i] / n
    z_star22 = var(data_clamp) + w2[i] / n
    
    m_star5[i] = z_star12
    s_star5[i] = sqrt(max(z_star22,0))
  }
  s_star5_sorted = sort(s_star5)
  new_quantile = mean(s_star5_sorted <= sa_0)
  theta_upper = quantile(s_star, new_quantile)
  
  # automatic percentile CI for std
  intervalS5 = c(theta_lower, theta_upper)
  if(intervalS5[1]<=true_normal_sd & intervalS5[2]>=true_normal_sd){
    coverageS5 = 1
  }
  else{coverageS5 = 0}
  widthS5 = intervalS5[2]-intervalS5[1]
  
  # automatic percentile type of estimate for std (0% CI)
  sa_0 = quantile(s_star, 0.5)
  w1 = rnorm(n=R,m=0,s=(b-a)/ep)
  w2 = rnorm(n=R,m=0,s=(b-a)^2/ep)
  for(i in 1:R){
    data = generate_X_func(u[i,],m=z_m,s=sa_0)
    data_clamp = clamp(data)
    z_star12 = mean(data_clamp) + w1[i] / n
    z_star22 = var(data_clamp) + w2[i] / n
    
    m_star5[i] = z_star12
    s_star5[i] = sqrt(max(z_star22,0))
  }
  sa_0_prime = quantile(s_star5, 0.5)
  for(i in 1:R){
    data = generate_X_func(u[i,],m=z_m,s=sa_0_prime)
    data_clamp = clamp(data)
    z_star12 = mean(data_clamp) + w1[i] / n
    z_star22 = var(data_clamp) + w2[i] / n
    
    m_star5[i] = z_star12
    s_star5[i] = sqrt(max(z_star22,0))
  }
  s_star5_sorted = sort(s_star5)
  new_quantile = mean(s_star5_sorted <= sa_0)
  corrected_s = quantile(s_star, new_quantile)
  
  # automatic percentile type of estimate
  estimate_list5 <- c(corrected_m, corrected_s)

  ############ 
  # all coverages, widths, and estimates
  c(coverage, coverage2, coverage3, coverage4, coverage5,
    coverageS, coverageS2, coverageS3, coverageS4, coverageS5,
    width, width2, width3, width4, width5, 
    widthS, widthS2, widthS3, widthS4, widthS5, 
    estimate_list, estimate_list2, estimate_list3, estimate_list4, estimate_list5)
}

num_est = 5
coverage = mean(results[,1])
coverage2 = mean(results[,2]) 
coverage3 = mean(results[,3])
coverage4 = mean(results[,4])
coverage5 = mean(results[,5])

coverageS = mean(results[,num_est+1])
coverageS2 = mean(results[,num_est+2]) 
coverageS3 = mean(results[,num_est+3])
coverageS4 = mean(results[,num_est+4])
coverageS5 = mean(results[,num_est+5])

width = results[,2*num_est+1]
width2 = results[,2*num_est+2]
width3 = results[,2*num_est+3]
width4 = results[,2*num_est+4]
width5 = results[,2*num_est+5]

widthS = results[,3*num_est+1]
widthS2 = results[,3*num_est+2]
widthS3 = results[,3*num_est+3]
widthS4 = results[,3*num_est+4]
widthS5 = results[,3*num_est+5]

estimate_list = results[,c(4*num_est+1, 4*num_est+2)]
estimate_list2 = results[,c(4*num_est+3, 4*num_est+4)]
estimate_list3 = results[,c(4*num_est+5, 4*num_est+6)]
estimate_list4 = results[,c(4*num_est+7, 4*num_est+8)]
estimate_list5 = results[,c(4*num_est+9, 4*num_est+10)]

### Parametric Bootstrap (percentile) for mu (mean)
coverage_se <- sqrt(coverage*(1-coverage)/nSIM)
c(mean(width), sqrt(var(width)/nSIM))

### Parametric Bootstrap (percentile) for sigma (standard deviation)
coverageS_se <- sqrt(coverageS*(1-coverageS)/nSIM)
c(mean(widthS), sqrt(var(widthS)/nSIM))


### Parametric Bootstrap (simplified t) for mu (mean)
c(coverage2, sqrt(coverage2*(1-coverage2)/nSIM))
c(mean(width2), sqrt(var(width2)/nSIM))

### Parametric Bootstrap (simplified t) for sigma (standard deviation)
c(coverageS2, sqrt(coverageS2*(1-coverageS2)/nSIM))
c(mean(widthS2), sqrt(var(widthS2)/nSIM))


### Parametric Bootstrap (ferrando) for mu (mean)
c(coverage3, sqrt(coverage3*(1-coverage3)/nSIM))
c(mean(width3), sqrt(var(width3)/nSIM))

### Parametric Bootstrap (ferrando) for sigma (standard deviation)
c(coverageS3, sqrt(coverageS3*(1-coverageS3)/nSIM))
c(mean(widthS3), sqrt(var(widthS3)/nSIM))


### Parametric Bootstrap (Efron BC) for mu (mean)
c(coverage4, sqrt(coverage4*(1-coverage4)/nSIM))
c(mean(width4), sqrt(var(width4)/nSIM))

### Parametric Bootstrap (Efron BC) for sigma (standard deviation)
c(coverageS4, sqrt(coverageS4*(1-coverageS4)/nSIM))
c(mean(widthS4), sqrt(var(widthS4)/nSIM))

### Parametric Bootstrap (automatic percentile) for mu (mean)
c(coverage5, sqrt(coverage5*(1-coverage5)/nSIM))
c(mean(width5), sqrt(var(width5)/nSIM))

### Parametric Bootstrap (automatic percentile) for sigma (standard deviation)
c(coverageS5, sqrt(coverageS5*(1-coverageS5)/nSIM))
c(mean(widthS5), sqrt(var(widthS5)))

method_stat <- data.frame(coverage_mu=c(coverage, coverage2, coverage3, coverage4, coverage5),
                          coverage_mu_se=sapply(c(coverage, coverage2, coverage3, coverage4, coverage5), function(x) sqrt(x*(1-x)/nSIM)),
                          coverage_sigma=c(coverageS, coverageS2, coverageS3, coverageS4, coverageS5),
                          coverage_sigma_se=sapply(c(coverageS, coverageS2, coverageS3, coverageS4, coverageS5), function(x) sqrt(x*(1-x)/nSIM)),
                          width_mu=c(mean(width), mean(width2), mean(width3), mean(width4), mean(width5)),
                          width_mu_se=c(sd(width), sd(width2), sd(width3), sd(width4), sd(width5))/sqrt(nSIM),
                          width_sigma=c(mean(widthS), mean(widthS2), mean(widthS3), mean(widthS4), mean(widthS5)),
                          width_sigma_se=c(sd(widthS), sd(widthS2), sd(widthS3), sd(widthS4), sd(widthS5))/sqrt(nSIM)
                          )
method_names <- c('naive','simplified t','Ferrando et al', 'Efron\'s BC', 'automatic percentile')
row.names(method_stat) <- method_names
write.table(method_stat, "results/paramboot_normal_comparison_results.csv", sep=",", row.names=TRUE, col.names=TRUE)

method_names <- c('naive','simplified t', 'Efron\'s BC', 'automatic percentile')
estimate_lists <- data.frame(cbind(rbind(estimate_list, estimate_list2, estimate_list4, estimate_list5),
                                   rep(method_names, each=nSIM)))
colnames(estimate_lists) <- c('mean', 'std', 'method')
write.table(estimate_lists, "results/paramboot_normal_comparison_estimates.csv", sep=",", row.names=FALSE, col.names=TRUE)



stopCluster(cl)
